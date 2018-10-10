//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4bool PrimaryGeneratorAction::isSpectrum=false;
double* PrimaryGeneratorAction::ws_array;//=new double[0];
int* PrimaryGeneratorAction::ibin_array;//=new int[0];
double* PrimaryGeneratorAction::en_array;
int PrimaryGeneratorAction::ndata;
G4double PrimaryGeneratorAction::incidentEnergy=0.100*MeV;
double PrimaryGeneratorAction::kv=0.0;
double PrimaryGeneratorAction::filtration=0.0;
double PrimaryGeneratorAction::anodAngle=0.0;
double PrimaryGeneratorAction::ripple=0.0;
std::string PrimaryGeneratorAction::anodMaterial="";


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();  
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(50.*MeV);
  G4double position = -0.5*(Detector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));
  
  rndmFlag = "off";
  isSpectrum=true;//to be handled via run.mac

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PrimaryGeneratorAction::IsSpectrum()
{
	return isSpectrum;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double x0 = -0.5*(Detector->GetWorldSizeX());
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (Detector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
      z0 = (Detector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
     } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  
  //=====================
  if (isSpectrum){
	  incidentEnergy = aliasSample();
	  particleGun->SetParticleEnergy(incidentEnergy);
  } else{
	  incidentEnergy = particleGun->GetParticleEnergy();
  }
  //==========================
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PrimaryGeneratorAction::aliasSample(){
	   //(int nsbin, double[] xs_array,// )
	   // double[] ws_array, int[] ibin_array) {
		// "===============================================================
		// "
		// " samples from an alias table which must have been prepared
		// " using prepare_alias_table
		// "
		// "===============================================================

		// ;Copyright NRC;
		// implicit none;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL xs_array(0:nsbin),ws_array(nsbin);
	    double alias_sample = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		double aj = 0.0;
		int j = 0;
		int nsbin=ndata;
		v1 = G4UniformRand();//random01();
		v2 = G4UniformRand();//random01();
		aj = 1.0 + v1 * nsbin;
		//Double dbl = new Double(aj);
		//j = dbl.intValue();// minim1 maxim nsbin!!
		j=aj;
		if (j > nsbin)
			j = nsbin; // " this happens only if $RANDOMSET produces
						// " numbers in (0,1]--------> is not the DEFAULT
						// case!!!
		aj = aj - j;
		if (aj > ws_array[j - 1]) {
			j = ibin_array[j - 1];
		}
		alias_sample = (1.0 - v2) * en_array[j - 1] + v2 * en_array[j];// ok, xs=0 biased

		alias_sample=alias_sample*MeV;

		return alias_sample;
}

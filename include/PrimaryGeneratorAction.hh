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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  virtual ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}

  G4ParticleGun* GetParticleGun() const { return particleGun; }
  //=====================
  static void SetIfSpectrum(G4String val){
		if (val=="yes")
		 PrimaryGeneratorAction::isSpectrum=true;
		else
		 PrimaryGeneratorAction::isSpectrum=false;
	};
  static G4bool IsSpectrum();
  static G4bool isSpectrum;
  static double* ws_array;
  static int* ibin_array;
  static double* en_array;
  static int ndata;
  
  static void SetSpectrumAliasSamplingData(double* ws_a, int* ibin_a, double*en_a, int ndat){
		ws_array=ws_a;
		ibin_array=ibin_a;
		en_array=en_a;
		ndata=ndat;
  }
  static G4double incidentEnergy;
  static double kv;
  static void SetKv(std::string val){std::istringstream buffer(val); buffer>>kv;};
  
  static double filtration;
  static void SetFiltration(std::string val){std::istringstream buffer(val); buffer>>filtration;};
  
  static double anodAngle;
  static void SetAnodAngle(std::string val){std::istringstream buffer(val); buffer>>anodAngle;};
  
  static double ripple;
  static void SetRipple(std::string val){std::istringstream buffer(val); buffer>>ripple;};
	
  static std::string anodMaterial;
  static void SetAnodMaterial(std::string val){anodMaterial=val;};

private:
  G4ParticleGun*           particleGun;  //pointer a to G4  class
  DetectorConstruction*    Detector;     //pointer to the geometry
    
  PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  G4String                   rndmFlag;     //flag for a rndm impact point

  G4double aliasSample();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



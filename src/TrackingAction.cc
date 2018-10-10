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
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "RunAction.hh"
//#include "HistoManager.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()//RunAction* RuAct)
//:fRunAction(RuAct),NbOfGamma(1)
{ 
	fRunAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
	eventaction = (EventAction*)
                G4RunManager::GetRunManager()->GetUserEventAction(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{ 
  //fEdepCavity = 0.;  
 if(aTrack->GetDefinition()!= G4Gamma::Gamma()){//charged particle

	 G4ThreeVector position = aTrack->GetPosition(); 
     G4double vx = position.x();//where we are

	 G4double kev=aTrack->GetKineticEnergy();///kinetic energy of particle;
	 G4int pID = aTrack->GetParentID();//ID of parent particle

	 if (pID==1){//it is a secondary particle
		 //parentID=1 comes from primary gamma=>save kinetic energy for KERMA!!!

		 //G4cout << G4endl<<"Number of e: "<<NbOfGamma<<"; energy [keV]: "<<kev
		 //<<"; x [cm]: "<<vx<<"; parentId "<<pID<< G4endl;
		 //===================
		 //SCORE IT for EACH EVENT (statistic reasons)!!!!!!
		 eventaction->AddKinetic(kev,vx);
		 //=======================

		//NbOfGamma++;
	 }
	 
	 
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*)
{
  //sum energy in cavity
  //
  //if (fEdepCavity > 0.) {
    //fRunAction->AddEdepCavity(fEdepCavity);
    //fHistoManager->FillHisto(11,fEdepCavity);
  //}  
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


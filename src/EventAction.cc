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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{
	fDetector=(DetectorConstruction*)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  eventMessenger = new EventActionMessenger(this);
  printModulo = 100;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
  delete [] dAbsArray;//// When done, free memory pointed to by kAbsArray.
  delete [] dGapArray;
  //delete [] dAbsArray2;//NOT HERE, AFTER ALL EVENT ARE PROCESSED=>in end of run!!
  //delete [] dGapArray2;

  delete [] EkinAbsArray;
  delete [] EkinGapArray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb
		<<" ;Incident energy: "<<G4BestUnit(PrimaryGeneratorAction::incidentEnergy,"Energy")<< G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }
 
 // initialisation per event
 EnergyAbs = EnergyGap = 0.;
 TrackLAbs = TrackLGap = 0.;

 EkinGap = 0.0;
 EkinAbs = 0.0;

 G4int n = fDetector->GetNbOfLayers();
 NreplicaEach=n;
 absorberThickness=fDetector->GetAbsorberThickness();
 gapThickness=fDetector->GetGapThickness();
 boundaryAG=(-gapThickness+absorberThickness)*NreplicaEach/2.0;

 dAbsArray = new G4double[n];
 dGapArray = new G4double[n];
 //dAbsArray2 = new G4double[n];
 //dGapArray2 = new G4double[n];
 EkinAbsArray = new G4double[n];
 EkinGapArray = new G4double[n];
 for (int i=0; i<n; i++){
	  dAbsArray[i]=0.0;
	  dGapArray[i]=0.0;
	  //dAbsArray2[i]=0.0;
	  //dGapArray2[i]=0.0;

	  EkinAbsArray[i]=0.0;
	  EkinGapArray[i]=0.0;
 }
}

void EventAction::AddKinetic(G4double de, G4double xposition) {
	 
	 if (xposition<=boundaryAG){
		  //we are in absorber
		EkinAbs += de;
		//=============now the layers=================
		G4double lowerBoundary=-NreplicaEach*(gapThickness+absorberThickness)/2.0;//start
		G4double higherBoundary=lowerBoundary+absorberThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				EkinAbsArray[0]+=de;		
		}

		for (int i=1;i<NreplicaEach;i++){
			lowerBoundary=lowerBoundary+absorberThickness;
			higherBoundary=lowerBoundary+absorberThickness;
			if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				EkinAbsArray[i]+=de;
				break;//optimization
			}
		}
		
		//============================================
	  } else {
		//We are in gap
		EkinGap += de;

		G4double lowerBoundary=boundaryAG;//start
		G4double higherBoundary=lowerBoundary+gapThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				EkinGapArray[0]+=de;
		}

		for (int i=1;i<NreplicaEach;i++){
			lowerBoundary=lowerBoundary+gapThickness;
			higherBoundary=lowerBoundary+gapThickness;
			if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				EkinGapArray[i]+=de;
				break;//optimization
			}
		}
	  }
	  
}

void EventAction::AddAbs(G4double de, G4double dl, G4double xposition) {
	EnergyAbs += de;
	TrackLAbs += dl;

	G4double absorberThickness=fDetector->GetAbsorberThickness();
	G4double gapThickness=fDetector->GetGapThickness();
	G4int NreplicaEach=fDetector->GetNbOfLayers();

	//=============now the layers=================
	G4double lowerBoundary=-NreplicaEach*(gapThickness+absorberThickness)/2.0;//start
	G4double higherBoundary=lowerBoundary+absorberThickness;
	if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
			dAbsArray[0]+=de;
			//dAbsArray2[0]+=de*de;
	}
		
	for (int i=1;i<NreplicaEach;i++){
		lowerBoundary=lowerBoundary+absorberThickness;
		higherBoundary=lowerBoundary+absorberThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
			dAbsArray[i]+=de;
			//dAbsArray2[i]+=de*de;
			break;//optimization
		}
	}
}

void EventAction::AddGap(G4double de, G4double dl, G4double xposition) {
	EnergyGap += de;
	TrackLGap += dl;

	G4double absorberThickness=fDetector->GetAbsorberThickness();
	G4double gapThickness=fDetector->GetGapThickness();
	G4int NreplicaEach=fDetector->GetNbOfLayers();
	G4double boundaryAGcm=(-gapThickness+absorberThickness)*NreplicaEach/2.0;///cm;

	//the layers=================
	G4double lowerBoundary=boundaryAGcm;//start
	G4double higherBoundary=lowerBoundary+gapThickness;
	if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
			dGapArray[0]+=de;
			//dGapArray2[0]+=de*de;
	}

	for (int i=1;i<NreplicaEach;i++){
		lowerBoundary=lowerBoundary+gapThickness;
		higherBoundary=lowerBoundary+gapThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
			dGapArray[i]+=de;
			//dGapArray2[i]+=de*de;
			break;//optimization
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  //
  runAct->fillPerEvent(EnergyAbs, EnergyGap, TrackLAbs, TrackLGap, dAbsArray, dGapArray);//,dAbsArray2, dGapArray2);
  runAct->fillEkinPerEvent(EkinAbs, EkinGap, EkinAbsArray, EkinGapArray);
  
  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;        

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLGap,"Length")
       << G4endl;
          
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

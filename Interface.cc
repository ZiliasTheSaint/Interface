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
//// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "TrackingAction.hh"//@@@@@@@@@@@@@@@@@@@
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  
  DetectorConstruction* det;
  PhysicsList* phys;

  runManager->SetUserInitialization(det  = new DetectorConstruction);//new DetectorConstruction);
  //
  runManager->SetUserInitialization(phys = new PhysicsList(det));//new PhysicsList);
    
  #ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif
  // Set user action classes
  //RunAction* run        = new RunAction();  
  //TrackingAction* track = new TrackingAction();//(run); 

  runManager->SetUserAction(new PrimaryGeneratorAction);
  //
  runManager->SetUserAction(new RunAction);//run);//new RunAction);//
  //runManager->SetUserAction(new TrackingAction);//track);//new RunAction);//
  runManager->SetUserAction(new EventAction);//AFTER RUN bcause from here runAction is called!!!
  runManager->SetUserAction(new TrackingAction);//AFTER EVENT OTHERWISE ERROR AT RUNTIME//track);//new RunAction);//  //
  runManager->SetUserAction(new SteppingAction);//AFTER EVENT
  
  // Initialize G4 kernel
  //runManager->Initialize();
  //INITIALIZATION IS CALLED FROM run.mac AFTER PHANTOM IS SET!!!OR NOT=>updateGeometry
  //The vis.mac for vizualization is called from run.mac also!!!
  
  // Initialize G4 kernel  
  //runManager->Initialize();
  


  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {        
      #ifdef G4UI_USE
             G4cout << " UI session starts ..." << G4endl;
             G4UIExecutive* ui = new G4UIExecutive(argc, argv);
             UImanager->ApplyCommand("/control/execute runInterf.mac");     
             ui->SessionStart();
             delete ui;
	 #endif
      //////////////////////////////////#######################
      
//#endif
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

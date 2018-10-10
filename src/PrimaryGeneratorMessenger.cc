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

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                          PrimaryGeneratorAction* Gun)
:Action(Gun)
{
  gunDir = new G4UIdirectory("/N03/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  RndmCmd = new G4UIcmdWithAString("/N03/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //=============================================
  isSpectrumCmd = new G4UIcmdWithAString("/xfield/isSpectrum?",this);
  isSpectrumCmd->SetGuidance("Set yes if we have an energy spectrum, and no otherwise.");
  isSpectrumCmd->SetGuidance("  Choice : yes(default), no");
  isSpectrumCmd->SetParameterName("choice",true);
  isSpectrumCmd->SetDefaultValue("yes");
  isSpectrumCmd->SetCandidates("yes no");
  isSpectrumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  kvCmd = new G4UIcmdWithAString("/xfield/kVp",this);
  kvCmd->SetGuidance("Set Xray tube kilovoltage.");  
  kvCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  filtrationCmd = new G4UIcmdWithAString("/xfield/filtration",this);
  filtrationCmd->SetGuidance("Set Xray tube filtration [mmAl].");  
  filtrationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  anodMaterialCmd = new G4UIcmdWithAString("/xfield/anodMaterial",this);
  anodMaterialCmd->SetGuidance("Set Xray tube anodMaterial.");  
  anodMaterialCmd->SetCandidates("W MO RH");
  anodMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  anodAngleCmd = new G4UIcmdWithAString("/xfield/anodAngle",this);
  anodAngleCmd->SetGuidance("Set Xray tube anode angle [deg].");  
  anodAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rippleCmd = new G4UIcmdWithAString("/xfield/ripple",this);
  rippleCmd->SetGuidance("Set Xray tube waveform ripple [%].");  
  rippleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete gunDir;

  delete isSpectrumCmd;
  delete kvCmd;
  delete filtrationCmd;
  delete anodAngleCmd;
  delete rippleCmd;
  delete anodMaterialCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == RndmCmd )
   { Action->SetRndmFlag(newValue);}

  if( command == isSpectrumCmd )
   { Action->SetIfSpectrum(newValue);
   }

	if( command == kvCmd )
   { Action->SetKv(newValue);
   }

	if( command == filtrationCmd )
   { Action->SetFiltration(newValue);
   }

	if( command == anodMaterialCmd )
   { Action->SetAnodMaterial(newValue);
   }

	if( command == anodAngleCmd )
   { Action->SetAnodAngle(newValue);
   }

	if( command == rippleCmd )
   { Action->SetRipple(newValue);
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


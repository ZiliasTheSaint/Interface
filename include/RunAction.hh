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

#ifndef RunAction_h
#define RunAction_h 1

#include "ProcessesCount.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

//class DetectorConstruction;
//class PrimaryGeneratorAction;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class XRayBuilder;
class RunAction : public G4UserRunAction
{
public:
  RunAction();
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
  void fillPerEvent(G4double, G4double, G4double, G4double, G4double*, G4double*);//,G4double*, G4double*); 
  
  //========================================
  void fillEkinPerEvent(G4double, G4double, G4double*, G4double*); 
  //void AddEkin(G4double, G4double);
  
  G4double basicAnalyze(G4double x, G4double x2, G4int n){
	  G4double temp=x;G4double temp2=x2;G4int nevents=n;
	  temp=temp/nevents;temp2=temp2/nevents;
	  temp2 = temp2 - temp*temp;
	  if (nevents>1) temp2=temp2/(nevents-1.0);
	  if (temp2 >0.)
	   temp2 = sqrt(temp2); 
	  //else temp2 = 99.99;//never!
	  temp2=std::abs(temp2);
	  //=========percent
	  if (temp!=0.0){
	   temp2 = std::min(100.0*temp2/temp,99.9);
      } //else temp2 = 99.9;//no score means no score not necessarly an error!
	  return temp2;
  }

private:
	DetectorConstruction*   fDetector;//pointer to geometry
    PrimaryGeneratorAction* fKinematic;//pointer to source
    ProcessesCount*         fProcCounter;

  G4double sumEAbs, sum2EAbs;
  G4double sumEGap, sum2EGap;
    
  G4double sumLAbs, sum2LAbs;
  G4double sumLGap, sum2LGap;    

   G4double                fEkinAbs, fEkinAbs2;
   G4long                  fNbEventAbs; 

   G4double                fEkinGap, fEkinGap2;
   G4long                  fNbEventGap; 

   G4double* kAbsArray;//kerma (energy) in Absorber
   G4double* kGapArray;
   G4double* ddAbsArray;//dose (energy) in Absorber
   G4double* ddGapArray;

   G4double* kAbsArray2;//kerma (energy) in Absorber
   G4double* kGapArray2;
   G4double* ddAbsArray2;//dose (energy) in Absorber
   G4double* ddGapArray2;

   void prepareAliasSampling(int nsbin, std::vector<double> fs_array);
	void prepareAliasSampling(int nsbin, double* fs_array);
	int ndata;
	double* enV_array;
	double* ws_array;
	int* ibin_array;
	double photonsPerDAP;
	double kermaPer_mas_at75cm;
	bool dummy;
	XRayBuilder* xraybuilder;

	double linInt(double x1, double y1, double x2, double y2,double x);
	int findNearestDataIndexFromArray(double* a, int arraySize, double value, bool lowerThanValue);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


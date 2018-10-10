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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
//#include "XRayBuilder.hh"
#include"XRAYBuilder.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
	// add new units for dose 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);   

	fDetector=(DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();

	fKinematic=(PrimaryGeneratorAction*)
                G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
	delete [] kAbsArray;//// When done, free memory pointed to by kAbsArray.
	delete [] kGapArray;

	delete [] ddAbsArray;//// When done, free memory pointed to by kAbsArray.
	delete [] ddGapArray;

	delete [] kAbsArray2;//// When done, free memory pointed to by kAbsArray.
	delete [] kGapArray2;

	delete [] ddAbsArray2;//// When done, free memory pointed to by kAbsArray.
	delete [] ddGapArray2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //initialize cumulative quantities
  //
  sumEAbs = sum2EAbs =sumEGap = sum2EGap = 0.;
  sumLAbs = sum2LAbs =sumLGap = sum2LGap = 0.; 

  fEkinGap = 0.0;fEkinGap2=0.0;
  fEkinAbs = 0.0;fEkinAbs2=0.0;
  fNbEventAbs=0;
  fNbEventGap=0;

  G4int n = fDetector->GetNbOfLayers();
  kAbsArray = new G4double[n];
  kGapArray = new G4double[n];
  ddAbsArray = new G4double[n];
  ddGapArray = new G4double[n];

  kAbsArray2 = new G4double[n];
  kGapArray2 = new G4double[n];
  ddAbsArray2 = new G4double[n];
  ddGapArray2 = new G4double[n];
  for (int i=0; i<n; i++){
	  kAbsArray[i]=0.0;
	  kGapArray[i]=0.0;
	  ddAbsArray[i]=0.0;
	  ddGapArray[i]=0.0;

	  kAbsArray2[i]=0.0;
	  kGapArray2[i]=0.0;
	  ddAbsArray2[i]=0.0;
	  ddGapArray2[i]=0.0;
  }

  if(PrimaryGeneratorAction::IsSpectrum()){
	xraybuilder = new XRayBuilder(
	PrimaryGeneratorAction::kv,//80.0
	PrimaryGeneratorAction::filtration,//2.5
	PrimaryGeneratorAction::anodAngle,//17,"W"
	PrimaryGeneratorAction::ripple,//0
	PrimaryGeneratorAction::anodMaterial//"W"
	);
	
	enV_array=xraybuilder->GetEnergyArray();
	ndata = xraybuilder->GetBins();
	double* prob=xraybuilder->GetProbArray();
	photonsPerDAP=xraybuilder->Get_photons_per_mm2_per_uGy();
	kermaPer_mas_at75cm=xraybuilder->GetAirKerma_per_mAs_at75cm();
	prepareAliasSampling(ndata,prob);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EAbs, G4double EGap,
                                  G4double LAbs, G4double LGap, G4double* dAbsArray, G4double* dGapArray)//, G4double* dAbsArray2, G4double* dGapArray2)
{
  //accumulate statistic
  //
  sumEAbs += EAbs;  sum2EAbs += EAbs*EAbs;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLAbs += LAbs;  sum2LAbs += LAbs*LAbs;
  sumLGap += LGap;  sum2LGap += LGap*LGap;  

  G4int n = fDetector->GetNbOfLayers();
  for (int i=0; i<n; i++){
	  ddAbsArray[i]=ddAbsArray[i]+dAbsArray[i];
	  ddGapArray[i]=ddGapArray[i]+dGapArray[i];

	  ddAbsArray2[i]=ddAbsArray2[i]+dAbsArray[i]*dAbsArray[i];//dAbsArray2[i];
	  ddGapArray2[i]=ddGapArray2[i]+dGapArray[i]*dGapArray[i];//dGapArray2[i];
 }

}

void RunAction::fillEkinPerEvent(G4double kinAbs, G4double kinGap, G4double* kinAbsArray, G4double* kinGapArray) { 
	
	G4int NreplicaEach=fDetector->GetNbOfLayers();
	
	fEkinAbs=fEkinAbs+kinAbs;fEkinAbs2=fEkinAbs2+kinAbs*kinAbs;
	fEkinGap=fEkinGap+kinGap;fEkinGap2=fEkinGap2+kinGap*kinGap;
	
	for (int i=0;i<NreplicaEach;i++){
		kAbsArray[i]=kAbsArray[i]+kinAbsArray[i];
		kAbsArray2[i]=kAbsArray2[i]+kinAbsArray[i]*kinAbsArray[i];
		
		kGapArray[i]=kGapArray[i]+kinGapArray[i];
		kGapArray2[i]=kGapArray2[i]+kinGapArray[i]*kinGapArray[i];
	}
}

/*void RunAction::AddEkin(G4double de, G4double xposition) { 
	  //decide wheter in gap or absorber
	//Left entrance is at x = -N(g+a)/2. Boundary is at x = -N(g+a)/2+Na = 1/2N(-g+a)
	//=====================================
	  G4double absorberThickness=fDetector->GetAbsorberThickness();
	  G4double gapThickness=fDetector->GetGapThickness();
	  G4int NreplicaEach=fDetector->GetNbOfLayers();

	  G4double boundaryAGcm=(-gapThickness+absorberThickness)*NreplicaEach/2.0;///cm;

	  if (xposition<=boundaryAGcm){
		  //we are in absorber
		fEkinAbs += de;
		fEkinAbs2 += de*de;
		fNbEventAbs++;  

		//=============now the layers=================
		G4double lowerBoundary=-NreplicaEach*(gapThickness+absorberThickness)/2.0;//start
		G4double higherBoundary=lowerBoundary+absorberThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				kAbsArray[0]+=de;
				kAbsArray2[0]+=de*de;
		}

		//lowerBoundary=-NreplicaEach*(gapThickness+absorberThickness)/2.0+absorberThickness;//start
		//higherBoundary=lowerBoundary+absorberThickness;
		//if ((lowerBoundary<=xposition) && (xposition<higherBoundary)){
		//		kAbsArray[1]+=de;
		//}
		for (int i=1;i<NreplicaEach;i++){
			lowerBoundary=lowerBoundary+absorberThickness;
			higherBoundary=lowerBoundary+absorberThickness;
			if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				kAbsArray[i]+=de;
				kAbsArray2[i]+=de*de;
				break;//optimization
			}
		}
		
		//============================================
	  } else {

		fEkinGap += de;
		fEkinGap2 += de*de;
		fNbEventGap++;  

		G4double lowerBoundary=boundaryAGcm;//start
		G4double higherBoundary=lowerBoundary+gapThickness;
		if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				kGapArray[0]+=de;
				kGapArray2[0]+=de*de;
		}

		for (int i=1;i<NreplicaEach;i++){
			lowerBoundary=lowerBoundary+gapThickness;
			higherBoundary=lowerBoundary+gapThickness;
			if ((lowerBoundary<xposition) && (xposition<=higherBoundary)){
				kGapArray[i]+=de;
				kGapArray2[i]+=de*de;
				break;//optimization
			}
		}
	  }
  };*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4double massAbsorberLayer = fDetector->GetAbsorber_logical()->GetMass();
  G4double massGapLayer = fDetector->GetGap_logical()->GetMass();
  //compute statistics: mean and rms
  //
  G4double err = basicAnalyze(sumEAbs, sum2EAbs, NbOfEvents);
  sumEAbs /= NbOfEvents; sum2EAbs /= NbOfEvents;
  G4double rmsEAbs = err;//sum2EAbs - sumEAbs*sumEAbs;
  //if (rmsEAbs >0.) rmsEAbs = std::sqrt(rmsEAbs); else rmsEAbs = 0.;
  
  err = basicAnalyze(sumEGap, sum2EGap, NbOfEvents);
  sumEGap /= NbOfEvents; sum2EGap /= NbOfEvents;
  G4double rmsEGap = err;//sum2EGap - sumEGap*sumEGap;
  //if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;
  
  err = basicAnalyze(sumLAbs, sum2LAbs, NbOfEvents);
  sumLAbs /= NbOfEvents; sum2LAbs /= NbOfEvents;
  G4double rmsLAbs = err;//sum2LAbs - sumLAbs*sumLAbs;
  //if (rmsLAbs >0.) rmsLAbs = std::sqrt(rmsLAbs); else rmsLAbs = 0.;
  
  err = basicAnalyze(sumLGap, sum2LGap, NbOfEvents);
  sumLGap /= NbOfEvents; sum2LGap /= NbOfEvents;
  G4double rmsLGap = err;//sum2LGap - sumLGap*sumLGap;
  //if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  G4int n = fDetector->GetNbOfLayers();
  G4double doseAbsorber=sumEAbs/(massAbsorberLayer*n);
  G4double doseGap=sumEGap/(massGapLayer*n);
  //print
  //
  
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Absorber : " << G4BestUnit(sumEAbs,"Energy")<<" ; dose : "<<G4BestUnit(doseAbsorber,"Dose")
     << " +- "                          << rmsEAbs<<" %"//G4BestUnit(rmsEAbs,"Energy")  
     << "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")<<" ; dose : "<<G4BestUnit(doseGap,"Dose")
     << " +- "                          << rmsEGap<<" %"//G4BestUnit(rmsEGap,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in Absorber : " << G4BestUnit(sumLAbs,"Length")
     << " +- "                               << rmsLAbs<<" %"//G4BestUnit(rmsLAbs,"Length")  
     << "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     << " +- "                               << rmsLGap<<" %"//G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;

  G4double test=0.0;//ok
  G4double doseInLayer=0.0;
  for (int i=0; i<n; i++){
	     err = basicAnalyze(ddAbsArray[i], ddAbsArray2[i], NbOfEvents);	

		ddAbsArray[i]= ddAbsArray[i]/NbOfEvents;
		doseInLayer=ddAbsArray[i]/(massAbsorberLayer);
	test=test+ddAbsArray[i];
		G4cout
		 << "\n Energy deposited in Absorber layer: " <<i<<"; "<< G4BestUnit(ddAbsArray[i],"Energy")<<" ; dose: "<<G4BestUnit(doseInLayer,"Dose")//<<" test "<<G4BestUnit(test,"Energy")
		  <<" +- "<<err<<" %"<< G4endl;
	}

  test=0.0;//some remains!!!!
   for (int i=0; i<n; i++){
	    err = basicAnalyze(ddGapArray[i], ddGapArray2[i], NbOfEvents);	

		ddGapArray[i]= ddGapArray[i]/NbOfEvents;
		doseInLayer=ddGapArray[i]/(massGapLayer);
    test=test+ddGapArray[i];
		G4cout
		 << "\n Energy deposited in Gap layer: " <<i<<"; "<< G4BestUnit(ddGapArray[i],"Energy")<<" ; dose: "<<G4BestUnit(doseInLayer,"Dose")//<<" test "<<G4BestUnit(test,"Energy")
		  <<" +- "<<err<<" %"<< G4endl;
	}
  //===========
  fDetector->PrintCalorParameters();
  //==================================
  err = basicAnalyze(fEkinAbs, fEkinAbs2, NbOfEvents);//fNbEventAbs);//NbOfEvents);//
  fEkinAbs /= NbOfEvents;//fNbEventAbs;//NbOfEvents;//
  fEkinAbs2 /= NbOfEvents;
  G4double rmsfEkinAbs = err;//fEkinAbs2 - fEkinAbs*fEkinAbs;
  //if (rmsfEkinAbs >0.) rmsfEkinAbs = std::sqrt(rmsfEkinAbs); else rmsfEkinAbs = 0.;

  err = basicAnalyze(fEkinGap, fEkinGap2, NbOfEvents);
  fEkinGap /= NbOfEvents; fEkinGap2 /= NbOfEvents;
  G4double rmsfEkinGap = err;//fEkinGap2 - fEkinGap*fEkinGap;
  //if (rmsfEkinGap >0.) rmsfEkinGap = std::sqrt(rmsfEkinGap); else rmsfEkinGap = 0.;

  G4double kermaAbsorber=fEkinAbs/(massAbsorberLayer*n);
  G4double kermaGap=fEkinGap/(massGapLayer*n);

    G4cout
     << "\n Energy released to charged particles in Absorber : " << G4BestUnit(fEkinAbs,"Energy")<<" ; kerma= "<< G4BestUnit(kermaAbsorber,"Dose")
     << " +- "                               << rmsfEkinAbs<<" %" //<<" from Nevents: "<<NbOfEvents//G4BestUnit(rmsfEkinAbs,"Energy")  
     << "\n Energy relesed to charged particles in Gap      : " << G4BestUnit(fEkinGap,"Energy")<<" ; kerma= "<< G4BestUnit(kermaGap,"Dose")
     << " +- "                               << rmsfEkinGap<<" %"//G4BestUnit(rmsfEkinGap,"Energy")
     << "\n------------------------------------------------------------\n"
     << G4endl;

	
	 test=0.0;//ok
	for (int i=0; i<n; i++){
	    err = basicAnalyze(kAbsArray[i], kAbsArray2[i], NbOfEvents);		

		kAbsArray[i]= kAbsArray[i]/NbOfEvents;
		kermaAbsorber=kAbsArray[i]/(massAbsorberLayer);
	test=test+kAbsArray[i];
		G4cout
		 << "\n Energy released in Absorber layer: " <<i<<"; "<< G4BestUnit(kAbsArray[i],"Energy")
		 <<" ; kerma: "<<G4BestUnit(kermaAbsorber,"Dose")//<<" test "<<G4BestUnit(test,"Energy")
		 <<" +- "<<err<<" %" << G4endl;
	}

	 test=0.0;//ok
	for (int i=0; i<n; i++){
		err = basicAnalyze(kGapArray[i], kGapArray2[i], NbOfEvents);	

		kGapArray[i]= kGapArray[i]/NbOfEvents;
		kermaGap=kGapArray[i]/(massGapLayer);
	test=test+kGapArray[i];
		G4cout
		 << "\n Energy released in Gap layer: " <<i<<"; "<< G4BestUnit(kGapArray[i],"Energy")<<" ; kerma: "<<G4BestUnit(kermaGap,"Dose")//<<" test "<<G4BestUnit(test,"Energy")
		 <<" +- "<<err<<" %" << G4endl;
	}

	//Dosimeteres are usualy calibrated in terms of kerma in air due to its direct
	//relation with ionization energy in air Wair=13.6eV. Number of e-ion times W gives
	//direct kerma in the dosimeter volume.
	//Kerma reffers to sum of kinetic energy of secondaries charged particles emmitted by 
	//primary photons. It has two parts: 
	//1. collision kerma which reffers to absorbed dose 
	//which is guaranteed to occur somewhere in target media (infinite size). 
	//2. radiative kerma which reffers to radiation energy emmited by secondaries
	//during variouse processes such as bremmstrahlung.
	//If target material has low Z, such as water, tissue and so on, brems is about 0.3 %
	//so radiative kerma can be neglected.

	//Absorbed dose in a particular location is not necessarly equal with collision kerma.
	//Only at electron equilibrium collision kerma equals absorbed dose.
	//It foloows that:
	//Small ionization chamber filled with air does not work in regime of electron equilibrium
	//because electron track in air is larger then the chamber size. Therefore, one cannot
	//compute dose in medium from dose in air in cavity by multiplying with the ratio of
	//mass-energy absorption coefficients for media and air because this ratio is equal with
	//collision kerma ratio which is not equal with dose ratio.
	//Also, one cannot compute dose in medium from dose in cavity by
	//multiplying with the stopping power ratio of electrons in medium and air because dose
	//in air is undefined.
	//Moreover, when primary radiation comes from a spectrum, the analytical approach is cumbersome.

	//what is always true is the statement that the reading is proportional with absorbed dose (and kerma).
	//if not, the dosimeter is not a reliable one.
	//SO, one can inferr the corection factor for dose in media versus reading of dosimeter:
	//simdose is internal and not actual dose (variance reduction, etc.). 
	//Reading/simdose in gap layer 0 (closer to inteface where dosimeter is placed)=Dosemedium/simdose in last absorber layer
	//notice dimdose in this simulation in Absorber layers is ~ constant
	//so: Dosemedium = reading x ConvFactor. ConvFactor=SimDoseAbsLastLayer/SimDoseGapFirstLayer.

	double conversionFactor=(ddAbsArray[n-1]/massAbsorberLayer)/(ddGapArray[0]/massGapLayer);

	//G4cout
		// << "\n Abs(n-1) = "<<ddAbsArray[n-1]/massAbsorberLayer<<"; gap(0)= "<<ddGapArray[0]/massGapLayer<<"\n";
	G4cout
		 << "\n Dose_in_medium = Dosimeter_Reading x "<<conversionFactor<<"\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::prepareAliasSampling(int nsbin, std::vector<double> fs_array){
	// "================double[] ws_array, int[] ibin_array)=================
		// "
		// " inputs: nsbin: number of bins in the histogram
		// " fs_array: bin probabilities
		// "
		// " Note that we don't need the bin limits at this point, they
		// " are needed for the actual sampling (in alias_sample)
		// "
		// " outputs: ws_array, ibin_array: alias table ready for sampling
		// "
		// "====================================================================

		// ;Copyright NRC;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL fs_array(nsbin),ws_array(nsbin);
	    ws_array=new double[nsbin];
		ibin_array=new int[nsbin];

		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;

		bool exit1b = false;
		bool exit2b = false;

		sum = 0;
		for (i = 1; i <= nsbin; i++) {
			if (fs_array[i - 1] < 1.e-30) {
				fs_array[i - 1] = 1.e-30;
			}
			ws_array[i - 1] = -fs_array[i - 1];
			ibin_array[i - 1] = 1;
			sum = sum + fs_array[i - 1];
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			exit1b = false;
			exit2b = false;
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (abs(ws_array[j_h - 1]) > sum) {
						// GOTO :AT_EXIT1:;
						exit1b = true;
						break;
					}
				}
			}
			if (!exit1b)
				j_h = nsbin;
			// :AT_EXIT1:

			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (abs(ws_array[j_l - 1]) < sum) {
						// GOTO :AT_EXIT2:;
						exit2b = true;
						break;
					}
				}
			}
			if (!exit2b)
				j_l = nsbin;
			// :AT_EXIT2:

			aux = sum - abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;

			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1;
			}

		}

		PrimaryGeneratorAction::SetSpectrumAliasSamplingData(ws_array,ibin_array, enV_array, ndata);
}

void RunAction::prepareAliasSampling(int nsbin, double* fs_array){
	// "================double[] ws_array, int[] ibin_array)=================
		// "
		// " inputs: nsbin: number of bins in the histogram
		// " fs_array: bin probabilities
		// "
		// " Note that we don't need the bin limits at this point, they
		// " are needed for the actual sampling (in alias_sample)
		// "
		// " outputs: ws_array, ibin_array: alias table ready for sampling
		// "
		// "====================================================================

		// ;Copyright NRC;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL fs_array(nsbin),ws_array(nsbin);
	    ws_array=new double[nsbin];
		ibin_array=new int[nsbin];

		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;

		bool exit1b = false;
		bool exit2b = false;

		sum = 0;
		for (i = 1; i <= nsbin; i++) {
			if (fs_array[i - 1] < 1.e-30) {
				fs_array[i - 1] = 1.e-30;
			}
			ws_array[i - 1] = -fs_array[i - 1];
			ibin_array[i - 1] = 1;
			sum = sum + fs_array[i - 1];
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			exit1b = false;
			exit2b = false;
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (abs(ws_array[j_h - 1]) > sum) {
						// GOTO :AT_EXIT1:;
						exit1b = true;
						break;
					}
				}
			}
			if (!exit1b)
				j_h = nsbin;
			// :AT_EXIT1:

			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (abs(ws_array[j_l - 1]) < sum) {
						// GOTO :AT_EXIT2:;
						exit2b = true;
						break;
					}
				}
			}
			if (!exit2b)
				j_l = nsbin;
			// :AT_EXIT2:

			aux = sum - abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;

			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1;
			}

		}

		PrimaryGeneratorAction::SetSpectrumAliasSamplingData(ws_array,ibin_array, enV_array, ndata);
}

double RunAction::linInt(double x1, double y1, double x2, double y2,double x){
	double result = -1.0;
	double mn[2];
    // insucces
	mn[0] = -1.0;// m
	mn[1] = -1.0;// n
	double num = x1 - x2;
	if (num != 0.0) {
		mn[0] = (y1 - y2) / num;
		mn[1] = (x1 * y2 - y1 * x2) / num;
		result = mn[0] * x + mn[1];
	}
	return result;
}

int RunAction::findNearestDataIndexFromArray(double* a, int arraySize, double value, bool lowerThanValue){
	bool b = true;
	int ip = 0;

	if (arraySize > 1) {
		while (b) {
			if (lowerThanValue) {
                if ((a[ip] <= value) && (a[ip + 1] > value)) {
					break;
				}
			} else {
				if (ip > 0)
					if ((a[ip] >= value) && (a[ip - 1] < value)) {
						break;
					}
			}

			ip++;
			if (ip == arraySize - 1) {
				b = false;
				break;
			}
        }
			//nearestposition = ip;// ----------------
        return ip;//a[ip];
    } else {
			//nearestposition = 0;// ---------------
        return 0;//a[0];
    }
}
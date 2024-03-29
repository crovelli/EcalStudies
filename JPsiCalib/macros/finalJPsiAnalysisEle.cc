#include <TStyle.h>
#include <TCanvas.h>
#include "TLine.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string.h>

#include <math.h>
#include <map>
#include <vector>

#include "finalJPsiAnalysisEle.hh"

#define HLTCUT 
#define PTcut 4

using namespace std;

finalJPsiAnalysisEle::finalJPsiAnalysisEle(TTree *tree) 
  : JPsiTreeBase(tree) {
}

finalJPsiAnalysisEle::~finalJPsiAnalysisEle(){ } 

void finalJPsiAnalysisEle::Loop(int theSample) {

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // booking histos
  bookHistos();

  // loop over events
  std::cout << "Number of entries = " << nentries << std::endl;


  // counters
  long AllIncluded         = 0;
  long totalEvents         = 0;
  long totalReco           = 0;
  long totalRecoGt4        = 0;
  long totalRecoGt4Charge  = 0;
  long totalTrigger        = 0;
  long totalIdentified     = 0;
  long numberOfPairsOk     = 0;
  long numbersOfInvMassOk  = 0;
  long totalIdentifiedMore = 0;
  long goodGene            = 0;
  long okEvent_McTruth     = 0;
  long okEvent_Highest     = 0;

  bool passStep[7];    


  // to choose the best pair
  int okHighest_many      = 0;
  int okBestS9S25_many    = 0;
  int okBestSEE_many      = 0;
  int okBestDeta_many     = 0;
  int wrongHighest_many   = 0;
  int wrongBestS9S25_many = 0;
  int wrongBestSEE_many   = 0;
  int wrongBestDeta_many  = 0;

  int data = 0; 
  
  for (int jentry=0; jentry<nentries; jentry++) {

    for (int i=0;i<7;i++) {
      passStep[i]= false;
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry ==1) {
      std::cout << ">>> Signal is " << signal << std::endl;
      if (pthat == -1) data = 1;
      std::cout << ">>> Data is " << data << std::endl;
    }

    if (jentry%100000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    AllIncluded++;

    if ( !signal && isJPsi ) continue;  // important! remove J/psi events in background samples

    totalEvents++;

    // counters for this entry
    int totEleGt4      = 0;
    int totEleGt4plus  = 0;
    int totEleGt4minus = 0;

    // ---------------------------------------------------
    // 1) plots over MC truth
    TVector3 trueEle_3P, truePos_3P;
    TLorentzVector trueEle_4P, truePos_4P;
    float maxPtGen = 0.;
    float minPtGen = 0.;
    
    HepHisto_ptHatAll->Fill(pthat); 
    /// if it is signal sample
    if (signal && !data) {
      
      HepHisto_size->Fill(numberOfGenerated);
    
      if(numberOfGenerated!=2) continue;

      if (chargeGenEle[0]==-1 && chargeGenEle[1]==1) {
	truePos_3P.SetXYZ (pxGenEle[1], pyGenEle[1], pzGenEle[1]);
	trueEle_3P.SetXYZ (pxGenEle[0], pyGenEle[0], pzGenEle[0]);
	truePos_4P.SetXYZT(pxGenEle[1], pyGenEle[1], pzGenEle[1], eneGenEle[1]); 
	trueEle_4P.SetXYZT(pxGenEle[0], pyGenEle[0], pzGenEle[0], eneGenEle[0]); 
      }
      
      if (chargeGenEle[1]==-1 && chargeGenEle[0]==1) {
	trueEle_3P.SetXYZ (pxGenEle[1], pyGenEle[1], pzGenEle[1]);
	truePos_3P.SetXYZ (pxGenEle[0], pyGenEle[0], pzGenEle[0]);
	trueEle_4P.SetXYZT(pxGenEle[1], pyGenEle[1], pzGenEle[1], eneGenEle[1]); 
	truePos_4P.SetXYZT(pxGenEle[0], pyGenEle[0], pzGenEle[0], eneGenEle[0]); 
      }    
      
      if ( trueEle_3P.Perp()>truePos_3P.Perp() ){ maxPtGen=trueEle_3P.Perp(); minPtGen=truePos_3P.Perp(); }
      if ( trueEle_3P.Perp()<truePos_3P.Perp() ){ maxPtGen=truePos_3P.Perp(); minPtGen=trueEle_3P.Perp(); }    
      
      float hep_deltaR = trueEle_3P.DeltaR(truePos_3P);
      float hep_mee = (trueEle_4P + truePos_4P).M();
      
      // filling histos
      HepHisto_eta->Fill(trueEle_3P.Eta());
      HepHisto_eta->Fill(truePos_3P.Eta());
      HepHisto_phi->Fill(trueEle_3P.Phi());
      HepHisto_phi->Fill(truePos_3P.Phi());
      HepHisto_deltaR->Fill(hep_deltaR);
      HepHisto_invMass->Fill(hep_mee);
      HepHisto_maxPt->Fill(maxPtGen);
      HepHisto_minPt->Fill(minPtGen);
      HepHisto_PtVsEta->Fill(trueEle_3P.Eta(),trueEle_3P.Perp());
      HepHisto_PtVsEta->Fill(truePos_3P.Eta(),truePos_3P.Perp());
    }

    // if it is backgroung sample or  number of generated == 2
    if(!signal || numberOfGenerated==2 || data) {
      goodGene++;

#ifdef HLTCUT      // selecting HLTbit 
      if (hltJpsi ) {
#endif

	totalTrigger++;
	//	cout<< "siamo qui?? " << totalTrigger << endl;
	
	passStep[0] = true;
	ScHistoEle_size -> Fill(numberOfElectrons);
	
	// few numbers: all electrons
	if (numberOfElectrons>1) {
	  totalReco++;
	  passStep[1] = true;
	}

	// few numbers: all electrons, Et > PTcut
	TLorentzVector tmp;
	std::vector<TLorentzVector> this4Pe, this4Pp;	
	for(int theEle=0; theEle<numberOfElectrons; theEle++) { 
	  if (etRecoEle[theEle]> PTcut) totEleGt4++;
	  if (etRecoEle[theEle]> PTcut && chargeRecoEle[theEle]>0) {
	    totEleGt4plus++;
	    tmp.SetPxPyPzE(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle], eneRecoEle[theEle]);
	    this4Pp.push_back(tmp);
	  }
	  if (etRecoEle[theEle]> PTcut && chargeRecoEle[theEle]<0) {
	    totEleGt4minus++;
	    tmp.SetPxPyPzE(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle], eneRecoEle[theEle]);
	    this4Pe.push_back(tmp);
	  }
	}
	/// Fill invMassplot with combinatorial mee
	for(int l=0;l<this4Pe.size();l++) {
	  for(int m=0;m<this4Pp.size();m++) {
	    //cout << l << " "  << m << this4Pe.size() << " " << this4Pp.size() <<  endl;
	    float mee_comb = (this4Pe[l] + this4Pp[m]).M();
	    ScHisto_invMassHighestEt[3]->Fill(mee_comb);
	  }
	}

	ScHistoGt4_size -> Fill(totEleGt4);
	if(totEleGt4plus>=1)  ScHistoGt4plus_size  -> Fill(totEleGt4plus);
	if(totEleGt4minus>=1) ScHistoGt4minus_size -> Fill(totEleGt4minus);
	if(totEleGt4>1) {
	  totalRecoGt4++;
	  passStep[2] = true;
	}
	if(totEleGt4plus>=1 && totEleGt4minus>=1) {
	  totalRecoGt4Charge++;
	  passStep[3] = true;
	}



	// --------------------------------------------------------------------------------
	// 2) loop to decide which criteria to apply to select the electron if more than 1 is reconstructed	
	TVector3 closerToMc3p,  highestEt_3p, bestS9S25_3p, bestSEE_3p, bestDeta_3p;
	TVector3 closerToMc3e,  highestEt_3e, bestS9S25_3e, bestSEE_3e, bestDeta_3e;
	TLorentzVector closerToMc4p,  highestEt_4p, bestS9S25_4p, bestSEE_4p, bestDeta_4p; 
	TLorentzVector closerToMc4e,  highestEt_4e, bestS9S25_4e, bestSEE_4e, bestDeta_4e; 

	float deltaRposi =  9999.; 
	float deltaRelec =  9999.; 
	float highestEtp = -9999.;
	float highestEte = -9999.;
	float bestS9S25p =  9999.;
	float bestS9S25e =  9999.;
	float bestSEEp   =  9999.;
	float bestSEEe   =  9999.; 
	float bestDetae  =  9999.;
	float bestDetap  =  9999.; 
	
	int theDeltaRp    = -1;
	int theDeltaRe    = -1;
	int theHighestEtp = -1;
	int theHighestEte = -1;
	int theBestS9S25p = -1;
	int theBestS9S25e = -1;
	int theBestSEEp   = -1; 
	int theBestSEEe   = -1; 
	int theBestDetap  = -1;
	int theBestDetae  = -1;
	
	float highestEt_S9S25p  = 999.; 
	float highestEt_S9S25e  = 999.; 
	float highestEt_SEEp    = 999.;
	float highestEt_SEEe    = 999.;
	float highestEt_HoEp    = 999.;
	float highestEt_HoEe    = 999.;
	float highestEt_dEtap   = 999.;
	float highestEt_dEtae   = 999.;
	float highestEt_dPhip   = 999.;
	float highestEt_dPhie   = 999.;
	float highestEt_TIso03p = 999.;
        float highestEt_TIso03e = 999.;
	
	int numberOfEleOkPlus   = 0;
	int numberOfEleOkMinus  = 0;


	for(int theEle=0; theEle<numberOfElectrons; theEle++) { 
	  
	  if (totEleGt4plus>=1 && totEleGt4plus>=1) {
	    
	    // skipping problematic electrons
	    if (pxRecoEle[theEle]>-700) {
	      
	      // skipping electrons below threshold
	      if (etRecoEle[theEle]> PTcut) {

		// searching for some criteria to select the two 'best' electrons (among those with Et>=4) ( analysis carried on using ECAL )	      
		TVector3 this3P; 
		TLorentzVector this4P;
		this3P.SetXYZ(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle]);
		this4P.SetPxPyPzE(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle], eneRecoEle[theEle]);
		
		// further selections
		bool isGood   = true;
		bool isBarrel = true;
		bool isPF = false;

		fillHistoCutsVariables( theEle );

		if(isPFlowRecoEle[theEle]) isPF = true;

		if(fabs(this3P.Eta())>1.479) isBarrel = false;


		if (isBarrel){		
		  if( fabs(dEtaAtVtxRecoEle[theEle]) > 0.010 ) isGood = false;
		  if( fabs(dPhiAtVtxRecoEle[theEle]) > 0.06 )   isGood = false;
		  if( fabs(EoverPRecoEle[theEle] -1) > 0.2 )   isGood = false;
		  ///if( HoverERecoEle[theEle]>0.003 )                    isGood = false;
		  if( sigmaIetaIetaRecoEle[theEle]>0.018 )             isGood = false;
		  if( dr03EcalSumEtRecoEle[theEle]>5) isGood = false;
		}
		
		if (!isBarrel){
		  if( fabs(dEtaAtVtxRecoEle[theEle])>0.012 ) isGood = false;
		  if( fabs(dPhiAtVtxRecoEle[theEle])>0.06 )   isGood = false;
		  if( fabs((EoverPRecoEle[theEle] -1))> 0.4 )   isGood = false;
		  //if( HoverERecoEle[theEle]>0.003 )          isGood = false;
		  //if( sigmaEtaEtaRecoEle[theEle]>0.06 )      isGood = false;
		  if( dr03EcalSumEtRecoEle[theEle]>11 )       isGood = false;
		}
		
		if (!isGood) continue;
		if (chargeRecoEle[theEle]>0) numberOfEleOkPlus++;
		if (chargeRecoEle[theEle]<0) numberOfEleOkMinus++;
		
		
		// search for the two reco electrons closer to the generated electrons in case of signal
		if (signal && !data) {
		  float deltaRmc=0;
		  if(chargeRecoEle[theEle]>0) deltaRmc = this3P.DeltaR(truePos_3P);    
		  if(chargeRecoEle[theEle]<0) deltaRmc = this3P.DeltaR(trueEle_3P);    
		  
		  if((chargeRecoEle[theEle]>0) && (deltaRmc<deltaRposi)) { 
		    theDeltaRp=theEle;
		    deltaRposi=deltaRmc; 
		    closerToMc3p=this3P; 
		    closerToMc4p=this4P; 
		  }
		  
		  if((chargeRecoEle[theEle]<0) && (deltaRmc<deltaRelec)) { 
		    theDeltaRe=theEle;
		    deltaRelec=deltaRmc; 
		    closerToMc3e=this3P; 
		    closerToMc4e=this4P; 
		  }
		}
		
		
		// search for highest Et electrons which charge hp = 1
		if ( (chargeRecoEle[theEle]>0) && (etRecoEle[theEle]>=highestEtp) ){ 
		  theHighestEtp        = theEle;
		  highestEtp           = etRecoEle[theEle];
		  highestEt_3p         = this3P;
		  highestEt_4p         = this4P;
		  highestEt_S9S25p     = e9e25RecoEle[theEle];
		  highestEt_SEEp       = sigmaEtaEtaRecoEle[theEle];
		  highestEt_HoEp       = HoverERecoEle[theEle];
		  highestEt_dEtap      = dEtaAtVtxRecoEle[theEle];
		  highestEt_dPhip      = dPhiAtVtxRecoEle[theEle];
		  highestEt_TIso03p    = dr03TkSumPtRecoEle[theEle];
		}
		
		// search for highest Et electrons which charge hp = -1
		if ( (chargeRecoEle[theEle]<0) && (etRecoEle[theEle]>=highestEte) ){
		  theHighestEte        = theEle;                    
		  highestEte           = etRecoEle[theEle];         
		  highestEt_3e         = this3P;                    
		  highestEt_4e         = this4P;                    
		  highestEt_S9S25e     = e9e25RecoEle[theEle];      
		  highestEt_SEEe       = sigmaEtaEtaRecoEle[theEle];
		  highestEt_HoEe       = HoverERecoEle[theEle];     
		  highestEt_dEtae      = dEtaAtVtxRecoEle[theEle];  
		  highestEt_dPhie      = dPhiAtVtxRecoEle[theEle];  
		  highestEt_TIso03e    = dr03TkSumPtRecoEle[theEle];
		}
		
		
		// search for best E9/E25 electrons with charge Hp = +1
		float s9s25diff = fabs(e9e25RecoEle[theEle] - 1.);
		if ( (chargeRecoEle[theEle]>0) && (s9s25diff<=bestS9S25p) ){
		  theBestS9S25p = theEle;
		  bestS9S25p    = s9s25diff;
		  bestS9S25_3p  = this3P;
		  bestS9S25_4p  = this4P;
		}
		
		// search for best E9/E25 electrons with charge Hp = -1
		if ( (chargeRecoEle[theEle]<0) && (s9s25diff<=bestS9S25e) ){
		  theBestS9S25e = theEle;
		  bestS9S25e    = s9s25diff;
		  bestS9S25_3e  = this3P;
		  bestS9S25_4e  = this4P;
		}
		
		
		// search for best sigmaEtaEta electrons with charge Hp = +1
		float sEEdiff = fabs(sigmaEtaEtaRecoEle[theEle]);
		if ( (chargeRecoEle[theEle]>0) && (sEEdiff<=bestSEEp) ){
		  theBestSEEp = theEle;
		  bestSEEp    = sEEdiff;
		  bestSEE_3p  = this3P;
		  bestSEE_4p  = this4P;
		}
		
		// search for best E9/E25 electrons with charge Hp = -1
		if ( (chargeRecoEle[theEle]<0) && (sEEdiff<=bestSEEe) ){
		  theBestSEEe = theEle;
		  bestSEEe    = sEEdiff;
		  bestSEE_3e  = this3P;
		  bestSEE_4e  = this4P;
		}
		
		// search for best dEta electrons with charge Hp = +1	      
		float absDeta = fabs(dEtaAtVtxRecoEle[theEle]);
		if ( (chargeRecoEle[theEle]>0) && (absDeta<=bestDetap) ){
		  theBestDetap = theEle;
		  bestDetap    = absDeta;
		  bestDeta_3p  = this3P;
		  bestDeta_4p  = this4P;
		}
		
		// search for best dEta electrons with charge Hp = +1	      
		if ( (chargeRecoEle[theEle]<0) && (absDeta<=bestDetae) ){
		  theBestDetae = theEle;
		  bestDetae    = absDeta;
		  bestDeta_3e  = this3P;
		  bestDeta_4e  = this4P;
		}
		
	      } // ok Et
	    }  // ok matching with track
	    
	  } // two matched electons
	  
	} // loop over electron collection


	
	if ( numberOfEleOkPlus>0 && numberOfEleOkMinus>0 ) {
	  totalIdentified++;	  
	  passStep[4] = true;
	}
	if ((numberOfEleOkPlus>0 && numberOfEleOkMinus>1 ) || (numberOfEleOkPlus>1 && numberOfEleOkMinus>0)) totalIdentifiedMore++;
	
	
	// to choose, signal only
	if (signal && !data) {
	  
	  if (numberOfEleOkPlus>=1 && numberOfEleOkMinus>=1) {
	    
	    if ( theBestDetap==-1 || theBestDetae==-1 ) cout << "non deve succede mai" << endl;
	    
	    float deltaR_WrtMcEle_closerToMc = closerToMc3e.DeltaR(trueEle_3P);
	    float deltaR_WrtMcPos_closerToMc = closerToMc3p.DeltaR(truePos_3P);
	    float deltaR_WrtMcEle_highestEt  = highestEt_3e.DeltaR(trueEle_3P);
	    float deltaR_WrtMcPos_highestEt  = highestEt_3p.DeltaR(truePos_3P);
	    float deltaR_WrtMcEle_bestS9S25  = bestS9S25_3e.DeltaR(trueEle_3P);
	    float deltaR_WrtMcPos_bestS9S25  = bestS9S25_3p.DeltaR(truePos_3P);
	    float deltaR_WrtMcEle_bestSEE    = bestSEE_3e.DeltaR(trueEle_3P);
	    float deltaR_WrtMcPos_bestSEE    = bestSEE_3p.DeltaR(truePos_3P);
	    float deltaR_WrtMcEle_bestDeta   = bestDeta_3e.DeltaR(trueEle_3P);
	    float deltaR_WrtMcPos_bestDeta   = bestDeta_3p.DeltaR(truePos_3P);
	  
	    // how many times the electron/positron pair we use is correct?
	    if (deltaR_WrtMcEle_closerToMc < 0.5 && deltaR_WrtMcPos_closerToMc<0.5) okEvent_McTruth++; 
	    if (deltaR_WrtMcEle_highestEt < 0.5 && deltaR_WrtMcPos_highestEt<0.5)   okEvent_Highest++; 
	  
	    if (numberOfEleOkPlus>1) {
	      ScHisto_deltaRWrtMc_CloserToMc_more1 -> Fill(deltaR_WrtMcPos_closerToMc);
	      ScHisto_deltaRWrtMc_HighestEt_more1  -> Fill(deltaR_WrtMcPos_highestEt);
	      ScHisto_deltaRWrtMc_BestS9S25_more1  -> Fill(deltaR_WrtMcPos_bestS9S25);
	      ScHisto_deltaRWrtMc_BestSEE_more1    -> Fill(deltaR_WrtMcPos_bestSEE);
	    
	      // count how many times we get the correct one
	      if (deltaRposi<0.5){
		if (theDeltaRp==theHighestEtp) okHighest_many++;
		else wrongHighest_many++;
		if (theDeltaRp==theBestS9S25p) okBestS9S25_many++;
		else wrongBestS9S25_many++;
		if (theDeltaRp==theBestSEEp)   okBestSEE_many++;
		else wrongBestSEE_many++;
		if (theDeltaRp==theBestDetap)  okBestDeta_many++;
		else wrongBestDeta_many++;
	      }
	    }
	  
	    if (numberOfEleOkMinus>1) {
	      ScHisto_deltaRWrtMc_CloserToMc_more1 -> Fill(deltaR_WrtMcEle_closerToMc);
	      ScHisto_deltaRWrtMc_HighestEt_more1  -> Fill(deltaR_WrtMcEle_highestEt);
	      ScHisto_deltaRWrtMc_BestS9S25_more1  -> Fill(deltaR_WrtMcEle_bestS9S25);
	      ScHisto_deltaRWrtMc_BestSEE_more1    -> Fill(deltaR_WrtMcEle_bestSEE);
	    
	      // count how many times we get the correct one
	      if (deltaRelec<0.5){
		if (theDeltaRe==theHighestEte) okHighest_many++;
		else wrongHighest_many++;
		if (theDeltaRe==theBestS9S25e) okBestS9S25_many++;
		else wrongBestS9S25_many++;
		if (theDeltaRe==theBestSEEe)   okBestSEE_many++;
		else wrongBestSEE_many++;
		if (theDeltaRe==theBestDetae)  okBestDeta_many++;
		else wrongBestDeta_many++;
	      }
	    }
	  }  // ok matched   
	} // signal
      
      
	// to compare the distributions: use the two highest Et electrons
	if (numberOfEleOkPlus>=1 && numberOfEleOkMinus>=1) {  

	  float mee_highestEt = (highestEt_4p + highestEt_4e).M();
	  ScHisto_invMassHighestEt[4]->Fill(mee_highestEt);

	
	  // tracker isolation added here
	  float deltaR_highestEt = highestEt_3p.DeltaR(highestEt_3e);	
	  if (deltaR_highestEt < 0.3) {
	    highestEt_TIso03p -= highestEt_4e.Et();
	    if (highestEt_TIso03p < 0.) highestEt_TIso03p = 0.;
	    highestEt_TIso03e -= highestEt_4p.Et();
	    if (highestEt_TIso03e < 0.) highestEt_TIso03e = 0.;
	  }
	
	  EleHisto_dr03tkcorr->Fill(highestEt_TIso03e);
	  EleHisto_dr03tkcorr->Fill(highestEt_TIso03p);

	  // BEST PAIR CUTS  it was 2.8 !!
	  if (highestEt_TIso03p < 3 && highestEt_TIso03e < 3 ) {
	    numberOfPairsOk++;
	    passStep[5] = true;
	    
	    float maxEt_highestEt = 0.;
	    float minEt_highestEt = 0.;
	    if (highestEt_4p.Et()>highestEt_4e.Et()){ 
	      maxEt_highestEt=highestEt_4p.Et(); 
	      minEt_highestEt=highestEt_4e.Et(); 
	    }
	    if (highestEt_4p.Et()<highestEt_4e.Et()){ 
	      maxEt_highestEt=highestEt_4e.Et(); 
	      minEt_highestEt=highestEt_4p.Et(); 
	    }
	    
	    float deltaR_highestEt = highestEt_3p.DeltaR(highestEt_3e);
	    
	    //  float mee_highestEt = (highestEt_4p + highestEt_4e).M();
	    ScHisto_invMassHighestEt[5]->Fill(mee_highestEt);
	    
	    if (mee_highestEt<3.5 && mee_highestEt>2.5) {
	      numbersOfInvMassOk++;
	      passStep[6] = true;

	      ScHisto_maxEtHighestEt->Fill(maxEt_highestEt);
	      ScHisto_minEtHighestEt->Fill(minEt_highestEt);
	      ScHisto_deltaRHighestEt->Fill(deltaR_highestEt);	
	      ScHisto_etaHighestEt->Fill(highestEt_3p.Eta());
	      ScHisto_etaHighestEt->Fill(highestEt_3e.Eta());
	      ScHisto_phiHighestEt->Fill(highestEt_3p.Phi());
	      ScHisto_phiHighestEt->Fill(highestEt_3e.Phi());

	    }

	    if (mee_highestEt>3.5 && mee_highestEt<5) {

	      ScHisto_maxEtBG ->Fill(maxEt_highestEt);
	      ScHisto_minEtBG ->Fill(minEt_highestEt);
	      ScHisto_deltaRBG->Fill(deltaR_highestEt);	
	      ScHisto_etaBG   ->Fill(highestEt_3p.Eta());
	      ScHisto_etaBG   ->Fill(highestEt_3e.Eta());
	      ScHisto_phiBG   ->Fill(highestEt_3p.Phi());
	      ScHisto_phiBG   ->Fill(highestEt_3e.Phi());

	    }

	    if (mee_highestEt>5 && mee_highestEt<6.5) {

	      ScHisto_maxEtBG2 ->Fill(maxEt_highestEt);
	      ScHisto_minEtBG2 ->Fill(minEt_highestEt);
	      ScHisto_deltaRBG2->Fill(deltaR_highestEt);	
	      ScHisto_etaBG2   ->Fill(highestEt_3p.Eta());
	      ScHisto_etaBG2   ->Fill(highestEt_3e.Eta());
	      ScHisto_phiBG2   ->Fill(highestEt_3p.Phi());
	      ScHisto_phiBG2   ->Fill(highestEt_3e.Phi());

	    }

	    // to study EB/ EE depencency: both in barrel only
            if ( fabs(highestEt_3p.Eta())<1.479 && fabs(highestEt_3e.Eta())<1.479 ){
              ScHisto_deltaRHighestEt_EB  -> Fill(deltaR_highestEt);
              ScHisto_invMassHighestEt_EB -> Fill(mee_highestEt);
            }
	    
            // both in endcap only
            if ( fabs(highestEt_3p.Eta())>1.479 && fabs(highestEt_3e.Eta())>1.479 ){
              ScHisto_deltaRHighestEt_EE  -> Fill(deltaR_highestEt);
              ScHisto_invMassHighestEt_EE -> Fill(mee_highestEt);
            }
	    
            // one and one
            if (  (fabs(highestEt_3p.Eta())<1.479 && fabs(highestEt_3e.Eta())>1.479) || (fabs(highestEt_3e.Eta())<1.479 && fabs(highestEt_3p.Eta())>1.479) ) {
              ScHisto_deltaRHighestEt_EBEE       -> Fill(deltaR_highestEt);
              ScHisto_invMassHighestEt_EBEE      -> Fill(mee_highestEt);
            }
	    
	    
	    // study of possible biases for signal
	    if ( (signal || data) && (mee_highestEt<3.5 && mee_highestEt>2.5)) {    
	      float jpsi_ene = (highestEt_4p + highestEt_4e).E();	  
	      float jpsi_et  = (highestEt_4p + highestEt_4e).Et();	  
	      float jpsi_eta = (highestEt_4p + highestEt_4e).Eta();	  
	      float jpsi_phi = (highestEt_4p + highestEt_4e).Phi();	  
	      float deltaR   =  highestEt_4p.DeltaR(highestEt_4e);
	      
	      ScHisto_JeneVsJeta_highestEt      -> Fill(jpsi_eta, jpsi_ene);
	      ScHisto_JetVsJeta_highestEt       -> Fill(jpsi_eta, jpsi_et);
	      ScHisto_InvMassVsJeta_highestEt   -> Fill(jpsi_eta, mee_highestEt);                
	      ScHisto_InvMassVsJphi_highestEt   -> Fill(jpsi_phi, mee_highestEt);                
	      ScHisto_InvMassVsJene_highestEt   -> Fill(jpsi_ene, mee_highestEt);        
	      ScHisto_InvMassVsJet_highestEt    -> Fill(jpsi_et,  mee_highestEt);        
	      ScHisto_InvMassVsDeltaR_highestEt -> Fill(deltaR,   mee_highestEt);        
	      if (fabs(jpsi_eta)<1.5) ScHisto_InvMassVsJeneEB_highestEt -> Fill(jpsi_ene, mee_highestEt);        
	      if (fabs(jpsi_eta)>1.5) ScHisto_InvMassVsJeneEE_highestEt -> Fill(jpsi_ene, mee_highestEt);        
	      if (fabs(jpsi_eta)<1.5) ScHisto_InvMassVsJetEB_highestEt  -> Fill(jpsi_et,  mee_highestEt);        
	      if (fabs(jpsi_eta)>1.5) ScHisto_InvMassVsJetEE_highestEt  -> Fill(jpsi_et,  mee_highestEt);        
	    }
	  
	  } // ok best pair with Tracker isol
	} // ok matched
      
	// plots with all combinatorics
	if ( numberOfEleOkPlus>=1 && numberOfEleOkMinus>=1 && numberOfPairsOk>=1 ){ 
	  for(int theEle1=0; theEle1<numberOfElectrons; theEle1++) { 
	    if (pxRecoEle[theEle1]<-700) continue;
	    if (etRecoEle[theEle1]< PTcut)     continue;
	  
	    for(int theEle2=(theEle1+1); theEle2<numberOfElectrons; theEle2++) { 
	      if (pxRecoEle[theEle2]<-700) continue;
	      if (etRecoEle[theEle2]< PTcut)     continue;
	    
	      TLorentzVector tlvTheEle1, tlvTheEle2;
	      TVector3 tv3TheEle1, tv3TheEle2;
	      tlvTheEle1.SetPxPyPzE(pxRecoEle[theEle1], pyRecoEle[theEle1], pzRecoEle[theEle1], eneRecoEle[theEle1]);
	      tlvTheEle2.SetPxPyPzE(pxRecoEle[theEle2], pyRecoEle[theEle2], pzRecoEle[theEle2], eneRecoEle[theEle2]);
	      tv3TheEle1.SetXYZ (pxRecoEle[theEle1], pyRecoEle[theEle1], pzRecoEle[theEle1]);
	      tv3TheEle2.SetXYZ (pxRecoEle[theEle2], pyRecoEle[theEle2], pzRecoEle[theEle2]);
	    
	      float deltaR = tv3TheEle1.DeltaR(tv3TheEle2);
	      float mee    = (tlvTheEle1 + tlvTheEle2).M();
	      ScHisto_deltaRComb->Fill(deltaR);
	      ScHisto_invMassComb->Fill(mee);
	      ScHisto_deltaRVsInvMassComb->Fill(mee,deltaR);
	    }
	  }
	}
#ifdef HLTCUT      // selecting HLTbit 
      } // ok HLT
#endif
    } // ok generated   

    /// filling Histograms per steps
    for (int i=0;i<7;i++) {
      if( passStep[i] ) HepHisto_ptHat[i]->Fill(pthat);
    } 


  } // loop over entries
  
  
  cout << endl;
  cout << "to choose the best criterium" << endl;
  cout << "highest ET: ok    = " << okHighest_many   << ", wrong = " << wrongHighest_many   << endl;
  cout << "best S9/S25: ok   = " << okBestS9S25_many << ", wrong = " << wrongBestS9S25_many << endl;
  cout << "best SEE: ok      = " << okBestSEE_many   << ", wrong = " << wrongBestSEE_many   << endl;
  cout << "best deltaEta: ok = " << okBestDeta_many  << ", wrong = " << wrongBestDeta_many  << endl;
  cout << endl;
  

  saveHistos();

  //  drawPlots();
  float effHLT    = (double) totalTrigger/goodGene;
  float erreffHLT = sqrt(effHLT *(1-effHLT)/goodGene);
  float effSelected    = (double) totalIdentified/totalTrigger;
  float erreffSelected = sqrt(effSelected *(1-effSelected)/totalTrigger);
  float effSelIsol    = (double) numberOfPairsOk/totalTrigger;
  float erreffSelIsol = sqrt(effSelIsol *(1-effSelIsol)/totalTrigger);

  // summary
  cout << "total number of events = "              << totalEvents     << " ALL " << AllIncluded << endl;
  cout << "good number of gene MC = "              << goodGene        << endl;
  cout << "total number of events passing HLT = "  << totalTrigger    << ";  eff-> " << effHLT <<"+-" << erreffHLT << endl;
  cout << "total number of reco = "                << totalReco       << endl;
  cout << "total number of reco & Et>" << PTcut << "= "          << totalRecoGt4    << endl;
  cout << "total number of reco & Et>" << PTcut << " & charge requirement = " << totalRecoGt4Charge << endl;
  cout << "total number of reco & Et>" << PTcut << " & id & isol = " << totalIdentified <<  ";  eff-> " << effSelected <<"+-" << erreffSelected << endl;
  cout << "total number of reco & ET>" << PTcut << " & id & full isol (tracker also) = " << numberOfPairsOk << ";  eff-> " << effSelIsol <<"+-" << erreffSelIsol << endl;
  cout << "total number of reco & ET>" << PTcut << " & id & full isol (tracker also) + invMass in 2.5 - 3.7 = " << numbersOfInvMassOk << endl;
  cout << endl;
  cout << "total number of reco & matched & id, more than 2 = " << totalIdentifiedMore << endl;
  cout << endl;
  cout << "total number of events with closer to MC match MC within dR=0.5 = " << okEvent_McTruth << endl; 
  cout << "total number of events with highest Et match MC within dR=0.5 = "   << okEvent_Highest << endl; 

  
  // normalizing to cross sections
  float filterEff_prod;
  long crossSection; // in picobarn
  if (theSample==1) { filterEff_prod = 0.00285;  crossSection =    13420000; }   // Jpsi7teV          13420000      0.00285 
  if (theSample==2) { filterEff_prod = 0.00046;  crossSection =   235500000;}    // BCtoE20to30       0.2355 mb     0.00046
  if (theSample==3) { filterEff_prod = 0.00234;  crossSection =    59300000;}    // BCtoE30to80       0.0593 mb     0.00234
  if (theSample==4) { filterEff_prod = 0.0104 ;  crossSection =      906000;}    // BCtoE80to170      0.906e-3 mb   0.0104
  if (theSample==5) { filterEff_prod = 0.0073 ;  crossSection =   235500000;}    // EMenrich20to30    0.2355 mb     0.0073
  if (theSample==6) { filterEff_prod = 0.059  ;  crossSection =    59300000;}    // EMenrich30to80    0.0593 mb     0.059
  if (theSample==7) { filterEff_prod = 0.148  ;  crossSection =      906000;}    // EMenrich80to170   0.906e-3 mb   0.148
  if (theSample==8) { filterEff_prod = 2.3    ;    crossSection = 208000000;}    // !!doubleem 6-20 20.8 mb 	    0.023   
  if (theSample==9) { filterEff_prod = 0.235  ;    crossSection = 297000000;}    // doubleem >20    0.297 mb       
  if (theSample==10) { filterEff_prod = 0.0094  ;  crossSection = 95540;    }    // DYee 1-10        95540.       0.0094         
  if (theSample==11) { filterEff_prod =0.074   ;  crossSection = 484400000;    }    // ppEleX   2409996 events 2417242 raw 48.44 mb  0.00074
  
 
  long numEvents;
  numEvents = goodGene; 
  float kineEff_reco          = (double) totalReco/numEvents;
  float kineEff_recoGt4       = (double) totalRecoGt4/numEvents;
  float kineEff_recoGt4charge = (double) totalRecoGt4Charge/numEvents;
  float kineEff_id            = (double) totalIdentified/numEvents;
  float kineEff_isol          = (double) numberOfPairsOk/numEvents;
  float kineEff_invMass       = (double) numbersOfInvMassOk/numEvents;

  float exp_hlt            = effHLT                * filterEff_prod*crossSection;
  float exp_reco           = kineEff_reco          * filterEff_prod*crossSection;
  float exp_recoGt4        = kineEff_recoGt4       * filterEff_prod*crossSection;
  float exp_recoGt4charge  = kineEff_recoGt4charge * filterEff_prod*crossSection;
  float exp_id             = kineEff_id            * filterEff_prod*crossSection;
  float exp_isol           = kineEff_isol          * filterEff_prod*crossSection;
  float exp_okMass         = kineEff_invMass       * filterEff_prod*crossSection;
  
  cout << "number of expected events in 10 pb-1: "            << endl;
  cout << "after HLT                  : " << 10.*exp_hlt      << endl; 
  cout << "after reco                 : " << 10.*exp_reco     << endl; 
  cout << "after reco + Et>4          : " << 10.*exp_recoGt4  << endl; 
  cout << "after reco + Et>4 + charge : " << 10.*exp_recoGt4charge << endl; 
  cout << "after id + ecal isol       : " << 10.*exp_id       << endl; 
  cout << "after isolation            : " << 10.*exp_isol     << endl; 
  cout << "in the inv mass region     : " << 10.*exp_okMass   << endl; 
    

} // end of program



void finalJPsiAnalysisEle::fillHistoCutsVariables(  int theEle ){

  EleHisto_eta->Fill(etaRecoEle[theEle]);
  EleHisto_phi->Fill(phiRecoEle[theEle]);
  EleHisto_detaVtx->Fill(dEtaAtVtxRecoEle[theEle]);   
  EleHisto_dphiVtx->Fill(dPhiAtVtxRecoEle[theEle]);   
  EleHisto_HoE->Fill(HoverERecoEle[theEle]); 
  EleHisto_EoP->Fill(EoverPRecoEle[theEle]);       
  EleHisto_fbrem->Fill(fBremRecoEle[theEle]);       
  EleHisto_sigmaIeta->Fill(sigmaIetaIetaRecoEle[theEle]);
  EleHisto_sigmaEta->Fill(sigmaEtaEtaRecoEle[theEle]);
  EleHisto_dr03ecal->Fill(dr03EcalSumEtRecoEle[theEle]);
  EleHisto_dr03tk->Fill(dr03TkSumPtRecoEle[theEle]);
  EleHisto_dr03hcal1->Fill(dr03Hcal1SumEtRecoEle[theEle]);
  EleHisto_dr03hcal2->Fill(dr03Hcal2SumEtRecoEle[theEle]);
  EleHisto_momErr->Fill(momentumErrorRecoEle[theEle]);
  EleHisto_isPF->Fill(isPFlowRecoEle[theEle]);
  EleHisto_isEcal->Fill(isEcalDrivenRecoEle[theEle]);
  EleHisto_mva->Fill(pFlowMvaRecoEle[theEle]);
}

      
void finalJPsiAnalysisEle::bookHistos() {

  ScHistoEle_size      = new TH1F("ScHistoEle_size",      "Num of electrons", 10, 0,10);      
  ScHistoGt4_size      = new TH1F("ScHistoGt4_size",      "Num of electrons with Et>4", 10, 0,10);
  ScHistoGt4plus_size  = new TH1F("ScHistoGt4plus_size",  "Num of positrons (+) with Et>4", 10, 0,10);
  ScHistoGt4minus_size = new TH1F("ScHistoGt4minus_size", "Num of electrons (-) with Et>4", 10, 0,10);

  HepHisto_size    = new TH1F("HepHisto_size",    "Num of GenElectrons ",  10, 0,  10 );
  HepHisto_ptHatAll = new TH1F("HepHisto_ptHatAll",   "Pt hat distributions",  800, 0., 170.);

  for (int i=0;i<7;i++) {
    ostringstream name,legend;
    name << "HepHisto_ptHat" << i + 1;
    legend << "PtHat - step" << i + 1;
    string strname = name.str();
    string strleg = legend.str();
    TH1F *ScHisto_ptHattmp = new TH1F( strname.c_str(),  strleg.c_str(),  800, 0., 170.);
    HepHisto_ptHat.push_back(ScHisto_ptHattmp);
  }

  HepHisto_deltaR  = new TH1F("HepHisto_deltaR",  "GenElectrons deltaR ", 200, 0.,  2.);
  HepHisto_maxPt   = new TH1F("HepHisto_maxPt",   "GenElectrons max pt ", 200, 0., 20.);
  HepHisto_minPt   = new TH1F("HepHisto_minPt",   "GenElectrons min pt ", 200, 0., 20.);
  HepHisto_eta     = new TH1F("HepHisto_eta",     "GenElectrons eta    ", 100, -3.,+3.);
  HepHisto_phi     = new TH1F("HepHisto_phi",     "GenElectrons phi    ", 100, -3.15,3.15);
  HepHisto_invMass = new TH1F("HepHisto_invMass", "GenElectrons invariant mass", 150,0.,6.); 
  HepHisto_PtVsEta = new TH2F("HepHisto_PtVsEta", "GenElectrons Pt vs Eta", 100, -3.,+3., 200, 0.,20.);

  EleHisto_eta       = new  TH1F("EleHisto_eta      ", " Eta recoEle", 200  , -3.,  3.);      
  EleHisto_phi       = new  TH1F("EleHisto_phi      ", " Phi recoEle", 200, -3.15,3.15);      
  EleHisto_detaVtx   = new  TH1F("EleHisto_detaVtx  ", " deta AtVtx ", 200,  -0.4, 0.4);  
  EleHisto_dphiVtx   = new  TH1F("EleHisto_dphiVtx  ", " dphy AtVtx ", 200,   -1.,  1.);  
  EleHisto_HoE       = new  TH1F("EleHisto_HoE      ", " HoverE Ele ", 200,    0., 0.1);      
  EleHisto_EoP       = new  TH1F("EleHisto_EoP      ", " EoverP Ele ", 200,    0.,  5.);      
  EleHisto_fbrem     = new  TH1F("EleHisto_fbrem    ", " fBrem Ele  ", 200,  -30.,  2.);   
  EleHisto_sigmaIeta = new  TH1F("EleHisto_sigmaIeta", " sigmaIetaIeta", 300,    0., 1.5);
  EleHisto_sigmaEta  = new  TH1F("EleHisto_sigmaEta ", " sigmaEtaEta", 300,    0., 0.5); 
  EleHisto_dr03ecal  = new  TH1F("EleHisto_dr03ecal ", " dr03ecal   ", 200,   -5.,  60); 
  EleHisto_dr03tk    = new  TH1F("EleHisto_dr03tk   ", " dr03tk     ", 200,   -5.,  60);   
  EleHisto_dr03tkcorr= new  TH1F("EleHisto_dr03tkCorr", " dr03tk corr ", 200,   -5.,  60);   
  EleHisto_dr03hcal1 = new  TH1F("EleHisto_dr03hcal1", " dr03hcal1  ", 200,   -5.,  30);
  EleHisto_dr03hcal2 = new  TH1F("EleHisto_dr03hcal2", " dr03hcal2  ", 200,   -5.,  30);
  EleHisto_momErr    = new  TH1F("EleHisto_momErr   ", " momentum error", 500,    0., 1100.);   
  EleHisto_isPF      = new  TH1F("EleHisto_isPF     ", " isPFlow Ele",   2,    0.,   2);     
  EleHisto_isEcal    = new  TH1F("EleHisto_isEcal   ", " isEcalDriven",  2,    0.,   2);   
  EleHisto_mva       = new  TH1F("EleHisto_mva      ", " MvaPFlow   ",   100,    -1.2,   1.2);      
				     

  ScHisto_deltaRWrtMc_CloserToMc_more1      = new TH1F("ScHisto_deltaRWrtMc_CloserToMc_more1", "deltaR with Mc - closer to MC",     100, 0.,5.);
  ScHisto_deltaRWrtMc_HighestEt_more1       = new TH1F("ScHisto_deltaRWrtMc_HighestEt_more1",  "deltaR with Mc - highest Et",       100, 0.,5.);
  ScHisto_deltaRWrtMc_BestS9S25_more1       = new TH1F("ScHisto_deltaRWrtMc_BestS9S25_more1",  "deltaR with Mc - best S9/S25",      100, 0.,5.);
  ScHisto_deltaRWrtMc_BestSEE_more1         = new TH1F("ScHisto_deltaRWrtMc_BestSEE_more1",    "deltaR with Mc - best SigmaEtaEta", 100, 0.,5.);

  ScHisto_deltaRComb            = new TH1F("ScHisto_deltaRComb",            "Sc deltaR",          100, 0.,3.);   
  ScHisto_invMassComb           = new TH1F("ScHisto_invMassComb",           "Sc invariant mass",  150, 0.,6.);   
  ScHisto_deltaRVsInvMassComb   = new TH2F("ScHisto_deltaRVsInvMassComb",   "Sc deltaR vs Invariant Mass", 150, 0.,6., 100, 0.,3.);

  ScHisto_maxEtHighestEt        = new TH1F("ScHisto_maxEtHighestEt",        "Sc max Et 2.5/3.5",        200, 0.,20.);
  ScHisto_minEtHighestEt        = new TH1F("ScHisto_minEtHighestEt",        "Sc min Et 2.5/3.5",        200, 0.,20.);
  ScHisto_deltaRHighestEt       = new TH1F("ScHisto_deltaRHighestEt",       "Sc deltaR 2.5/3.5",        200, 0.,3.2);   
  ScHisto_etaHighestEt          = new TH1F("ScHisto_etaHighestEt",          "Sc eta 2.5/3.5 ",          200, -3.,+3.);
  ScHisto_phiHighestEt          = new TH1F("ScHisto_phiHighestEt",          "Sc phi 2.5/3.5",           200, -3.15,+3.15);

  ScHisto_maxEtBG               = new TH1F("ScHisto_maxEtBG",        "Sc max Et 3.5/5",           200, 0.,20.);
  ScHisto_minEtBG               = new TH1F("ScHisto_minEtBG",        "Sc min Et 3.5/5",           200, 0.,20.); 
  ScHisto_deltaRBG	        = new TH1F("ScHisto_deltaRBG",       "Sc deltaR 3.5/5",           200, 0.,3.2); 
  ScHisto_etaBG                 = new TH1F("ScHisto_etaBG",          "Sc eta 3.5/5",              200, -3.,+3.); 
  ScHisto_phiBG                 = new TH1F("ScHisto_phiBG",          "Sc phi 3.5/5",              200, -3.15,+3.15);

  ScHisto_maxEtBG2               = new TH1F("ScHisto_maxEtBG2",        "Sc max Et 5/6.5",          200, 0.,20.);
  ScHisto_minEtBG2               = new TH1F("ScHisto_minEtBG2",        "Sc min Et 5/6.5",          200, 0.,20.); 
  ScHisto_deltaRBG2	        = new TH1F("ScHisto_deltaRBG2",        "Sc deltaR 5/6.5",          200, 0.,3.2); 
  ScHisto_etaBG2                 = new TH1F("ScHisto_etaBG2",          "Sc eta 5/6.5",             200, -3.,+3.); 
  ScHisto_phiBG2                 = new TH1F("ScHisto_phiBG2",          "Sc phi 5/6.5",             200, -3.15,+3.15);


  //  ScHisto_invMassHighestEt      = new TH1F("ScHisto_invMassHighestEt",      "Sc invariant mass",  150, 0.,6.);   

  for (int i=0;i<6;i++) {
    ostringstream name;
    name << "ScHisto_invMassHighestEt" << i+1 ;
    string strname = name.str();
    TH1F *ScHistotmp = new TH1F( strname.c_str(),    "Sc invariant mass",  150, 0.,6.); 
    ScHisto_invMassHighestEt.push_back(ScHistotmp);
  }

  ScHisto_s9s25HighestEt        = new TH1F("ScHisto_s9s25HighestEt",        "S9/S25",             100, 0.,1.);   
  ScHisto_sEEHighestEt          = new TH1F("ScHisto_sEEHighestEt",          "#sigma_{#eta #eta}", 100, 0.,0.1);  
  ScHisto_hoeHighestEt          = new TH1F("ScHisto_hoeHighestEt",          "H/E",                100, -0.1, 1.);
  ScHisto_dEtaTrHighestEt       = new TH1F("ScHisto_dEtaTrHighestEt",       "#Delta #eta",        100, 0., 0.3);
  ScHisto_dPhiTrHighestEt       = new TH1F("ScHisto_dPhiTrHighestEt",       "#Delta #phi",        100, 0., 0.3);

  ScHisto_deltaRHighestEt_EB    = new TH1F("ScHisto_deltaRHighestEt_EB",       "Sc deltaR",          150, 0.,5.);
  ScHisto_invMassHighestEt_EB   = new TH1F("ScHisto_invMassHighestEt_EB",      "Sc invariant mass",  150, 0.,6.);
  ScHisto_deltaRHighestEt_EE    = new TH1F("ScHisto_deltaRHighestEt_EE",       "Sc deltaR",          150, 0.,5.);
  ScHisto_invMassHighestEt_EE   = new TH1F("ScHisto_invMassHighestEt_EE",      "Sc invariant mass",  150, 0.,6.);
  ScHisto_deltaRHighestEt_EBEE  = new TH1F("ScHisto_deltaRHighestEt_EBEE",       "Sc deltaR",          150, 0.,5.);
  ScHisto_invMassHighestEt_EBEE = new TH1F("ScHisto_invMassHighestEt_EBEE",      "Sc invariant mass",  150, 0.,6.);



  ScHisto_JeneVsJeta_highestEt       = new TH2F("ScHisto_JeneVsJeta_highestEt",      "ScHisto_JeneVsJeta_highestEt",      100,  0., 2.5,    50, 5., 100.);
  ScHisto_JetVsJeta_highestEt        = new TH2F("ScHisto_JetVsJeta_highestEt",       "ScHisto_JetVsJeta_highestEt",       100,  0., 2.5,    50, 5., 30.);
  ScHisto_InvMassVsJene_highestEt    = new TH2F("ScHisto_InvMassVsJene_highestEt",   "ScHisto_InvMassVsJene_highestEt",    50,  5., 100.,   75, 0., 6.);
  ScHisto_InvMassVsJet_highestEt     = new TH2F("ScHisto_InvMassVsJet_highestEt",    "ScHisto_InvMassVsJet_highestEt",     50,  5.,  30.,   75, 0., 6.);
  ScHisto_InvMassVsJeneEB_highestEt  = new TH2F("ScHisto_InvMassVsJeneEB_highestEt", "ScHisto_InvMassVsJeneEB_highestEt",  50,  5., 100.,   75, 0., 6.);
  ScHisto_InvMassVsJeneEE_highestEt  = new TH2F("ScHisto_InvMassVsJeneEE_highestEt", "ScHisto_InvMassVsJeneEE_highestEt",  50,  5., 100.,   75, 0., 6.);
  ScHisto_InvMassVsJetEB_highestEt   = new TH2F("ScHisto_InvMassVsJetEB_highestEt",  "ScHisto_InvMassVsJetEB_highestEt",   50,  5.,  30.,   75, 0., 6.);
  ScHisto_InvMassVsJetEE_highestEt   = new TH2F("ScHisto_InvMassVsJetEE_highestEt",  "ScHisto_InvMassVsJetEE_highestEt",   50,  5.,  30.,   75, 0., 6.);
  ScHisto_InvMassVsJeta_highestEt    = new TH2F("ScHisto_InvMassVsJeta_highestEt",   "ScHisto_InvMassVsJeta_highestEt",    50, -2.5,  2.5,  75, 0., 6.);
  ScHisto_InvMassVsJphi_highestEt    = new TH2F("ScHisto_InvMassVsJphi_highestEt",   "ScHisto_InvMassVsJphi_highestEt",    50, -3.14, 3.14, 75, 0., 6.);
  ScHisto_InvMassVsDeltaR_highestEt  = new TH2F("ScHisto_InvMassVsDeltaR_highestEt", "ScHisto_InvMassVsDeltaR_highestEt",  100,  0.,   3.14,   75, 0., 6.);
}


void finalJPsiAnalysisEle::drawPlots() {

  TCanvas c; c.cd();

  ScHistoEle_size       -> Draw();  c.Print("AllEle.eps");
  ScHistoGt4_size       -> Draw();  c.Print("EleGt4.eps");
  ScHistoGt4plus_size   -> Draw();  c.Print("EleGt4plus.eps");
  ScHistoGt4minus_size  -> Draw();  c.Print("EleGt4minus.eps");

  c.SetLogy();
  HepHisto_size    -> Draw();  c.Print("HepSize.eps");

  c.SetLogy(0);
  HepHisto_deltaR  -> Draw();  c.Print("HepDeltaR.eps");
  HepHisto_maxPt   -> Draw();  c.Print("HepMaxPt.eps");
  HepHisto_minPt   -> Draw();  c.Print("HepMinPt.eps");
  HepHisto_eta     -> Draw();  c.Print("HepEta.eps");
  HepHisto_phi     -> Draw();  c.Print("HepPhi.eps");
  HepHisto_invMass -> Draw();  c.Print("HepInvMass.eps");

  if (signal) {
    c.SetLogy();
    ScHisto_deltaRWrtMc_CloserToMc_more1 -> Draw();  c.Print("dRwithMC_closer_many.eps");
    ScHisto_deltaRWrtMc_HighestEt_more1  -> Draw();  c.Print("dRwithMC_highestEt_many.eps");
    ScHisto_deltaRWrtMc_BestS9S25_more1  -> Draw();  c.Print("dRwithMC_bestS9S25_many.eps");
    ScHisto_deltaRWrtMc_BestSEE_more1    -> Draw();  c.Print("dRwithMC_bestSEE_many.eps");
    
    // superimposed, to understand
    ScHisto_deltaRWrtMc_CloserToMc_more1 -> SetLineColor(1); ScHisto_deltaRWrtMc_CloserToMc_more1 -> Draw();  
    ScHisto_deltaRWrtMc_HighestEt_more1  -> SetLineColor(2); ScHisto_deltaRWrtMc_HighestEt_more1  -> Draw("same");
    ScHisto_deltaRWrtMc_BestS9S25_more1  -> SetLineColor(3); ScHisto_deltaRWrtMc_BestS9S25_more1  -> Draw("same");  
    ScHisto_deltaRWrtMc_BestSEE_more1    -> SetLineColor(4); ScHisto_deltaRWrtMc_BestSEE_more1    -> Draw("same");  
    c.Print("dRwithMC_all_many.eps");
    
    HepHisto_PtVsEta  -> Draw("colz");
    c.Print("henPtVsEta.eps");
  }

  c.SetLogy(0);
  ScHisto_deltaRComb            -> Draw();        c.Print("scDeltaR_combinatorics.eps");           c.Print("scDeltaR_combinatorics.root");
  ScHisto_invMassComb           -> Draw();        c.Print("scInvMass_combinatorics.eps");          c.Print("scInvMass_combinatorics.root");
  ScHisto_deltaRVsInvMassComb   -> Draw("colz");  c.Print("scDeltaRVsInvMass_combinatorics.eps");  c.Print("scDeltaRVsInvMass_combinatorics.root"); 

  c.SetLogy(0);
  ScHisto_etaHighestEt          -> Draw();  c.Print("scEta.eps");
  ScHisto_phiHighestEt          -> Draw();  c.Print("scPhi.eps");
  ScHisto_maxEtHighestEt        -> Draw();  c.Print("scMaxEt.eps");
  ScHisto_minEtHighestEt        -> Draw();  c.Print("scMinEt.eps");
  ScHisto_deltaRHighestEt       -> Draw();  c.Print("scDeltaR.eps");
  //ScHisto_invMassHighestEt      -> Draw();  c.Print("scInvMass.eps");
  ScHisto_s9s25HighestEt        -> Draw();  c.Print("scS9S25.eps");
  ScHisto_sEEHighestEt          -> Draw();  c.Print("scSee.eps");
  c.SetLogy(1);
  ScHisto_hoeHighestEt          -> Draw();  c.Print("scHoE.eps");
  //c.SetLogy(0);
  ScHisto_dEtaTrHighestEt       -> Draw();  c.Print("scDetaTr.eps");
  ScHisto_dPhiTrHighestEt       -> Draw();  c.Print("scDphiTr.eps");


  // to study the invariant mass
  if (signal) {
    c.SetLogy(0);
    ScHisto_InvMassVsJene_highestEt   -> FitSlicesY();
    ScHisto_InvMassVsJet_highestEt    -> FitSlicesY();
    ScHisto_InvMassVsJeneEB_highestEt -> FitSlicesY();
    ScHisto_InvMassVsJeneEE_highestEt -> FitSlicesY();
    ScHisto_InvMassVsJetEB_highestEt  -> FitSlicesY();
    ScHisto_InvMassVsJetEE_highestEt  -> FitSlicesY();
    ScHisto_InvMassVsJeta_highestEt   -> FitSlicesY();
    ScHisto_InvMassVsJphi_highestEt   -> FitSlicesY();
    ScHisto_InvMassVsDeltaR_highestEt -> FitSlicesY();

    TH1F* GHisto_InvMassVsJene_highestEt   = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJene_highestEt_1");
    TH1F* GHisto_InvMassVsJet_highestEt    = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJet_highestEt_1");
    TH1F* GHisto_InvMassVsJeneEB_highestEt = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJeneEB_highestEt_1");
    TH1F* GHisto_InvMassVsJeneEE_highestEt = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJeneEE_highestEt_1");
    TH1F* GHisto_InvMassVsJetEB_highestEt  = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJetEB_highestEt_1");
    TH1F* GHisto_InvMassVsJetEE_highestEt  = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJetEE_highestEt_1");
    TH1F* GHisto_InvMassVsJeta_highestEt   = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJeta_highestEt_1");
    TH1F* GHisto_InvMassVsJphi_highestEt   = (TH1F*)gDirectory->Get("ScHisto_InvMassVsJphi_highestEt_1");
    TH1F* GHisto_InvMassVsDeltaR_highestEt = (TH1F*)gDirectory->Get("ScHisto_InvMassVsDeltaR_highestEt_1");

    GHisto_InvMassVsJene_highestEt   -> Draw();  c.Print("GHisto_InvMassVsJene_highestEt.eps");   c.Print("GHisto_InvMassVsJene_highestEt.root");
    GHisto_InvMassVsJet_highestEt    -> Draw();  c.Print("GHisto_InvMassVsJet_highestEt.eps");    c.Print("GHisto_InvMassVsJet_highestEt.root");
    GHisto_InvMassVsJeneEB_highestEt -> Draw();  c.Print("GHisto_InvMassVsJeneEB_highestEt.eps"); c.Print("GHisto_InvMassVsJeneEB_highestEt.root");
    GHisto_InvMassVsJeneEE_highestEt -> Draw();  c.Print("GHisto_InvMassVsJeneEE_highestEt.eps"); c.Print("GHisto_InvMassVsJeneEE_highestEt.root");
    GHisto_InvMassVsJetEB_highestEt  -> Draw();  c.Print("GHisto_InvMassVsJetEB_highestEt.eps");  c.Print("GHisto_InvMassVsJetEB_highestEt.root");
    GHisto_InvMassVsJetEE_highestEt  -> Draw();  c.Print("GHisto_InvMassVsJetEE_highestEt.eps");  c.Print("GHisto_InvMassVsJetEE_highestEt.root");
    GHisto_InvMassVsJeta_highestEt   -> Draw();  c.Print("GHisto_InvMassVsJeta_highestEt.eps");   c.Print("GHisto_InvMassVsJeta_highestEt.root");
    GHisto_InvMassVsJphi_highestEt   -> Draw();  c.Print("GHisto_InvMassVsJphi_highestEt.eps");   c.Print("GHisto_InvMassVsJphi_highestEt.root");
    GHisto_InvMassVsDeltaR_highestEt -> Draw();  c.Print("GHisto_InvMassVsDeltaR_highestEt.eps"); c.Print("GHisto_InvMassVsDeltaR_highestEt.root");

    // 
    ScHisto_JeneVsJeta_highestEt -> FitSlicesY();
    ScHisto_JetVsJeta_highestEt  -> FitSlicesY();
    TH1F* GScHisto_JeneVsJeta_highestEt = (TH1F*)gDirectory->Get("ScHisto_JeneVsJeta_highestEt_1");
    TH1F* GScHisto_JetVsJeta_highestEt  = (TH1F*)gDirectory->Get("ScHisto_JetVsJeta_highestEt_1");
    GScHisto_JeneVsJeta_highestEt -> Draw();   c.Print("GHisto_JeneVsJeta_highestEt.eps");   c.Print("GHisto_JeneVsJeta_highestEt.root");
    GScHisto_JetVsJeta_highestEt  -> Draw();   c.Print("GHisto_JetVsJeta_highestEt.eps");    c.Print("GHisto_JetVsJeta_highestEt.root");    
    }
}

void finalJPsiAnalysisEle::saveHistos() {

  TFile fOut("JPsi_histo.root", "RECREATE");
  fOut.cd();
  
  ScHistoEle_size       -> Write();
  ScHistoGt4_size       -> Write();
  ScHistoGt4plus_size   -> Write();
  ScHistoGt4minus_size  -> Write();

  HepHisto_ptHatAll   -> Write();
  for (int i=0;i<7;i++)   {
    HepHisto_ptHat[i] -> Write();
  }

  if (signal) {

    HepHisto_size   -> Write();
    HepHisto_deltaR   -> Write();
    HepHisto_maxPt    -> Write();
    HepHisto_minPt    -> Write();
    HepHisto_eta      -> Write();
    HepHisto_phi      -> Write();
    HepHisto_invMass  -> Write();
    HepHisto_PtVsEta  -> Write();
  }

  EleHisto_eta        -> Write();
  EleHisto_phi        -> Write();
  EleHisto_detaVtx    -> Write();
  EleHisto_dphiVtx    -> Write();
  EleHisto_HoE        -> Write();
  EleHisto_EoP        -> Write();
  EleHisto_fbrem      -> Write();
  EleHisto_sigmaIeta  -> Write();
  EleHisto_sigmaEta   -> Write();
  EleHisto_dr03ecal   -> Write();
  EleHisto_dr03tk     -> Write();
  EleHisto_dr03tkcorr -> Write();
  EleHisto_dr03hcal1  -> Write();
  EleHisto_dr03hcal2  -> Write();
  EleHisto_momErr     -> Write();
  EleHisto_isPF       -> Write();
  EleHisto_isEcal     -> Write();
  EleHisto_mva        -> Write();
		      
  ScHisto_deltaRComb        -> Write();
  ScHisto_invMassComb       -> Write();

  ScHisto_etaHighestEt      -> Write();
  ScHisto_phiHighestEt      -> Write();
  ScHisto_maxEtHighestEt    -> Write();
  ScHisto_minEtHighestEt    -> Write();
  ScHisto_deltaRHighestEt   -> Write();
  ScHisto_maxEtBG           -> Write();
  ScHisto_minEtBG           -> Write(); 
  ScHisto_deltaRBG          -> Write();
  ScHisto_etaBG             -> Write(); 
  ScHisto_phiBG             -> Write();
  ScHisto_maxEtBG2           -> Write();
  ScHisto_minEtBG2           -> Write(); 
  ScHisto_deltaRBG2          -> Write();
  ScHisto_etaBG2             -> Write(); 
  ScHisto_phiBG2             -> Write();

  for (int i=0;i<6;i++)   {
    ScHisto_invMassHighestEt[i]-> Write();
  }
  ScHisto_s9s25HighestEt    -> Write();
  ScHisto_sEEHighestEt      -> Write();
  ScHisto_hoeHighestEt      -> Write();
  ScHisto_dEtaTrHighestEt   -> Write();
  ScHisto_dPhiTrHighestEt   -> Write();

  ScHisto_deltaRHighestEt_EB    -> Write();
  ScHisto_invMassHighestEt_EB   -> Write();
  
  ScHisto_deltaRHighestEt_EE    -> Write();
  ScHisto_invMassHighestEt_EE   -> Write();
  
  ScHisto_deltaRHighestEt_EBEE    -> Write();
  ScHisto_invMassHighestEt_EBEE   -> Write();
  ScHisto_InvMassVsDeltaR_highestEt   -> Write();


  ScHisto_InvMassVsJene_highestEt  -> Write();         
  ScHisto_InvMassVsJeta_highestEt  -> Write();         
  ScHisto_InvMassVsJphi_highestEt  -> Write();         
  ScHisto_InvMassVsDeltaR_highestEt   -> Write();

}


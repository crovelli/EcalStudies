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

#include "finalJPsiAnalysis.hh"

using namespace std;

finalJPsiAnalysis::finalJPsiAnalysis(TTree *tree) 
  : JPsiTreeBase(tree) {
}

finalJPsiAnalysis::~finalJPsiAnalysis(){ } 

void finalJPsiAnalysis::Loop() {

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // booking histos
  bookHistos();

  // loop over events
  std::cout << "Number of entries = " << nentries << std::endl;

  // counters
  int totalEvents      = 0;
  int totalReco        = 0;
  int totalTrigger     = 0;
  int totalMatched     = 0;
  int totalMatchedMore = 0;
  int totalIdentified  = 0;
  int totalIdentifiedMore  = 0;
  int goodGene         = 0;
  int okEvent_McTruth  = 0;
  int okEvent_Highest  = 0;

  // to choose the best pair
  int okHighest_many      = 0;
  int okBestS9S25_many    = 0;
  int okBestSEE_many      = 0;
  int okBestDeta_many     = 0;
  int wrongHighest_many   = 0;
  int wrongBestS9S25_many = 0;
  int wrongBestSEE_many   = 0;
  int wrongBestDeta_many  = 0;

  for (int jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    totalEvents++;

    // ---------------------------------------------------
    // 1) plots over MC truth
    TVector3 trueEle_3P, truePos_3P;
    TLorentzVector trueEle_4P, truePos_4P;
    float maxPtGen = 0.;
    float minPtGen = 0.;
    
    if (signal) {
      
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




    // only for good number of generated
    if(!signal || numberOfGenerated==2) {
      goodGene++;
      
      // HLT
      if (hlt29) {
	totalTrigger++;
	

	// few numbers
	if (numberOfScGt4>1) totalReco++;
	if (numberOfScMatchedPlus>0 && numberOfScMatchedMinus>0) totalMatched++;
	if ((numberOfScMatchedMinus>0 && numberOfScMatchedPlus>1) || (numberOfScMatchedMinus>1 && numberOfScMatchedPlus>0)) totalMatchedMore++;
	ScHistoGt4_size -> Fill(numberOfScGt4);
	if(numberOfScGt4>1) ScHistoMatchedPlusSide_size  -> Fill(numberOfScMatchedPlus);
	if(numberOfScGt4>1) ScHistoMatchedMinusSide_size -> Fill(numberOfScMatchedMinus);
	
	
	// --------------------------------------------------------------------------------
	// 2) loop to decide which criteria to apply to select sc if more than 1 is reconstructed	
	TVector3 closerToMc3p,  highestEt_3p, bestS9S25_3p, bestSEE_3p, bestDeta_3p;
	TVector3 closerToMc3e,  highestEt_3e, bestS9S25_3e, bestSEE_3e, bestDeta_3e;
	TLorentzVector closerToMc4p,  highestEt_4p, bestS9S25_4p, bestSEE_4p, bestDeta_4p; 
	TLorentzVector closerToMc4e,  highestEt_4e, bestS9S25_4e, bestSEE_4e, bestDeta_4e; 
	//
	TVector3 highestEt_pv_3p, highestEt_tracker_3p;
	TVector3 highestEt_pv_3e, highestEt_tracker_3e;
	TLorentzVector highestEt_pv_4p, highestEt_tracker_4p; 
	TLorentzVector highestEt_pv_4e, highestEt_tracker_4e; 

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
	
	float highestEt_S9S25p = 999.; 
	float highestEt_S9S25e = 999.; 
	float highestEt_SEEp   = 999.;
	float highestEt_SEEe   = 999.;
	float highestEt_HoEp   = 999.;
	float highestEt_HoEe   = 999.;
	float highestEt_dEtap  = 999.;
	float highestEt_dEtae  = 999.;
	float highestEt_dPhip  = 999.;
	float highestEt_dPhie  = 999.;
	
	int numberOfScOkPlus  = 0;
	int numberOfScOkMinus = 0;
	for(int theSc=0; theSc<numberOfScGt4; theSc++) { 
	  
	  if (numberOfScMatchedPlus>=1 && numberOfScMatchedMinus>=1) {
	    
	    // skipping superclusters without matching to tracker
	    if (xOrRecoSc[theSc]>-700) {
	      
	      // searching for some criteria to select the two 'best' electrons (among those with Et>=4) ( analysis carried on using ECAL )	      
	      TVector3 this3P; 
	      TLorentzVector this4P;
	      this3P.SetXYZ(xOrRecoSc[theSc], yOrRecoSc[theSc], zOrRecoSc[theSc]);
	      this4P.SetPxPyPzE(xOrRecoSc[theSc], yOrRecoSc[theSc], zOrRecoSc[theSc], eneRecoSc[theSc]);
	      
	      // for further studies on the invariant mass: ECAL wrt PV
	      TVector3 this3P_pv; 
	      TLorentzVector this4P_pv;
	      this3P_pv.SetXYZ(xPVRecoSc[theSc], yPVRecoSc[theSc], zPVRecoSc[theSc]);
	      this4P_pv.SetPxPyPzE(xPVRecoSc[theSc], yPVRecoSc[theSc], zPVRecoSc[theSc], eneRecoSc[theSc]);

	      TVector3 this3P_tracker; 
	      TLorentzVector this4P_tracker;
	      this3P_tracker.SetXYZ(pxFromTrackerSc[theSc], pyFromTrackerSc[theSc], pzFromTrackerSc[theSc]);
	      this4P_tracker.SetPxPyPzE(pxFromTrackerSc[theSc], pyFromTrackerSc[theSc], pzFromTrackerSc[theSc], this3P_tracker.Mag());

	      
	      
	      // --------------------------------------------
	      // further selections
	      bool isGood   = true;
	      bool isBarrel = true;
	      if(fabs(this3P.Eta())>1.479) isBarrel = false;

	      if (isBarrel){		
		if(fabs(dEtaWithTrackerRecoSc[theSc])>0.009) isGood = false;
		if(fabs(dPhiWithTrackerRecoSc[theSc])>0.09)  isGood = false;
		if(HoverERecoSc[theSc]>0.115)                isGood = false;
		if(sigmaEtaEtaRecoSc[theSc]>0.0140)          isGood = false;
	      }
	      
	      if (!isBarrel){
		if(fabs(dEtaWithTrackerRecoSc[theSc])>0.0105) isGood = false;
		if(fabs(dPhiWithTrackerRecoSc[theSc])>0.092)  isGood = false;
		if(HoverERecoSc[theSc]>0.150)                 isGood = false;
		if(sigmaEtaEtaRecoSc[theSc]>0.0275)           isGood = false;
		}

	      if (!isGood) continue;
	      if (chargeSc[theSc]>0) numberOfScOkPlus++;
	      if (chargeSc[theSc]<0) numberOfScOkMinus++;
	      // --------------------------------------------
	      
	      
	      // search for the two reco superclusters closer to the generated electrons in case of signal
	      if (signal) {
		float deltaRmc=0;
		if(chargeSc[theSc]>0) deltaRmc = this3P.DeltaR(truePos_3P);    
		if(chargeSc[theSc]<0) deltaRmc = this3P.DeltaR(trueEle_3P);    
		
		if((chargeSc[theSc]>0) && (deltaRmc<deltaRposi)) { 
		  theDeltaRp=theSc;
		  deltaRposi=deltaRmc; 
		  closerToMc3p=this3P; 
		  closerToMc4p=this4P; 
		}
		
		if((chargeSc[theSc]<0) && (deltaRmc<deltaRelec)) { 
		  theDeltaRe=theSc;
		  deltaRelec=deltaRmc; 
		  closerToMc3e=this3P; 
		  closerToMc4e=this4P; 
		}
	      }
	      
	      
	      // search for highest Et superclusters which charge hp = 1
	      if ( (chargeSc[theSc]>0) && (etRecoSc[theSc]>=highestEtp) ){ 
		theHighestEtp        = theSc;
		highestEtp           = etRecoSc[theSc];
		highestEt_3p         = this3P;
		highestEt_4p         = this4P;
		highestEt_pv_3p      = this3P_pv;
		highestEt_pv_4p      = this4P_pv;
		highestEt_tracker_3p = this3P_tracker;
		highestEt_tracker_4p = this4P_tracker;
		highestEt_S9S25p     = e9e25RecoSc[theSc];
		highestEt_SEEp       = sigmaEtaEtaRecoSc[theSc];
		highestEt_HoEp       = HoverERecoSc[theSc];
		highestEt_dEtap      = dEtaWithTrackerRecoSc[theSc];
		highestEt_dPhip      = dPhiWithTrackerRecoSc[theSc];
	      }
	      
	      // search for highest Et superclusters which charge hp = -1
	      if ( (chargeSc[theSc]<0) && (etRecoSc[theSc]>=highestEte) ){
		theHighestEte        = theSc;
		highestEte           = etRecoSc[theSc];
		highestEt_3e         = this3P;
		highestEt_4e         = this4P;
		highestEt_pv_3e      = this3P_pv;
		highestEt_pv_4e      = this4P_pv;
		highestEt_tracker_3e = this3P_tracker;
		highestEt_tracker_4e = this4P_tracker;
		highestEt_S9S25e     = e9e25RecoSc[theSc];
		highestEt_SEEe       = sigmaEtaEtaRecoSc[theSc];
		highestEt_HoEe       = HoverERecoSc[theSc];
		highestEt_dEtae      = dEtaWithTrackerRecoSc[theSc];
		highestEt_dPhie      = dPhiWithTrackerRecoSc[theSc];
	      }
	      
	      
	      // search for best E9/E25 electrons with charge Hp = +1
	      float s9s25diff = fabs(e9e25RecoSc[theSc] - 1.);
	      if ( (chargeSc[theSc]>0) && (s9s25diff<=bestS9S25p) ){
		theBestS9S25p = theSc;
		bestS9S25p    = s9s25diff;
		bestS9S25_3p  = this3P;
		bestS9S25_4p  = this4P;
	      }
	      
	      // search for best E9/E25 electrons with charge Hp = -1
	      if ( (chargeSc[theSc]<0) && (s9s25diff<=bestS9S25e) ){
		theBestS9S25e = theSc;
		bestS9S25e    = s9s25diff;
		bestS9S25_3e  = this3P;
		bestS9S25_4e  = this4P;
	      }
	      
	      
	      // search for best sigmaEtaEta electrons with charge Hp = +1
	      float sEEdiff = fabs(sigmaEtaEtaRecoSc[theSc]);
	      if ( (chargeSc[theSc]>0) && (sEEdiff<=bestSEEp) ){
		theBestSEEp = theSc;
		bestSEEp    = sEEdiff;
		bestSEE_3p  = this3P;
		bestSEE_4p  = this4P;
	      }
	      
	      // search for best E9/E25 electrons with charge Hp = -1
	      if ( (chargeSc[theSc]<0) && (sEEdiff<=bestSEEe) ){
		theBestSEEe = theSc;
		bestSEEe    = sEEdiff;
		bestSEE_3e  = this3P;
		bestSEE_4e  = this4P;
	      }
	      
	      // search for best dEta electrons with charge Hp = +1	      
	      float absDeta = fabs(dEtaWithTrackerRecoSc[theSc]);
	      if ( (chargeSc[theSc]>0) && (absDeta<=bestDetap) ){
		theBestDetap = theSc;
		bestDetap    = absDeta;
		bestDeta_3p  = this3P;
		bestDeta_4p  = this4P;
	      }
	      
	      // search for best dEta electrons with charge Hp = +1	      
	      if ( (chargeSc[theSc]<0) && (absDeta<=bestDetae) ){
		theBestDetae = theSc;
		bestDetae    = absDeta;
		bestDeta_3e  = this3P;
		bestDeta_4e  = this4P;
	      }

	    }  // ok matching with track
	    
	  } // two matched electons
	  
	} // loop over superclusters
	

	if (numberOfScOkPlus>0 && numberOfScOkMinus>0) totalIdentified++;
	if ((numberOfScOkPlus>0 && numberOfScOkMinus>1) || (numberOfScOkPlus>1 && numberOfScOkMinus>0)) totalIdentifiedMore++;
	
	
	// to choose, signal only
	if (signal) {
	  
	  if (numberOfScOkPlus>=1 && numberOfScOkMinus>=1) {
	    
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

	    if (numberOfScOkPlus>1) {
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
	    
	    if (numberOfScOkMinus>1) {
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
	
      
	// plots with all combinatorics
	if ( (numberOfScMatchedMinus+numberOfScMatchedPlus)>1 ){ 
	  for(int theSc1=0; theSc1<numberOfScGt4; theSc1++) { 
	    if (xOrRecoSc[theSc1]<-700) continue;
	    for(int theSc2=(theSc1+1); theSc2<numberOfScGt4; theSc2++) { 
	      if (xOrRecoSc[theSc2]<-700) continue;

		TLorentzVector tlvTheSc1, tlvTheSc2;
		TVector3 tv3TheSc1, tv3TheSc2;
		tlvTheSc1.SetPxPyPzE(xOrRecoSc[theSc1], yOrRecoSc[theSc1], zOrRecoSc[theSc1], eneRecoSc[theSc1]);
		tlvTheSc2.SetPxPyPzE(xOrRecoSc[theSc2], yOrRecoSc[theSc2], zOrRecoSc[theSc2], eneRecoSc[theSc2]);
		tv3TheSc1.SetXYZ (xOrRecoSc[theSc1], yOrRecoSc[theSc1], zOrRecoSc[theSc1]);
		tv3TheSc2.SetXYZ (xOrRecoSc[theSc2], yOrRecoSc[theSc2], zOrRecoSc[theSc2]);
		
		float deltaR = tv3TheSc1.DeltaR(tv3TheSc2);
		float mee    = (tlvTheSc1 + tlvTheSc2).M();
		ScHisto_deltaRComb->Fill(deltaR);
		ScHisto_invMassComb->Fill(mee);
		ScHisto_deltaRVsInvMassComb->Fill(mee,deltaR);
	    }}
	}
	

	// to compare the distributions: use the two highest Et electrons
	if (numberOfScMatchedPlus>=1 && numberOfScMatchedMinus>=1 && numberOfScOkPlus>=1 && numberOfScOkMinus>=1) {
	  
	  ScHisto_etaHighestEt->Fill(highestEt_3p.Eta());
	  ScHisto_etaHighestEt->Fill(highestEt_3e.Eta());
	  ScHisto_phiHighestEt->Fill(highestEt_3p.Phi());
	  ScHisto_phiHighestEt->Fill(highestEt_3e.Phi());
	  
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
	  ScHisto_maxEtHighestEt->Fill(maxEt_highestEt);
	  ScHisto_minEtHighestEt->Fill(minEt_highestEt);
	  
	  float deltaR_highestEt = highestEt_3p.DeltaR(highestEt_3e);
	  ScHisto_deltaRHighestEt->Fill(deltaR_highestEt);	
	  
	  float mee_highestEt = (highestEt_4p + highestEt_4e).M();
	  ScHisto_invMassHighestEt->Fill(mee_highestEt);

	  /*
	    ScHisto_s9s25HighestEt->Fill(highestEt_S9S25p);
	    ScHisto_s9s25HighestEt->Fill(highestEt_S9S25e);
	    ScHisto_sEEHighestEt->Fill(highestEt_SEEp);
	    ScHisto_sEEHighestEt->Fill(highestEt_SEEe);
	    ScHisto_hoeHighestEt->Fill(highestEt_HoEp);
	    ScHisto_hoeHighestEt->Fill(highestEt_HoEe);	    
	    ScHisto_dEtaTrHighestEt -> Fill(highestEt_dEtap);
	    ScHisto_dEtaTrHighestEt -> Fill(highestEt_dEtae);
	    ScHisto_dPhiTrHighestEt -> Fill(highestEt_dPhip);
	    ScHisto_dPhiTrHighestEt -> Fill(highestEt_dPhie);
	  */
	  
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
	  
	  // further invariant mass studies
	  float mee_highestEt_pv = (highestEt_pv_4p + highestEt_pv_4e).M();
	  ScHisto_pv_invMassHighestEt -> Fill(mee_highestEt_pv);	  	  
	  if ( fabs(highestEt_pv_3p.Eta())<1.479 && fabs(highestEt_pv_3e.Eta())<1.479 ) ScHisto_pv_invMassHighestEt_EB -> Fill(mee_highestEt_pv);
	  if ( fabs(highestEt_pv_3p.Eta())>1.479 && fabs(highestEt_pv_3e.Eta())>1.479 ) ScHisto_pv_invMassHighestEt_EE -> Fill(mee_highestEt_pv);
	  if ( (fabs(highestEt_pv_3p.Eta())<1.479 && fabs(highestEt_pv_3e.Eta())>1.479) || (fabs(highestEt_pv_3e.Eta())<1.479 && fabs(highestEt_pv_3p.Eta())>1.479)) 
	    ScHisto_pv_invMassHighestEt_EBEE -> Fill(mee_highestEt_pv);
	  
	  float mee_highestEt_tracker = (highestEt_tracker_4p + highestEt_tracker_4e).M();
	  ScHisto_tracker_invMassHighestEt -> Fill(mee_highestEt_tracker);	  	  
	  if ( fabs(highestEt_tracker_3p.Eta())<1.479 && fabs(highestEt_tracker_3e.Eta())<1.479 ) ScHisto_tracker_invMassHighestEt_EB -> Fill(mee_highestEt_tracker);
	  if ( fabs(highestEt_tracker_3p.Eta())>1.479 && fabs(highestEt_tracker_3e.Eta())>1.479 ) ScHisto_tracker_invMassHighestEt_EE -> Fill(mee_highestEt_tracker);
	  if ( (fabs(highestEt_tracker_3p.Eta())<1.479 && fabs(highestEt_tracker_3e.Eta())>1.479) || (fabs(highestEt_tracker_3e.Eta())<1.479 && fabs(highestEt_tracker_3p.Eta())>1.479)) 
	    ScHisto_tracker_invMassHighestEt_EBEE -> Fill(mee_highestEt_tracker);


	  

	  // study of possible biases
	  if (signal) {    
	    float jpsi_ene    = (highestEt_4p + highestEt_4e).E();	  // chiara, check se funziona cosi'
	    float jpsi_eta    = (highestEt_4p + highestEt_4e).Eta();	  
	    float jpsi_phi    = (highestEt_4p + highestEt_4e).Phi();	  
	    //float cosTheta    = cos(highestEt_3p.Theta(highestEt_3e)); 
	    //float cosThetaGen = cos(trueEle_3P.Theta(truePos_3P));
	    
	    ScHisto_InvMassVsJene_highestEt -> Fill(jpsi_ene, mee_highestEt);        
	    ScHisto_InvMassVsJeta_highestEt -> Fill(jpsi_eta, mee_highestEt);                
	    ScHisto_InvMassVsJphi_highestEt -> Fill(jpsi_phi, mee_highestEt);                
	    // ScHisto_thRecOverThTrueVsJene_highestEt -> Fill(cosTheta/cosThetaGen, jpsi_ene);
	    if (highestEt_4e.E() > highestEt_4p.E()) { 
	      ScHisto_EoEtrueVsErec_highestEt1 -> Fill(highestEt_4e.E(), highestEt_4e.E()/maxPtGen);
	      ScHisto_EoEtrueVsErec_highestEt2 -> Fill(highestEt_4p.E(), highestEt_4p.E()/minPtGen);
	    } else {
	      ScHisto_EoEtrueVsErec_highestEt1 -> Fill(highestEt_4p.E(), highestEt_4p.E()/maxPtGen);
	      ScHisto_EoEtrueVsErec_highestEt2 -> Fill(highestEt_4e.E(), highestEt_4e.E()/minPtGen);
	    }
	  }
	  

	  /*
	    // correlations                             
	    ScHisto_deltaRVsDeltaPhi_HighestEt  -> Fill(highestEt_dPhip, deltaR_highestEt);
	    ScHisto_deltaRVsDeltaPhi_HighestEt  -> Fill(highestEt_dPhie, deltaR_highestEt);
	    ScHisto_deltaRVsDeltaEta_HighestEt  -> Fill(highestEt_dEtap, deltaR_highestEt);
	    ScHisto_deltaRVsDeltaEta_HighestEt  -> Fill(highestEt_dEtae, deltaR_highestEt);
	    ScHisto_deltaRVsHoE_HighestEt       -> Fill(highestEt_HoEp,  deltaR_highestEt);
	    ScHisto_deltaRVsHoE_HighestEt       -> Fill(highestEt_HoEe,  deltaR_highestEt);
	    ScHisto_deltaRVsSee_HighestEt       -> Fill(highestEt_SEEp,  deltaR_highestEt);    
	    ScHisto_deltaRVsSee_HighestEt       -> Fill(highestEt_SEEe,  deltaR_highestEt);    
	    ScHisto_deltaRVsInvMass_HighestEt   -> Fill(mee_highestEt,   deltaR_highestEt);    
	  */
	  

	  
	} // ok matched
      } // ok HLT
      
    } // ok generated
    

    
    
  } // loop over entries

  cout << endl;
  cout << "to choose the best criterium" << endl;
  cout << "highest ET: ok    = " << okHighest_many   << ", wrong = " << wrongHighest_many   << endl;
  cout << "best S9/S25: ok   = " << okBestS9S25_many << ", wrong = " << wrongBestS9S25_many << endl;
  cout << "best SEE: ok      = " << okBestSEE_many   << ", wrong = " << wrongBestSEE_many   << endl;
  cout << "best deltaEta: ok = " << okBestDeta_many  << ", wrong = " << wrongBestDeta_many  << endl;
  cout << endl;
  

  saveHistos();

  drawPlots();


  // summary
  cout << "total number of events = "              << totalEvents     << endl;
  cout << "good number of gene MC = "              << goodGene        << endl;
  cout << "total number of events passing HLT = "  << totalTrigger    << endl;
  cout << "total number of reco = "                << totalReco       << endl;
  cout << "total number of reco & matched = "      << totalMatched    << endl;
  cout << "total number of reco & matched & id = " << totalIdentified << endl;
  cout << "total number of reco & matched & id, more than 2 = " << totalIdentifiedMore << endl;
  cout << "total number of reco & matched, more than 2 sc = " << totalMatchedMore << endl;
  cout << endl;
  cout << "total number of events with closer to MC match MC within dR=0.5 = " << okEvent_McTruth << endl; 
  cout << "total number of events with highest Et match MC within dR=0.5 = "   << okEvent_Highest << endl; 

} // end of program


void finalJPsiAnalysis::bookHistos() {

  ScHistoGt4_size              = new TH1F("ScHistoGt4_size",  "Num of superclusters with Et > 4", 10, 0,10);
  ScHistoMatchedPlusSide_size  = new TH1F("ScHistoMatchedPlusSide_size", "Num of superclusters with Et>1 matching a track with charge plus",  10, 0,10);
  ScHistoMatchedMinusSide_size = new TH1F("ScHistoMatchedMinusSide_size","Num of superclusters with Et>1 matching a track with charge minus", 10, 0,10);

  HepHisto_size    = new TH1F("HepHisto_size",    "Num of GenElectrons", 10, 0,10);
  HepHisto_deltaR  = new TH1F("HepHisto_deltaR",  "GenElectrons deltaR", 100, 0.,1.);
  HepHisto_maxPt   = new TH1F("HepHisto_maxPt",   "GenElectrons max pt", 200, 0.,20.);
  HepHisto_minPt   = new TH1F("HepHisto_minPt",   "GenElectrons min pt", 200, 0.,20.);
  HepHisto_eta     = new TH1F("HepHisto_eta",     "GenElectrons eta", 100, -3.,+3.);
  HepHisto_phi     = new TH1F("HepHisto_phi",     "GenElectrons phi", 100, -3.15,3.15);
  HepHisto_invMass = new TH1F("HepHisto_invMass", "GenElectrons invariant mass", 150,0.,6.); 
  HepHisto_PtVsEta = new TH2F("HepHisto_PtVsEta", "GenElectrons Pt vs Eta", 100, -3.,+3., 200, 0.,20.);

  ScHisto_deltaRWrtMc_CloserToMc_more1      = new TH1F("ScHisto_deltaRWrtMc_CloserToMc_more1", "deltaR with Mc - closer to MC",     100, 0.,5.);
  ScHisto_deltaRWrtMc_HighestEt_more1       = new TH1F("ScHisto_deltaRWrtMc_HighestEt_more1",  "deltaR with Mc - highest Et",       100, 0.,5.);
  ScHisto_deltaRWrtMc_BestS9S25_more1       = new TH1F("ScHisto_deltaRWrtMc_BestS9S25_more1",  "deltaR with Mc - best S9/S25",      100, 0.,5.);
  ScHisto_deltaRWrtMc_BestSEE_more1         = new TH1F("ScHisto_deltaRWrtMc_BestSEE_more1",    "deltaR with Mc - best SigmaEtaEta", 100, 0.,5.);

  ScHisto_deltaRComb            = new TH1F("ScHisto_deltaRComb",            "Sc deltaR",          100, 0.,3.);   
  ScHisto_invMassComb           = new TH1F("ScHisto_invMassComb",           "Sc invariant mass",  150, 0.,6.);   
  ScHisto_deltaRVsInvMassComb   = new TH2F("ScHisto_deltaRVsInvMassComb",   "Sc deltaR vs Invariant Mass", 150, 0.,6., 100, 0.,3.);

  ScHisto_etaHighestEt          = new TH1F("ScHisto_etaHighestEt",          "Sc eta",             100, -3.,+3.);
  ScHisto_phiHighestEt          = new TH1F("ScHisto_phiHighestEt",          "Sc phi",             100, -3.15,+3.15);
  ScHisto_maxEtHighestEt        = new TH1F("ScHisto_maxEtHighestEt",        "Sc max Et",           50, 0.,20.);
  ScHisto_minEtHighestEt        = new TH1F("ScHisto_minEtHighestEt",        "Sc min Et",           50, 0.,20.);
  ScHisto_deltaRHighestEt       = new TH1F("ScHisto_deltaRHighestEt",       "Sc deltaR",          100, 0.,5.);   
  ScHisto_invMassHighestEt      = new TH1F("ScHisto_invMassHighestEt",      "Sc invariant mass",  150, 0.,6.);   
  ScHisto_s9s25HighestEt        = new TH1F("ScHisto_s9s25HighestEt",        "S9/S25",             100, 0.,1.);   
  ScHisto_sEEHighestEt          = new TH1F("ScHisto_sEEHighestEt",          "#sigma_{#eta #eta}", 100, 0.,0.1);  
  ScHisto_hoeHighestEt          = new TH1F("ScHisto_hoeHighestEt",          "H/E",                100, -0.1, 1.);
  ScHisto_dEtaTrHighestEt       = new TH1F("ScHisto_dEtaTrHighestEt",       "#Delta #eta",        100, 0., 0.3);
  ScHisto_dPhiTrHighestEt       = new TH1F("ScHisto_dPhiTrHighestEt",       "#Delta #phi",        100, 0., 0.3);

  ScHisto_deltaRHighestEt_EB    = new TH1F("ScHisto_deltaRHighestEt_EB",       "Sc deltaR",          100, 0.,5.);   
  ScHisto_invMassHighestEt_EB   = new TH1F("ScHisto_invMassHighestEt_EB",      "Sc invariant mass",  150, 0.,6.);   
  ScHisto_deltaRHighestEt_EE    = new TH1F("ScHisto_deltaRHighestEt_EE",       "Sc deltaR",          100, 0.,5.);
  ScHisto_invMassHighestEt_EE   = new TH1F("ScHisto_invMassHighestEt_EE",      "Sc invariant mass",  150, 0.,6.);
  ScHisto_deltaRHighestEt_EBEE  = new TH1F("ScHisto_deltaRHighestEt_EBEE",       "Sc deltaR",          100, 0.,5.);
  ScHisto_invMassHighestEt_EBEE = new TH1F("ScHisto_invMassHighestEt_EBEE",      "Sc invariant mass",  150, 0.,6.);

  ScHisto_pv_invMassHighestEt           = new TH1F("ScHisto_pv_invMassHighestEt",           "Sc invariant mass",  150, 0.,6.);   
  ScHisto_pv_invMassHighestEt_EB        = new TH1F("ScHisto_pv_invMassHighestEt_EB",        "Sc invariant mass",  150, 0.,6.);   
  ScHisto_pv_invMassHighestEt_EE        = new TH1F("ScHisto_pv_invMassHighestEt_EE",        "Sc invariant mass",  150, 0.,6.);
  ScHisto_pv_invMassHighestEt_EBEE      = new TH1F("ScHisto_pv_invMassHighestEt_EBEE",      "Sc invariant mass",  150, 0.,6.);
  ScHisto_tracker_invMassHighestEt      = new TH1F("ScHisto_tracker_invMassHighestEt",      "Sc invariant mass",  150, 0.,6.);   
  ScHisto_tracker_invMassHighestEt_EB   = new TH1F("ScHisto_tracker_invMassHighestEt_EB",   "Sc invariant mass",  150, 0.,6.);   
  ScHisto_tracker_invMassHighestEt_EE   = new TH1F("ScHisto_tracker_invMassHighestEt_EE",   "Sc invariant mass",  150, 0.,6.);
  ScHisto_tracker_invMassHighestEt_EBEE = new TH1F("ScHisto_tracker_invMassHighestEt_EBEE", "Sc invariant mass",  150, 0.,6.);

  ScHisto_deltaRVsDeltaPhi_HighestEt = new TH2F("ScHisto_deltaRVsDeltaPhi_HighestEt", "sc, e+e-deltaR vs dphi",                100, 0., 0.05,  100, 0.,1.);
  ScHisto_deltaRVsDeltaEta_HighestEt = new TH2F("ScHisto_deltaRVsDeltaEta_HighestEt", "sc, e+e-deltaR vs deta",                100, 0., 0.05,  100, 0.,1.);
  ScHisto_deltaRVsHoE_HighestEt      = new TH2F("ScHisto_deltaRVsHoE_HighestEt",      "sc, e+e-deltaR vs H/E",                 100, 0., 0.001, 100, 0.,1.);
  ScHisto_deltaRVsSee_HighestEt      = new TH2F("ScHisto_deltaRVsSee_HighestEt",      "sc, e+e-deltaR vs #sigma_{#eta #eta}",  100, 0., 0.04,  100, 0.,1.);
  ScHisto_deltaRVsInvMass_HighestEt  = new TH2F("ScHisto_deltaRVsInvMass_HighestEt",  "sc, e+e-deltaR vs invariant mass",      150, 0., 6.,    100, 0.,1.);

  ScHisto_InvMassVsJene_highestEt         = new TH2F("ScHisto_InvMassVsJene_highestEt",         "ScHisto_InvMassVsJene_highestEt",         100, 0.,  30., 150, 0., 6.);
  ScHisto_InvMassVsJeta_highestEt         = new TH2F("ScHisto_InvMassVsJeta_highestEt",         "ScHisto_InvMassVsJeta_highestEt",         100, 0., 2.5,  150, 0., 6.);
  ScHisto_InvMassVsJphi_highestEt         = new TH2F("ScHisto_InvMassVsJphi_highestEt",         "ScHisto_InvMassVsJphi_highestEt",         100, 0., 3.14, 150, 0., 6.);
  ScHisto_thRecOverThTrueVsJene_highestEt = new TH2F("ScHisto_thRecOverThTrueVsJene_highestEt", "ScHisto_thRecOverThTrueVsJene_highestEt", 100, 0.,  1.,  150, 0., 6.);
  ScHisto_EoEtrueVsErec_highestEt1        = new TH2F("ScHisto_EoEtrueVsErec_highestEt1",        "ScHisto_EoEtrueVsErec_highestEt1",        100, 0.,  30., 100, 0., 1.);
  ScHisto_EoEtrueVsErec_highestEt2        = new TH2F("ScHisto_EoEtrueVsErec_highestEt2",        "ScHisto_EoEtrueVsErec_highestEt2",        100, 0.,  30., 100, 0., 1.);
}

void finalJPsiAnalysis::drawPlots() {

  TCanvas c; c.cd();

  ScHistoGt4_size              -> Draw();  c.Print("ScGt4.eps");
  ScHistoMatchedPlusSide_size  -> Draw();  c.Print("ScMatchedPlus.eps");
  ScHistoMatchedMinusSide_size -> Draw();  c.Print("ScMatchedMinus.eps");

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


  // to study biases
  if (signal) {
    c.SetLogy(0);
    TProfile *ScP_InvMassVsJene_highestEt  = ScHisto_InvMassVsJene_highestEt  -> ProfileX();
    TProfile *ScP_InvMassVsJeta_highestEt  = ScHisto_InvMassVsJeta_highestEt  -> ProfileX();
    TProfile *ScP_InvMassVsJphi_highestEt  = ScHisto_InvMassVsJphi_highestEt  -> ProfileX();
    TProfile *ScP_EoEtrueVsErec_highestEt1 = ScHisto_EoEtrueVsErec_highestEt1 -> ProfileX();
    TProfile *ScP_EoEtrueVsErec_highestEt2 = ScHisto_EoEtrueVsErec_highestEt2 -> ProfileX();

    ScHisto_InvMassVsJene_highestEt  -> Draw();  c.Print("H_InvMassVsJpsiEnergy.eps");
    ScP_InvMassVsJene_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiEnergy.eps");
    ScP_InvMassVsJene_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiEnergy.root");
    ScHisto_InvMassVsJeta_highestEt  -> Draw();  c.Print("H_InvMassVsJpsiEta.eps");
    ScP_InvMassVsJeta_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiEta.eps");
    ScP_InvMassVsJeta_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiEta.root");
    ScHisto_InvMassVsJphi_highestEt  -> Draw();  c.Print("H_InvMassVsJpsiPhi.eps");
    ScP_InvMassVsJphi_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiPhi.eps");
    ScP_InvMassVsJphi_highestEt      -> Draw();  c.Print("P_InvMassVsJpsiPhi.root");
    ScHisto_EoEtrueVsErec_highestEt1 -> Draw();  c.Print("H_EoEtrueVsErec1.eps");
    ScP_EoEtrueVsErec_highestEt1     -> Draw();  c.Print("P_EoEtrueVsErec1.eps");
    ScP_EoEtrueVsErec_highestEt1     -> Draw();  c.Print("P_EoEtrueVsErec1.root");
    ScHisto_EoEtrueVsErec_highestEt2 -> Draw();  c.Print("H_EoEtrueVsErec2.eps");
    ScP_EoEtrueVsErec_highestEt2     -> Draw();  c.Print("P_EoEtrueVsErec2.eps");
    ScP_EoEtrueVsErec_highestEt2     -> Draw();  c.Print("P_EoEtrueVsErec2.root");
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
  ScHisto_invMassHighestEt      -> Draw();  c.Print("scInvMass.eps");
  ScHisto_s9s25HighestEt        -> Draw();  c.Print("scS9S25.eps");
  ScHisto_sEEHighestEt          -> Draw();  c.Print("scSee.eps");
  c.SetLogy(1);
  ScHisto_hoeHighestEt          -> Draw();  c.Print("scHoE.eps");
  //c.SetLogy(0);
  ScHisto_dEtaTrHighestEt       -> Draw();  c.Print("scDetaTr.eps");
  ScHisto_dPhiTrHighestEt       -> Draw();  c.Print("scDphiTr.eps");

  c.SetLogy(0);
  ScHisto_deltaRHighestEt_EB    -> Draw();  c.Print("scDeltaR_EB.eps");
  ScHisto_invMassHighestEt_EB   -> Draw();  c.Print("scInvMass_EB.eps");

  c.SetLogy(0);
  ScHisto_deltaRHighestEt_EE    -> Draw();  c.Print("scDeltaR_EE.eps");
  ScHisto_invMassHighestEt_EE   -> Draw();  c.Print("scInvMass_EE.eps");

  c.SetLogy(0);
  ScHisto_deltaRHighestEt_EBEE   -> Draw();  c.Print("scDeltaR_EBEE.eps");
  ScHisto_invMassHighestEt_EBEE  -> Draw();  c.Print("scInvMass_EBEE.eps");

  c.SetLogy(0);
  ScHisto_pv_invMassHighestEt           -> Draw();  c.Print("scInvMass_pv.eps");          
  ScHisto_pv_invMassHighestEt_EB        -> Draw();  c.Print("scInvMass_EB_pv.eps");          
  ScHisto_pv_invMassHighestEt_EE        -> Draw();  c.Print("scInvMass_EE_pv.eps");          
  ScHisto_pv_invMassHighestEt_EBEE      -> Draw();  c.Print("scInvMass_EBEE_pv.eps");           
  ScHisto_tracker_invMassHighestEt      -> Draw();  c.Print("scInvMass_tracker.eps");              
  ScHisto_tracker_invMassHighestEt_EB   -> Draw();  c.Print("scInvMass_EB_tracker.eps");              
  ScHisto_tracker_invMassHighestEt_EE   -> Draw();  c.Print("scInvMass_EE_tracker.eps");                 
  ScHisto_tracker_invMassHighestEt_EBEE -> Draw();  c.Print("scInvMass_EBEE_tracker.eps");               

  c.SetLogy(0);
  ScHisto_deltaRVsDeltaPhi_HighestEt -> Draw("colz");  c.Print("deltaRVsDeltaPhi.eps");
  ScHisto_deltaRVsDeltaEta_HighestEt -> Draw("colz");  c.Print("deltaRVsDeltaEta.eps");
  ScHisto_deltaRVsHoE_HighestEt      -> Draw("colz");  c.Print("deltaRVsHoE.eps");
  ScHisto_deltaRVsSee_HighestEt      -> Draw("colz");  c.Print("deltaRVsSee.eps"); 
  ScHisto_deltaRVsInvMass_HighestEt  -> Draw("colz");  c.Print("deltaRVsInvMass.eps");   
}

void finalJPsiAnalysis::saveHistos() {

  TFile fOut("Outfile_histo.root", "RECREATE");
  fOut.cd();
  
  ScHistoGt4_size              -> Write();
  ScHistoMatchedPlusSide_size  -> Write();
  ScHistoMatchedMinusSide_size -> Write();

  if (signal) {
    /*
    HepHisto_deltaR   -> Write();
    HepHisto_maxPt    -> Write();
    HepHisto_minPt    -> Write();
    HepHisto_eta      -> Write();
    HepHisto_phi      -> Write();
    HepHisto_invMass  -> Write();
    HepHisto_PtVsEta  -> Write();
    */
  }

  /*
  ScHisto_deltaRComb            -> Write();
  ScHisto_invMassComb           -> Write();

  ScHisto_etaHighestEt          -> Write();
  ScHisto_phiHighestEt          -> Write();
  ScHisto_maxEtHighestEt        -> Write();
  ScHisto_minEtHighestEt        -> Write();
  ScHisto_deltaRHighestEt       -> Write();
  ScHisto_invMassHighestEt      -> Write();
  ScHisto_s9s25HighestEt        -> Write();
  ScHisto_sEEHighestEt          -> Write();
  ScHisto_hoeHighestEt          -> Write();
  ScHisto_dEtaTrHighestEt       -> Write();
  ScHisto_dPhiTrHighestEt       -> Write();

  ScHisto_deltaRHighestEt_EB    -> Write();   
  ScHisto_invMassHighestEt_EB   -> Write();   
  ScHisto_s9s25HighestEt_EB     -> Write();   
  ScHisto_sEEHighestEt_EB       -> Write();   
  ScHisto_hoeHighestEt_EB       -> Write();   
  ScHisto_dEtaTrHighestEt_EB    -> Write();   
  ScHisto_dPhiTrHighestEt_EB    -> Write();   

  ScHisto_deltaRHighestEt_EE    -> Write();      
  ScHisto_invMassHighestEt_EE   -> Write();      
  ScHisto_s9s25HighestEt_EE     -> Write();      
  ScHisto_sEEHighestEt_EE       -> Write();      
  ScHisto_hoeHighestEt_EE       -> Write();      
  ScHisto_dEtaTrHighestEt_EE    -> Write();      
  ScHisto_dPhiTrHighestEt_EE    -> Write();      

  ScHisto_deltaRHighestEt_EBEE    -> Write();    
  ScHisto_invMassHighestEt_EBEE   -> Write();    
  ScHisto_s9s25HighestEt_EBEE     -> Write();    
  ScHisto_sEEHighestEt_EBEE       -> Write();    
  ScHisto_hoeHighestEt_EBEE       -> Write();    
  ScHisto_dEtaTrHighestEt_EBEE    -> Write();    
  ScHisto_dPhiTrHighestEt_EBEE    -> Write();    

  ScHisto_InvMassVsJene_highestEt  -> Write();         
  ScHisto_InvMassVsJeta_highestEt  -> Write();         
  ScHisto_InvMassVsJphi_highestEt  -> Write();         
  ScHisto_thRecOverThTrueVsJene_highestEt -> Write();         
  ScHisto_EoEtrueVsErec_highestEt1 -> Write();        
  ScHisto_EoEtrueVsErec_highestEt2 -> Write();        
  */
}


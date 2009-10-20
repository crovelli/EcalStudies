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
  float totalEvents         = 0.;
  float totalReco           = 0.;
  float totalRecoGt4        = 0.;
  float totalRecoGt4Charge  = 0.;
  float totalTrigger        = 0.;
  float totalIdentified     = 0.;
  float numberOfPairsOk     = 0.;
  float numbersOfInvMassOk  = 0.;
  float totalIdentifiedMore = 0.;
  float goodGene            = 0.;
  float okEvent_McTruth     = 0.;
  float okEvent_Highest     = 0.;

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
    if (jentry ==1) std::cout << ">>> Signal is " << signal << std::endl;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    if (theSample==2) signal = 0;  // not-prompt J/psi have the wrong signal values ...
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
    
    /// if it is signal sample
    if (signal) {
      
      HepHisto_size->Fill(numberOfGenerated);
      HepHisto_ptHat->Fill(pthat); 

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
    if(!signal || numberOfGenerated==2) {
      goodGene++;
      
      // selecting HLTbit 
      if (hltJpsi ) {
	totalTrigger++;
	
	// few numbers: all electrons
	if (numberOfElectrons>1) totalReco++;
	ScHistoEle_size -> Fill(numberOfElectrons);
	
	// few numbers: all electrons, Et > 4
	for(int theEle=0; theEle<numberOfElectrons; theEle++) { 
	  if (etRecoEle[theEle]>4) totEleGt4++;
	  if (etRecoEle[theEle]>4 && chargeRecoEle[theEle]>0) totEleGt4plus++;
	  if (etRecoEle[theEle]>4 && chargeRecoEle[theEle]<0) totEleGt4minus++;
	}
	ScHistoGt4_size -> Fill(totEleGt4);
	if(totEleGt4plus>=1)  ScHistoGt4plus_size  -> Fill(totEleGt4plus);
	if(totEleGt4minus>=1) ScHistoGt4minus_size -> Fill(totEleGt4minus);
	if(totEleGt4>1) totalRecoGt4++;
	if(totEleGt4plus>=1 && totEleGt4minus>=1) totalRecoGt4Charge++;
	
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
	      if (etRecoEle[theEle]>4) {


		// searching for some criteria to select the two 'best' electrons (among those with Et>=4) ( analysis carried on using ECAL )	      
		TVector3 this3P; 
		TLorentzVector this4P;
		this3P.SetXYZ(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle]);
		this4P.SetPxPyPzE(pxRecoEle[theEle], pyRecoEle[theEle], pzRecoEle[theEle], eneRecoEle[theEle]);
		
		// further selections
		bool isGood   = true;
		bool isBarrel = true;
		bool isPF = false;

		if(isPFlowRecoEle[theEle]) isPF = true;
		if(fabs(this3P.Eta())>1.479) isBarrel = false;
		

		fillHistoCutsVariables( theEle );


		if (isBarrel){		
		  if( fabs(dEtaAtVtxRecoEle[theEle]) > 0.012 ) isGood = false;
		  if( fabs(dPhiAtVtxRecoEle[theEle]) > 0.1 )   isGood = false;
		  if( HoverERecoEle[theEle]>0.003 )                    isGood = false;
		  if( sigmaIetaIetaRecoEle[theEle]>0.022 )             isGood = false;
		  if( dr03EcalSumEtRecoEle[theEle]>2 ) isGood = false;
		}
		
		if (!isBarrel){
		  if( fabs(dEtaAtVtxRecoEle[theEle])>0.012 ) isGood = false;
		  if( fabs(dPhiAtVtxRecoEle[theEle])>0.1 )   isGood = false;
		  if( HoverERecoEle[theEle]>0.003 )          isGood = false;
		  if( sigmaEtaEtaRecoEle[theEle]>0.06 )      isGood = false;
		  if( dr03EcalSumEtRecoEle[theEle]>2 )       isGood = false;
		}
		
		if (!isGood) continue;
		if (chargeRecoEle[theEle]>0) numberOfEleOkPlus++;
		if (chargeRecoEle[theEle]<0) numberOfEleOkMinus++;
		
		
		// search for the two reco electrons closer to the generated electrons in case of signal
		if (signal) {
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

	
	if (numberOfEleOkPlus>0 && numberOfEleOkMinus>0) totalIdentified++;
	if ((numberOfEleOkPlus>0 && numberOfEleOkMinus>1) || (numberOfEleOkPlus>1 && numberOfEleOkMinus>0)) totalIdentifiedMore++;
	
	
	// to choose, signal only
	if (signal) {
	  
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
	
	  // tracker isolation added here
	  float deltaR_highestEt = highestEt_3p.DeltaR(highestEt_3e);	
	  if (deltaR_highestEt < 0.3) {
	    highestEt_TIso03p -= (highestEt_4e.Et()/highestEt_4p.Et());
	    if (highestEt_TIso03p < 0.) highestEt_TIso03p = 0.;
	    highestEt_TIso03e -= (highestEt_4p.Et()/highestEt_4e.Et());
	    if (highestEt_TIso03e < 0.) highestEt_TIso03e = 0.;
	  }
	
	  // BEST PAIR CUTS
	  if (highestEt_TIso03p*highestEt_4p.Et() < 2.8 && highestEt_TIso03e*highestEt_4e.Et() < 2.8) numberOfPairsOk++;
	
	  if (highestEt_TIso03p*highestEt_4p.Et() < 2.8 && highestEt_TIso03e*highestEt_4e.Et() < 2.8) {
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
	    if (mee_highestEt<4. && mee_highestEt>2.5) numbersOfInvMassOk++;
	  
	  
	    // study of possible biases for signal
	    if (signal) {    
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
	  
	  } // ok best pair
	} // ok matched
      
	// plots with all combinatorics
	if ( numberOfEleOkPlus>=1 && numberOfEleOkMinus>=1 && numberOfPairsOk>=1 ){ 
	  for(int theEle1=0; theEle1<numberOfElectrons; theEle1++) { 
	    if (pxRecoEle[theEle1]<-700) continue;
	    if (etRecoEle[theEle1]<4)     continue;
	  
	    for(int theEle2=(theEle1+1); theEle2<numberOfElectrons; theEle2++) { 
	      if (pxRecoEle[theEle2]<-700) continue;
	      if (etRecoEle[theEle2]<4)     continue;
	    
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

  //  drawPlots();


  // summary
  cout << "total number of events = "              << totalEvents     << endl;
  cout << "good number of gene MC = "              << goodGene        << endl;
  cout << "total number of events passing HLT = "  << totalTrigger    << endl;
  cout << "total number of reco = "                << totalReco       << endl;
  cout << "total number of reco & Et>4= "          << totalRecoGt4    << endl;
  cout << "total number of reco & Et>4 & charge requirement = " << totalRecoGt4Charge << endl;
  cout << "total number of reco & Et>4 & id & isol = " << totalIdentified << endl;
  cout << "total number of reco & ET>4 & id & full isol (tracker also) = " << numberOfPairsOk << endl;
  cout << "total number of reco & ET>4 & id & full isol (tracker also) + invMass in 2.5 - 4 = " << numbersOfInvMassOk << endl;
  cout << endl;
  cout << "total number of reco & matched & id, more than 2 = " << totalIdentifiedMore << endl;
  cout << endl;
  cout << "total number of events with closer to MC match MC within dR=0.5 = " << okEvent_McTruth << endl; 
  cout << "total number of events with highest Et match MC within dR=0.5 = "   << okEvent_Highest << endl; 

  
  // normalizing to cross sections
  float hltEff_prod; 
  float filterEff_prod;
  float crossSection;
  if (theSample==1) { hltEff_prod=1.;       filterEff_prod = 1.;      crossSection = 17000.;       }    // hltEff_prod = j/psi only 
  if (theSample==2) { hltEff_prod=1.;       filterEff_prod = 1.;      crossSection = 2500;         }
  if (theSample==3) { hltEff_prod=0.0345;   filterEff_prod = 0.00048; crossSection = 400000000.;   } 
  if (theSample==4) { hltEff_prod=0.1119;   filterEff_prod = 0.0024;  crossSection = 100000000.;   } 
  if (theSample==5) { hltEff_prod=0.2051;   filterEff_prod = 0.012;   crossSection = 1900000.;     } 
  if (theSample==6) { hltEff_prod=1.;       filterEff_prod = 0.0080;  crossSection = 400000000.;   }
  if (theSample==7) { hltEff_prod=0.0684;   filterEff_prod = 0.047;   crossSection = 100000000.;   }  
  if (theSample==8) { hltEff_prod=0.12728;  filterEff_prod = 0.15;    crossSection = 1900000.;     }
  if (theSample==9) { hltEff_prod=0.000079; filterEff_prod = 1.;      crossSection = 75280000000.; }

  float numEvents;
  if(theSample > 2) numEvents = totalTrigger/hltEff_prod; 
  if(theSample <=2) numEvents = goodGene;
  float kineEff_hlt           = totalTrigger/numEvents;
  float kineEff_reco          = totalReco/numEvents;
  float kineEff_recoGt4       = totalRecoGt4/numEvents;
  float kineEff_recoGt4charge = totalRecoGt4Charge/numEvents;
  float kineEff_id            = totalIdentified/numEvents;
  float kineEff_isol          = numberOfPairsOk/numEvents;
  float kineEff_invMass       = numbersOfInvMassOk/numEvents;

  float exp_hlt            = kineEff_hlt           * filterEff_prod*crossSection;
  float exp_reco           = kineEff_reco          * filterEff_prod*crossSection;
  float exp_recoGt4        = kineEff_recoGt4       * filterEff_prod*crossSection;
  float exp_recoGt4charge  = kineEff_recoGt4charge * filterEff_prod*crossSection;
  float exp_id             = kineEff_id            * filterEff_prod*crossSection;
  float exp_isol           = kineEff_isol          * filterEff_prod*crossSection;
  float exp_okMass         = kineEff_invMass       * filterEff_prod * crossSection;
  
  cout << "number of expected events in 10 pb-1: "            << endl;
  cout << "after HLT                  : " << 10.*exp_hlt      << endl; 
  cout << "after reco                 : " << 10.*exp_reco     << endl; 
  cout << "after reco + Et>4          : " << 10.*exp_recoGt4  << endl; 
  cout << "after reco + Et>4 + charge : " << 10.*exp_recoGt4charge << endl; 
  cout << "after id                   : " << 10.*exp_id       << endl; 
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
  HepHisto_ptHat   = new TH1F("HepHisto_ptHat",   "Pt hat distributions", 800, 0, 200.);
  HepHisto_deltaR  = new TH1F("HepHisto_deltaR",  "GenElectrons deltaR ", 100, 0.,  1.);
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
  EleHisto_HoE       = new  TH1F("EleHisto_HoE      ", " HoverE Ele ", 200,    0., 0.3);      
  EleHisto_EoP       = new  TH1F("EleHisto_EoP      ", " EoverP Ele ", 200,    0.,  5.);      
  EleHisto_fbrem     = new  TH1F("EleHisto_fbrem    ", " fBrem Ele  ", 200,  -30.,  2.);   
  EleHisto_sigmaIeta = new  TH1F("EleHisto_sigmaIeta", " sigmaIetaIeta", 300,    0., 1.5);
  EleHisto_sigmaEta  = new  TH1F("EleHisto_sigmaEta ", " sigmaEtaEta", 300,    0., 0.5); 
  EleHisto_dr03ecal  = new  TH1F("EleHisto_dr03ecal ", " dr03ecal   ", 200,   -5.,  60); 
  EleHisto_dr03tk    = new  TH1F("EleHisto_dr03tk   ", " dr03tk     ", 200,   -5.,  60);   
  EleHisto_dr03hcal1 = new  TH1F("EleHisto_dr03hcal1", " dr03hcal1  ", 200,   -5.,  30);
  EleHisto_dr03hcal2 = new  TH1F("EleHisto_dr03hcal2", " dr03hcal2  ", 200,   -5.,  30);
  EleHisto_momErr    = new  TH1F("EleHisto_momErr   ", " momentum error", 500,    0., 10000.);   
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
  ScHisto_InvMassVsDeltaR_highestEt  = new TH2F("ScHisto_InvMassVsDeltaR_highestEt", "ScHisto_InvMassVsDeltaR_highestEt",  50,  0.,   1.,   75, 0., 6.);
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
  ScHisto_invMassHighestEt      -> Draw();  c.Print("scInvMass.eps");
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
  ScHistoGt4plus_size   -> Write();
  ScHistoGt4minus_size  -> Write();

  ScHisto_invMassHighestEt      -> Write();

  if (signal) {

    HepHisto_size   -> Write();
    HepHisto_ptHat   -> Write();
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
  ScHisto_invMassHighestEt  -> Write();
  ScHisto_s9s25HighestEt    -> Write();
  ScHisto_sEEHighestEt      -> Write();
  ScHisto_hoeHighestEt      -> Write();
  ScHisto_dEtaTrHighestEt   -> Write();
  ScHisto_dPhiTrHighestEt   -> Write();


  ScHisto_InvMassVsJene_highestEt  -> Write();         
  ScHisto_InvMassVsJeta_highestEt  -> Write();         
  ScHisto_InvMassVsJphi_highestEt  -> Write();         

}


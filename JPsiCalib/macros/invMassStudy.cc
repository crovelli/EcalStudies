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

#include "invMassStudy.hh"

using namespace std;

invMassStudy::invMassStudy(TTree *tree) 
  : JPsiTreeBase(tree) {
}

invMassStudy::~invMassStudy(){ } 

void invMassStudy::Loop() {

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
  int totalIdentified  = 0;
  int goodGene         = 0;

  for (int jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // signal = 0;
    
    totalEvents++;

    // only for good number of generated
    if(!signal || numberOfGenerated==2) {
      goodGene++;
      
      // HLT
      if (hlt29) {
	totalTrigger++;
	
	// few numbers
	if (numberOfScGt4>1) totalReco++;
	if (numberOfScMatchedPlus>0 && numberOfScMatchedMinus>0) totalMatched++;
	
	TVector3 highestEt_3p, highestEt_pv_3p, highestEt_tracker_3p;
	TVector3 highestEt_3e, highestEt_pv_3e, highestEt_tracker_3e;
	TLorentzVector highestEt_4p, highestEt_pv_4p, highestEt_tracker_4p; 
	TLorentzVector highestEt_4e, highestEt_pv_4e, highestEt_tracker_4e; 
	
	float highestEtp  = -9999.;
	float highestEte  = -9999.;
	int theHighestEtp = -1;
	int theHighestEte = -1;
	
	int numberOfScOkPlus  = 0;
	int numberOfScOkMinus = 0;
	for(int theSc=0; theSc<numberOfScGt4; theSc++) { 
	  
	  if (numberOfScMatchedPlus>=1 && numberOfScMatchedMinus>=1) {
	    
	    // skipping superclusters without matching to tracker
	    if (xOrRecoSc[theSc]>-700) {
	      
	      TVector3 this3P; 
	      TLorentzVector this4P;
	      this3P.SetXYZ(xOrRecoSc[theSc], yOrRecoSc[theSc], zOrRecoSc[theSc]);
	      this4P.SetPxPyPzE(xOrRecoSc[theSc], yOrRecoSc[theSc], zOrRecoSc[theSc], eneRecoSc[theSc]);

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
		/* old
		if(fabs(dEtaWithTrackerRecoSc[theSc])>0.009) isGood = false;
		if(fabs(dPhiWithTrackerRecoSc[theSc])>0.09)  isGood = false;
		if(HoverERecoSc[theSc]>0.115)                isGood = false;
		if(sigmaEtaEtaRecoSc[theSc]>0.0140)          isGood = false;
		*/

		// optimized
		if(fabs(dEtaWithTrackerRecoSc[theSc])>0.0065)                 isGood = false;
		if(fabs(dPhiWithTrackerRecoSc[theSc])>0.045)                  isGood = false;
		if(HoverERecoSc[theSc]>0.032)                                 isGood = false;
		if(sigmaEtaEtaRecoSc[theSc]>0.0115)                           isGood = false;		
		if(emIsolRecoSc_03[theSc]*this4P.E()*sin(this4P.Theta())>18.) isGood = false;		
		if(trkIsolRecoSc_03[theSc]*this4P.E()*sin(this4P.Theta())>3.) isGood = false;
		if(hadIsolRecoSc_03[theSc]>0.45)                              isGood = false;
	      }
	      
	      if (!isBarrel){
		/* old
		   if(fabs(dEtaWithTrackerRecoSc[theSc])>0.0105) isGood = false;
		   if(fabs(dPhiWithTrackerRecoSc[theSc])>0.092)  isGood = false;
		   if(HoverERecoSc[theSc]>0.150)                 isGood = false;
		   if(sigmaEtaEtaRecoSc[theSc]>0.0275)           isGood = false;
		*/
		
		// optimized
		if(fabs(dEtaWithTrackerRecoSc[theSc])>0.0075)                 isGood = false;
		if(fabs(dPhiWithTrackerRecoSc[theSc])>0.045)                  isGood = false;
		if(HoverERecoSc[theSc]>0.032)                                 isGood = false;
		if(sigmaEtaEtaRecoSc[theSc]>0.029)                            isGood = false; 		
		if(emIsolRecoSc_03[theSc]*this4P.E()*sin(this4P.Theta())>18.) isGood = false;		
		if(trkIsolRecoSc_03[theSc]*this4P.E()*sin(this4P.Theta())>3.) isGood = false;
		if(hadIsolRecoSc_03[theSc]>0.45)                              isGood = false;
	      }

	      if (!isGood) continue;
	      if (chargeSc[theSc]>0) numberOfScOkPlus++;
	      if (chargeSc[theSc]<0) numberOfScOkMinus++;
	      // ---------------------------------------------
	      

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
	      }

	    }  // ok matching with track
	    
	  } // two matched electons
	  
	} // loop over superclusters
	
	if (numberOfScOkPlus>0 && numberOfScOkMinus>0) totalIdentified++;
	
	
	// to compare the distributions: use the two highest Et electrons
	if (numberOfScMatchedPlus>=1 && numberOfScMatchedMinus>=1 && numberOfScOkPlus>=1 && numberOfScOkMinus>=1) {
	  
	  float mee_highestEt         = (highestEt_4p + highestEt_4e).M();
	  float mee_highestEt_pv      = (highestEt_pv_4p + highestEt_pv_4e).M();
	  float mee_highestEt_tracker = (highestEt_tracker_4p + highestEt_tracker_4e).M();

	  ScHisto_invMassHighestEt->Fill(mee_highestEt);
	  ScHisto_pv_invMassHighestEt -> Fill(mee_highestEt_pv);	  	  
	  ScHisto_tracker_invMassHighestEt -> Fill(mee_highestEt_tracker);	  	  
	  
	  // to study EB/ EE depencency: both in barrel only
	  if ( fabs(highestEt_3p.Eta())        <1.479 && fabs(highestEt_3e.Eta())        <1.479 ) ScHisto_invMassHighestEt_EB         -> Fill(mee_highestEt);
	  if ( fabs(highestEt_pv_3p.Eta())     <1.479 && fabs(highestEt_pv_3e.Eta())     <1.479 ) ScHisto_pv_invMassHighestEt_EB      -> Fill(mee_highestEt_pv);
	  if ( fabs(highestEt_tracker_3p.Eta())<1.479 && fabs(highestEt_tracker_3e.Eta())<1.479 ) ScHisto_tracker_invMassHighestEt_EB -> Fill(mee_highestEt_tracker);
	  
	  // both in endcap only
	  if ( fabs(highestEt_3p.Eta())         >1.479 && fabs(highestEt_3e.Eta())        >1.479 ) ScHisto_invMassHighestEt_EE         -> Fill(mee_highestEt);                
	  if ( fabs(highestEt_pv_3p.Eta())      >1.479 && fabs(highestEt_pv_3e.Eta())     >1.479 ) ScHisto_pv_invMassHighestEt_EE      -> Fill(mee_highestEt_pv);
	  if ( fabs(highestEt_tracker_3p.Eta()) >1.479 && fabs(highestEt_tracker_3e.Eta())>1.479 ) ScHisto_tracker_invMassHighestEt_EE -> Fill(mee_highestEt_tracker);	    
	  
	  // one and one
	  if (  (fabs(highestEt_3p.Eta())<1.479 && fabs(highestEt_3e.Eta())>1.479) || (fabs(highestEt_3e.Eta())<1.479 && fabs(highestEt_3p.Eta())>1.479) ) 
	    ScHisto_invMassHighestEt_EBEE -> Fill(mee_highestEt);
	  if ( (fabs(highestEt_pv_3p.Eta())<1.479 && fabs(highestEt_pv_3e.Eta())>1.479) || (fabs(highestEt_pv_3e.Eta())<1.479 && fabs(highestEt_pv_3p.Eta())>1.479)) 
	    ScHisto_pv_invMassHighestEt_EBEE -> Fill(mee_highestEt_pv);
	  if ( (fabs(highestEt_tracker_3p.Eta())<1.479 && fabs(highestEt_tracker_3e.Eta())>1.479) || (fabs(highestEt_tracker_3e.Eta())<1.479 && fabs(highestEt_tracker_3p.Eta())>1.479)) 
	    ScHisto_tracker_invMassHighestEt_EBEE -> Fill(mee_highestEt_tracker);


	  // study of possible biases
	  if (signal) {    
	    float jpsi_ene    = (highestEt_4p + highestEt_4e).E();	  
	    float jpsi_et     = (highestEt_4p + highestEt_4e).Et();	  
	    float jpsi_eta    = (highestEt_4p + highestEt_4e).Eta();	  
	    float jpsi_phi    = (highestEt_4p + highestEt_4e).Phi();	  
	    float deltaR      = highestEt_4p.DeltaR(highestEt_4e);

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
	  
	} // ok matched
      } // ok HLT
      
    } // ok generated
    

    
    
  } // loop over entries


  saveHistos();

  drawPlots();


  // summary
  cout << "total number of events = "              << totalEvents     << endl;
  cout << "good number of gene MC = "              << goodGene        << endl;
  cout << "total number of events passing HLT = "  << totalTrigger    << endl;
  cout << "total number of reco = "                << totalReco       << endl;
  cout << "total number of reco & matched = "      << totalMatched    << endl;
  cout << "total number of reco & matched & id = " << totalIdentified << endl;
  cout << endl;

} // end of program


void invMassStudy::bookHistos() {

  ScHisto_invMassHighestEt                = new TH1F("ScHisto_invMassHighestEt",           "Sc invariant mass",  75, 0.,6.);   
  ScHisto_invMassHighestEt_EB             = new TH1F("ScHisto_invMassHighestEt_EB",        "Sc invariant mass",  75, 0.,6.);   
  ScHisto_invMassHighestEt_EE             = new TH1F("ScHisto_invMassHighestEt_EE",        "Sc invariant mass",  75, 0.,6.);
  ScHisto_invMassHighestEt_EBEE           = new TH1F("ScHisto_invMassHighestEt_EBEE",      "Sc invariant mass",  75, 0.,6.);

  ScHisto_pv_invMassHighestEt             = new TH1F("ScHisto_pv_invMassHighestEt",           "Sc invariant mass",  75, 0.,6.);   
  ScHisto_pv_invMassHighestEt_EB          = new TH1F("ScHisto_pv_invMassHighestEt_EB",        "Sc invariant mass",  75, 0.,6.);   
  ScHisto_pv_invMassHighestEt_EE          = new TH1F("ScHisto_pv_invMassHighestEt_EE",        "Sc invariant mass",  75, 0.,6.);
  ScHisto_pv_invMassHighestEt_EBEE        = new TH1F("ScHisto_pv_invMassHighestEt_EBEE",      "Sc invariant mass",  75, 0.,6.);

  ScHisto_tracker_invMassHighestEt        = new TH1F("ScHisto_tracker_invMassHighestEt",      "Sc invariant mass",  75, 0.,6.);   
  ScHisto_tracker_invMassHighestEt_EB     = new TH1F("ScHisto_tracker_invMassHighestEt_EB",   "Sc invariant mass",  75, 0.,6.);   
  ScHisto_tracker_invMassHighestEt_EE     = new TH1F("ScHisto_tracker_invMassHighestEt_EE",   "Sc invariant mass",  75, 0.,6.);
  ScHisto_tracker_invMassHighestEt_EBEE   = new TH1F("ScHisto_tracker_invMassHighestEt_EBEE", "Sc invariant mass",  75, 0.,6.);

  ScHisto_JeneVsJeta_highestEt            = new TH2F("ScHisto_JeneVsJeta_highestEt",            "ScHisto_JeneVsJeta_highestEt",            100, 0., 2.5, 50, 5., 100.);
  ScHisto_JetVsJeta_highestEt             = new TH2F("ScHisto_JetVsJeta_highestEt",             "ScHisto_JetVsJeta_highestEt",             100, 0., 2.5, 50, 5., 30.);
  ScHisto_InvMassVsJene_highestEt         = new TH2F("ScHisto_InvMassVsJene_highestEt",         "ScHisto_InvMassVsJene_highestEt",          50, 5., 100., 75, 0., 6.);
  ScHisto_InvMassVsJet_highestEt          = new TH2F("ScHisto_InvMassVsJet_highestEt",          "ScHisto_InvMassVsJet_highestEt",           50, 5.,  30., 75, 0., 6.);
  ScHisto_InvMassVsJeneEB_highestEt       = new TH2F("ScHisto_InvMassVsJeneEB_highestEt",       "ScHisto_InvMassVsJeneEB_highestEt",        50, 5., 100., 75, 0., 6.);
  ScHisto_InvMassVsJeneEE_highestEt       = new TH2F("ScHisto_InvMassVsJeneEE_highestEt",       "ScHisto_InvMassVsJeneEE_highestEt",        50, 5., 100., 75, 0., 6.);
  ScHisto_InvMassVsJetEB_highestEt        = new TH2F("ScHisto_InvMassVsJetEB_highestEt",        "ScHisto_InvMassVsJetEB_highestEt",         50, 5.,  30., 75, 0., 6.);
  ScHisto_InvMassVsJetEE_highestEt        = new TH2F("ScHisto_InvMassVsJetEE_highestEt",        "ScHisto_InvMassVsJetEE_highestEt",         50, 5.,  30., 75, 0., 6.);
  ScHisto_InvMassVsJeta_highestEt         = new TH2F("ScHisto_InvMassVsJeta_highestEt",         "ScHisto_InvMassVsJeta_highestEt",          50, -2.5, 2.5,  75, 0., 6.);
  ScHisto_InvMassVsJphi_highestEt         = new TH2F("ScHisto_InvMassVsJphi_highestEt",         "ScHisto_InvMassVsJphi_highestEt",          50, -3.14,3.14, 75, 0., 6.);
  ScHisto_InvMassVsDeltaR_highestEt       = new TH2F("ScHisto_InvMassVsDeltaR_highestEt",       "ScHisto_InvMassVsDeltaR_highestEt",        50, 0., 1., 75, 0., 6.);
}

void invMassStudy::drawPlots() {

  TCanvas c; c.cd();

  if (signal) {
    c.SetLogy(0);

    // 
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
		

  c.SetLogy(0);
  ScHisto_invMassHighestEt      -> Draw();  c.Print("scInvMass.eps");
  ScHisto_invMassHighestEt_EB   -> Draw();  c.Print("scInvMass_EB.eps");
  ScHisto_invMassHighestEt_EE   -> Draw();  c.Print("scInvMass_EE.eps");
  ScHisto_invMassHighestEt_EBEE -> Draw();  c.Print("scInvMass_EBEE.eps");

  c.SetLogy(0);
  ScHisto_pv_invMassHighestEt           -> Draw();  c.Print("scInvMass_pv.eps");          
  ScHisto_pv_invMassHighestEt_EB        -> Draw();  c.Print("scInvMass_EB_pv.eps");          
  ScHisto_pv_invMassHighestEt_EE        -> Draw();  c.Print("scInvMass_EE_pv.eps");          
  ScHisto_pv_invMassHighestEt_EBEE      -> Draw();  c.Print("scInvMass_EBEE_pv.eps");           

  c.SetLogy(0);
  ScHisto_tracker_invMassHighestEt      -> Draw();  c.Print("scInvMass_tracker.eps");              
  ScHisto_tracker_invMassHighestEt_EB   -> Draw();  c.Print("scInvMass_EB_tracker.eps");              
  ScHisto_tracker_invMassHighestEt_EE   -> Draw();  c.Print("scInvMass_EE_tracker.eps");                 
  ScHisto_tracker_invMassHighestEt_EBEE -> Draw();  c.Print("scInvMass_EBEE_tracker.eps");               
}

void invMassStudy::saveHistos() {

  TFile fOut("Outfile_histo.root", "RECREATE");
  fOut.cd();

  ScHisto_invMassHighestEt      -> Write();
  ScHisto_invMassHighestEt_EB   -> Write();
  ScHisto_invMassHighestEt_EE   -> Write();
  ScHisto_invMassHighestEt_EBEE -> Write();

  ScHisto_pv_invMassHighestEt      -> Write();
  ScHisto_pv_invMassHighestEt_EB   -> Write();
  ScHisto_pv_invMassHighestEt_EE   -> Write();
  ScHisto_pv_invMassHighestEt_EBEE -> Write();

  ScHisto_tracker_invMassHighestEt      -> Write();
  ScHisto_tracker_invMassHighestEt_EB   -> Write();
  ScHisto_tracker_invMassHighestEt_EE   -> Write();
  ScHisto_tracker_invMassHighestEt_EBEE -> Write();
}


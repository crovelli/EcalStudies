#ifndef finalJPsiAnalysis_h
#define finalJPsiAnalysis_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


class finalJPsiAnalysis : public JPsiTreeBase{
public:
  
  finalJPsiAnalysis(TTree *tree=0);
  virtual ~finalJPsiAnalysis();
  void Loop();

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  
  // histos
  TH1F *ScHistoGt4_size;
  TH1F *ScHistoMatchedPlusSide_size; 
  TH1F *ScHistoMatchedMinusSide_size;

  TH1F *HepHisto_size;
  TH1F *HepHisto_deltaR;
  TH1F *HepHisto_maxPt;
  TH1F *HepHisto_minPt;
  TH1F *HepHisto_eta;
  TH1F *HepHisto_phi;
  TH1F *HepHisto_invMass;
  TH2F *HepHisto_PtVsEta;

  TH1F *ScHisto_deltaRWrtMc_CloserToMc;  
  TH1F *ScHisto_deltaRWrtMc_HighestEt; 
  TH1F *ScHisto_deltaRWrtMc_BestS9S25; 
  TH1F *ScHisto_deltaRWrtMc_BestSEE;
  TH1F *ScHisto_deltaRWrtMc_CloserToMc_more1;
  TH1F *ScHisto_deltaRWrtMc_HighestEt_more1;
  TH1F *ScHisto_deltaRWrtMc_BestS9S25_more1;
  TH1F *ScHisto_deltaRWrtMc_BestSEE_more1;
  TH1F *ScHisto_deltaRWrtMc_CloserToMc_more1_zoom;
  TH1F *ScHisto_deltaRWrtMc_HighestEt_more1_zoom;
  TH1F *ScHisto_deltaRWrtMc_BestS9S25_more1_zoom;
  TH1F *ScHisto_deltaRWrtMc_BestSEE_more1_zoom;

  TH1F *ScHisto_deltaRComb;       
  TH1F *ScHisto_invMassComb;
  TH2F *ScHisto_deltaRVsInvMassComb;

  TH1F *ScHisto_etaHighestEt;          
  TH1F *ScHisto_phiHighestEt;          
  TH1F *ScHisto_maxEtHighestEt;        
  TH1F *ScHisto_minEtHighestEt;        
  TH1F *ScHisto_deltaRHighestEt;       
  TH1F *ScHisto_invMassHighestEt;      
  TH1F *ScHisto_s9s25HighestEt;        
  TH1F *ScHisto_sEEHighestEt;          
  TH1F *ScHisto_hcalIsolHighestEt;     
  TH1F *ScHisto_hcalIsolHighestEtZoom; 
  TH1F *ScHisto_hoeHighestEt;          
  TH1F *ScHisto_dEtaTrHighestEt;
  TH1F *ScHisto_dPhiTrHighestEt;

  TH1F *ScHisto_deltaRHighestEt_EB;       
  TH1F *ScHisto_invMassHighestEt_EB;      
  TH1F *ScHisto_s9s25HighestEt_EB;        
  TH1F *ScHisto_sEEHighestEt_EB;          
  TH1F *ScHisto_hcalIsolHighestEt_EB;     
  TH1F *ScHisto_hoeHighestEt_EB;          
  TH1F *ScHisto_dEtaTrHighestEt_EB;
  TH1F *ScHisto_dPhiTrHighestEt_EB;

  TH1F *ScHisto_pv_invMassHighestEt;
  TH1F *ScHisto_pv_invMassHighestEt_EB;
  TH1F *ScHisto_pv_invMassHighestEt_EE;
  TH1F *ScHisto_pv_invMassHighestEt_EBEE;
  TH1F *ScHisto_tracker_invMassHighestEt;
  TH1F *ScHisto_tracker_invMassHighestEt_EB;
  TH1F *ScHisto_tracker_invMassHighestEt_EE;
  TH1F *ScHisto_tracker_invMassHighestEt_EBEE;

  TH1F *ScHisto_deltaRHighestEt_EE;       
  TH1F *ScHisto_invMassHighestEt_EE;      
  TH1F *ScHisto_s9s25HighestEt_EE;        
  TH1F *ScHisto_sEEHighestEt_EE;          
  TH1F *ScHisto_hcalIsolHighestEt_EE;     
  TH1F *ScHisto_hoeHighestEt_EE;          
  TH1F *ScHisto_dEtaTrHighestEt_EE;
  TH1F *ScHisto_dPhiTrHighestEt_EE;

  TH1F *ScHisto_deltaRHighestEt_EBEE;       
  TH1F *ScHisto_invMassHighestEt_EBEE;      
  TH1F *ScHisto_s9s25HighestEt_EBEE;        
  TH1F *ScHisto_sEEHighestEt_EBEE;          
  TH1F *ScHisto_hcalIsolHighestEt_EBEE;     
  TH1F *ScHisto_hoeHighestEt_EBEE;
  TH1F *ScHisto_dEtaTrHighestEt_EBEE;
  TH1F *ScHisto_dPhiTrHighestEt_EBEE;

  TH2F *ScHisto_deltaRVsDeltaPhi_HighestEt;
  TH2F *ScHisto_deltaRVsDeltaEta_HighestEt;
  TH2F *ScHisto_deltaRVsHoE_HighestEt;     
  TH2F *ScHisto_deltaRVsSee_HighestEt;    
  TH2F *ScHisto_deltaRVsInvMass_HighestEt; 

  TH2F *ScHisto_InvMassVsJene_highestEt;
  TH2F *ScHisto_InvMassVsJeta_highestEt;
  TH2F *ScHisto_InvMassVsJphi_highestEt;
  TH2F *ScHisto_thRecOverThTrueVsJene_highestEt;
  TH2F *ScHisto_EoEtrueVsErec_highestEt1;
  TH2F *ScHisto_EoEtrueVsErec_highestEt2;
};
#endif


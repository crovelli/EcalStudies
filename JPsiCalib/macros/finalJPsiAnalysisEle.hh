#ifndef finalJPsiAnalysisEle_h
#define finalJPsiAnalysisEle_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


class finalJPsiAnalysisEle : public JPsiTreeBase{
public:
  
  finalJPsiAnalysisEle(TTree *tree=0);
  virtual ~finalJPsiAnalysisEle();
  void Loop(int theSample);

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  void fillHistoCutsVariables( int Eleindex);

  // histos
  TH1F *ScHistoEle_size;
  TH1F *ScHistoGt4_size;
  TH1F *ScHistoGt4plus_size;
  TH1F *ScHistoGt4minus_size;

  TH1F *HepHisto_size;
  TH1F *HepHisto_ptHat;
  TH1F *HepHisto_deltaR;
  TH1F *HepHisto_maxPt;
  TH1F *HepHisto_minPt;
  TH1F *HepHisto_eta;
  TH1F *HepHisto_phi;
  TH1F *HepHisto_invMass;
  TH2F *HepHisto_PtVsEta;

  TH1F *EleHisto_eta;      
  TH1F *EleHisto_phi;      
  TH1F *EleHisto_detaVtx;  
  TH1F *EleHisto_dphiVtx;  
  TH1F *EleHisto_HoE;      
  TH1F *EleHisto_EoP;      
  TH1F *EleHisto_fbrem;   
  TH1F *EleHisto_sigmaIeta;
  TH1F *EleHisto_sigmaEta; 
  TH1F *EleHisto_dr03ecal; 
  TH1F *EleHisto_dr03tk;   
  TH1F *EleHisto_dr03tkcorr;   
  TH1F *EleHisto_dr03hcal1;
  TH1F *EleHisto_dr03hcal2;
  TH1F *EleHisto_momErr;   
  TH1F *EleHisto_isPF;     
  TH1F *EleHisto_isEcal;   
  TH1F *EleHisto_mva;      


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

  TH1F *ScHisto_deltaRHighestEt_EE;  
  TH1F *ScHisto_invMassHighestEt_EE;
  
  TH1F *ScHisto_deltaRHighestEt_EBEE;
  TH1F *ScHisto_invMassHighestEt_EBEE;
  

  TH2F *ScHisto_InvMassVsJene_highestEt;         
  TH2F *ScHisto_InvMassVsJet_highestEt;          
  TH2F *ScHisto_InvMassVsJeneEB_highestEt;         
  TH2F *ScHisto_InvMassVsJeneEE_highestEt;         
  TH2F *ScHisto_InvMassVsJetEB_highestEt;          
  TH2F *ScHisto_InvMassVsJetEE_highestEt;          
  TH2F *ScHisto_InvMassVsJeta_highestEt;         
  TH2F *ScHisto_InvMassVsJphi_highestEt;         
  TH2F *ScHisto_JeneVsJeta_highestEt;
  TH2F *ScHisto_JetVsJeta_highestEt;
  TH2F *ScHisto_InvMassVsDeltaR_highestEt;
};
#endif


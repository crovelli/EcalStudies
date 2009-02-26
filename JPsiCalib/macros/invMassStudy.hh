#ifndef invMassStudy_h
#define invMassStudy_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


class invMassStudy : public JPsiTreeBase{
public:
  
  invMassStudy(TTree *tree=0);
  virtual ~invMassStudy();
  void Loop();

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  
  // histos
  TH1F *ScHisto_invMassHighestEt;      
  TH1F *ScHisto_invMassHighestEt_EB;      
  TH1F *ScHisto_invMassHighestEt_EE;      
  TH1F *ScHisto_invMassHighestEt_EBEE;      

  TH1F *ScHisto_pv_invMassHighestEt;
  TH1F *ScHisto_pv_invMassHighestEt_EB;
  TH1F *ScHisto_pv_invMassHighestEt_EE;
  TH1F *ScHisto_pv_invMassHighestEt_EBEE;

  TH1F *ScHisto_tracker_invMassHighestEt;
  TH1F *ScHisto_tracker_invMassHighestEt_EB;
  TH1F *ScHisto_tracker_invMassHighestEt_EE;
  TH1F *ScHisto_tracker_invMassHighestEt_EBEE;

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


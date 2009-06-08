#ifndef finalLowPtAnalysis_h
#define finalLowPtAnalysis_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "LowPtTreeBase.h"


class finalLowPtAnalysis : public LowPtTreeBase{
public:
  
  finalLowPtAnalysis(TTree *tree=0);
  virtual ~finalLowPtAnalysis();
  void Loop();

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  
  // histos
  TH1F *ScHistoEle_size;
  TH1F *ScHistoGt4_size;
  TH1F *ScHistoGt4plus_size;
  TH1F *ScHistoGt4minus_size;

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

  TH1F *ScHisto_eta;          
  std::vector<TH1F*> ScHisto_etaFilt;

  TH1F *ScHisto_et;          
  std::vector<TH1F*> ScHisto_etFilt;

  TH1F *ScHisto_charge;          
  std::vector<TH1F*> ScHisto_chargeFilt;


  TH2F *Histo_PassVsEffl; 
  TH2F *Histo_PassVsEffh; 
  TH2F *Histo_PassVs2eEffl; 
  TH2F *Histo_PassVs2eEffh; 

};
#endif


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

  TH1F *HepHisto_size;
  TH1F *HepHisto_ptHat;
  std::vector<TH1F*> HepHisto_ptHatFilt;

  TH1F *HepHisto_MCEledR;
  TH1F *HepHisto_PiD;
  TH1F *HepHisto_eta;

  TH1F *HepHisto_PiDStable;
  std::vector<TH1F*> HepHisto_PiDStableF;

  TH1F *HepHisto_ptHatReco;
  std::vector<TH1F*> HepHisto_ptHatRecoFilt;

  TH1F *ScHisto_deltaR;
  std::vector<TH1F*> ScHisto_deltaRFilt;

  TH1F *ScHisto_invMass; 

  TH1F *ScHisto_eta;          
  std::vector<TH1F*> ScHisto_etaFilt;

  TH1F *ScHisto_phi;          
  std::vector<TH1F*> ScHisto_phiFilt;

  TH1F *ScHisto_et;          
  std::vector<TH1F*> ScHisto_etFilt;

  TH1F *ScHisto_charge;          
  std::vector<TH1F*> ScHisto_chargeFilt;

  TH1F *ScHisto_detaTk;          
  std::vector<TH1F*> ScHisto_detaTkFilt;

  TH1F *ScHisto_dphiTk;          
  std::vector<TH1F*> ScHisto_dphiTkFilt;

  TH1F *ScHisto_HoE;          
  std::vector<TH1F*> ScHisto_HoEFilt;

  TH1F *ScHisto_EoP;          
  std::vector<TH1F*> ScHisto_EoPFilt;


  TH2F *Histo_PassVsEffl; 
  TH2F *Histo_PassVsEffh; 
  TH2F *Histo_PassVs2eEffl; 
  TH2F *Histo_PassVs2eEffh; 

};
#endif


#ifndef CALIBRATION_ECALCALIBALGOS_JPSICALIBRATION
#define CALIBRATION_ECALCALIBALGOS_JPSICALIBRATION

// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Calibration/Tools/interface/ZIterativeAlgorithmWithFit.h"
#include "Calibration/Tools/interface/CalibElectron.h"
#include "Calibration/EcalCalibAlgos/interface/JPsiPlots.h"
#include "Calibration/EcalCalibAlgos/interface/ZeeRescaleFactorPlots.h"
#include "Calibration/Tools/interface/calibXMLwriter.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include<vector>
#include<string>

// class declaration
//

class JPsiCalibration : public edm::ESProducerLooper {
   
 public:
  
  /// Constructor
  JPsiCalibration( const edm::ParameterSet& iConfig );
  
  /// Destructor
  ~JPsiCalibration();
  
  /// Dummy implementation (job done in duringLoop)
  virtual void produce(edm::Event&, const edm::EventSetup&) {};
  
  /// Called at beginning of job
  virtual void beginOfJob(const edm::EventSetup&);
  
  /// Called at end of job
  virtual void endOfJob();
  
  /// Called at beginning of loop
  virtual void startingNewLoop( unsigned int iLoop );
  
  /// Called at end of loop
  virtual Status endOfLoop( const edm::EventSetup&, unsigned int iLoop );

  /// Called at each event
  virtual Status duringLoop( const edm::Event&, const edm::EventSetup& );
  
  /// Produce Ecal interCalibrations
  virtual boost::shared_ptr<EcalIntercalibConstants> produceEcalIntercalibConstants( const EcalIntercalibConstantsRcd& iRecord );
  
 private:
  
  double fEtaBarrelBad(double scEta) const;
  double fEtaBarrelGood(double scEta) const;
  double fEtaEndcapBad(double scEta) const;
  double fEtaEndcapGood(double scEta) const;

  int ringNumberCorrector(int k);
  double getEtaCorrection(const reco::PixelMatchGsfElectron*);
  
  double getE5x5(std::vector<DetId> mySCRecHits, const EBRecHitCollection* ebhits, const EERecHitCollection* eehits);
  void fillEleInfo(std::vector<HepMC::GenParticle*>& a, std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>& b);
  void fillMCInfo(HepMC::GenParticle* mcele);

  void fillMCmap(const std::vector<const reco::PixelMatchGsfElectron*>* electronCollection, const std::vector<HepMC::GenParticle*>& mcEle,std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>& myMCmap);
  
  float EvalDPhi(float Phi,float Phi_ref);
  float EvalDR(float Eta,float Eta_ref,float Phi,float Phi_ref);


  void bookHistograms();

  void resetVariables();

  void resetHistograms();

  void printStatistics();

  std::pair<DetId, double> getHottestDetId(std::vector<DetId> mySCRecHits, const EBRecHitCollection* ebhits , const EERecHitCollection* eehits);

  bool xtalIsOnModuleBorder( EBDetId myEBDetId );

  float computeCoefficientDistanceAtIteration( float v1[250], float v2[250], int size);



  // ----------member data ---------------------------

  TTree* myTree;

  edm::ParameterSet theParameterSet;
  
  std::string outputFileName_;
  std::string ZCalib_InvMass_;
  std::string mcProducer_;
  std::string calibMode_;

  std::string rechitProducerEB_;
  std::string rechitCollectionEB_;
  std::string rechitProducerEE_;
  std::string rechitCollectionEE_;
  std::string scProducerEB_;
  std::string scCollectionEB_;
  std::string scProducerEE_;
  std::string scCollectionEE_;
  std::string electronProducer_;
  std::string electronCollection_;
  edm::InputTag tracksCollection_;
  edm::InputTag calotowersCollection_;

  std::string barrelfile_, endcapfile_;
  TFile* outputFile_;

  int passedDeta, passedDphi;
  int passedHoE,  passedSee;

  double resonanceMass_, resonancePdgId_;

  unsigned int etaBins_, etBins_;  
  double etaMin_, etaMax_, etMin_, etMax_;
  double P0_,  P1_,  P2_;
  double P0e_, P1e_, P2e_;
  double mass;
  float massCorr4tree, massCorrDiff4tree;
  float massCorr4treeBB, massCorr4treeEB, massCorr4treeEE;
  
  int read_events;
  int loopFlag_;
  
  float calibCoeff[nMaxChannels];
  float NewCalibCoeff[nMaxChannels];
  float calibCoeffError[nMaxChannels];
  float initCalibCoeff[nMaxChannels];

  boost::shared_ptr<EcalIntercalibConstants> ical;
  
  ZIterativeAlgorithmWithFit* theAlgorithm_;
  
  JPsiPlots* myJPsiPlots_;
  ZeeRescaleFactorPlots* myJPsiRescaleFactorPlots_;


  TH1F *h1_eventsBeforeHLTSelection_,     *h1_eventsAfterHLTSelection_;
  TH1F *h1_eventsBeforeMCSelection_,      *h1_eventsAfterMCSelection_;
  TH1F *h1_eventsBefore2GsfSelection_,    *h1_eventsAfter2GsfSelection_;
  TH1F *h1_eventsBeforeEtSelection_,      *h1_eventsAfterEtSelection_;
  TH1F *h1_eventsBeforeEleIDSelection_,   *h1_eventsAfterEleIDSelection_;
  TH1F *h1_eventsBeforeEleIsolSelection_, *h1_eventsAfterEleIsolSelection_;
  TH1F *h1_eventsBeforeTrackerIsolSelection_, *h1_eventsAfterTrackerIsolSelection_;
  TH1F *h1_eventsBeforeBorderSelection_,  *h1_eventsAfterBorderSelection_;
  TH1F *h1_eventsBeforeInvMassSelection,  *h1_eventsAfterInvMassSelection;
  TH1F *h1_ZCandMult_;
  TH2F *h2_occupancy_vsEtaPhi;
  TH1F *h1_seedOverSC_Barrel_,   *h1_seedOverSC_Endcap_,    *h1_preshowerOverSC_;
  TH1F *h1_electronCosTheta_TK_, *h1_electronCosTheta_SC_,  *h1_electronCosTheta_SC_TK_;
  TH2F *h_ESCEtrue_vs_ESC;
  TH2F *h_ESCEtrue_vs_ESC_barrel, *h_ESCEtrue_vs_ESC_endcap;

  TH1F *h1_Theta12Resolution_TK_BB_;
  TH1F *h1_Theta12Resolution_SC_BB_;
  TH1F *h1_Theta12Resolution_TK_EE_;
  TH1F *h1_Theta12Resolution_SC_EE_;
  TH1F *h1_Theta12_TK_BB_;
  TH1F *h1_Theta12_TK_EE_;
  TH1F *h1_Theta12_SC_BB_;
  TH1F *h1_Theta12_SC_EE_;

  TH1F *h1_sEEeb_before_, *h1_sEEee_before_;
  TH1F *h1_sEEeb_after_,  *h1_sEEee_after_;

  TH1F *h1_zMassResol_, *h1_eleEtaResol_, *h1_elePhiResol_;
  TH1F *h1_zEtaResol_,  *h1_zPhiResol_;

  TH1F *h1_massResolutionTerm_total_BB_;
  TH1F *h1_massResolutionTerm_angle_BB_;
  TH1F *h1_massResolutionTerm_energy_BB_;
  TH1F *h1_massResolutionTerm_angle_ratio_BB_;
  TH1F *h1_massResolutionTerm_total_EE_;  
  TH1F *h1_massResolutionTerm_angle_EE_;
  TH1F *h1_massResolutionTerm_energy_EE_;
  TH1F *h1_massResolutionTerm_angle_ratio_EE_;
  
  TH1F *h1_reco_ZMass_,     *h1_reco_ZMassBB_,     *h1_reco_ZMassEB_,     *h1_reco_ZMassEE_;     
  TH1F *h1_reco_ZMassCorr_, *h1_reco_ZMassCorrBB_, *h1_reco_ZMassCorrEB_, *h1_reco_ZMassCorrEE_; 

  TH2F *h2_xtalMiscalibCoeffBarrel_;
  TH2F *h2_xtalMiscalibCoeffEndcapMinus_;
  TH2F *h2_xtalMiscalibCoeffEndcapPlus_;

  TH2F *h_ESCEtrueVsEta_[25];
  TH1F *h_ESCEtrue_[25];
  TH1F *h_PTres_[25];
  TH1F *h_PTres_Endcap_[25];
  TH2F *h2_chi2_[25];
  TH2F *h2_iterations_[25];
  TH2F *h_ESCcorrEtrueVsEta_[25];
  TH1F *h_ESCcorrEtrue_[25];
  TH2F *h2_xtalRecalibCoeffBarrel_[25]; 
  TH2F *h2_xtalRecalibCoeffEndcapMinus_[25];
  TH2F *h2_xtalRecalibCoeffEndcapPlus_[25];
  TH1F *h_ESCEtrue_array_barrel[25];
  TH1F *h_ESCEtrue_array_endcap[25];
  TH1F *h_eleEffEta_[2];
  TH1F *h_eleEffPhi_[2];
  TH1F *h_eleEffPt_[2]; 
  TH1F *h1_efficiencySummary_;

  TH1F *h1_occupancyVsEta_;
  TH1F *h1_weightSumMeanBarrel_;
  TH1F *h1_weightSumMeanEndcap_;  
  TH1F *h1_occupancy_, *h1_occupancyBarrel_, *h1_occupancyEndcap_;

  TH2F *h2_coeffVsEta_, *h2_coeffVsEtaGrouped_;
  TH2F *h2_zMassVsLoop_;
  TH2F *h2_zMassDiffVsLoop_;
  TH2F *h2_zWidthVsLoop_;
  TH2F *h2_coeffVsLoop_;
  TH2F *h2_residualSigma_;
  TH2F *h2_miscalRecal_;
  TH2F *h2_miscalRecalEB_;
  TH2F *h2_miscalRecalEE_;
  
  TH1F *h1_mc_, *h1_mcEB_, *h1_mcEE_;
  TH1F *h1_mcParz_[25], *h1_mcEBParz_[25], *h1_mcEEParz_[25];


  Int_t BBZN,EBZN,EEZN;
  Int_t NEVT, MCZBB, MCZEB, MCZEE;      
  
  unsigned int theMaxLoops;
  bool wantEtaCorrection_;
  
  bool useHoESelection_;
  bool useSigmaEtaEtaSelection_;
  bool useDeltaPhiSelection_;
  bool useDeltaEtaSelection_;
  bool useEmIsolSelection_;
  bool useHadIsolSelection_;
  bool useTrackerIsolSelection_;
  bool useInvMassSelection_;
  double detaEBCut_;
  double detaEECut_;               
  double dphiEBCut_;               
  double dphiEECut_;               
  double sigmaEtaEtaEBCut_;        
  double sigmaEtaEtaEECut_;        
  double HoeEBCut_;                
  double HoeEECut_;                
  double minInvMassCut_;
  double maxInvMassCut_;
  double etCut_;
  double emIsolCut_, hcalIsolCut_, trackerIsolCut_;

  unsigned int electronSelection_; 

  double loopArray[50];
  double sigmaArray[50];
  double sigmaErrorArray[50];
  double coefficientDistanceAtIteration[50];

  int BARREL_ELECTRONS_BEFORE_BORDER_CUT;
  int BARREL_ELECTRONS_AFTER_BORDER_CUT;
  int TOTAL_ELECTRONS_IN_BARREL;
  int TOTAL_ELECTRONS_IN_ENDCAP;
  int CRACK_ELECTRONS_IN_BARREL;
  int CRACK_ELECTRONS_IN_ENDCAP;

  double eventWeightForHistos;
  
  std::string weightName_;

  edm::InputTag hlTriggerResults_;
  edm::TriggerNames triggerNames_;  
  int intHlt29;

  unsigned int  nEvents_;           
  bool init_;                          
};
#endif

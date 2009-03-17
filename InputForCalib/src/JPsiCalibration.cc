#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/CachedProducts.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Framework/interface/TriggerNames.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include <FWCore/Framework/interface/Selector.h>

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "Calibration/Tools/interface/calibXMLwriter.h"

#include "Calibration/Tools/interface/CalibrationCluster.h"
#include "Calibration/Tools/interface/HouseholderDecomposition.h"
#include "Calibration/Tools/interface/MinL3Algorithm.h"
#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"
#include "Calibration/Tools/interface/EcalIndexingTools.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "Calibration/EcalAlCaRecoProducers/interface/AlCaPhiSymRecHitsProducer.h"
#include "Calibration/EcalCalibAlgos/interface/JPsiCalibration.h"
#include "Calibration/EcalCalibAlgos/interface/ZeeKinematicTools.h"

#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalBarrel.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalEndcap.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibMapEcal.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include "HLTrigger/HLTanalyzers/interface/HLTrigReport.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "EcalStudies/JPsiCalib/interface/TrackerIsolation.h"		
#include "EcalStudies/JPsiCalib/interface/CalotowerIsolation.h"

#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <utility>
#include <map>
#include <fstream>

using namespace edm;
using namespace std;
using namespace reco;

JPsiCalibration::JPsiCalibration(const edm::ParameterSet& iConfig) {
  
  std::cout<<"[JPsiCalibration] Starting the ctor"<<std::endl;
  
  // --------------------    reading from config file -----------------------------------------------
  theParameterSet=iConfig;

  // input collections
  rechitProducerEB_     = iConfig.getParameter<std::string>("rechitProducerEB");
  rechitCollectionEB_   = iConfig.getParameter<std::string>("rechitCollectionEB");
  rechitProducerEE_     = iConfig.getParameter<std::string>("rechitProducerEE");
  rechitCollectionEE_   = iConfig.getParameter<std::string>("rechitCollectionEE");
  scProducerEB_         = iConfig.getParameter<std::string>("scProducerEB");
  scCollectionEB_       = iConfig.getParameter<std::string>("scCollectionEB");
  scProducerEE_         = iConfig.getParameter<std::string>("scProducerEE");
  scCollectionEE_       = iConfig.getParameter<std::string>("scCollectionEE");
  electronProducer_     = iConfig.getParameter<std::string>("electronProducer");
  electronCollection_   = iConfig.getParameter<std::string>("electronCollection");
  tracksCollection_     = iConfig.getParameter<InputTag>("tracksCollection");
  calotowersCollection_ = iConfig.getParameter<InputTag>("calotowersCollection");
  hlTriggerResults_     = iConfig.getParameter<edm::InputTag>("HLTriggerResults");
  mcProducer_           = iConfig.getUntrackedParameter<std::string>("mcProducer","");

  // algo parameters
  resonanceMass_  = iConfig.getUntrackedParameter<double>("resonanceMass", 3.096);
  resonancePdgId_ = iConfig.getUntrackedParameter<double>("resonancePdgId", 443);  
  theMaxLoops     = iConfig.getUntrackedParameter<unsigned int>("maxLoops",0);
  calibMode_      = iConfig.getUntrackedParameter<std::string>("ZCalib_CalibType");
  ZCalib_InvMass_ = iConfig.getUntrackedParameter<std::string>("ZCalib_InvMass");
  barrelfile_     = iConfig.getUntrackedParameter<std::string>("initialMiscalibrationBarrel","");
  endcapfile_     = iConfig.getUntrackedParameter<std::string>("initialMiscalibrationEndcap","");
  
  // corrections
  wantEtaCorrection_ = iConfig.getUntrackedParameter<bool>("wantEtaCorrection",true);   
  P0_  = iConfig.getUntrackedParameter<double>("p0_EB", 1.);   
  P1_  = iConfig.getUntrackedParameter<double>("p2_EB", 0.);   
  P2_  = iConfig.getUntrackedParameter<double>("p4_EB", 0.);   
  P0e_ = iConfig.getUntrackedParameter<double>("p0_EE", 1.);   
  P1e_ = iConfig.getUntrackedParameter<double>("p2_EE", 0.);   
  P2e_ = iConfig.getUntrackedParameter<double>("p4_EE", 0.);   

  // for calbration
  etaBins_ = iConfig.getUntrackedParameter<unsigned int>("etaBins", 10);   
  etBins_  = iConfig.getUntrackedParameter<unsigned int>("etBins", 10);   
  etaMin_  = iConfig.getUntrackedParameter<double>("etaMin", 0.);   
  etMin_   = iConfig.getUntrackedParameter<double>("etMin", 0.);   
  etaMax_  = iConfig.getUntrackedParameter<double>("etaMax", 3.);   
  etMax_   = iConfig.getUntrackedParameter<double>("etMax", 100.);   
  
  // electron combination selection 
  electronSelection_ = iConfig.getUntrackedParameter<unsigned int> ("electronSelection",0); 

  // ewk selection 
  etCut_                   = iConfig.getUntrackedParameter<double>("etCut", 4.);
  useHoESelection_         = iConfig.getUntrackedParameter<bool>  ("useHoESelection",true);
  useSigmaEtaEtaSelection_ = iConfig.getUntrackedParameter<bool>  ("useSigmaEtaEtaSelection",true);
  useDeltaPhiSelection_    = iConfig.getUntrackedParameter<bool>  ("useDeltaPhiSelection",true);
  useDeltaEtaSelection_    = iConfig.getUntrackedParameter<bool>  ("useDeltaEtaSelection",true);
  useEmIsolSelection_      = iConfig.getUntrackedParameter<bool>  ("useEmIsolSelection", true);
  useHadIsolSelection_     = iConfig.getUntrackedParameter<bool>  ("useHadIsolSelection",true);
  useTrackerIsolSelection_ = iConfig.getUntrackedParameter<bool>  ("useTrackerIsolSelection",true);
  useInvMassSelection_     = iConfig.getUntrackedParameter<bool>  ("useInvMassSelection",true);   
  detaEBCut_               = iConfig.getUntrackedParameter<double>("dEta_barrel", 0.006);
  detaEECut_               = iConfig.getUntrackedParameter<double>("dEta_endcap", 0.0075);
  dphiEBCut_               = iConfig.getUntrackedParameter<double>("dPhi_barrel", 0.048);
  dphiEECut_               = iConfig.getUntrackedParameter<double>("dPhi_endcap", 0.048);
  sigmaEtaEtaEBCut_        = iConfig.getUntrackedParameter<double>("sigmaEtaEta_barrel", 0.0115);
  sigmaEtaEtaEECut_        = iConfig.getUntrackedParameter<double>("sigmaEtaEta_endcap", 0.03);  
  HoeEBCut_                = iConfig.getUntrackedParameter<double>("HoE_barrel", 0.041);
  HoeEECut_                = iConfig.getUntrackedParameter<double>("HoE_endcap", 0.041);
  emIsolCut_               = iConfig.getUntrackedParameter<double>("emIsol", 17.);
  hcalIsolCut_             = iConfig.getUntrackedParameter<double>("hcalIsol", 1.5);
  trackerIsolCut_          = iConfig.getUntrackedParameter<double>("trackerIsol", 2.1);
  minInvMassCut_           = iConfig.getUntrackedParameter<double>("minInvMassCut", 70.);   
  maxInvMassCut_           = iConfig.getUntrackedParameter<double>("maxInvMassCut", 110.);   

  

  // --------------------    initializations  -----------------------------------------------

  // output file / tree
  outputFileName_  = iConfig.getParameter<std::string>("outputFile");
  outputFile_ = TFile::Open(outputFileName_.c_str(),"RECREATE"); 
  myTree = new TTree("myTree","myTree");
  myTree->Branch("zMassCorr",    &massCorr4tree,    "massCorr/F");
  myTree->Branch("zMassCorrBB",  &massCorr4treeBB,  "massCorrBB/F");
  myTree->Branch("zMassCorrEB",  &massCorr4treeEB,  "massCorrEB/F");
  myTree->Branch("zMassCorrEE",  &massCorr4treeEE,  "massCorrEE/F");
  myTree->Branch("zMassCorrDiff",&massCorrDiff4tree,"massDiff/F");
  
  // ECAL topology  
  EcalIndexingTools* myIndexTool=0;
  myIndexTool = EcalIndexingTools::getInstance();  
  myIndexTool -> setBinRange( etaBins_, etaMin_, etaMax_, etBins_, etMin_, etMax_ );
  
  // creating the algorithm
  theAlgorithm_ = new ZIterativeAlgorithmWithFit(iConfig);
  
  // tell the framework what data is being produced
  setWhatProduced (this, &JPsiCalibration::produceEcalIntercalibConstants ) ;
  findingRecord<EcalIntercalibConstantsRcd> () ;

  for(int ii=0; ii<50; ii++){
    coefficientDistanceAtIteration[ii] = -1.;
    loopArray[ii]  = -1.;
    sigmaArray[ii] = -1.;
    sigmaErrorArray[ii] = -1.;
  }
  
  std::cout<<"[JPsiCalibration] Done with the ctor"<<std::endl;
}


JPsiCalibration::~JPsiCalibration() { }

boost::shared_ptr<EcalIntercalibConstants> JPsiCalibration::produceEcalIntercalibConstants( const EcalIntercalibConstantsRcd& iRecord ) {

  std::cout << "@SUB=JPsiCalibration::produceEcalIntercalibConstants" << std::endl;
  return ical;
}

void JPsiCalibration::beginOfJob( const edm::EventSetup& iSetup ) {
  
  std::cout<<"[JPsiCalibration] Entering beginOfJob"<<std::endl;

  loopFlag_ = 0;

  // counters 
  passedDeta = 0;
  passedDphi = 0;
  passedSee  = 0;
  passedHoE  = 0;
  
  // geometry initialization
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);     
  EcalRingCalibrationTools::setCaloGeometry(&(*pG));  
     
  // histos and plots
  myJPsiPlots_ = new JPsiPlots("jpsiPlots.root");  
  outputFile_ -> cd();
  bookHistograms();
  
  // read miscalibration map if requested
  CaloMiscalibMapEcal* miscalibMap=0;
  if (!barrelfile_.empty() || !barrelfile_.empty()) {
    miscalibMap = new CaloMiscalibMapEcal();
    miscalibMap -> prefillMap();
  }
  
  if(!barrelfile_.empty()) {
    MiscalibReaderFromXMLEcalBarrel barrelreader_(*miscalibMap);
    barrelreader_.parseXMLMiscalibFile(barrelfile_);
    std::cout<<"[JPsiCalibration::beginOfJob] Parsed EB miscal file"<<std::endl;
  }
  
  if(!endcapfile_.empty()) {
    MiscalibReaderFromXMLEcalEndcap endcapreader_(*miscalibMap);
    endcapreader_.parseXMLMiscalibFile(endcapfile_);
    std::cout<<"[JPsiCalibration::beginOfJob] Parsed EE miscal file"<<std::endl;
  }
  
  // set miscalibration  
  for(int k = 0; k < theAlgorithm_->getNumberOfChannels(); k++) {
    calibCoeff[k]      = 1.;
    calibCoeffError[k] = 0.;
     
    std::vector<DetId> ringIds;
    if(calibMode_ == "RING")   ringIds = EcalRingCalibrationTools::getDetIdsInRing(k);
    if(calibMode_ == "MODULE") ringIds = EcalRingCalibrationTools::getDetIdsInModule(k);
    if(calibMode_ == "ABS_SCALE" || calibMode_ == "ETA_ET_MODE" ) ringIds = EcalRingCalibrationTools::getDetIdsInECAL();
	      
    if (miscalibMap) {
      initCalibCoeff[k]=0.;	      
      for (unsigned int iid=0; iid<ringIds.size();++iid) {
	float miscalib=* (miscalibMap->get().getMap().find(ringIds[iid])  );
	initCalibCoeff[k]+=miscalib;
      }
      initCalibCoeff[k]/=(float)ringIds.size();
      std::cout << "[JPsiCalibration] Set miscalib.for channel "  << k << " " << initCalibCoeff[k] << " " << ringIds.size() << std::endl;
    }
    else {
      std::cout << "[JPsiCalibration] Using 1. as calibration coefficient ... " << endl;
      initCalibCoeff[k]=1.;
    }
  }

  ical = boost::shared_ptr<EcalIntercalibConstants>( new EcalIntercalibConstants() );
  for(int k=0; k<theAlgorithm_->getNumberOfChannels(); k++) {

    std::vector<DetId> ringIds;
    if(calibMode_ == "RING")   ringIds = EcalRingCalibrationTools::getDetIdsInRing(k);
    if(calibMode_ == "MODULE") ringIds = EcalRingCalibrationTools::getDetIdsInModule(k);
    if(calibMode_ == "ABS_SCALE" || calibMode_ == "ETA_ET_MODE") ringIds = EcalRingCalibrationTools::getDetIdsInECAL();
	
    for (unsigned int iid=0; iid<ringIds.size();++iid){
	
      if(ringIds[iid].subdetId() == EcalBarrel){
	EBDetId myEBDetId(ringIds[iid]);  
	h2_xtalMiscalibCoeffBarrel_->SetBinContent( myEBDetId.ieta() + 86, myEBDetId.iphi(), * (miscalibMap->get().getMap().find(ringIds[iid]) ) ); // fill TH2 with miscalibCoeff
      }
      
      if(ringIds[iid].subdetId() == EcalEndcap){
	EEDetId myEEDetId(ringIds[iid]);
	if(myEEDetId.zside() < 0)
	  h2_xtalMiscalibCoeffEndcapMinus_->SetBinContent( myEEDetId.ix(), myEEDetId.iy(), * ( miscalibMap->get().getMap().find(ringIds[iid]) ) );  // fill TH2 with miscalibCoeff
	if(myEEDetId.zside() > 0)
	  h2_xtalMiscalibCoeffEndcapPlus_->SetBinContent( myEEDetId.ix(), myEEDetId.iy(), * (miscalibMap->get().getMap().find(ringIds[iid]) ) );    // fill TH2 with miscalibCoeff
      }
      
      ical->setValue( ringIds[iid], *(miscalibMap->get().getMap().find(ringIds[iid])  ) );
    }

  
  read_events = 0;
  init_ = false;
  }
  std::cout<<"[JPsiCalibration] Done with beginOfJob"<<std::endl;
}


void JPsiCalibration::endOfJob() {
  
  printStatistics();

  // writing out calibration coefficients  
  if(calibMode_ != "ETA_ET_MODE"){
    
    calibXMLwriter* barrelWriter = new calibXMLwriter(EcalBarrel);
    for(int ieta=-EBDetId::MAX_IETA; ieta<=EBDetId::MAX_IETA ;++ieta) {
      if(ieta==0) continue;
      for(int iphi=EBDetId::MIN_IPHI; iphi<=EBDetId::MAX_IPHI; ++iphi) {
	if (EBDetId::validDetId(ieta,iphi)) {
	  EBDetId ebid(ieta,iphi);
	  barrelWriter->writeLine(ebid,* (ical->getMap().find(ebid.rawId()) ));
	}
      }
    }
    
    calibXMLwriter* endcapWriter = new calibXMLwriter(EcalEndcap);
    for(int iX=EEDetId::IX_MIN; iX<=EEDetId::IX_MAX ;++iX) {
      for(int iY=EEDetId::IY_MIN; iY<=EEDetId::IY_MAX; ++iY) {
	if (EEDetId::validDetId(iX,iY,1)) {	  
	  EEDetId eeid(iX,iY,1);
	  endcapWriter->writeLine(eeid,*(ical->getMap().find(eeid.rawId())  ) );
	}
	if (EEDetId::validDetId(iX,iY,-1)) {
	  EEDetId eeid(iX,iY,-1);
	  endcapWriter->writeLine(eeid, *(ical->getMap().find(eeid.rawId())) );
	}
      }
    }
  } 


  // saving histos
  std::cout<<"Writing  histos..."<<std::endl;
  outputFile_->cd();
  
  h1_eventsBeforeHLTSelection_     -> Write();
  h1_eventsAfterHLTSelection_      -> Write();
  h1_eventsBefore2GsfSelection_    -> Write();
  h1_eventsAfter2GsfSelection_     -> Write();  
  h1_eventsBeforeEtSelection_      -> Write();
  h1_eventsAfterEtSelection_       -> Write();
  h1_eventsBeforeEleIDSelection_   -> Write();
  h1_eventsAfterEleIDSelection_    -> Write();
  h1_eventsBeforeEleIsolSelection_ -> Write();
  h1_eventsAfterEleIsolSelection_  -> Write();
  h1_eventsBeforeTrackerIsolSelection_ -> Write();
  h1_eventsAfterTrackerIsolSelection_  -> Write();
  h1_eventsBeforeBorderSelection_  -> Write();
  h1_eventsAfterBorderSelection_   -> Write();
  h1_eventsBeforeInvMassSelection  -> Write();
  h1_eventsAfterInvMassSelection   -> Write();

  h1_ZCandMult_->Write(); 

  h2_occupancy_vsEtaPhi       -> Write();
  h1_seedOverSC_Barrel_       -> Write();
  h1_seedOverSC_Endcap_       -> Write();  
  h1_preshowerOverSC_         -> Write();
  h1_electronCosTheta_SC_     -> Write();
  h1_electronCosTheta_TK_     -> Write();
  h1_electronCosTheta_SC_TK_  -> Write();

  h_ESCEtrue_vs_ESC           -> Write();

  h1_Theta12Resolution_TK_BB_ -> Write();
  h1_Theta12Resolution_SC_BB_ -> Write();
  h1_Theta12Resolution_TK_EE_ -> Write();
  h1_Theta12Resolution_SC_EE_ -> Write();
  h1_Theta12_TK_BB_           -> Write();
  h1_Theta12_SC_BB_           -> Write();
  h1_Theta12_TK_EE_           -> Write();
  h1_Theta12_SC_EE_           -> Write();
  h1_zMassResol_              -> Write();
  h1_eleEtaResol_             -> Write();
  h1_elePhiResol_             -> Write(); 

  h1_reco_ZMass_              -> Write();
  h1_reco_ZMassBB_            -> Write();
  h1_reco_ZMassEB_            -> Write();
  h1_reco_ZMassEE_            -> Write();
  h1_reco_ZMassCorr_          -> Write();  
  h1_reco_ZMassCorrBB_        -> Write();
  h1_reco_ZMassCorrEB_        -> Write();
  h1_reco_ZMassCorrEE_        -> Write();

  h1_massResolutionTerm_total_BB_       -> Write();
  h1_massResolutionTerm_angle_BB_       -> Write();
  h1_massResolutionTerm_energy_BB_      -> Write();
  h1_massResolutionTerm_angle_ratio_BB_ -> Write();  
  h1_massResolutionTerm_total_EE_       -> Write();
  h1_massResolutionTerm_angle_EE_       -> Write();
  h1_massResolutionTerm_energy_EE_      -> Write();
  h1_massResolutionTerm_angle_ratio_EE_ -> Write();

  h2_xtalMiscalibCoeffBarrel_      -> Write();
  h2_xtalMiscalibCoeffEndcapMinus_ -> Write();
  h2_xtalMiscalibCoeffEndcapPlus_  -> Write();

  h1_sEEeb_before_ -> Write();
  h1_sEEee_before_ -> Write();
  h1_sEEeb_after_  -> Write();
  h1_sEEee_after_  -> Write();
  
  for(int i =0; i<25; i++){
    h_ESCEtrue_array_barrel[i]->Write();
    h_ESCEtrue_array_endcap[i]->Write();    
    if( i < theMaxLoops ){      
      h_PTres_[i]             -> Write();
      h_PTres_Endcap_[i]      -> Write();
      h_ESCEtrueVsEta_[i]     -> Write();
      h_ESCEtrue_[i]          -> Write();
      h_ESCcorrEtrueVsEta_[i] -> Write();
      h_ESCcorrEtrue_[i]      -> Write();
      h2_chi2_[i]             -> Write();
      h2_iterations_[i]       -> Write();
    }}

  for(int i =0; i<2; i++){
    h_eleEffEta_[i] -> Write();
    h_eleEffPhi_[i] -> Write();
    h_eleEffPt_[i]  -> Write();
  }
 
  h1_efficiencySummary_ ->Fill("All",          h1_eventsBeforeHLTSelection_    -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("HLT",          h1_eventsAfterHLTSelection_     -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("2 GSF",        h1_eventsAfter2GsfSelection_    -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("E_{T}",        h1_eventsAfterEtSelection_      -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("eleID",        h1_eventsAfterEleIDSelection_   -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("isol",         h1_eventsAfterEleIsolSelection_ -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("tracker isol", h1_eventsAfterTrackerIsolSelection_ -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("Borders",      h1_eventsAfterBorderSelection_  -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );
  h1_efficiencySummary_ ->Fill("M_{ee}",       h1_eventsAfterInvMassSelection  -> GetEntries() / h1_eventsBeforeHLTSelection_->GetEntries() );

  // computing coefficients
  int j = 0;
  int flag=0;  
  Double_t mean[25] = {0.};
  Double_t num[25]  = {0.};
  Double_t meanErr[25] = {0.};
  Float_t rms[25]   = {0.};
  Float_t tempRms[10][25];
  for(int ia = 0; ia<10; ia++){
    for(int ib = 0; ib<25; ib++){ tempRms[ia][ib] = 0.; }}
      
  int aa = 0;  
  for( int k = 0; k<theAlgorithm_->getNumberOfChannels(); k++ ) {
    
    bool isNearCrack = false;
    
    if( calibMode_ == "RING"){
      
      isNearCrack = ( abs( ringNumberCorrector(k) ) == 1 || abs( ringNumberCorrector(k) ) == 25 ||
		      abs( ringNumberCorrector(k) ) == 26 || abs( ringNumberCorrector(k) ) == 45 ||
		      abs( ringNumberCorrector(k) ) == 46 || abs( ringNumberCorrector(k) ) == 65 ||
		      abs( ringNumberCorrector(k) ) == 66 || abs( ringNumberCorrector(k) ) == 85 ||
		      abs( ringNumberCorrector(k) ) == 86 || abs( ringNumberCorrector(k) ) == 124 );
    }
      
    if(k<85){ 	  
      if((k+1)%5!=0){
	      
	if(!isNearCrack){
	  num[j] += 2.;
	  mean[j]+=calibCoeff[k];
	  mean[j]+=calibCoeff[169 - k];
	  meanErr[j]+= 1./ pow ( calibCoeffError[k], 2 );
	  meanErr[j]+= 1./ pow ( calibCoeffError[169 - k], 2);

	  tempRms[aa][j]+=calibCoeff[k];
	  aa++;
	  tempRms[aa][j]+=calibCoeff[169 - k];
	  aa++;
	}
      }
      else {
	if(!isNearCrack){
	  num[j] += 2.;
	  mean[j]+=calibCoeff[k];
	  mean[j]+=calibCoeff[169 - k];
	  meanErr[j]+= 1./ pow ( calibCoeffError[k], 2 );
	  meanErr[j]+= 1./ pow ( calibCoeffError[169 - k], 2);
	  
	  tempRms[aa][j]+=calibCoeff[k];
	  aa++;
	  tempRms[aa][j]+=calibCoeff[169 - k];
	  aa++;
	}
	j++;
	aa = 0;
      }
    }

    // EE begin           
    if(k>=170 && k<=204){	
      if(flag<4){	
	mean[j]+=calibCoeff[k]/10.;                          // make groups of 5 Xtals in #eta
	mean[j]+=calibCoeff[k+39]/10.;
	meanErr[j]+= calibCoeffError[k]/30.;
	meanErr[j]+= calibCoeffError[k + 39]/30.;
	
	tempRms[aa][j]+=calibCoeff[k];
	aa++;
	tempRms[aa][j]+=calibCoeff[k + 39];
	aa++;
	
	flag++;
      }
      else if(flag==4){
	mean[j]+=calibCoeff[k]/10.;                       // make groups of 5 Xtals in #eta
	mean[j]+=calibCoeff[k+39]/10.;
	meanErr[j]+= calibCoeffError[k]/30.;
	meanErr[j]+= calibCoeffError[k + 39]/30.;
	
	tempRms[aa][j]+=calibCoeff[k];
	aa++;
	tempRms[aa][j]+=calibCoeff[k + 39];
	aa++;
	
	flag=0;
	j++;
	aa = 0;
      }
    }
    if(k>=205 && k<=208){
      mean[j]+=calibCoeff[k]/8.;
      mean[j]+=calibCoeff[k+39]/8.;      
      meanErr[j]+= calibCoeffError[k]/30.;
      meanErr[j]+= calibCoeffError[k + 39]/30.;
            
      tempRms[aa][j]+=calibCoeff[k];
      aa++;
      tempRms[aa][j]+=calibCoeff[k + 39];
      aa++;
    } //EE end
          
    if(!isNearCrack){
      h2_coeffVsEta_  -> Fill( ringNumberCorrector(k), calibCoeff[k] );
      h2_miscalRecal_ -> Fill( initCalibCoeff[k], 1./calibCoeff[k] );
      h1_mc_          -> Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
	
      if(k<170){
	h2_miscalRecalEB_ -> Fill( initCalibCoeff[k], 1./calibCoeff[k] );
	h1_mcEB_          -> Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
      }
      
      if(k>=170){
	h2_miscalRecalEE_ -> Fill( initCalibCoeff[k], 1./calibCoeff[k] );
	h1_mcEE_          -> Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
      }    
    }
  }
  
  for(int ic = 0; ic< 17; ic++){
    mean[ic] = mean[ic] / num[ic];        // find mean of recalib coeff on group of rings
    meanErr[ic] = 1. / sqrt(meanErr[ic]); // find mean of recalib coeff on group of rings
  }

  // build array of RMS
  for(int ic = 0; ic< 25; ic++){
    for(int id = 0; id< 10; id++){
      if(tempRms[id][ic] > 0.) rms[ic] += (tempRms[id][ic] - mean[j])*(tempRms[id][ic] - mean[j]);
    }
    rms[ic]/= 10.;//this is approximate
    rms[ic] = sqrt(rms[ic]);
  }
  
  Double_t xtalEta[25] = {1.4425, 1.3567,1.2711,1.1855,
			  1.10,1.01,0.92,0.83,
			  0.7468,0.6612,0.5756,0.4897,0.3985,0.3117,0.2250,0.1384,0.0487,
			  1.546, 1.651, 1.771, 1.908, 2.071, 2.267, 2.516, 2.8};
  
  Double_t zero[25] = {0.026};   //interval/sqrt(12)

  for(int j = 0; j <25; j++) h2_coeffVsEtaGrouped_->Fill( xtalEta[j],mean[j]);

  TProfile *px = h2_coeffVsEta_->ProfileX("coeffVsEtaProfile");
  px -> SetXTitle("Eta channel");
  px -> SetYTitle("recalibCoeff");
  px -> Write();

  h2_coeffVsEta_        -> Write();
  h2_coeffVsEtaGrouped_ -> Write();
  h2_zMassVsLoop_       -> Write();
  h2_zMassDiffVsLoop_   -> Write();
  h2_zWidthVsLoop_      -> Write();
  h2_coeffVsLoop_       -> Write();
  h2_miscalRecal_       -> Write();
  h1_mc_                -> Write();
  h2_miscalRecalEB_     -> Write();
  h1_mcEB_              -> Write();
  h2_miscalRecalEE_     -> Write();
  h1_mcEE_              -> Write();  
  h2_residualSigma_     -> Write();
  
  const ZIterativeAlgorithmWithFit::ZIterativeAlgorithmWithFitPlots* algoHistos=theAlgorithm_->getHistos();

  
  // write algorithm plots
  for (int iIteration=0;iIteration<theAlgorithm_->getNumberOfIterations();iIteration++)
    for (int iChannel=0;iChannel<theAlgorithm_->getNumberOfChannels();iChannel++) {      
      if( iChannel%10 == 0 ){
	algoHistos->weightedRescaleFactor[iIteration][iChannel]->Write();
	algoHistos->unweightedRescaleFactor[iIteration][iChannel]->Write();
	algoHistos->weight[iIteration][iChannel]->Write();
      }
    }
  
  double weightSumMeanBarrel = 0.;
  double weightSumMeanEndcap = 0.;
  
  for (int iIteration=0;iIteration<theAlgorithm_->getNumberOfIterations();iIteration++)
    for (int iChannel=0;iChannel<theAlgorithm_->getNumberOfChannels();iChannel++) {

      if( iIteration==(theAlgorithm_->getNumberOfIterations()-1) ){
	
	if(iChannel < 170)  weightSumMeanBarrel += algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral()/170.; 
	if(iChannel >= 170) weightSumMeanEndcap += algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral()/78.; 
	
	h1_occupancyVsEta_ -> Fill((Double_t)ringNumberCorrector(iChannel), algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral() );
	h1_occupancy_      -> Fill( algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral() );
	  
	if(iChannel < 170)  h1_occupancyBarrel_->Fill( algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral() );
	if(iChannel >= 170) h1_occupancyEndcap_->Fill( algoHistos->weightedRescaleFactor[iIteration][iChannel]->Integral() );
      }
    }
  
  h1_weightSumMeanBarrel_ ->Fill(weightSumMeanBarrel);
  h1_weightSumMeanEndcap_ ->Fill(weightSumMeanEndcap);
  std::cout<<"Weight sum mean on channels in Barrel is :"<<weightSumMeanBarrel<<std::endl;
  std::cout<<"Weight sum mean on channels in Endcap is :"<<weightSumMeanEndcap<<std::endl;  
  h1_weightSumMeanBarrel_ -> Write();
  h1_weightSumMeanEndcap_ -> Write();
  h1_occupancyVsEta_      -> Write();
  h1_occupancy_           -> Write();
  h1_occupancyBarrel_     -> Write();
  h1_occupancyEndcap_     -> Write();
  
  myTree->Write();

  TGraphErrors* graph = new TGraphErrors(25,xtalEta,mean,zero,meanErr);
  graph->Draw("APL");
  graph->Write();
  
  double zero50[50] = { 0. };

  TGraphErrors* residualSigmaGraph = new TGraphErrors(50,loopArray,sigmaArray,zero50,sigmaErrorArray);
  residualSigmaGraph->SetName("residualSigmaGraph");
  residualSigmaGraph->Draw("APL");
  residualSigmaGraph->Write();

  TGraphErrors* coefficientDistanceAtIterationGraph = new TGraphErrors(50,loopArray,coefficientDistanceAtIteration,zero50,zero50);
  coefficientDistanceAtIterationGraph->SetName("coefficientDistanceAtIterationGraph");
  coefficientDistanceAtIterationGraph->Draw("APL");
  coefficientDistanceAtIterationGraph->Write();

  Float_t noError[250] = {0.};

  Float_t ringInd[250];
  for(int i =0; i<250; i++)
    ringInd[i]=ringNumberCorrector(i);

  TGraphErrors* graphCoeff = new TGraphErrors(theAlgorithm_->getNumberOfChannels(),ringInd,calibCoeff,noError,calibCoeffError);
  graphCoeff->SetName("graphCoeff");
  graphCoeff->Draw("APL");
  graphCoeff->Write();

  outputFile_->Close();
  
  myJPsiPlots_ ->writeEleHistograms();
  myJPsiPlots_ ->writeEleHistogramsAS();
  myJPsiPlots_ ->writeMCEleHistograms();
  myJPsiPlots_ ->writeJPsiHistograms();
  myJPsiPlots_ ->writeMCJPsiHistograms();
}


// Called at each event
edm::EDLooper::Status JPsiCalibration::duringLoop( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  
  // -----------------------------------------
  // MC infos 
  std::vector<HepMC::GenParticle*> mcEle;
  float myGenJPsiMass(-1);      
  
  if (!mcProducer_.empty()) {
    
    Handle< HepMCProduct > hepProd ;
    iEvent.getByLabel( mcProducer_.c_str(), hepProd ) ;
    const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();
    
    if (loopFlag_ == 0) myJPsiPlots_ -> fillJPsiMCInfo( & (*myGenEvent), resonancePdgId_ );
    if (loopFlag_ == 0) myJPsiPlots_ -> fillEleMCInfo( & (*myGenEvent) );                      
    
    for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
      if ( ((*p)->pdg_id()==resonancePdgId_) && ((*p)->status()==2) ) myGenJPsiMass = (*p)->momentum().m();                      // j/psi infos
    }
    
    for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
      if (  abs((*p)->pdg_id())==11 ) mcEle.push_back( (*p) );                                                                   // e+ e- from j/psi
    }

    if(mcEle.size()==2 && fabs(mcEle[0]->momentum().eta())<2.4 &&  fabs(mcEle[1]->momentum().eta())<2.4 ){   // counting electrons in the acceptance
      NEVT++;      
      if(  fabs(mcEle[0]->momentum().eta())<1.479 && fabs(mcEle[1]->momentum().eta())<1.479 ) MCZBB++;
      if( (fabs(mcEle[0]->momentum().eta())>1.479 && fabs(mcEle[1]->momentum().eta())<1.479) || (fabs(mcEle[0]->momentum().eta())<1.479 && fabs(mcEle[1]->momentum().eta())>1.479) ) MCZEB++;
      if(  fabs(mcEle[0]->momentum().eta())>1.479 && fabs(mcEle[1]->momentum().eta())>1.479 ) MCZEE++;
    }
  }    

  
  // -----------------------------------------
  // HLT selection
  Handle<TriggerResults> HLTR;
  try { iEvent.getByLabel(hlTriggerResults_, HLTR); }
  catch(cms::Exception& ex){ LogError("ProblemHltTriggserResults") << "Trigger results: " << hlTriggerResults_ << " not found"; }
  if (!HLTR.isValid()) throw cms::Exception("ProductNotValid") << "TriggerResults product not valid";
  
  bool hlt29accept = false;
  if (HLTR.isValid()){
    triggerNames_.init(*HLTR);
    for ( unsigned int iHlt=0; iHlt < HLTR->size(); iHlt++ ) {
      if (triggerNames_.triggerName(iHlt)=="HLT_DoubleEM8e29_Jpsi" && HLTR->accept(iHlt)==1) hlt29accept = true;
    }
  }
  if ( hlt29accept) intHlt29 = 1;
  if (!hlt29accept) intHlt29 = 0;
  
  h1_eventsBeforeHLTSelection_ -> Fill(1);
  if (intHlt29==0) return kContinue;           
  h1_eventsAfterHLTSelection_  -> Fill(1);
  
  

  // -----------------------------------------  
  // taking all the needed collections   
  const CaloTopology *theTopology;         // calo topology 
  edm::ESHandle<CaloTopology> topology;
  iSetup.get<CaloTopologyRecord>().get(topology);
  if (topology.isValid()) theTopology = topology.product();
  
  const CaloGeometry *theGeometry;        // calo geometry
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);
  if(geometry.isValid()) theGeometry = geometry.product();

  Handle<reco::GsfElectronCollection> pElectrons;                // electrons
  iEvent.getByLabel(electronProducer_, electronCollection_, pElectrons);
  const reco::PixelMatchGsfElectronCollection* electronCollection = pElectrons.product();

  Handle<reco::SuperClusterCollection> superClustersEB;          // superclusters in EB
  iEvent.getByLabel(scProducerEB_, scCollectionEB_, superClustersEB);
  const reco::SuperClusterCollection* scCollectionEB = superClustersEB.product();
  
  Handle<reco::SuperClusterCollection> superClustersEE;         // superclusters in EE
  iEvent.getByLabel(scProducerEE_, scCollectionEE_, superClustersEE);
  const reco::SuperClusterCollection* scCollectionEE = superClustersEE.product();

  Handle<EBRecHitCollection> pRechitsEB;                        // EBRecHits (from recalib collection)
  iEvent.getByLabel( rechitProducerEB_, rechitCollectionEB_, pRechitsEB);
  const EBRecHitCollection* rechitsEB = pRechitsEB.product(); 

  Handle<EERecHitCollection> pRechitsEE;                        // EERecHits (from recalib collection)
  iEvent.getByLabel( rechitProducerEE_, rechitCollectionEE_, pRechitsEE);
  const EERecHitCollection* rechitsEE = pRechitsEE.product(); 
  
  const TrackCollection *theTracks;                             // all tracker tracks, for isolation 
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(tracksCollection_, tracks); 
  if(tracks.isValid()) theTracks = tracks.product();

  const CaloTowerCollection *theCaloTowers;                     // calotowers, for HCAL/ECAL based isolation  chiara: prendile miscalibrate
  Handle<CaloTowerCollection> calotowers;
  iEvent.getByLabel(calotowersCollection_, calotowers); 
  if (calotowers.isValid()) theCaloTowers = calotowers.product();


  // basic histos
  if (loopFlag_ == 0) myJPsiPlots_->fillEleInfo(electronCollection);   
  

  // -----------------------------------------    
  // skipping events with problems 
  if (!rechitsEB && !rechitsEE)                       return kContinue;
  if ( rechitsEB->size()==0 && rechitsEE->size()==0 ) return kContinue;
  if (!electronCollection)                            return kContinue;
  
  // skipping events not passing a very basic selection: at least 2 electrons
  h1_eventsBefore2GsfSelection_  -> Fill(1);
  if( (scCollectionEB->size()+scCollectionEE->size()) < 2) return kContinue;
  if ( electronCollection->size() < 2)                     return kContinue;
  h1_eventsAfter2GsfSelection_ -> Fill(1);
  



  
  
  // ---------------------------------------------------------------
  // start here with the analysis
  read_events++;
  
  // filling new ElectronCollection with new SC ref and calibElectron container
  std::vector<calib::CalibElectron> calibElectronsAll;
  for(unsigned int e_it=0 ; e_it!=electronCollection->size(); e_it++) {
    calibElectronsAll.push_back(calib::CalibElectron(&((*electronCollection)[e_it]),rechitsEB,rechitsEE));
  }
  if (calibElectronsAll.size() < 2) return kContinue;   // non dovrebbe succedere mai, no?
  


  // ---------------------------------------------------------------
  // at least one electron and one positron above threshold
  std::vector<calib::CalibElectron> calibElectronsGtThr;     // electrons above Et threshold
  std::vector<calib::CalibElectron> calibPositronsGtThr;     // positrons above Et threshold
  for(unsigned int e_it=0 ; e_it < calibElectronsAll.size(); e_it++){
    
    float theta = 2. * atan( exp(- calibElectronsAll[e_it].getRecoElectron()->superCluster()->eta()) );
    float et    = (calibElectronsAll[e_it].getRecoElectron()->superCluster()->energy()) * sin( theta);
    if (et>etCut_) {
      if (calibElectronsAll[e_it].getRecoElectron()->charge()>0) calibPositronsGtThr.push_back(calib::CalibElectron(calibElectronsAll[e_it].getRecoElectron(),rechitsEB,rechitsEE));
      if (calibElectronsAll[e_it].getRecoElectron()->charge()<0) calibElectronsGtThr.push_back(calib::CalibElectron(calibElectronsAll[e_it].getRecoElectron(),rechitsEB,rechitsEE));      
    }
  }

  h1_eventsBeforeEtSelection_ -> Fill(1);
  if (calibElectronsGtThr.size() < 1) return kContinue;   
  if (calibPositronsGtThr.size() < 1) return kContinue;   
  h1_eventsAfterEtSelection_  -> Fill(1);



  // ---------------------------------------------------------------
  // electron id requirements
  std::vector<calib::CalibElectron> calibElectronsID;     // electrons above Et threshold and identified
  std::vector<calib::CalibElectron> calibPositronsID;     // positrons above Et threshold and identified

  double detaCut_(-1), dphiCut_(-1);
  double sigmaEtaEtaCut_(-1), HoeCut_(-1);

  for(unsigned int e_it=0 ; e_it < calibElectronsGtThr.size(); e_it++){
    
    // cluster shape variables
    BasicClusterRef theSeedRef = calibElectronsGtThr[e_it].getRecoElectron()->superCluster()->seed();
    const EcalRecHitCollection *rechits = 0;
    float seedEta = theSeedRef->position().eta();      
    if( fabs(seedEta) < 1.479 ) rechits = rechitsEB;
    else rechits = rechitsEE;     
    std::vector<float> covMatrix = EcalClusterTools::covariances( *theSeedRef, rechits, theTopology, theGeometry );
    float sEtaEta = sqrt(covMatrix[0]); 
    if( fabs(seedEta) >= 1.479 ) sEtaEta = sEtaEta - 0.02*(fabs(seedEta) - 2.3);
    
    int theEleClass = calibElectronsGtThr[e_it].getRecoElectron()->classification();
    if(theEleClass < 100) { 
      detaCut_        = detaEBCut_;
      dphiCut_        = dphiEBCut_;
      sigmaEtaEtaCut_ = sigmaEtaEtaEBCut_;
      HoeCut_         = HoeEBCut_;
    }
    else {
      detaCut_        = detaEECut_;
      dphiCut_        = dphiEECut_;
      sigmaEtaEtaCut_ = sigmaEtaEtaEECut_;
      HoeCut_         = HoeEECut_;
    }

    bool DeltaEtaIn = ( calibElectronsGtThr[e_it].getRecoElectron()->deltaEtaSuperClusterTrackAtVtx() < detaCut_);
    bool DeltaPhiIn = ( calibElectronsGtThr[e_it].getRecoElectron()->deltaPhiSuperClusterTrackAtVtx() < dphiCut_);
    bool HoE        = ( calibElectronsGtThr[e_it].getRecoElectron()->hadronicOverEm() < HoeCut_);
    bool sigmaEE    = ( sEtaEta < sigmaEtaEtaCut_ );

    bool isGood = true;
    if( useDeltaEtaSelection_ && !DeltaEtaIn ) isGood = false;
    if (isGood) passedDeta++;   // to debug
    if( useDeltaPhiSelection_ && !DeltaPhiIn ) isGood = false;
    if (isGood) passedDphi++;
    if( useHoESelection_ && !HoE )             isGood = false;
    if (isGood) passedHoE++;
    if( useSigmaEtaEtaSelection_ && !sigmaEE ) isGood = false; 
    if (isGood) passedSee++;

    // to debug
    if (theEleClass < 100) h1_sEEeb_before_ -> Fill(sEtaEta);
    if (theEleClass >=100) h1_sEEee_before_ -> Fill(sEtaEta);
    if (theEleClass < 100 && useSigmaEtaEtaSelection_ && sigmaEE) h1_sEEeb_after_ -> Fill(sEtaEta);
    if (theEleClass >=100 && useSigmaEtaEtaSelection_ && sigmaEE) h1_sEEee_after_ -> Fill(sEtaEta);
    
    if (isGood) calibElectronsID.push_back(calib::CalibElectron(calibElectronsGtThr[e_it].getRecoElectron(),rechitsEB,rechitsEE));
  }

  for(unsigned int e_it=0 ; e_it < calibPositronsGtThr.size(); e_it++){

    // cluster shape variables
    BasicClusterRef theSeedRef = calibPositronsGtThr[e_it].getRecoElectron()->superCluster()->seed();
    const EcalRecHitCollection *rechits = 0;
    float seedEta = theSeedRef->position().eta();      
    if( fabs(seedEta) < 1.479 ) rechits = rechitsEB;
    else rechits = rechitsEE;     
    std::vector<float> covMatrix = EcalClusterTools::covariances( *theSeedRef, rechits, theTopology, theGeometry );
    float sEtaEta = sqrt(covMatrix[0]); 
    if( fabs(seedEta) >= 1.479 ) sEtaEta = sEtaEta - 0.02*(fabs(seedEta) - 2.3);
    
    int theEleClass = calibPositronsGtThr[e_it].getRecoElectron()->classification();
    if(theEleClass < 100) { 
      detaCut_        = detaEBCut_;
      dphiCut_        = dphiEBCut_;
      sigmaEtaEtaCut_ = sigmaEtaEtaEBCut_;
      HoeCut_         = HoeEBCut_;
    }
    else {
      detaCut_        = detaEECut_;
      dphiCut_        = dphiEECut_;
      sigmaEtaEtaCut_ = sigmaEtaEtaEECut_;
      HoeCut_         = HoeEECut_;
    }

    bool DeltaEtaIn = ( calibPositronsGtThr[e_it].getRecoElectron()->deltaEtaSuperClusterTrackAtVtx() < detaCut_);
    bool DeltaPhiIn = ( calibPositronsGtThr[e_it].getRecoElectron()->deltaPhiSuperClusterTrackAtVtx() < dphiCut_);
    bool HoE        = ( calibPositronsGtThr[e_it].getRecoElectron()->hadronicOverEm() < HoeCut_);
    bool sigmaEE    = ( sEtaEta < sigmaEtaEtaCut_ );    

    bool isGood = true;
    if( useDeltaEtaSelection_ && !DeltaEtaIn ) isGood = false;
    if (isGood) passedDeta++;   // to debug
    if( useDeltaPhiSelection_ && !DeltaPhiIn ) isGood = false;
    if (isGood) passedDphi++;
    if( useHoESelection_ && !HoE )             isGood = false;
    if (isGood) passedHoE++;
    if( useSigmaEtaEtaSelection_ && !sigmaEE ) isGood = false; 
    if (isGood) passedSee++;
    
    if (isGood) calibPositronsID.push_back(calib::CalibElectron(calibPositronsGtThr[e_it].getRecoElectron(),rechitsEB,rechitsEE));
  }

  h1_eventsBeforeEleIDSelection_ -> Fill(1);
  if (calibElectronsID.size() < 1) return kContinue;   
  if (calibPositronsID.size() < 1) return kContinue;   
  h1_eventsAfterEleIDSelection_  -> Fill(1);



  // ---------------------------------------------------------------
  // electron isolation requirements
  std::vector<calib::CalibElectron> calibElectronsIsol;     // electrons above Et threshold and identified and isolated
  std::vector<calib::CalibElectron> calibPositronsIsol;     // positrons above Et threshold and identified and isolated
  
  for(unsigned int e_it=0 ; e_it < calibElectronsID.size(); e_it++){
    
    // calo towers based isolation in ECAL and HCAL   
    float theta  = 2. * atan( exp(- calibElectronsID[e_it].getRecoElectron()->superCluster()->eta()) );
    float energy = calibElectronsID[e_it].getRecoElectron()->superCluster()->energy();

    SuperClusterRef theScRef = calibElectronsID[e_it].getRecoElectron()->superCluster();
    CalotowerIsolation calotowerIsolation( &(*theScRef), theCaloTowers);
    calotowerIsolation.setIntRadius(0.1);
    calotowerIsolation.setExtRadius(0.3);	
    float sumHadEt03 = calotowerIsolation.getEtHcal(true);
    float sumEmEt03  = calotowerIsolation.getEtEcal(true);	
    
    bool emIsol   = ( sumEmEt03*energy*sin(theta) < emIsolCut_);
    bool hcalIsol = ( sumHadEt03*energy*sin(theta) < hcalIsolCut_);
    
    bool isGood = true;
    if ( useEmIsolSelection_  && !emIsol )   isGood = false;
    if ( useHadIsolSelection_ && !hcalIsol ) isGood = false;
    if (isGood) calibElectronsIsol.push_back(calib::CalibElectron(calibElectronsID[e_it].getRecoElectron(),rechitsEB,rechitsEE));
  }

  for(unsigned int e_it=0 ; e_it < calibPositronsID.size(); e_it++){
    
    // calo towers based isolation in ECAL and HCAL   
    float theta  = 2. * atan( exp(- calibPositronsID[e_it].getRecoElectron()->superCluster()->eta()) );
    float energy = calibPositronsID[e_it].getRecoElectron()->superCluster()->energy();

    SuperClusterRef theScRef = calibPositronsID[e_it].getRecoElectron()->superCluster();
    CalotowerIsolation calotowerIsolation( &(*theScRef), theCaloTowers);
    calotowerIsolation.setIntRadius(0.1);
    calotowerIsolation.setExtRadius(0.3);	
    float sumHadEt03 = calotowerIsolation.getEtHcal(true);
    float sumEmEt03  = calotowerIsolation.getEtEcal(true);	
    
    bool emIsol   = ( sumEmEt03*energy*sin(theta) < emIsolCut_);
    bool hcalIsol = ( sumHadEt03*energy*sin(theta) < hcalIsolCut_);
    
    bool isGood = true;
    if ( useEmIsolSelection_  && !emIsol )   isGood = false;
    if ( useHadIsolSelection_ && !hcalIsol ) isGood = false;
    if (isGood) calibPositronsIsol.push_back(calib::CalibElectron(calibPositronsID[e_it].getRecoElectron(),rechitsEB,rechitsEE));
  }

  h1_eventsBeforeEleIsolSelection_ -> Fill(1);
  if (calibElectronsIsol.size() < 1) return kContinue;   
  if (calibPositronsIsol.size() < 1) return kContinue;   
  h1_eventsAfterEleIsolSelection_  -> Fill(1);



  // ---------------------------------------------------------------
  // after the full selection (a part for tracker isolation, which needs a choice) we look for the two highest ET electrons (with charge requirement)
  unsigned int theHighestEtp = 999;
  unsigned int theHighestEte = 999;
  float highestEtp  = -1.;
  float highestEte  = -1.;
  
  for(unsigned int e_it=0 ; e_it < calibElectronsIsol.size(); e_it++){
    float theta = 2. * atan( exp(- calibElectronsIsol[e_it].getRecoElectron()->superCluster()->eta()) );
    float et    = (calibElectronsIsol[e_it].getRecoElectron()->superCluster()->energy()) * sin( theta);
    if (et>highestEte) {
      highestEte    = et;
      theHighestEte = e_it;
    }}

  for(unsigned int e_it=0 ; e_it < calibPositronsIsol.size(); e_it++){  
    float theta = 2. * atan( exp(- calibPositronsIsol[e_it].getRecoElectron()->superCluster()->eta()) );
    float et    = (calibPositronsIsol[e_it].getRecoElectron()->superCluster()->energy()) * sin( theta);
    if (et>highestEtp) {
      highestEtp    = et;
      theHighestEtp = e_it;
    }}
  


  // -------------------------------------
  // tracker isolation. Here maybe eff loss, but anyway.... (to be fixed)
  TVector3 candEle, candPos;
  candEle.SetXYZ(calibElectronsIsol[theHighestEte].getRecoElectron()->px(), calibElectronsIsol[theHighestEte].getRecoElectron()->py(), calibElectronsIsol[theHighestEte].getRecoElectron()->pz());
  candPos.SetXYZ(calibPositronsIsol[theHighestEtp].getRecoElectron()->px(), calibPositronsIsol[theHighestEtp].getRecoElectron()->py(), calibPositronsIsol[theHighestEtp].getRecoElectron()->pz());   
  float deltaRBetweenEles = candEle.DeltaR(candPos);

  GsfTrackRef theGsfEleTrack = calibElectronsIsol[theHighestEte].getRecoElectron()->gsfTrack();
  GsfTrackRef theGsfPosTrack = calibPositronsIsol[theHighestEtp].getRecoElectron()->gsfTrack();
  TrackerIsolation trackIsolationEle( &(*theGsfEleTrack), theTracks );
  TrackerIsolation trackIsolationPos( &(*theGsfPosTrack), theTracks );
  trackIsolationEle.setIntRadius(0.02);
  trackIsolationPos.setIntRadius(0.02);
  trackIsolationEle.setExtRadius(0.3);
  trackIsolationPos.setExtRadius(0.3);
  float trackerIsolEleRel = trackIsolationEle.getPtTracks(true);
  float trackerIsolPosRel = trackIsolationPos.getPtTracks(true);

  float thetaEle  = 2. * atan( exp(- calibElectronsIsol[theHighestEte].getRecoElectron()->superCluster()->eta()) );
  float thetaPos  = 2. * atan( exp(- calibPositronsIsol[theHighestEtp].getRecoElectron()->superCluster()->eta()) );
  float etEle     = (calibElectronsIsol[theHighestEte].getRecoElectron()->superCluster()->energy()) * sin( thetaEle);
  float etPos     = (calibPositronsIsol[theHighestEtp].getRecoElectron()->superCluster()->energy()) * sin( thetaPos);
  float energyEle = calibElectronsIsol[theHighestEte].getRecoElectron()->superCluster()->energy();
  float energyPos = calibPositronsIsol[theHighestEtp].getRecoElectron()->superCluster()->energy();

  if (deltaRBetweenEles < 0.3) {
    trackerIsolPosRel -= (etEle/etPos);
    if (trackerIsolPosRel < 0.) trackerIsolPosRel = 0;
    trackerIsolEleRel -= (etPos/etEle);
    if (trackerIsolEleRel < 0.) trackerIsolEleRel = 0;
  }
  
  float trackerIsolEle = ( trackerIsolEleRel*energyEle*sin(thetaEle) < trackerIsolCut_);
  float trackerIsolPos = ( trackerIsolPosRel*energyPos*sin(thetaPos) < trackerIsolCut_);
  h1_eventsBeforeTrackerIsolSelection_ -> Fill(1);
  if ( useTrackerIsolSelection_  && !trackerIsolEle ) return kContinue;   
  if ( useTrackerIsolSelection_  && !trackerIsolPos ) return kContinue;   
  h1_eventsAfterTrackerIsolSelection_ -> Fill(1);



  // ---------------------------------------------------------------
  // combinatory for j/psi mass: now the J/psi from the highest et candidates is chosen
  std::vector< pair<calib::CalibElectron*,calib::CalibElectron*> > zeeCandidates;   // all possible combinations
  int  myBestZ = -1;
  mass = -1.;
  
  for(unsigned int e_it = 0 ; e_it < calibElectronsIsol.size(); e_it++){
    for(unsigned int p_it = 0; p_it < calibPositronsIsol.size() ; p_it++) {
  
      mass = ZeeKinematicTools::calculateZMass_withTK(std::pair<calib::CalibElectron*,calib::CalibElectron*>(&(calibElectronsIsol[e_it]),&(calibPositronsIsol[p_it])));
      if (mass<0) continue;
      
      zeeCandidates.push_back(std::pair<calib::CalibElectron*,calib::CalibElectron*>(&(calibElectronsIsol[e_it]),&(calibPositronsIsol[p_it])));
      if ( (e_it==theHighestEte) && (p_it==theHighestEtp) ) myBestZ=zeeCandidates.size()-1;        
    }}

  h1_ZCandMult_->Fill(zeeCandidates.size());  
  if(zeeCandidates.size()==0 || myBestZ==-1 ) return kContinue;   // non dovrebbe succedere mai, credo
  
  

  // ---------------------------------------------------------------
  // electron class to separate barrel and encap	
  int class1 = zeeCandidates[myBestZ].first->getRecoElectron()->classification();
  int class2 = zeeCandidates[myBestZ].second->getRecoElectron()->classification();
  
  // some statistics
  if( class1 < 100 )  TOTAL_ELECTRONS_IN_BARREL++;    
  if( class1 >= 100 ) TOTAL_ELECTRONS_IN_ENDCAP++;
  if( class2 < 100 )  TOTAL_ELECTRONS_IN_BARREL++;
  if( class2 >= 100 ) TOTAL_ELECTRONS_IN_ENDCAP++;
  if( class1==40)     CRACK_ELECTRONS_IN_BARREL++;
  if( class1==140)    CRACK_ELECTRONS_IN_ENDCAP++;
  if( class2==40)     CRACK_ELECTRONS_IN_BARREL++;
  if( class2==140)    CRACK_ELECTRONS_IN_ENDCAP++;



  // ---------------------------------------------------------------
  // exclude electrons with the hottest crystal at border of modules and electrons in cracks
  // exclude electrons in barrel or endcap if required
  h1_eventsBeforeBorderSelection_->Fill(1);

  DetId firstElehottestDetId  = getHottestDetId( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->seed()->getHitsByDetId() , rechitsEB,  rechitsEE).first;
  DetId secondElehottestDetId = getHottestDetId( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->seed()->getHitsByDetId(), rechitsEB,  rechitsEE).first;
  
  double firstEle_E1    = getHottestDetId( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->seed()->getHitsByDetId(), rechitsEB, rechitsEE).second;
  double secondEle_E1   = getHottestDetId( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->seed()->getHitsByDetId(), rechitsEB, rechitsEE).second;
  double firstEle_e5x5  = getE5x5( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->seed()->getHitsByDetId(), rechitsEB, rechitsEE);
  double secondEle_e5x5 = getE5x5( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->seed()->getHitsByDetId(), rechitsEB, rechitsEE);
 
  bool firstElectronIsOnModuleBorder(false);
  bool secondElectronIsOnModuleBorder(false);

  if(class1<100){
    if( firstElehottestDetId.subdetId()==EcalBarrel) firstElectronIsOnModuleBorder = xtalIsOnModuleBorder( firstElehottestDetId  );
    BARREL_ELECTRONS_BEFORE_BORDER_CUT++;
    if( firstElehottestDetId.subdetId()==EcalBarrel && !firstElectronIsOnModuleBorder ) BARREL_ELECTRONS_AFTER_BORDER_CUT++;
  }
  
  if(class2<100){  
    if( secondElehottestDetId.subdetId()==EcalBarrel) secondElectronIsOnModuleBorder = xtalIsOnModuleBorder( secondElehottestDetId  );
    BARREL_ELECTRONS_BEFORE_BORDER_CUT++;
    if( secondElehottestDetId.subdetId()==EcalBarrel && !secondElectronIsOnModuleBorder ) BARREL_ELECTRONS_AFTER_BORDER_CUT++;
  }
  
  if(class1<100){
    if ( firstElehottestDetId.subdetId()==EcalBarrel &&  firstElectronIsOnModuleBorder ) return kContinue;
  }
  
  if(class2<100){
    if ( secondElehottestDetId.subdetId() == EcalBarrel &&  secondElectronIsOnModuleBorder ) return kContinue;
  }
  
  // to decide the ECAL region
  bool selectionBool=false;  
  if(electronSelection_==0)selectionBool=( // no crack
					  zeeCandidates[myBestZ].first->getRecoElectron()->classification() !=  40 && 
					  zeeCandidates[myBestZ].first->getRecoElectron()->classification() != 140 && 
					  zeeCandidates[myBestZ].second->getRecoElectron()->classification()!=  40 && 
					  zeeCandidates[myBestZ].second->getRecoElectron()->classification()!= 140);

  if(electronSelection_==1)selectionBool=( myBestZ != -1 &&  // all barrel electrons but not crack
					   zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && 
					   zeeCandidates[myBestZ].second->getRecoElectron()->classification()< 100 && 
					   zeeCandidates[myBestZ].first->getRecoElectron()->classification()!= 40 &&
					   zeeCandidates[myBestZ].second->getRecoElectron()->classification()!= 40);

  if(electronSelection_==2)selectionBool=( myBestZ != -1 &&  // all endcap electrons but not crack
					   zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 && 
					   zeeCandidates[myBestZ].second->getRecoElectron()->classification()>= 100 && 
					   zeeCandidates[myBestZ].first->getRecoElectron()->classification()!= 140 &&
					   zeeCandidates[myBestZ].second->getRecoElectron()->classification()!= 140);
  
  if(electronSelection_==3)selectionBool=( myBestZ != -1 &&   // only barrel/barrel or endcap/endcap
					   !(zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && 
					     zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100) &&
					   !(zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 &&
					     zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100) );

  if(!selectionBool) return kContinue;
  h1_eventsAfterBorderSelection_->Fill(1);

  // counting
  if( class1<100  && class2<100 )  BBZN++;
  if( class1>=100 && class2>=100 ) EEZN++;
  if((class1<100  && class2>=100) || (class2<100 && class1>=100)) EBZN++;


  


  // ---------------------------------------------------------------
  // possible cut on the j/psi invariant mass
  h1_eventsBeforeInvMassSelection->Fill(1);
  bool invMassBool = useInvMassSelection_ ? ( (mass >= minInvMassCut_) && (mass <= maxInvMassCut_) ) : true;
  if( useInvMassSelection_ && !invMassBool)  return kContinue;  
  h1_eventsAfterInvMassSelection->Fill(1);

  float eta1 = zeeCandidates[myBestZ].first->getRecoElectron()->eta();
  float phi1 = zeeCandidates[myBestZ].first->getRecoElectron()->phi();
  float eta2 = zeeCandidates[myBestZ].second->getRecoElectron()->eta();  
  float phi2 = zeeCandidates[myBestZ].second->getRecoElectron()->phi();  
  
  // taking energy corrections
  float ele1EnergyCorrection(1.);
  float ele2EnergyCorrection(1.);
  if(wantEtaCorrection_){
    ele1EnergyCorrection=1./getEtaCorrection(zeeCandidates[myBestZ].first->getRecoElectron());
    ele2EnergyCorrection=1./getEtaCorrection(zeeCandidates[myBestZ].second->getRecoElectron());
  }




  // several distributions after the full selection
  h2_occupancy_vsEtaPhi ->Fill( eta1 , phi1);
  h2_occupancy_vsEtaPhi ->Fill( eta2 , phi2);
  
  // E1/Esupercluster histos
  if( zeeCandidates[myBestZ].first->getRecoElectron()->classification() < 100 )
    h1_seedOverSC_Barrel_->Fill( firstEle_E1 / zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->rawEnergy() );
  
  if( zeeCandidates[myBestZ].first->getRecoElectron()->classification() >= 100){
    h1_preshowerOverSC_->Fill( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->preshowerEnergy() / zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->rawEnergy() );
    h1_seedOverSC_Endcap_->Fill( firstEle_E1 / zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->energy() );
  }

  if( zeeCandidates[myBestZ].second->getRecoElectron()->classification() < 100 )
    h1_seedOverSC_Barrel_->Fill( secondEle_E1 / zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->rawEnergy() );

  if( zeeCandidates[myBestZ].second->getRecoElectron()->classification() >= 100){
    h1_preshowerOverSC_->Fill( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->preshowerEnergy() / zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->energy() );
    h1_seedOverSC_Endcap_->Fill( secondEle_E1 / zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->rawEnergy() );
  }
    
  if (invMassBool && selectionBool) {
    h1_electronCosTheta_SC_    -> Fill( ZeeKinematicTools::cosThetaElectrons_SC(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) );
    h1_electronCosTheta_TK_    -> Fill( ZeeKinematicTools::cosThetaElectrons_TK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) );
    h1_electronCosTheta_SC_TK_ -> Fill( ZeeKinematicTools::cosThetaElectrons_SC(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection)/ZeeKinematicTools::cosThetaElectrons_TK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) - 1.);
  }
  
  // mc - reco association and comparison plots
  if (!mcProducer_.empty()) {

    std::vector<const reco::PixelMatchGsfElectron*> dauElectronCollection;
    dauElectronCollection.push_back(zeeCandidates[myBestZ].first ->getRecoElectron() );
    dauElectronCollection.push_back(zeeCandidates[myBestZ].second->getRecoElectron() );
   
    std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*> myMCmap;	    
    fillMCmap(&dauElectronCollection,mcEle,myMCmap);
    fillEleInfo(mcEle,myMCmap);
    HepMC::GenParticle* first_mc_ele  = 0;
    HepMC::GenParticle* second_mc_ele = 0;

    for(int i=0; i<mcEle.size(); i++){
      std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>::const_iterator mIter = myMCmap.find(mcEle[i]);
      if (mIter == myMCmap.end()) continue;
      if((*mIter).second) {
	const reco::PixelMatchGsfElectron* myEle=(*mIter).second;
	if(     myEle->superCluster()->energy() == zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->energy())  first_mc_ele = mcEle[i];
	else if(myEle->superCluster()->energy() == zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->energy()) second_mc_ele = mcEle[i];
      }
    }	

    h1_zMassResol_ -> Fill(mass-myGenJPsiMass);
    
    if( first_mc_ele && second_mc_ele){
      
      float theta12_mc = acos( ZeeKinematicTools::cosThetaMCElectrons( first_mc_ele, second_mc_ele )  );
      float theta12_TK = acos( ZeeKinematicTools::cosThetaElectrons_TK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) );
      float theta12_SC = acos(ZeeKinematicTools::cosThetaElectrons_SC(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) );
	    
      float massResolutionTerm_E1 = ( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->energy() - first_mc_ele->momentum().e() ) / first_mc_ele->momentum().e();
      float massResolutionTerm_E2 = ( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->energy() - second_mc_ele->momentum().e() ) / second_mc_ele->momentum().e();
      
      int energyBin_1 = zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->energy() / 20.;
      int energyBin_2 = zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->energy() / 20.;
 
      h_ESCEtrue_vs_ESC -> Fill( zeeCandidates[myBestZ].first->getRecoElectron()->superCluster()->energy() ,  fabs(massResolutionTerm_E1)  );
      h_ESCEtrue_vs_ESC -> Fill( zeeCandidates[myBestZ].second->getRecoElectron()->superCluster()->energy() , fabs(massResolutionTerm_E2)  );      
      if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 ) h_ESCEtrue_array_barrel[energyBin_1]->Fill( massResolutionTerm_E1);
	
      if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100  ) h_ESCEtrue_array_endcap[energyBin_1]->Fill(massResolutionTerm_E1);
      if(zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100  ) h_ESCEtrue_array_barrel[energyBin_2]->Fill( massResolutionTerm_E2);
      if(zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100 ) h_ESCEtrue_array_endcap[energyBin_2]->Fill( massResolutionTerm_E2);
	
      float massResolutionTerm_E     = sqrt( pow(massResolutionTerm_E1,2) + pow(massResolutionTerm_E2,2));
      float massResolutionTerm_angle = ( theta12_TK - theta12_mc )/tan( 0.5 * theta12_mc ) ; 
      float massResolutionTerm_total = sqrt( pow(massResolutionTerm_E1,2) + pow(massResolutionTerm_E2,2) + pow(massResolutionTerm_angle,2) );
      float massResolutionTerm_angle_ratio = massResolutionTerm_angle / massResolutionTerm_total;

      if(zeeCandidates[myBestZ].first->getRecoElectron()->classification() < 100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification() < 100){
	h1_Theta12Resolution_TK_BB_           -> Fill( theta12_TK - theta12_mc);
	h1_Theta12Resolution_SC_BB_           -> Fill( theta12_SC - theta12_mc);	      
	h1_Theta12_TK_BB_                     -> Fill( theta12_TK, theta12_mc  );
	h1_Theta12_SC_BB_                     -> Fill( theta12_SC, theta12_mc );
	h1_massResolutionTerm_total_BB_       -> Fill(massResolutionTerm_total);
	h1_massResolutionTerm_angle_BB_       -> Fill(massResolutionTerm_angle);
	h1_massResolutionTerm_energy_BB_      -> Fill(massResolutionTerm_E);	
	h1_massResolutionTerm_angle_ratio_BB_ -> Fill( massResolutionTerm_angle_ratio );
      }
      
      if(zeeCandidates[myBestZ].first->getRecoElectron()->classification() >= 100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification() >= 100){
	h1_Theta12_TK_EE_                     -> Fill( theta12_TK, theta12_mc  );
	h1_Theta12_SC_EE_                     -> Fill( theta12_SC, theta12_mc );
	h1_Theta12Resolution_TK_EE_           -> Fill( theta12_TK - theta12_mc);
	h1_Theta12Resolution_SC_EE_           -> Fill( theta12_SC - theta12_mc);
	h1_massResolutionTerm_total_EE_       -> Fill(massResolutionTerm_total);
	h1_massResolutionTerm_angle_EE_       -> Fill(massResolutionTerm_angle);
	h1_massResolutionTerm_energy_EE_      -> Fill(massResolutionTerm_E);
	h1_massResolutionTerm_angle_ratio_EE_ -> Fill( massResolutionTerm_angle_ratio );
      }
    }
  }


  // put f(eta) in the algorithm
  theAlgorithm_->getEventWeight(1);
  
  double E5x5_correctionFactor = 1.;
  if( ZCalib_InvMass_ == "E5x5TRMass" || ZCalib_InvMass_ == "E5x5Mass"  )  // with or wo tracker, correct to move to 5x5 instead of supercluster   
    E5x5_correctionFactor = sqrt(zeeCandidates[myBestZ].first->getRecoElectron()->caloEnergy() * zeeCandidates[myBestZ].second->getRecoElectron()->caloEnergy() / (firstEle_e5x5*secondEle_e5x5) );
  
  theAlgorithm_->addEvent(zeeCandidates[myBestZ].first, zeeCandidates[myBestZ].second, 
			  resonanceMass_ * E5x5_correctionFactor * sqrt(ele1EnergyCorrection*ele2EnergyCorrection) );  // here we add this Z to the algo
  
  
  // final distributions, Z mass
  h1_reco_ZMass_     -> Fill(ZeeKinematicTools::calculateZMass_withTK(zeeCandidates[myBestZ]));
  h1_reco_ZMassCorr_ -> Fill(ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection));

  if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100 ) {
    h1_reco_ZMassBB_     -> Fill(ZeeKinematicTools::calculateZMass_withTK(zeeCandidates[myBestZ])); 
    h1_reco_ZMassCorrBB_ -> Fill(ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection));
  }  

  if( ( zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100 ) ||
      ( zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100 ) ) {
    h1_reco_ZMassEB_     -> Fill(ZeeKinematicTools::calculateZMass_withTK(zeeCandidates[myBestZ])); 
    h1_reco_ZMassCorrEB_ -> Fill(ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection));
  }

  if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100 ) {
    h1_reco_ZMassEE_     -> Fill(ZeeKinematicTools::calculateZMass_withTK(zeeCandidates[myBestZ])); 
    h1_reco_ZMassCorrEE_->Fill(ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection));
  }


  
  // tree filling
  massCorr4tree     = ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection);
  massCorrDiff4tree = ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection) - myGenJPsiMass;

  if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100 ) {
    massCorr4treeBB = ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection);
    massCorr4treeEB = -100;
    massCorr4treeEE = -100;
  }

  if( ( zeeCandidates[myBestZ].first->getRecoElectron()->classification()<100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100 ) ||
      ( zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()<100 ) ) {
    massCorr4treeEB = ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection);
    massCorr4treeBB = -100;
    massCorr4treeEE = -100;
  }

  if(zeeCandidates[myBestZ].first->getRecoElectron()->classification()>=100 && zeeCandidates[myBestZ].second->getRecoElectron()->classification()>=100 ) {
    massCorr4treeEE = ZeeKinematicTools::calculateZMassWithCorrectedElectrons_withTK(zeeCandidates[myBestZ],ele1EnergyCorrection,ele2EnergyCorrection);
    massCorr4treeBB = -100;
    massCorr4treeEB = -100;
  }

  myTree->Fill();
  
  if (loopFlag_ == 0){
    myJPsiPlots_ -> fillEleInfoAS( zeeCandidates[myBestZ].first->getRecoElectron() );
    myJPsiPlots_ -> fillEleInfoAS( zeeCandidates[myBestZ].second->getRecoElectron() );
    myJPsiPlots_ -> fillJPsiInfo ( zeeCandidates[myBestZ] );    
  }
  
  return kContinue;
}


// called at beginning of loop
void JPsiCalibration::startingNewLoop ( unsigned int iLoop ) {
  
  theAlgorithm_->resetIteration();
  resetVariables();
  resetHistograms(); 
}


// called at end of loop
edm::EDLooper::Status JPsiCalibration::endOfLoop(const edm::EventSetup& iSetup, unsigned int iLoop) {

  cout << endl;  
  cout << "fine, eleID: " << endl;
  cout << "passedDeta = " << passedDeta << endl;
  cout << "passedDphi = " << passedDphi << endl;
  cout << "passedHoE  = " << passedHoE  << endl;
  cout << "passedSee  = " << passedSee  << endl;
  cout << endl;

  double par[3];
  double errpar[3];
  double zChi2;
  int zIters;
  
  // iterative gaussian fit to the reconstructed mass
  ZIterativeAlgorithmWithFit::gausfit(h1_reco_ZMass_,par,errpar,2.,2., &zChi2, &zIters );
  h2_zMassVsLoop_     -> Fill(loopFlag_,  par[1] );
  h2_zMassDiffVsLoop_ -> Fill(loopFlag_, (par[1]-resonanceMass_)/ resonanceMass_ );
  h2_zWidthVsLoop_    -> Fill(loopFlag_, par[2] );
  
  // j/psi peak fit
  theAlgorithm_->iterate();   // run the algorithm
  
  const std::vector<float>& optimizedCoefficients      = theAlgorithm_->getOptimizedCoefficients();
  const std::vector<float>& optimizedCoefficientsError = theAlgorithm_->getOptimizedCoefficientsError();
  const std::vector<float>& optimizedChi2              = theAlgorithm_->getOptimizedChiSquare();
  const std::vector<int>& optimizedIterations          = theAlgorithm_->getOptimizedIterations();

  // new calib coefficients
  for (unsigned int ieta=0;ieta<optimizedCoefficients.size();ieta++) {
    
    NewCalibCoeff[ieta] = calibCoeff[ieta] * optimizedCoefficients[ieta];
    h2_chi2_[loopFlag_]       -> Fill( ringNumberCorrector( ieta ), optimizedChi2[ieta] );
    h2_iterations_[loopFlag_] -> Fill( ringNumberCorrector( ieta ), optimizedIterations[ieta] );
  }
  
  coefficientDistanceAtIteration[loopFlag_]= computeCoefficientDistanceAtIteration(calibCoeff, NewCalibCoeff, optimizedCoefficients.size() );
  for (unsigned int ieta=0;ieta<optimizedCoefficients.size();ieta++){ 
    calibCoeff[ieta] *= optimizedCoefficients[ieta];
    calibCoeffError[ieta] = calibCoeff[ieta] * sqrt ( pow( optimizedCoefficientsError[ieta]/optimizedCoefficients[ieta], 2 ) + pow( calibCoeffError[ieta]/calibCoeff[ieta] , 2 )  );

    std::vector<DetId> ringIds;
    if(calibMode_ == "RING")   ringIds = EcalRingCalibrationTools::getDetIdsInRing(ieta);
    if(calibMode_ == "MODULE") ringIds = EcalRingCalibrationTools::getDetIdsInModule(ieta);
    if(calibMode_ == "ABS_SCALE" || calibMode_ == "ETA_ET_MODE" ) ringIds = EcalRingCalibrationTools::getDetIdsInECAL();
    
    for (unsigned int iid=0; iid<ringIds.size();++iid){
      
      if(ringIds[iid].subdetId() == EcalBarrel){
	EBDetId myEBDetId(ringIds[iid]);  
	h2_xtalRecalibCoeffBarrel_[loopFlag_]->SetBinContent( myEBDetId.ieta() + 86, myEBDetId.iphi(), 100 * (calibCoeff[ieta]*initCalibCoeff[ieta] - 1.) );
      }
      
      if(ringIds[iid].subdetId() == EcalEndcap){
	EEDetId myEEDetId(ringIds[iid]);
	if(myEEDetId.zside() < 0)
	  h2_xtalRecalibCoeffEndcapMinus_[loopFlag_]->SetBinContent( myEEDetId.ix(), myEEDetId.iy(), 100 * (calibCoeff[ieta]*initCalibCoeff[ieta] - 1.) );
	if(myEEDetId.zside() > 0)
	  h2_xtalRecalibCoeffEndcapPlus_[loopFlag_]->SetBinContent( myEEDetId.ix(), myEEDetId.iy(), 100 * (calibCoeff[ieta]*initCalibCoeff[ieta] - 1.) );
      }
      
      ical->setValue( ringIds[iid], *(ical->getMap().find(ringIds[iid])  ) * optimizedCoefficients[ieta] );
    }    
  }
  
  
  // dump residual miscalibration at each loop
  for( int k = 0; k<theAlgorithm_->getNumberOfChannels(); k++ ) {

    bool isNearCrack = ( abs( ringNumberCorrector(k) ) == 1 || abs( ringNumberCorrector(k) ) == 25 ||
			 abs( ringNumberCorrector(k) ) == 26 || abs( ringNumberCorrector(k) ) == 45 ||
			 abs( ringNumberCorrector(k) ) == 46 || abs( ringNumberCorrector(k) ) == 65 ||
			 abs( ringNumberCorrector(k) ) == 66 || abs( ringNumberCorrector(k) ) == 85 ||
			 abs( ringNumberCorrector(k) ) == 86 || abs( ringNumberCorrector(k) ) == 124
			 );    
    if(!isNearCrack){
      h1_mcParz_[iLoop]->Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
      if(k<170)  h1_mcEBParz_[iLoop]->Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
      if(k>=170) h1_mcEEParz_[iLoop]->Fill( initCalibCoeff[k]*calibCoeff[k] -1. );
    }
  }
  
  
  double parResidual[3];
  double errparResidual[3];
  double zResChi2;
  int zResIters;
  
  ZIterativeAlgorithmWithFit::gausfit(h1_mcParz_[iLoop],parResidual,errparResidual,3.,3., &zResChi2, &zResIters);
  h2_residualSigma_          -> Fill(loopFlag_ + 1,  parResidual[2]);
  loopArray[loopFlag_]       = loopFlag_ + 1;
  sigmaArray[loopFlag_]      = parResidual[2];
  sigmaErrorArray[loopFlag_] = errparResidual[2];

  outputFile_->cd();
  h1_mcParz_[iLoop]   -> Write();
  h1_mcEBParz_[iLoop] -> Write();
  h1_mcEEParz_[iLoop] -> Write();
  h2_xtalRecalibCoeffBarrel_[loopFlag_]      -> Write();
  h2_xtalRecalibCoeffEndcapPlus_[loopFlag_]  -> Write();
  h2_xtalRecalibCoeffEndcapMinus_[loopFlag_] -> Write();
  
  loopFlag_++;  
  if ( iLoop == theMaxLoops-1 || iLoop >= theMaxLoops ) return kStop;
  else return kContinue;
}


void JPsiCalibration::bookHistograms() {

  // to count
  h1_eventsBeforeHLTSelection_     =  new TH1F("h1_eventsBeforeHLTSelection",     "h1_eventsBeforeHLTSelection",     5,0,5); 
  h1_eventsAfterHLTSelection_      =  new TH1F("h1_eventsAfterHLTSelection",      "h1_eventsAfterHLTSelection",      5,0,5);
  h1_eventsBefore2GsfSelection_    =  new TH1F("h1_eventsBefore2GsfSelection",    "h1_eventsBefore2GsfSelection",    5,0,5); 
  h1_eventsAfter2GsfSelection_     =  new TH1F("h1_eventsAfter2GsfSelection",     "h1_eventsAfter2GsfSelection",     5,0,5);
  h1_eventsBeforeEtSelection_      =  new TH1F("h1_eventsBeforeEtSelection",      "h1_eventsBeforeEtSelection",      5,0,5); 
  h1_eventsAfterEtSelection_       =  new TH1F("h1_eventsAfterEtSelection",       "h1_eventsAfterEtSelection",       5,0,5); 
  h1_eventsBeforeEleIDSelection_   =  new TH1F("h1_eventsBeforeEleIDSelection",   "h1_eventsBeforeEleIDSelection",   5,0,5); 
  h1_eventsAfterEleIDSelection_    =  new TH1F("h1_eventsAfterEleIDSelection",    "h1_eventsAfterEleIDSelection",    5,0,5); 
  h1_eventsBeforeEleIsolSelection_ =  new TH1F("h1_eventsBeforeEleIsolSelection", "h1_eventsBeforeEleIsolSelection", 5,0,5); 
  h1_eventsAfterEleIsolSelection_  =  new TH1F("h1_eventsAfterEleIsolSelection",  "h1_eventsAfterEleIsolSelection",  5,0,5); 
  h1_eventsBeforeTrackerIsolSelection_ =  new TH1F("h1_eventsBeforeTrackerIsolSelection", "h1_eventsBeforeTrackerIsolSelection", 5,0,5); 
  h1_eventsAfterTrackerIsolSelection_  =  new TH1F("h1_eventsAfterTrackerIsolSelection",  "h1_eventsAfterTrackerIsolSelection",  5,0,5); 
  h1_eventsBeforeBorderSelection_  =  new TH1F("h1_eventsBeforeBorderSelection",  "h1_eventsBeforeBorderSelection",  5,0,5); 
  h1_eventsAfterBorderSelection_   =  new TH1F("h1_eventsAfterBorderSelection",   "h1_eventsAfterBorderSelection",   5,0,5);
  h1_eventsBeforeInvMassSelection  =  new TH1F("h1_eventsBeforeInvMassSelection", "h1_eventsBeforeInvMassSelection", 5,0,5);
  h1_eventsAfterInvMassSelection   =  new TH1F("h1_eventsAfterInvMassSelection",  "h1_eventsAfterInvMassSelection",  5,0,5);
  
  // z candidates
  h1_ZCandMult_ = new TH1F("ZCandMult","Multiplicity of Z candidates in one event",10,0.,10.);
  h1_ZCandMult_ -> SetXTitle("ZCandMult");

  // kinematics
  h2_occupancy_vsEtaPhi      = new TH2F("h2_occupancy_vsEtaPhi",  "h2_occupancy_vsEtaPhi",  200, -3, 3, 200, -3.15, 3.15);  
  h1_seedOverSC_Barrel_      = new TH1F("h1_seedOverSC_Barrel",   "h1_seedOverSC_Barrel",    50,  0.,2);
  h1_seedOverSC_Endcap_      = new TH1F("h1_seedOverSC_Endcap",   "h1_seedOverSC_Endcap",    50,  0.,2);
  h1_preshowerOverSC_        = new TH1F("h1_preshowerOverSC",     "h1_preshowerOverSC",     400,  0.,1.);  
  h1_electronCosTheta_TK_    = new TH1F("electronCosTheta_TK",    "electronCosTheta_TK",    100, -1., 1.);
  h1_electronCosTheta_SC_    = new TH1F("electronCosTheta_SC",    "electronCosTheta_SC",    100, -1, 1);
  h1_electronCosTheta_SC_TK_ = new TH1F("electronCosTheta_SC_TK", "electronCosTheta_SC_TK", 200, -0.1, 0.1);

  h1_electronCosTheta_TK_    -> SetXTitle("cos #theta_{12}");  
  h1_electronCosTheta_SC_    -> SetXTitle("cos #theta_{12}");
  h1_electronCosTheta_SC_TK_ -> SetXTitle("cos #theta_{12}^{SC}/ cos #theta_{12}^{TK} - 1");
  
  // comparison with MC
  h_ESCEtrue_vs_ESC          = new TH2F("h2_ESCEtrue_vs_ESC",        "h2_ESCEtrue_vs_ESC",        10, 0., 200., 200, 0., 0.1); 

  h1_Theta12Resolution_TK_BB_ = new TH1F("Theta12Resolution_TK_BB",   "Theta12Resolution_TK_BB",  200, -0.01, 0.01);
  h1_Theta12Resolution_TK_EE_ = new TH1F("Theta12Resolution_TK_EE",   "Theta12Resolution_TK_EE",  200, -0.01, 0.01);
  h1_Theta12Resolution_SC_BB_ = new TH1F("Theta12Resolution_SC_BB",   "Theta12Resolution_SC_BB",  200, -0.10, 0.10);
  h1_Theta12Resolution_SC_EE_ = new TH1F("Theta12Resolution_SC_EE",   "Theta12Resolution_SC_EE",  200, -0.10, 0.10);

  h1_Theta12_TK_BB_           = new TH1F("Theta12_TK_BB",    "Theta12_TK_BB",    200, -1.5, 1.5);
  h1_Theta12_SC_BB_           = new TH1F("Theta12_SC_BB",    "Theta12_SC_BB",    200, -1.5, 1.5);
  h1_Theta12_TK_EE_           = new TH1F("Theta12_TK_EE",    "Theta12_TK_EE",    200, -1.5, 1.5);
  h1_Theta12_SC_EE_           = new TH1F("Theta12_SC_EE",    "Theta12_SC_EE",    200, -1.5, 1.5);

  h1_zMassResol_              = new TH1F("zMassResol",  "zMassResol",  100, -5., 5.);
  h1_eleEtaResol_             = new TH1F("eleEtaResol", "eleEtaResol", 100, -0.01, 0.01);
  h1_elePhiResol_             = new TH1F("elePhiResol", "elePhiResol", 100, -0.01, 0.01);

  h1_massResolutionTerm_total_BB_       = new TH1F("h1_massResolutionTerm_total_BB",      "h1_massResolutionTerm_total_BB",      200,  0.,  0.5);
  h1_massResolutionTerm_total_EE_       = new TH1F("h1_massResolutionTerm_total_EE",      "h1_massResolutionTerm_total_EE",      200,  0.,  0.5);
  h1_massResolutionTerm_energy_BB_      = new TH1F("h1_massResolutionTerm_energy_BB",     "h1_massResolutionTerm_energy_BB",     200,  0.,  0.5);
  h1_massResolutionTerm_energy_EE_      = new TH1F("h1_massResolutionTerm_energy_EE",     "h1_massResolutionTerm_energy_EE",     200,  0.,  0.5);
  h1_massResolutionTerm_angle_BB_       = new TH1F("h1_massResolutionTerm_angle_BB",      "h1_massResolutionTerm_angle_BB",      200, -0.5, 0.5);
  h1_massResolutionTerm_angle_EE_       = new TH1F("h1_massResolutionTerm_angle_EE",      "h1_massResolutionTerm_angle_EE",      200, -0.5, 0.5);
  h1_massResolutionTerm_angle_ratio_BB_ = new TH1F("h1_massResolutionTerm_angle_ratio_BB","h1_massResolutionTerm_angle_ratio_BB",100, -1.,  1.);
  h1_massResolutionTerm_angle_ratio_EE_ = new TH1F("h1_massResolutionTerm_angle_ratio_EE","h1_massResolutionTerm_angle_ratio_EE",100, -1.,  1.);
  
  h1_Theta12Resolution_TK_BB_ -> SetXTitle("#theta_{12}^{TK} - #theta_{12}^{MC}"); 
  h1_Theta12Resolution_SC_BB_ -> SetXTitle("#theta_{12}^{SC} - #theta_{12}^{MC}");   
  h1_Theta12Resolution_TK_EE_ -> SetXTitle("#theta_{12}^{TK} - #theta_{12}^{MC}"); 
  h1_Theta12Resolution_SC_EE_ -> SetXTitle("#theta_{12}^{SC} - #theta_{12}^{MC}"); 
  h1_Theta12_TK_BB_ -> SetXTitle("#theta_{12}^{TK}"); 
  h1_Theta12_TK_BB_ -> SetYTitle("#theta_{12}^{MC}");
  h1_Theta12_SC_BB_ -> SetXTitle("#theta_{12}^{SC}"); 
  h1_Theta12_SC_BB_ -> SetYTitle("#theta_{12}^{MC}");
  h1_Theta12_TK_EE_ -> SetXTitle("#theta_{12}^{TK}"); 
  h1_Theta12_TK_EE_ -> SetYTitle("#theta_{12}^{MC}");
  h1_Theta12_SC_EE_ -> SetXTitle("#theta_{12}^{SC}"); 
  h1_Theta12_SC_EE_ -> SetYTitle("#theta_{12}^{MC}");
  h1_zMassResol_    -> SetXTitle("M_{Z, reco} - M_{Z, MC}");
  h1_eleEtaResol_   -> SetXTitle("#eta_{reco} - #eta_{MC}");
  h1_elePhiResol_   -> SetXTitle("#phi_{reco} - #phi_{MC}");

  // Z mass
  h1_reco_ZMass_       = new TH1F("reco_ZMass",      "Inv. mass of 2 reco Electrons",           80,0.,5.);
  h1_reco_ZMassBB_     = new TH1F("reco_ZMassBB",    "Inv. mass of 2 reco Electrons",           80,0.,5.);
  h1_reco_ZMassEB_     = new TH1F("reco_ZMassEB",    "Inv. mass of 2 reco Electrons",           80,0.,5.);
  h1_reco_ZMassEE_     = new TH1F("reco_ZMassEE",    "Inv. mass of 2 reco Electrons",           80,0.,5.);
  h1_reco_ZMassCorr_   = new TH1F("reco_ZMassCorr",  "Inv. mass of 2 corrected reco Electrons", 80,0.,5.);
  h1_reco_ZMassCorrBB_ = new TH1F("reco_ZMassCorrBB","Inv. mass of 2 corrected reco Electrons", 80,0.,5.);
  h1_reco_ZMassCorrEB_ = new TH1F("reco_ZMassCorrEB","Inv. mass of 2 corrected reco Electrons", 80,0.,5.);
  h1_reco_ZMassCorrEE_ = new TH1F("reco_ZMassCorrEE","Inv. mass of 2 corrected reco Electrons", 80,0.,5.);
  h1_reco_ZMass_       -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassBB_     -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassEB_     -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassEE_     -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassCorr_   -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassCorrBB_ -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassCorrEB_ -> SetXTitle("Mee (GeV/c^{2})");
  h1_reco_ZMassCorrEE_ -> SetXTitle("Mee (GeV/c^{2})");

  // to debug
  h1_sEEeb_before_ = new TH1F("h1_sEEeb_before_","sigma_EtaEta of SC, barrel", 100, 0., 0.05);
  h1_sEEee_before_ = new TH1F("h1_sEEee_before_","sigma_EtaEta of SC, endcap", 100, 0., 0.10);
  h1_sEEeb_after_  = new TH1F("h1_sEEeb_after_", "sigma_EtaEta of SC, barrel", 100, 0., 0.05);
  h1_sEEee_after_  = new TH1F("h1_sEEee_after_", "sigma_EtaEta of SC, endcap", 100, 0., 0.10);

  // coefficients
  h2_xtalMiscalibCoeffBarrel_      = new TH2F("h2_xtalMiscalibCoeffBarrel","h2_xtalMiscalibCoeffBarrel", 171, -85, 85, 360, 0, 360);
  h2_xtalMiscalibCoeffEndcapMinus_ = new TH2F("h2_xtalMiscalibCoeffEndcapMinus", "h2_xtalMiscalibCoeffEndcapMinus", 100, 0,100, 100, 0, 100);
  h2_xtalMiscalibCoeffEndcapPlus_  = new TH2F("h2_xtalMiscalibCoeffEndcapPlus", "h2_xtalMiscalibCoeffEndcapPlus", 100, 0,100, 100, 0, 100);
  h2_xtalMiscalibCoeffBarrel_      -> SetXTitle("ieta");
  h2_xtalMiscalibCoeffBarrel_      -> SetYTitle("iphi");
  h2_xtalMiscalibCoeffEndcapMinus_ -> SetXTitle("ix");
  h2_xtalMiscalibCoeffEndcapMinus_ -> SetYTitle("iy");
  h2_xtalMiscalibCoeffEndcapPlus_  -> SetXTitle("ix");
  h2_xtalMiscalibCoeffEndcapPlus_  -> SetYTitle("iy");

  // ele infos
  for (int i=0;i<25;i++) { 
    
    char histoName[50];

    sprintf(histoName,"h_ESCEtrueVsEta_%d",i);
    h_ESCEtrueVsEta_[i] = new TH2F(histoName,histoName, 150, 0., 2.7, 300,0.,1.5);
    
    sprintf(histoName,"h_ESCEtrue_%d",i);    
    h_ESCEtrue_[i] = new TH1F(histoName,histoName, 300,0.,1.5);
    
    sprintf(histoName,"h_PTres_%d",i);
    h_PTres_[i] = new TH1F(histoName,histoName, 300,0.,1.5);
    
    sprintf(histoName,"h_PTres_Endcap_%d",i);
    h_PTres_Endcap_[i] = new TH1F(histoName,histoName, 300,0.,1.5);
    
    sprintf(histoName,"h2_chi2_%d",i);
    h2_chi2_[i] = new TH2F(histoName,histoName, 1000,-150,150, 1000, -1, 5);
    
    sprintf(histoName,"h2_iterations_%d",i);
    h2_iterations_[i] = new TH2F(histoName,histoName, 1000,-150,150, 1000, -1, 15);
    
    sprintf(histoName,"h_ESCcorrEtrueVsEta_%d",i);    
    h_ESCcorrEtrueVsEta_[i] = new TH2F(histoName,histoName, 150, 0., 2.7, 300,0.,1.5);
    
    sprintf(histoName,"h_ESCcorrEtrue_%d",i);    
    h_ESCcorrEtrue_[i] = new TH1F(histoName,histoName, 300,0.,1.5);
    
    sprintf(histoName,"h2_xtalRecalibCoeffBarrel_%d",i);
    h2_xtalRecalibCoeffBarrel_[i] = new TH2F(histoName,histoName, 171, -85, 85, 360, 0, 360);    
    h2_xtalRecalibCoeffBarrel_[i] -> SetXTitle("ieta");
    h2_xtalRecalibCoeffBarrel_[i] -> SetYTitle("iphi");
    
    sprintf(histoName,"h2_xtalRecalibCoeffEndcapMinus_%d",i);
    h2_xtalRecalibCoeffEndcapMinus_[i] = new TH2F(histoName,histoName, 100, 0,100, 100, 0, 100);
    h2_xtalRecalibCoeffEndcapMinus_[i]->SetXTitle("ix");
    h2_xtalRecalibCoeffEndcapMinus_[i]->SetYTitle("iy");
    
    sprintf(histoName,"h2_xtalRecalibCoeffEndcapPlus_%d",i);
    h2_xtalRecalibCoeffEndcapPlus_[i] = new TH2F(histoName,histoName, 100, 0,100, 100, 0, 100);
    h2_xtalRecalibCoeffEndcapPlus_[i]->SetXTitle("ix");
    h2_xtalRecalibCoeffEndcapPlus_[i]->SetYTitle("iy");
    
    sprintf(histoName,"h_ESCEtrue_array_barrel_Ebin_%d",i);
    h_ESCEtrue_array_barrel[i] =  new TH1F(histoName,histoName, 100, -0.2, 0.2);
    
    sprintf(histoName,"h_ESCEtrue_array_endcap_Ebin_%d",i);
    h_ESCEtrue_array_endcap[i] =  new TH1F(histoName,histoName, 100, -0.2, 0.2);
  }                         

  // efficiencies
  for (int i=0;i<2;i++) { 
    char histoName[50];
    sprintf(histoName,"h_eleEffEta_%d",i);
    h_eleEffEta_[i] = new TH1F(histoName,histoName, 150, 0., 2.7);
    h_eleEffEta_[i] -> SetXTitle("|#eta|");    

    sprintf(histoName,"h_eleEffPhi_%d",i);
    h_eleEffPhi_[i] = new TH1F(histoName,histoName, 400, -4., 4.);
    h_eleEffPhi_[i] -> SetXTitle("Phi");
    
    sprintf(histoName,"h_eleEffPt_%d",i);
    h_eleEffPt_[i] = new TH1F(histoName,histoName, 200, 0., 200.);
    h_eleEffPt_[i]->SetXTitle("p_{T}(GeV/c)");
  }

  h1_efficiencySummary_ =  new TH1F("h1_efficiencySummary", "h1_efficiencySummary", 12, 0, 12 );
  h1_efficiencySummary_ ->SetYTitle("Efficiency");

  // algo results histos
  h1_occupancyVsEta_ = new TH1F("occupancyVsEta","occupancyVsEta",249, -124.5, 124.5);
  h1_occupancyVsEta_ -> SetYTitle("Weighted electron statistics");
  h1_occupancyVsEta_ -> SetXTitle("Eta channel");

  h1_weightSumMeanBarrel_ = new TH1F("weightSumMeanBarrel","weightSumMeanBarrel",10000, 0, 10000);
  h1_weightSumMeanEndcap_ = new TH1F("weightSumMeanEndcap","weightSumMeanEndcap",10000, 0, 10000);
  
  h1_occupancy_       = new TH1F("occupancy","occupancy",1000,0,10000);
  h1_occupancyBarrel_ = new TH1F("occupancyBarrel","occupancyBarrel",1000,0,10000);
  h1_occupancyEndcap_ = new TH1F("occupancyEndcap","occupancyEndcap",1000,0,10000);
  h1_occupancy_       -> SetXTitle("Weighted electron statistics");
  h1_occupancyBarrel_ -> SetXTitle("Weighted electron statistics");
  h1_occupancyEndcap_ -> SetXTitle("Weighted electron statistics");

  h2_coeffVsEta_        = new TH2F("h2_calibCoeffVsEta","       h2_calibCoeffVsEta",        249,-124,125, 200, 0., 2.);
  h2_coeffVsEtaGrouped_ = new TH2F("h2_calibCoeffVsEtaGrouped","h2_calibCoeffVsEtaGrouped", 200, 0.,  3., 200, 0.6, 1.4);
  h2_zMassVsLoop_       = new TH2F("h2_zMassVsLoop",           "h2_zMassVsLoop",           1000, 0,  40,   80, 0.,10.);
  h2_zMassDiffVsLoop_   = new TH2F("h2_zMassDiffVsLoop",       "h2_zMassDiffVsLoop",       1000, 0,  40,  100, -1., 1.);
  h2_zWidthVsLoop_      = new TH2F("h2_zWidthVsLoop",          "h2_zWidthVsLoop",          1000, 0,  40,  100, 0.,10.);
  h2_coeffVsLoop_       = new TH2F("h2_coeffVsLoop",           "h2_coeffVsLoop",           1000, 0,  40,  100, 0., 2.);
  h2_residualSigma_     = new TH2F("h2_residualSigma",         "h2_residualSigma",         1000, 0,  40,  100, 0., .5);
  h2_miscalRecal_       = new TH2F("h2_miscalRecal",           "h2_miscalRecal",            500, 0.,  2., 500, 0., 2.);
  h2_miscalRecalEB_     = new TH2F("h2_miscalRecalEB",         "h2_miscalRecalEB",          500, 0., 2., 500, 0., 2.);
  h2_miscalRecalEE_     = new TH2F("h2_miscalRecalEE",         "h2_miscalRecalEE",          500, 0., 2., 500, 0., 2.);
  h2_coeffVsEta_        -> SetXTitle("Eta channel");
  h2_coeffVsEta_        -> SetYTitle("recalibCoeff");
  h2_coeffVsEtaGrouped_ -> SetXTitle("|#eta|");
  h2_coeffVsEtaGrouped_ -> SetYTitle("recalibCoeff");
  h2_zMassDiffVsLoop_   -> SetXTitle("Iteration");
  h2_zMassDiffVsLoop_   -> SetYTitle("M_{Z, reco peak} - M_{Z, true}");
  h2_miscalRecal_       -> SetXTitle("initCalibCoeff");
  h2_miscalRecal_       -> SetYTitle("1/RecalibCoeff");
  h2_miscalRecalEB_     -> SetXTitle("initCalibCoeff");
  h2_miscalRecalEB_     -> SetYTitle("1/RecalibCoeff");
  h2_miscalRecalEE_     -> SetXTitle("initCalibCoeff");
  h2_miscalRecalEE_     -> SetYTitle("1/RecalibCoeff");

  // miscalibration
  h1_mc_   = new TH1F("h1_residualMiscalib",  "h1_residualMiscalib",   200, -0.2, 0.2);
  h1_mcEB_ = new TH1F("h1_residualMiscalibEB","h1_residualMiscalibEB", 200, -0.2, 0.2);
  h1_mcEE_ = new TH1F("h1_residualMiscalibEE","h1_residualMiscalibEE", 200, -0.2, 0.2);
  for (int i=0;i<25;i++){ 
    char histoName[50];
    sprintf(histoName,"h1_residualMiscalibParz_%d",i);
    h1_mcParz_[i]   = new TH1F(histoName,histoName, 200, -0.2, 0.2);
    sprintf(histoName,"h1_residualMiscalibEBParz_%d",i);
    h1_mcEBParz_[i] = new TH1F(histoName,histoName, 200, -0.2, 0.2);
    sprintf(histoName,"h1_residualMiscalibEEParz_%d",i);
    h1_mcEEParz_[i] = new TH1F(histoName,histoName, 200, -0.2, 0.2);
  }

  // others
  myJPsiPlots_ ->bookJPsiMCHistograms();  
  myJPsiPlots_ ->bookJPsiHistograms();  
  myJPsiPlots_ ->bookEleMCHistograms();	  
  myJPsiPlots_ ->bookEleHistogramsAS();	  
  myJPsiPlots_ ->bookEleHistograms();		
}

double JPsiCalibration::fEtaBarrelBad(double scEta) const{
  
  double p0 = P0_;
  double p1 = P1_;
  double p2 = P2_;
  
  double x  = (double) fabs(scEta);
  
  return  ( p0 + p1*x*x + p2*x*x*x*x );  
}
  
double JPsiCalibration::fEtaEndcapGood(double scEta) const{

  double p0 = P0e_;
  double p1 = P1e_;
  double p2 = P2e_;

  double x  = (double) fabs(scEta);

  return  ( p0 + p1*x*x + p2*x*x*x*x );  
}

double JPsiCalibration::fEtaEndcapBad(double scEta) const{
  
  double p0 = P0e_;
  double p1 = P1e_;
  double p2 = P2e_;

  double x  = (double) fabs(scEta);

 return  ( p0 + p1*x*x + p2*x*x*x*x );  
}
  
double JPsiCalibration::fEtaBarrelGood(double scEta) const{

  double p0 = P0_;
  double p1 = P1_;
  double p2 = P2_;

  double x  = (double) fabs(scEta);

 return  ( p0 + p1*x*x + p2*x*x*x*x );  
}

void JPsiCalibration::fillMCmap(const std::vector<const reco::PixelMatchGsfElectron*>* electronCollection,const std::vector<HepMC::GenParticle*>& mcEle,std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>& myMCmap) {

  for (unsigned int i=0;i<mcEle.size();i++) {

    float minDR=0.1;
    const reco::PixelMatchGsfElectron* myMatchEle=0;
    for (unsigned int j=0;j<electronCollection->size();j++) {

      float dr=EvalDR(mcEle[i]->momentum().pseudoRapidity(),(*(*electronCollection)[j]).eta(),mcEle[i]->momentum().phi(),(*(*electronCollection)[j]).phi());
      if (dr < minDR ) {
	myMatchEle = (*electronCollection)[j];
	minDR = dr;
      }}
    myMCmap.insert(std::pair<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>(mcEle[i],myMatchEle));  
  }
}
                                                                                                                             
float JPsiCalibration::EvalDR(float Eta,float Eta_ref,float Phi,float Phi_ref) {

  if (Phi<0) Phi = 2*TMath::Pi() + Phi;
  if (Phi_ref<0) Phi_ref = 2*TMath::Pi() + Phi_ref;
  float DPhi = Phi - Phi_ref ;
  if (fabs(DPhi)>TMath::Pi()) DPhi = 2*TMath::Pi() - fabs(DPhi);
  
  float DEta = Eta - Eta_ref ;
  float DR = sqrt( DEta*DEta + DPhi*DPhi );
  return DR;
}

float JPsiCalibration::EvalDPhi(float Phi,float Phi_ref){ 

  if (Phi<0) Phi = 2*TMath::Pi() + Phi;
  if (Phi_ref<0) Phi_ref = 2*TMath::Pi() + Phi_ref;
  return (Phi - Phi_ref);
}

void JPsiCalibration::fillEleInfo( std::vector<HepMC::GenParticle*>& mcEle, std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>& associationMap) {
  
  for (unsigned int i=0;i<mcEle.size();i++){

    h_eleEffEta_[0]->Fill(fabs(mcEle[i]->momentum().pseudoRapidity()));
    h_eleEffPhi_[0]->Fill(mcEle[i]->momentum().phi());
    h_eleEffPt_[0] ->Fill(mcEle[i]->momentum().perp());
    
    std::map<HepMC::GenParticle*,const reco::PixelMatchGsfElectron*>::const_iterator mIter = associationMap.find(mcEle[i]);
    if (mIter == associationMap.end() ) continue;
      
    if((*mIter).second) {

      const reco::PixelMatchGsfElectron* myEle=(*mIter).second;      
      h_eleEffEta_[1] -> Fill(fabs(mcEle[i]->momentum().pseudoRapidity()));
      h_eleEffPhi_[1] -> Fill(mcEle[i]->momentum().phi());
      h_eleEffPt_[1]  -> Fill(mcEle[i]->momentum().perp());
      h1_eleEtaResol_ -> Fill( myEle->eta() - mcEle[i]->momentum().eta() );
      h1_elePhiResol_ -> Fill( myEle->phi() - mcEle[i]->momentum().phi() );
      
      const reco::SuperCluster* mySC=&(*(myEle->superCluster()));
      
      if(myEle->classification() < 100 )  h_PTres_[loopFlag_]->Fill( myEle->gsfTrack()->ptMode() / mcEle[i]->momentum().perp() );
      if(myEle->classification() >= 100 ) h_PTres_Endcap_[loopFlag_]->Fill( myEle->gsfTrack()->ptMode() / mcEle[i]->momentum().perp() );
	  
      h_ESCEtrue_[loopFlag_]      -> Fill(mySC->energy()/mcEle[i]->momentum().e());
      h_ESCEtrueVsEta_[loopFlag_] -> Fill(fabs(mySC->position().eta()),mySC->energy()/mcEle[i]->momentum().e());
	  
      double corrSCenergy = ( mySC->energy() )*getEtaCorrection(myEle);
      h_ESCcorrEtrue_[loopFlag_]      -> Fill(corrSCenergy/mcEle[i]->momentum().e());
      h_ESCcorrEtrueVsEta_[loopFlag_] -> Fill(fabs(mySC->position().eta()),corrSCenergy/mcEle[i]->momentum().e());
      
      std::vector<DetId> mySCRecHits = mySC->seed()->getHitsByDetId();
    }
  }
}

int JPsiCalibration::ringNumberCorrector(int k) {

  int index=-999;
  
  if( calibMode_ == "RING"){
    if(k>=0 && k<=84)    index = k - 85;    
    if(k>=85 && k<=169)  index = k - 84;
    if(k>=170 && k<=208) index = -k + 84;
    if(k>=209 && k<=247) index = k - 123;    
  }
  
  else if( calibMode_ == "MODULE"){
    if(k>=0 && k<=71)   index = k - 72;
    if(k>=72 && k<=143) index = k - 71;
  }
  return index;
}


double JPsiCalibration::getEtaCorrection(const reco::PixelMatchGsfElectron* ele){

  double correction(1.);

  if(ele->classification() ==0 ||
     ele->classification() ==10 ||
     ele->classification() ==20)
    correction = fEtaBarrelGood(ele->superCluster()->eta());
                                                                                                                                               
  if(ele->classification() ==100 ||
     ele->classification() ==110 ||
     ele->classification() ==120)
    correction = fEtaEndcapGood(ele->superCluster()->eta());
                                                                                                                                               
  if(ele->classification() ==30 ||
     ele->classification() ==31 ||
     ele->classification() ==32 ||
     ele->classification() ==33 ||
     ele->classification() ==34)
    correction = fEtaBarrelBad(ele->superCluster()->eta());


  if(ele->classification() ==130 ||
     ele->classification() ==131 ||
     ele->classification() ==132 ||
     ele->classification() ==133 ||
     ele->classification() ==134)
    correction = fEtaEndcapBad(ele->superCluster()->eta());
 
  return correction;                                                                                                                                              
}

std::pair<DetId, double> JPsiCalibration::getHottestDetId(std::vector<DetId> mySCRecHits, const EBRecHitCollection* ebhits, const EERecHitCollection* eehits){
  
  double maxEnergy = -9999.;
  const EcalRecHit* hottestRecHit;
  
  std::pair<DetId, double> myPair (DetId(0), -9999.);

  for( std::vector<DetId>::const_iterator idIt=mySCRecHits.begin(); idIt != mySCRecHits.end(); idIt++){
    
    if (idIt->subdetId() == EcalBarrel ) {
      hottestRecHit = & (* ( ebhits->find(*idIt) ) );
      if( hottestRecHit == & (*( ebhits->end())) ) continue;
    }
    else if (idIt->subdetId() == EcalEndcap ) {
      hottestRecHit  = & (* ( eehits->find(*idIt) ) );
      if( hottestRecHit == & (*( eehits->end())) ) continue;
    }
    
    if(hottestRecHit && hottestRecHit->energy() > maxEnergy){
      maxEnergy = hottestRecHit->energy();    
      myPair.first  = hottestRecHit ->id();
      myPair.second = maxEnergy;
    }
    
  } //end loop to find hottest RecHit    
  
  return myPair;
  }

double JPsiCalibration::getE5x5(std::vector<DetId> mySCRecHits, const EBRecHitCollection* ebhits, const EERecHitCollection* eehits){
  
  DetId hottest = getHottestDetId(mySCRecHits, ebhits, eehits).first;
  
  double e5x5 = 0.;
  if(hottest.subdetId() == EcalBarrel ){
    EBDetId ebhottest(hottest);
    
    for(int i = -2; i <= 2; i++){
      for(int j = -2; j <= 2; j++){
	const EcalRecHit* rh = 0;       
	if(!EBDetId::validDetId(ebhottest.ieta() + i,ebhottest.iphi() + j ) ) continue;
	  
	EBDetId b(ebhottest.ieta() + i, ebhottest.iphi() + j );
	rh  = & (* ( ebhits->find( b ) ) );
	
	if(rh) e5x5 += rh->energy();
      }
    }
  }
  
  if(hottest.subdetId() == EcalEndcap ){
    EEDetId eehottest(hottest);
    
    for(int i = -2; i <= 2; i++){
      for(int j = -2; j <= 2; j++){
	const EcalRecHit* rh = 0;   
	
	if(!EEDetId::validDetId(eehottest.ix() + i, eehottest.iy() + j, eehottest.zside()  ) ) continue;
	
	EEDetId b(eehottest.ix() + i, eehottest.iy() + j, eehottest.zside()  );
	rh  = & (* ( eehits->find( b ) ) );
	
	if(rh) e5x5 += rh->energy();
      }
    }
  }
  
  return e5x5;
}


bool JPsiCalibration::xtalIsOnModuleBorder( EBDetId myEBDetId ){
  
  bool myBool(false); 
  
  short ieta = myEBDetId.ieta();
  short iphi = myEBDetId.iphi();
  
  myBool = ( abs( ieta )  == 1 || abs( ieta ) == 25
	     || abs( ieta )  ==26 || abs( ieta ) == 45
	     || abs( ieta )  ==46 || abs( ieta ) == 65
	     || abs( ieta )  ==66 || abs( ieta ) == 85 );
  
  for(int i = 0; i < 19; i++){
    if(iphi == ( 20*i + 1 ) || iphi == 20*i ) myBool = true;
  }
  
  return myBool;
}


float JPsiCalibration::computeCoefficientDistanceAtIteration( float v1[250], float v2[250], int size ){

  float dist(0.);
 
  for(int i =0; i < size; i++){
    
    bool isNearCrack = false;
    if( calibMode_ == "RING"){  //exclude non-calibrated rings from computation
      isNearCrack = ( abs( ringNumberCorrector(i) ) == 1 || abs( ringNumberCorrector(i) ) == 25 ||
		      abs( ringNumberCorrector(i) ) == 26 || abs( ringNumberCorrector(i) ) == 45 ||
		      abs( ringNumberCorrector(i) ) == 46 || abs( ringNumberCorrector(i) ) == 65 ||
		      abs( ringNumberCorrector(i) ) == 66 || abs( ringNumberCorrector(i) ) == 85 ||
		      abs( ringNumberCorrector(i) ) == 86 || abs( ringNumberCorrector(i) ) == 124 );
    }
    
    if(!isNearCrack) dist += pow( v1[i]-v2[i], 2 );
  }
  
  dist = sqrt(dist) / size;
  return dist;
}


void JPsiCalibration::resetVariables(){

 BBZN=0;
 EBZN=0;
 EEZN=0;

 TOTAL_ELECTRONS_IN_BARREL = 0;
 TOTAL_ELECTRONS_IN_ENDCAP = 0;
 CRACK_ELECTRONS_IN_BARREL = 0;
 CRACK_ELECTRONS_IN_ENDCAP = 0;


 BARREL_ELECTRONS_BEFORE_BORDER_CUT = 0;
 BARREL_ELECTRONS_AFTER_BORDER_CUT = 0;

 return;

}


void JPsiCalibration::resetHistograms(){

 h1_eventsBeforeHLTSelection_     -> Reset();
 h1_eventsAfterHLTSelection_      -> Reset();
 h1_eventsBefore2GsfSelection_    -> Reset();
 h1_eventsAfter2GsfSelection_     -> Reset();
 h1_eventsBeforeEtSelection_      -> Reset();
 h1_eventsAfterEtSelection_       -> Reset();
 h1_eventsBeforeEleIDSelection_   -> Reset();
 h1_eventsAfterEleIDSelection_    -> Reset();
 h1_eventsBeforeEleIsolSelection_ -> Reset();
 h1_eventsAfterEleIsolSelection_  -> Reset();
 h1_eventsBeforeTrackerIsolSelection_ -> Reset();
 h1_eventsAfterTrackerIsolSelection_  -> Reset();
 h1_eventsBeforeBorderSelection_  -> Reset();
 h1_eventsAfterBorderSelection_   -> Reset();
 h1_eventsBeforeInvMassSelection  -> Reset();
 h1_eventsAfterInvMassSelection   -> Reset();

 h1_ZCandMult_-> Reset();

 h2_occupancy_vsEtaPhi     -> Reset();
 h1_seedOverSC_Barrel_     -> Reset();
 h1_seedOverSC_Endcap_     -> Reset();
 h1_preshowerOverSC_       -> Reset();
 h1_electronCosTheta_TK_   -> Reset();
 h1_electronCosTheta_SC_   -> Reset();
 h1_electronCosTheta_SC_TK_-> Reset();

 h_ESCEtrue_vs_ESC         -> Reset();
 
 h1_Theta12Resolution_TK_BB_ -> Reset();
 h1_Theta12Resolution_SC_BB_ -> Reset();
 h1_Theta12Resolution_TK_EE_ -> Reset();
 h1_Theta12Resolution_SC_EE_ -> Reset();
 h1_Theta12_TK_BB_           -> Reset();
 h1_Theta12_SC_BB_           -> Reset();
 h1_Theta12_TK_EE_           -> Reset();
 h1_Theta12_SC_EE_           -> Reset();

 h1_eleEtaResol_ -> Reset();
 h1_elePhiResol_ -> Reset();
 h1_zMassResol_  -> Reset(); 

 h1_reco_ZMass_             -> Reset();
 h1_reco_ZMassBB_           -> Reset();
 h1_reco_ZMassEB_           -> Reset();
 h1_reco_ZMassEE_           -> Reset();
 h1_reco_ZMassCorr_         -> Reset();
 h1_reco_ZMassCorrBB_       -> Reset();
 h1_reco_ZMassCorrEB_       -> Reset();
 h1_reco_ZMassCorrEE_       -> Reset();

 h1_sEEeb_before_ -> Reset();
 h1_sEEee_before_ -> Reset();
 h1_sEEeb_after_  -> Reset();
 h1_sEEee_after_  -> Reset();

 h1_massResolutionTerm_total_BB_       -> Reset();
 h1_massResolutionTerm_angle_BB_       -> Reset();
 h1_massResolutionTerm_energy_BB_      -> Reset();
 h1_massResolutionTerm_angle_ratio_BB_ -> Reset();
 h1_massResolutionTerm_total_EE_       -> Reset();
 h1_massResolutionTerm_angle_EE_       -> Reset();
 h1_massResolutionTerm_energy_EE_      -> Reset();
 h1_massResolutionTerm_angle_ratio_EE_ -> Reset();

 for(int j = 0 ; j <25; j++){
   h_ESCEtrue_array_barrel[j] -> Reset();
   h_ESCEtrue_array_endcap[j] -> Reset();
  }

 for (int i=0;i<2;i++) {
   h_eleEffEta_[i] ->Reset();
   h_eleEffPhi_[i] ->Reset(); 
   h_eleEffPt_[i]  ->Reset();
 }
 
 h1_occupancyVsEta_  -> Reset();
 h1_occupancy_       -> Reset();
 h1_occupancyBarrel_ -> Reset();
 h1_occupancyEndcap_ -> Reset();

 return;
}


void JPsiCalibration::printStatistics(){


  std::cout<< "[ CHECK ON BARREL ELECTRON NUMBER ]"<<" first "<<BARREL_ELECTRONS_BEFORE_BORDER_CUT<<" second "<<TOTAL_ELECTRONS_IN_BARREL << std::endl;
    
  std::cout<< "[ EFFICIENCY OF THE BORDER SELECTION ]" << (float)BARREL_ELECTRONS_AFTER_BORDER_CUT / (float) BARREL_ELECTRONS_BEFORE_BORDER_CUT << endl;
  
  std::cout<< "[ EFFICIENCY OF THE CRACK SELECTION ] BARREL: " << (float)CRACK_ELECTRONS_IN_BARREL / (float) TOTAL_ELECTRONS_IN_BARREL << " ENDCAP: "<< (float)CRACK_ELECTRONS_IN_ENDCAP / (float) TOTAL_ELECTRONS_IN_ENDCAP << endl;
  
  
  ofstream fout("JPsiStatistics.txt");
  
  if(!fout) std::cout << "Cannot open output file.\n";

  fout<<"JPsiStatistics"<<std::endl;
  fout<<"\n"<<std::endl;

  fout<<"##########################GEN#########################"<<std::endl;
  fout<<"##################JPsi with Barrel-Barrel electrons: "<<(float)MCZBB/NEVT<<std::endl;
  fout<<"##################JPsi with Barrel-Endcap electrons: "<<(float)MCZEB/NEVT<<std::endl;
  fout<<"##################JPsi with Endcap-Endcap electrons: "<<(float)MCZEE/NEVT<<std::endl;

  fout.close();
}

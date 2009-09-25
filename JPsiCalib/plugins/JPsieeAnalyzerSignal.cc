// framework
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// root & others
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TFile.h>
#include <TTree.h>
#include "TMath.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include <map>
#include <cmath>
#include <string>

// my classes
#include "JPsieeAnalyzerSignal.h"

using namespace reco;
using namespace std;
using namespace edm;

JPsieeAnalyzerSignal::JPsieeAnalyzerSignal(const edm::ParameterSet& iConfig) { 
  
  fOutFileTreeName_        = iConfig.getUntrackedParameter<string>("fileTree");
  electronCollection_      = iConfig.getParameter<InputTag>("electronCollection");
  ecalRechitsCollectionEB_ = iConfig.getParameter<InputTag>("ecalrechitsCollectionEB");
  ecalRechitsCollectionEE_ = iConfig.getParameter<InputTag>("ecalrechitsCollectionEE");  
  triggerResults_	   = iConfig.getParameter<InputTag>("triggerResults");
  isSignal_                = iConfig.getUntrackedParameter<bool>("isSignal");
}


void JPsieeAnalyzerSignal::beginJob(const EventSetup&) { 

  OutputTree = new JPsiTree(fOutFileTreeName_.c_str()); 
} 

void JPsieeAnalyzerSignal::endJob() { 

  OutputTree->save(); 
}

void JPsieeAnalyzerSignal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // 1) HLT
  Handle<edm::TriggerResults> HLTR;  
  try { iEvent.getByLabel(triggerResults_, HLTR); }
  catch(cms::Exception& ex){ edm::LogError("ProblemHltTriggserResults") << "Trigger results: " << triggerResults_ << " not found"; }
  if (!HLTR.isValid()) throw cms::Exception("ProductNotValid") << "TriggerResults product not valid";
  
  bool hltJPsi    = false;
  bool hltUpsilon = false;
  bool hltBoth    = false;
  if (HLTR.isValid()){
    triggerNames_.init(*HLTR);
    for ( unsigned int iHlt=0; iHlt < HLTR->size(); iHlt++ ) {
      if (triggerNames_.triggerName(iHlt)=="HLT_DoublePhoton5_Jpsi_L1R"    && HLTR->accept(iHlt)==1) hltJPsi = true;
      if (triggerNames_.triggerName(iHlt)=="HLT_DoublePhoton5_Upsilon_L1R" && HLTR->accept(iHlt)==1) hltUpsilon = true;
      if (triggerNames_.triggerName(iHlt)=="HLT_DoublePhoton5_eeRes_L1R"   && HLTR->accept(iHlt)==1) hltBoth = true;
    }
  }
  if ( hltJPsi) intHltJPsi = 1;
  if (!hltJPsi) intHltJPsi = 0;
  if ( hltUpsilon) intHltUpsilon = 1;
  if (!hltUpsilon) intHltUpsilon = 0;
  if ( hltBoth) intHltBoth = 1;
  if (!hltBoth) intHltBoth = 0;

  
  // 2) calo topology
  const CaloTopology *theTopology;
  edm::ESHandle<CaloTopology> topology;
  iSetup.get<CaloTopologyRecord>().get(topology);
  if (topology.isValid()) theTopology = topology.product();
  
  // 3) generated electrons
  const HepMC::GenEvent *myGenEvent;
  Handle<edm::HepMCProduct> hepMC;
  iEvent.getByLabel("generator", hepMC);
  myGenEvent = hepMC->GetEvent();

  // 4) reconstructed electrons  
  Handle<GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons); 

  // 5) ECAL rechits to build cluster shapes
  const EcalRecHitCollection *rechitsEB;
  const EcalRecHitCollection *rechitsEE;
  Handle<EcalRecHitCollection> pRechitsEB;
  Handle<EcalRecHitCollection> pRechitsEE;
  iEvent.getByLabel(ecalRechitsCollectionEB_,pRechitsEB); 
  iEvent.getByLabel(ecalRechitsCollectionEE_,pRechitsEE); 
  if (pRechitsEB.isValid()) rechitsEB = pRechitsEB.product();
  if (pRechitsEE.isValid()) rechitsEE = pRechitsEE.product();
  
  // 6) electron identification (default cut based)
  eleIdResults_ = new eleIdContainer(4);
  iEvent.getByLabel( "eidLoose",       (*eleIdResults_)[0] ); 
  iEvent.getByLabel( "eidRobustLoose", (*eleIdResults_)[1] );
  iEvent.getByLabel( "eidRobustTight", (*eleIdResults_)[2] );
  iEvent.getByLabel( "eidTight",       (*eleIdResults_)[3] );
  
  // 7) pt hat
  Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByLabel("generator", genEventInfo);
  float ptHat = genEventInfo->qScale();
  

  // ---------------------------------------------------------------------
  // generator level informations
  
  // if running on background samples we skip events with true j/psi 
  bool isAJPsi = false;
  if( !isSignal_ ) {
    HepMC::GenEvent::particle_const_iterator mcIter;
    for ( mcIter=myGenEvent->particles_begin(); mcIter != myGenEvent->particles_end(); mcIter++ ) {
      if ( fabs((*mcIter)->pdg_id())==11) {
	HepMC::GenParticle* mother=0;
	if ( (*mcIter)->production_vertex() ) {
	  if ((*mcIter)->production_vertex()->particles_begin(HepMC::parents) != (*mcIter)->production_vertex()->particles_end(HepMC::parents)) { 
	    mother = *((*mcIter)->production_vertex()->particles_begin(HepMC::parents));
	  }}
	if ((mother != 0) && (mother->pdg_id() == 443)) isAJPsi = true; 
      }}
  }
  
  // if background sample and a j/psi is found, skip the event
  if( isSignal_ || !isAJPsi) {

    // ------------------------------------------------
    // run infos:
    int theRun         = iEvent.id().run();
    int theEvent       = iEvent.id().event();
    int theLumiSection = iEvent.luminosityBlock();
    OutputTree->fillRunInfos( theRun, theEvent, theLumiSection, ptHat );
    

    // counters for the tree
    int numberMcParticle  = 0;
    int numberOfElectrons = 0;
    

    // for signal only, MC info: electrons from j/psi
    int theMcPc = 0;
    if( isSignal_ ) {
      
      float pxGenEle[2];
      float pyGenEle[2];
      float pzGenEle[2];
      float eneGenEle[2];
      int chargeGenEle[2];
      for(int ii=0; ii<2; ii++) chargeGenEle[ii] = 0;
      
      HepMC::GenEvent::particle_const_iterator mcIter;
      for ( mcIter=myGenEvent->particles_begin(); mcIter != myGenEvent->particles_end(); mcIter++ ) {
	
	// select electrons
	if ( fabs((*mcIter)->pdg_id()) == 11) {
	  
	  // electrons from J/Psi
	  HepMC::GenParticle* mother=0;
	  if ( (*mcIter)->production_vertex() ) {
	    if ((*mcIter)->production_vertex()->particles_begin(HepMC::parents) != (*mcIter)->production_vertex()->particles_end(HepMC::parents)) { 
	      mother = *((*mcIter)->production_vertex()->particles_begin(HepMC::parents));
	    }}
	  
	  if ((mother != 0) && (mother->pdg_id() == 443)) {
	    
	    // filling the tree
	    pxGenEle[theMcPc]  = ((*mcIter)->momentum()).x();
	    pyGenEle[theMcPc]  = ((*mcIter)->momentum()).y();
	    pzGenEle[theMcPc]  = ((*mcIter)->momentum()).z();
	    eneGenEle[theMcPc] = ((*mcIter)->momentum()).e();	
	    if((*mcIter)->pdg_id()==11)  chargeGenEle[theMcPc] = -1;
	    if((*mcIter)->pdg_id()==-11) chargeGenEle[theMcPc] = +1;
	    if(chargeGenEle[theMcPc]) OutputTree->fillGenerated( chargeGenEle[theMcPc], pxGenEle[theMcPc], pyGenEle[theMcPc], pzGenEle[theMcPc], eneGenEle[theMcPc] );
	    theMcPc++;
	  }
	}
      }
          
    } // is signal
    
    // for the tree
    numberMcParticle = theMcPc;
    
    

    // reconstructed electrons
    int countMP = 0;
    numberOfElectrons = gsfElectrons->size();
    GsfElectronCollection::const_iterator eleIter;
    for (eleIter=gsfElectrons->begin(); eleIter != gsfElectrons->end() ; eleIter++) { 
      
      // reference to the electron
      const reco::GsfElectronRef eleRef(gsfElectrons,countMP);

      // taking the supercluster
      SuperClusterRef theScRef = eleIter->superCluster();

      // cluster shape variables
      const EcalRecHitCollection *rechits = 0;
      float seedEta = theScRef->seed()->position().eta();      
      if( fabs(seedEta) < 1.479 ) rechits = rechitsEB;
      else rechits = rechitsEE;     
      float e3x3   = EcalClusterTools::e3x3( *(theScRef->seed()), &(*rechits), theTopology );
      float e5x5   = EcalClusterTools::e5x5( *(theScRef->seed()), &(*rechits), theTopology );
      float s9s25  = e3x3/e5x5;
      float sigmaEtaEta   = eleIter->sigmaEtaEta();        
      float sigmaIetaIeta = eleIter->sigmaIetaIeta();    
      float hcalOverEcal  = eleIter->hcalOverEcal();

      // POG computed electronID
      const eleIdMap & eleIdLooseVal       = *( (*eleIdResults_)[0] );
      const eleIdMap & eleIdRobustLooseVal = *( (*eleIdResults_)[1] );
      const eleIdMap & eleIdRobustTightVal = *( (*eleIdResults_)[2] );
      const eleIdMap & eleIdTightVal       = *( (*eleIdResults_)[3] );
      int eleIdLoose    = 0;
      int eleIdRobLoose = 0;
      int eleIdRobTight = 0;
      int eleIdTight    = 0;
      if ( eleIdLooseVal[eleRef] )       eleIdLoose    = 1;
      if ( eleIdRobustLooseVal[eleRef] ) eleIdRobLoose = 1;
      if ( eleIdRobustTightVal[eleRef] ) eleIdRobTight = 1;
      if ( eleIdTightVal[eleRef] )       eleIdTight    = 1;
    
      // PF electron id
      float pfMva = eleIter->mva();

      // other useful quantities to study electron id
      float eSuperClusterOverP        = eleIter->eSuperClusterOverP();        
      float eSeedClusterOverPout      = eleIter->eSeedClusterOverPout();      
      float deltaEtaSuperClusterAtVtx = eleIter->deltaEtaSuperClusterTrackAtVtx(); 
      float deltaEtaSeedClusterAtCalo = eleIter->deltaEtaSeedClusterTrackAtCalo(); 
      float deltaPhiSuperClusterAtVtx = eleIter->deltaPhiSuperClusterTrackAtVtx(); 
      float deltaPhiSeedClusterAtCalo = eleIter->deltaPhiSeedClusterTrackAtCalo(); 
      float fbrem                     = eleIter->fbrem();

      // POG computed isolation
      float dr03TkSumPt              = eleIter->dr03TkSumPt();
      float dr04TkSumPt              = eleIter->dr04TkSumPt();
      float dr03EcalRecHitSumEt      = eleIter->dr03EcalRecHitSumEt();
      float dr04EcalRecHitSumEt      = eleIter->dr04EcalRecHitSumEt();
      float dr03HcalDepth1TowerSumEt = eleIter->dr03HcalDepth1TowerSumEt();
      float dr04HcalDepth1TowerSumEt = eleIter->dr04HcalDepth1TowerSumEt();
      float dr03HcalDepth2TowerSumEt = eleIter->dr03HcalDepth2TowerSumEt();
      float dr04HcalDepth2TowerSumEt = eleIter->dr04HcalDepth2TowerSumEt();

      // ECAL based quantities
      float rawSCenergy     = theScRef->rawEnergy();
      float rawESenergy     = theScRef->preshowerEnergy();      
      float ecalEnergy      = eleIter->ecalEnergy();          // ecal corrected energy (if !isEcalEnergyCorrected this value is identical to the supercluster energy
      float ecalEnergyError = eleIter->ecalEnergyError();     // error on correctedCaloEnergy

      // tracker based quantities
      float trackPx            = eleIter->trackMomentumAtVtx().x();
      float trackPy            = eleIter->trackMomentumAtVtx().y();
      float trackPz            = eleIter->trackMomentumAtVtx().z();
      float trackEta           = eleIter->trackMomentumAtVtx().eta();
      float trackPhi           = eleIter->trackMomentumAtVtx().phi();
      float trackMomentumError = eleIter->trackMomentumError();          // track momentum error from gsf fit

      // corrections
      bool isEcalEnergyCorrected = eleIter->isEcalEnergyCorrected();     // true if ecal energy has been corrected
      bool isMomentumCorrected   = eleIter->isMomentumCorrected();       // true if E-p combination has been applied (if not the electron momentum is the ecal corrected energy)
      int intIsEcalEnergyCorrected, intIsMomentumCorrected;
      if ( isEcalEnergyCorrected ) intIsEcalEnergyCorrected = 1;
      if (!isEcalEnergyCorrected ) intIsEcalEnergyCorrected = 0;
      if ( isMomentumCorrected)    intIsMomentumCorrected = 1;
      if (!isMomentumCorrected)    intIsMomentumCorrected = 0;
     
      // particle flow or standard electron
      bool isEcalDriven   = eleIter->isEcalDriven();
      bool isParticleFlow = eleIter->isTrackerDriven();
      int intIsEcalDriven, intIsParticleFlow;
      if ( isEcalDriven )   intIsEcalDriven = 1;
      if (!isEcalDriven )   intIsEcalDriven = 0;
      if ( isParticleFlow ) intIsParticleFlow = 1; 
      if (!isParticleFlow ) intIsParticleFlow = 0; 

      // general infos from candidate: should be the better estimate
      int eleClass    = eleIter->classification();
      int eleCharge   = eleIter->charge();
      float elePx     = eleIter->momentum().x();
      float elePy     = eleIter->momentum().y();
      float elePz     = eleIter->momentum().z();
      float eleEta    = eleIter->momentum().eta();
      float elePhi    = eleIter->momentum().phi();
      float eleEnergy = eleIter->energy();
      float eleEt     = eleIter->et();
      float electronMomentumError = eleIter->electronMomentumError(); // the final electron momentum error

      OutputTree->fillElectrons( e3x3, e5x5, s9s25, sigmaEtaEta, sigmaIetaIeta, hcalOverEcal, eleIdLoose, eleIdRobLoose, eleIdRobTight, eleIdTight, pfMva, eSuperClusterOverP, eSeedClusterOverPout, deltaEtaSuperClusterAtVtx, deltaEtaSeedClusterAtCalo, deltaPhiSuperClusterAtVtx, deltaPhiSeedClusterAtCalo, fbrem, dr03TkSumPt, dr04TkSumPt, dr03EcalRecHitSumEt, dr04EcalRecHitSumEt, dr03HcalDepth1TowerSumEt, dr04HcalDepth1TowerSumEt, dr03HcalDepth2TowerSumEt, dr04HcalDepth2TowerSumEt, rawSCenergy, rawESenergy, ecalEnergy, ecalEnergyError,trackPx, trackPy, trackPz, trackEta, trackPhi, trackMomentumError,intIsEcalEnergyCorrected, intIsMomentumCorrected,intIsEcalDriven, intIsParticleFlow, eleClass, eleCharge, elePx, elePy, elePz, eleEta, elePhi, eleEnergy, eleEt, electronMomentumError);				

      countMP++;

    } // loop over electrons

    
    // filling the tree: some summary numbers
    int signal = 0;
    if (isSignal_) signal = 1;
    OutputTree->fillGeneral(signal, numberMcParticle, numberOfElectrons, intHltJPsi, intHltUpsilon, intHltBoth);
    
    OutputTree->store();

  } // signal or background


}


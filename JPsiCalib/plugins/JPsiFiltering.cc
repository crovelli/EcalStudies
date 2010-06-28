// framework
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// root & others
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
#include "JPsiFiltering.h"

using namespace reco;
using namespace std;
using namespace edm;


JPsiFiltering::JPsiFiltering(const edm::ParameterSet& iConfig) { 
  
  fOutFileTreeName_    = iConfig.getUntrackedParameter<string>("fileTree");
  electronCollection_  = iConfig.getParameter<InputTag>("electronCollection");
  triggerResults_      = iConfig.getParameter<InputTag>("triggerResults");
}

void JPsiFiltering::beginJob(const EventSetup&) { 
  
  OutputTree = new FilteringTree(fOutFileTreeName_.c_str()); 
} 

void JPsiFiltering::endJob() { 
  
  OutputTree->save(); 
}

void JPsiFiltering::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // 1) generated electrons
  const HepMC::GenEvent *myGenEvent;
  Handle<edm::HepMCProduct> hepMC;
  iEvent.getByLabel("generator", hepMC);
  myGenEvent = hepMC->GetEvent();

  // 2) reconstructed electrons  
  Handle<GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons); 
  
  // 3) to get default cut based eleID
  eleIdResults_ = new eleIdContainer(4);
  iEvent.getByLabel( "eidLoose",       (*eleIdResults_)[0] ); 
  iEvent.getByLabel( "eidRobustLoose", (*eleIdResults_)[1] );
  iEvent.getByLabel( "eidRobustTight", (*eleIdResults_)[2] );
  iEvent.getByLabel( "eidTight",       (*eleIdResults_)[3] );

  // 4) pt hat
  Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByLabel("generator", genEventInfo);
  float ptHat = genEventInfo->qScale();

  // ---------------------------------------------------------------------
  // filters infos:
  Handle<edm::TriggerResults> FILTERR;  
  iEvent.getByLabel(triggerResults_, FILTERR);
  //  cout << "EEE triggerresults valid " << FILTERR.isValid() << endl;
  //cout << "EEE triggerresults failedtoget " << FILTERR.failedToGet() << endl;
   
  int filter1 = 0;
  int filter2 = 0;
  int filter3 = 0;
  int filter4 = 0;
  int filter5 = 0;
  int filter6 = 0;
  int filter7 = 0;
  int filter8 = 0;
  if (FILTERR.isValid()){
    triggerNames_.init(*FILTERR);
    //bool res = triggerNames_.init(*FILTERR);
    //   cout << "tirgger names is " << res << " " << FILTERR->size() << endl;
    for ( unsigned int iFilter=0; iFilter < FILTERR->size(); iFilter++ ) {
	//cout<<" &&&& filter " << triggerNames_.triggerName(iFilter) << " " <<  FILTERR->accept(iFilter) << endl;
      if (triggerNames_.triggerName(iFilter)=="p1" && FILTERR->accept(iFilter)==1) filter1 = 1;
      if (triggerNames_.triggerName(iFilter)=="p2" && FILTERR->accept(iFilter)==1) filter2 = 1;
      if (triggerNames_.triggerName(iFilter)=="p3" && FILTERR->accept(iFilter)==1) filter3 = 1;
      if (triggerNames_.triggerName(iFilter)=="p4" && FILTERR->accept(iFilter)==1) filter4 = 1;
      if (triggerNames_.triggerName(iFilter)=="p5" && FILTERR->accept(iFilter)==1) filter5 = 1;
      if (triggerNames_.triggerName(iFilter)=="p6" && FILTERR->accept(iFilter)==1) filter6 = 1;
      if (triggerNames_.triggerName(iFilter)=="p7" && FILTERR->accept(iFilter)==1) filter7 = 1;
      if (triggerNames_.triggerName(iFilter)=="p8" && FILTERR->accept(iFilter)==1) filter8 = 1;
    }
  }

  // ---------------------------------------------------------------------
  // run infos:
  int theRun         = iEvent.id().run();
  int theEvent       = iEvent.id().event();
  int theLumiSection = iEvent.luminosityBlock();
  OutputTree->fillRunInfos( theRun, theEvent, theLumiSection, ptHat );
  
  // counters for the tree
  int numberMcParticle  = 0;
  int numberOfElectrons = 0;
  
  // MC info
  int theMcPc        = 0;
  int theMother      = 0;  
  int theMotherIndex = 0;  
  HepMC::GenEvent::particle_const_iterator mcIter;
  HepMC::GenEvent::particle_const_iterator motherIter;

  for ( mcIter=myGenEvent->particles_begin(); mcIter != myGenEvent->particles_end(); mcIter++ ) {
    
    if (theMcPc>1499) continue; 
    
    // the mother of this particle
    HepMC::GenParticle* mother=0;
    if ( (*mcIter)->production_vertex() ) {
      if ((*mcIter)->production_vertex()->particles_begin(HepMC::parents) != (*mcIter)->production_vertex()->particles_end(HepMC::parents)) { 
	mother = *((*mcIter)->production_vertex()->particles_begin(HepMC::parents));
      }}
    
    // to find the mother index
    for ( motherIter=myGenEvent->particles_begin(); motherIter != myGenEvent->particles_end(); motherIter++ ) {
      if (theMother>5499) continue; 
      HepMC::GenParticle* mother2=0;
      if ( (*motherIter)->production_vertex() ) {
	if ((*motherIter)->production_vertex()->particles_begin(HepMC::parents) != (*motherIter)->production_vertex()->particles_end(HepMC::parents)) { 
	  mother2 = *((*motherIter)->production_vertex()->particles_begin(HepMC::parents));
	}}
      if (mother2==mother) { theMotherIndex = theMother; break; }
      theMother++;
    }
    
    // filling the tree    
    int mothId = 0;
    if (mother != 0) mothId = mother->pdg_id();
    if (mother == 0) mothId = -9999;
    
    OutputTree->fillGenerated( ((*mcIter)->momentum()).x(),
			       ((*mcIter)->momentum()).y(),
			       ((*mcIter)->momentum()).z(),
			       ((*mcIter)->momentum()).e(),
			       ((*mcIter)->momentum()).eta(),
			       ((*mcIter)->momentum()).phi(),
			       (*mcIter)->pdg_id(), 
			       (*mcIter)->status(), 
			       mothId, theMother );
    theMcPc++;
  }
  numberMcParticle = theMcPc;
  
  
  // reconstructed electrons
  int countMP = 0;
  numberOfElectrons = gsfElectrons->size();
  GsfElectronCollection::const_iterator eleIter;
  for (eleIter=gsfElectrons->begin(); eleIter != gsfElectrons->end() ; eleIter++) { 
    
    SuperClusterRef theScRef   = eleIter->get<SuperClusterRef>();
    const reco::GsfElectronRef eleRef(gsfElectrons,countMP);
    
    // electronID
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
    
    // POG computed isolation    
    float dr03TkSumPt              = eleIter->dr03TkSumPt();
    float dr04TkSumPt              = eleIter->dr04TkSumPt();
    float dr03EcalRecHitSumEt      = eleIter->dr03EcalRecHitSumEt();
    float dr04EcalRecHitSumEt      = eleIter->dr04EcalRecHitSumEt();
    float dr03HcalDepth1TowerSumEt = eleIter->dr03HcalDepth1TowerSumEt();
    float dr04HcalDepth1TowerSumEt = eleIter->dr04HcalDepth1TowerSumEt();
    float dr03HcalDepth2TowerSumEt = eleIter->dr03HcalDepth2TowerSumEt();
    float dr04HcalDepth2TowerSumEt = eleIter->dr04HcalDepth2TowerSumEt();


    // filling the tree with infos on all the reconstructed electrons
    int eleCharge   = eleIter->charge();
    float elePx     = eleIter->momentum().x();
    float elePy     = eleIter->momentum().y();
    float elePz     = eleIter->momentum().z();
    float eleEta    = eleIter->momentum().eta();
    float elePhi    = eleIter->momentum().phi();
    float eleEnergy = eleIter->energy();
    float eleEt     = eleIter->et();
    float HoE       = eleIter->hadronicOverEm();
    float deta      = eleIter->deltaEtaSuperClusterTrackAtVtx();
    float dphi      = eleIter->deltaPhiSuperClusterTrackAtVtx();
    float eleEoP    = eleIter->eSuperClusterOverP();
    OutputTree->fillElectrons( eleCharge, elePx, elePy, elePz, eleEta, elePhi, eleEnergy, eleEt, HoE, deta, dphi, eleEoP, eleIdLoose, eleIdRobLoose, eleIdRobTight, eleIdTight, dr03TkSumPt, dr04TkSumPt, dr03EcalRecHitSumEt, dr04EcalRecHitSumEt, dr03HcalDepth1TowerSumEt, dr04HcalDepth1TowerSumEt, dr03HcalDepth2TowerSumEt, dr04HcalDepth2TowerSumEt );
    
    countMP++;
    
  } // loop over electrons
  
  
  // filling the tree: some summary numbers
  OutputTree->fillGeneral(numberMcParticle, numberOfElectrons, filter1, filter2, filter3, filter4, filter5, filter6, filter7, filter8);
  
  OutputTree->store();
  
  // deleting
  delete eleIdResults_;
}


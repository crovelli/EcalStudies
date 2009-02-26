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
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"	
#include "DataFormats/VertexReco/interface/Vertex.h"	
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h" 
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/FTSFromVertexToPointFactory.h"
#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 

// root & others
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TFile.h>
#include <TH1.h>
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
#include "TrackerIsolation.h"		
#include "CalotowerIsolation.h"

using namespace reco;
using namespace std;
using namespace edm;


// class etGreater declaration (et sort)
class etGreater {
public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return ( (i.energy())*sin((i).position().theta()) >  ((j.energy())*sin((j).position().theta()) ));
  }
};




JPsieeAnalyzerSignal::JPsieeAnalyzerSignal(const edm::ParameterSet& iConfig) { 
  
  fOutFileTreeName_          = iConfig.getUntrackedParameter<string>("fileTree");
  electronCollection_        = iConfig.getParameter<InputTag>("electronCollection");
  superclusterCollectionEB_  = iConfig.getParameter<InputTag>("superclusterCollectionEB");
  superclusterCollectionEE_  = iConfig.getParameter<InputTag>("superclusterCollectionEE");
  ecalRechitsCollectionEB_   = iConfig.getParameter<InputTag>("ecalrechitsCollectionEB");
  ecalRechitsCollectionEE_   = iConfig.getParameter<InputTag>("ecalrechitsCollectionEE");  
  hcalRechitsCollectionHBHE_ = iConfig.getParameter<InputTag>("hcalrechitsCollectionHBHE");  
  tracksCollection_          = iConfig.getParameter<InputTag>("tracksCollection");
  calotowersCollection_      = iConfig.getParameter<InputTag>("calotowersCollection");
  vertices_		     = iConfig.getParameter<InputTag>("vertices"); 
  triggerResults_	     = iConfig.getParameter<InputTag>("triggerResults");
  isSignal_                  = iConfig.getUntrackedParameter<bool>("isSignal");
}


void JPsieeAnalyzerSignal::beginJob(const EventSetup&) { 

  OutputTree = new JPsiTree(fOutFileTreeName_.c_str()); 

  H_deltaR        = new TH1F("H_deltaR",       "H_deltaR",100,0.,1.);
  H_deltaRNotCorr = new TH1F("H_deltaRNotCorr","H_deltaRNotCorr",100,0.,1.);
} 

void JPsieeAnalyzerSignal::endJob() { 

  OutputTree->save(); 

  TFile fOut("outputHisto.root","RECREATE");
  H_deltaRNotCorr->Write();
  H_deltaR       ->Write();
  fOut.Close();
}

void JPsieeAnalyzerSignal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // 1) HLT
  Handle<edm::TriggerResults> HLTR;  
  try { iEvent.getByLabel(triggerResults_, HLTR); }
  catch(cms::Exception& ex){ edm::LogError("ProblemHltTriggserResults") << "Trigger results: " << triggerResults_ << " not found"; }
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
  

  // 2) vertex: if more than 1, take the highest pt; if no vertex found take the beam spot (0,0,0)   
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByType(theBeamSpot);
  edm::Handle<reco::VertexCollection> hVtx;
  iEvent.getByLabel(vertices_, hVtx);
  
  if ( hVtx->size()>0 ){ 
    float theMaxPt  = -999.;
    VertexCollection::const_iterator thisVertex;
    for(thisVertex = hVtx->begin(); thisVertex != hVtx->end(); ++thisVertex){      
      float SumPt = 0.0;
      if((*thisVertex).tracksSize() > 0){
	std::vector<TrackBaseRef >::const_iterator thisTrack;
	for( thisTrack=(*thisVertex).tracks_begin(); thisTrack!=(*thisVertex).tracks_end(); thisTrack++){
	  if((**thisTrack).charge()==-1 || (**thisTrack).charge()==1) SumPt += (**thisTrack).pt();
	}}
      if (SumPt>theMaxPt){ 
	theMaxPt = SumPt; 
	gpVertexPos  = GlobalPoint((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
	xyzVertexPos = math::XYZVector((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
      }}
  }
  else{
    gpVertexPos  = GlobalPoint(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
    xyzVertexPos = math::XYZVector(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  }  
  origin = GlobalPoint(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  // so: origin      -> beam spot
  //     gPVertexPos -> highest pt vertex; beam spot if not found 


  // magnetic field
  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  

  // calo topology
  const CaloTopology *theTopology;
  edm::ESHandle<CaloTopology> topology;
  iSetup.get<CaloTopologyRecord>().get(topology);
  if (topology.isValid()) theTopology = topology.product();
  

  // calo geometry
  const CaloGeometry *theGeometry;
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);
  if(geometry.isValid()) theGeometry = geometry.product();
  

  // generated electrons
  const HepMC::GenEvent *myGenEvent;
  // if (isSignal_){  // chiara, solo x minimum bias
  Handle<edm::HepMCProduct> hepMC;
  iEvent.getByLabel("source", hepMC);
  myGenEvent = hepMC->GetEvent();
  // }


  // reconstructed electrons  
  Handle<GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons); 
  
  // reconstructed superclusters
  Handle<SuperClusterCollection> recoClustersEB;
  Handle<SuperClusterCollection> recoClustersEE;
  iEvent.getByLabel(superclusterCollectionEB_,recoClustersEB); 
  iEvent.getByLabel ("multi5x5SuperClusters", "multi5x5EndcapSuperClusters", recoClustersEE); 
  

  // ECAL rechits to build cluster shapes
  const EcalRecHitCollection *rechitsEB;
  const EcalRecHitCollection *rechitsEE;
  Handle<EcalRecHitCollection> pRechitsEB;
  Handle<EcalRecHitCollection> pRechitsEE;
  iEvent.getByLabel(ecalRechitsCollectionEB_,pRechitsEB); 
  iEvent.getByLabel(ecalRechitsCollectionEE_,pRechitsEE); 
  if (pRechitsEB.isValid()) rechitsEB = pRechitsEB.product();
  if (pRechitsEE.isValid()) rechitsEE = pRechitsEE.product();
  
  // Barrel HCAL hits
  const HBHERecHitCollection *rechitsHBHE;
  Handle<HBHERecHitCollection> hitsHBHE;
  

  // HoE calculation
  mhbhe_=0;
  bool got = iEvent.getByLabel(hcalRechitsCollectionHBHE_,hitsHBHE);  
  if (got) rechitsHBHE = hitsHBHE.product();
  if (got) mhbhe_= new HBHERecHitMetaCollection(*hitsHBHE);


  // all tracker tracks to study curvature in B
  const TrackCollection *theTracks;
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(tracksCollection_, tracks); 
  if(tracks.isValid()) theTracks = tracks.product() ;

  // calotowers, to study HCAL baed isolation
  const CaloTowerCollection *theCaloTowers;
  Handle<CaloTowerCollection> calotowers;
  iEvent.getByLabel(calotowersCollection_, calotowers); 
  if (calotowers.isValid()) theCaloTowers = calotowers.product();

  // jurassic isolation for electrons
  eIsoFromDepsValueMap_ = new isoContainer(3);
  iEvent.getByLabel( "eleIsoFromDepsTk", (*eIsoFromDepsValueMap_)[0] );
  iEvent.getByLabel( "eleIsoFromDepsEcalFromHits", (*eIsoFromDepsValueMap_)[1] );
  iEvent.getByLabel( "eleIsoFromDepsHcalFromHits", (*eIsoFromDepsValueMap_)[2] );
  

  // ---------------------------------------------------------------------
  // generator level informations
  
  // if running on background samples we skip events with true j/psi 
  bool isAJPsi = false;
  if( !isSignal_ ) {
    HepMC::GenEvent::particle_const_iterator mcIter;
    for ( mcIter=myGenEvent->particles_begin(); mcIter != myGenEvent->particles_end(); mcIter++ ) {
      if ( fabs((*mcIter)->pdg_id()) == 11) {
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
    OutputTree->fillRunInfos( theRun, theEvent, theLumiSection );
    
    // counters for the tree
    int numberMcParticle     = 0;
    int allRecoScGt4         = 0;
    int scMatchedPlusNumber  = 0;
    int scMatchedMinusNumber = 0;
    int scMatchedAllNumber   = 0;
    int numberOfElectrons    = 0;
    
    // MC info: electrons from j/psi
    int theMcPc = 0;
    TVector3 trueEle_3P, truePos_3P;
    TLorentzVector trueEle_4P, truePos_4P;
    
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
    
    
    
    

    
    
    // ------------ analysis based on superclusters ----------------------------
    
    // filling an unique collection with all superclusters in EB & EE with Et > 4
    reco::SuperClusterCollection myRecoClusters;   
    SuperClusterCollection::const_iterator scIter;
    
    for (scIter = recoClustersEB->begin() ; scIter != recoClustersEB->end() ; scIter++) { 
      float Et = (scIter->energy())*sin(scIter->position().theta());
      if (Et>=4) myRecoClusters.push_back(*scIter); 
    }
    
    for (scIter = recoClustersEE->begin() ; scIter != recoClustersEE->end() ; scIter++) { 
      float Et = (scIter->energy())*sin(scIter->position().theta());
      if (Et>=4) myRecoClusters.push_back(*scIter); 
    }
    
    // putting them in order
    std::sort( myRecoClusters.begin(),  myRecoClusters.end(), etGreater());
    allRecoScGt4 = myRecoClusters.size();
    




    
    // loop over all superclusters with Et>4 and fill the tree for them
    for (scIter=myRecoClusters.begin(); scIter!=myRecoClusters.end(); scIter++ ) {
      
      float scEnergy  = scIter->energy();
      float scEt      = (scIter->energy())*sin(scIter->position().theta());
      float scEta     = scIter->position().eta();
      float etaScCorr = 0;

      const SuperCluster* scIterRef = &(*scIter);
      const GlobalPoint clusterPos(scIter->position().x(), scIter->position().y(), scIter->position().z());    
            
      // to find the curvature in the magnetic field: looking for the closer track
      int bestDr_q          = 9999; 
      float bestDr          = 9999.;
      float bestDrNotCorr   = 9999.;
      float bestDr_eop	    = 9999.;
      float dEtaForSelected = 9999.;
      float dPhiForSelected = 9999.;
      const reco::Track *bestDr_Track = 0; 

      TrackCollection::const_iterator trIter;      
      for (trIter=theTracks->begin(); trIter!=theTracks->end(); trIter++) {
	
	int trackQ          = trIter->charge();
	float trackPt       = trIter->p()*sin(trIter->theta());
	float ScEoverTrackP = scEnergy/trIter->p();
	
	if (trackPt > 1) {

	  /* old
	  // propagating the track to the supercluster position (preso da RecoEgamma/EgammaElectronAlgos/src/GsfElectronAlgo.cc)
	  ECALPositionCalculator posCalc;
	  float trackEta = posCalc.ecalEta(trIter->innerMomentum(),trIter->innerPosition());
	  float trackPhi = posCalc.ecalPhi(&(*theMagField),trIter->innerMomentum(),trIter->innerPosition(),trackQ);
	  float deltaEta = fabs(scIter->position().eta() - trackEta);   
	  float deltaPhi = fabs(scIter->position().phi() - trackPhi); 
	  if(deltaPhi>6.283185308) deltaPhi -= 6.283185308;
	  if(deltaPhi>3.141592654) deltaPhi = 6.283185308-deltaPhi;
	  float deltaR = sqrt (deltaEta*deltaEta + deltaPhi*deltaPhi);
	  */  // old 
	  

	  // preso da HLTrigger/Egamma/src/HLTElectronDetaDphiFilter.cc
	  const math::XYZVector trackMom = trIter->momentum();
	  math::XYZPoint SCcorrPosition(scIterRef->x()-gpVertexPos.x(), scIterRef->y()-gpVertexPos.y(), scIterRef->z()-gpVertexPos.z());
	  etaScCorr = SCcorrPosition.eta();                                                            // eta sc va corretto per il beam spot / vertex
	  float deltaEta  = fabs(etaScCorr - trIter->eta());                                           // eta traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
                                                                                                       // eta traccia non va propagato al calorimetro tanto non curva 
	  ECALPositionCalculator posCalc;
	  float phiTrCorr = posCalc.ecalPhi(&(*theMagField), trackMom, xyzVertexPos, trackQ);          // phi traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
	                                                                                               // ma phi traccia va propagato al calo e qui serve il constraint del vtx 
	  float deltaPhi  = fabs(scIter->phi() - phiTrCorr);
	  if(deltaPhi>6.283185308) deltaPhi -= 6.283185308;
	  if(deltaPhi>3.141592654) deltaPhi = 6.283185308-deltaPhi;
	  float deltaR = sqrt (deltaEta*deltaEta + deltaPhi*deltaPhi);

	  // to check the agreement between dR with/without correction for beam spot:
	  float deltaEtaNotCorr = fabs(scIter->eta() - trIter->eta());
	  float phiTrNotCorr    = posCalc.ecalPhi(&(*theMagField), trackMom, trIter->innerPosition(), trackQ); 
	  float deltaPhiNotCorr = fabs(scIter->phi() - phiTrNotCorr);
	  if(deltaPhiNotCorr>6.283185308) deltaPhiNotCorr -= 6.283185308;
	  if(deltaPhiNotCorr>3.141592654) deltaPhiNotCorr = 6.283185308-deltaPhiNotCorr;
	  float deltaRNotCorr = sqrt (deltaEtaNotCorr*deltaEtaNotCorr + deltaPhiNotCorr*deltaPhiNotCorr);	  
	  if (deltaRNotCorr < bestDrNotCorr) bestDrNotCorr = deltaRNotCorr;
	  // to check the agreement between dR with/without correction for beam spot; end


	  if (deltaR < bestDr){
	    bestDr          = deltaR;
	    bestDr_q        = trackQ;
	    dEtaForSelected = deltaEta;      // al calorimetro
	    dPhiForSelected = deltaPhi;      // al calorimetro
	    bestDr_eop 	    = ScEoverTrackP;
	    bestDr_Track    = &*trIter;	
	  }
	}
      }

      // to be used on signal only - chiara, per decidere il taglio  
      if (bestDr<1000)        H_deltaR        -> Fill(bestDr);   
      if (bestDrNotCorr<1000) H_deltaRNotCorr -> Fill(bestDrNotCorr);   
      
      bool isMatchingOk = false;
      if (bestDr < 0.3) isMatchingOk = true;
      
      if (!isMatchingOk) OutputTree->fillSuperclusters( -999., -999., -999, -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.);
      
      
      // if a track is found, full analysis
      float s9s25   = 0.;
      float sigmaEE = 0.;
      double HoE    = 0.;

      if (isMatchingOk) { 
	
	// chiara: io userei eta e phi della traccia e basta, ma non so imporre il passaggio per il vertice
	
	// correct phi of the SC using the curvature in the magnetic field --> phi info at vertex
	// qui propaghiamo la posizione del supercluster indietro al vertice, imponendo o meno il passaggio per il vertice/beam spot
	FreeTrajectoryState ftsPV = myftsPV(&(*theMagField), clusterPos, gpVertexPos, scEnergy, -bestDr_q);    
        FreeTrajectoryState ftsOr = myftsOr(&(*theMagField), clusterPos, origin, scEnergy, -bestDr_q);    
	

	// electron ID --------------------------------

	// looking to cluster shape variables
	BasicClusterRef theSeed = scIterRef->seed();
	const EcalRecHitCollection *rechits = 0;
	float seedEta = theSeed->position().eta();      
	if( fabs(seedEta) < 1.479 ) rechits = rechitsEB;
	else rechits = rechitsEE;     
	float e3x3   = EcalClusterTools::e3x3( *theSeed, &(*rechits), theTopology );
	float e5x5   = EcalClusterTools::e5x5( *theSeed, &(*rechits), theTopology );
	s9s25  = e3x3/e5x5;
	std::vector<float> covMatrix = EcalClusterTools::covariances( *theSeed, rechits, theTopology, theGeometry );
	sigmaEE = sqrt(covMatrix[0]); 
	if( fabs(seedEta) >= 1.479 ) sigmaEE = sigmaEE - 0.02*(fabs(seedEta) - 2.3);
	
	// H over E
	HoECalculator calc(theGeometry);
	HoE=calc(&(*scIterRef),mhbhe_);


	// isolation ---------------------------------
		
	/*
	// looking to HCAL deposits 
	hcalIsol=0.;
	HBHERecHitCollection::const_iterator hbheIter;
	for(hbheIter=rechitsHBHE->begin(); hbheIter!=rechitsHBHE->end(); ++hbheIter){  
	  double HcalHit_eta = theGeometry->getPosition(hbheIter->id()).eta();  
	  if(fabs(HcalHit_eta-(scIter->position().eta()))<0.2) {         
	    float HcalHit_pth=hbheIter->energy()*sin(2*atan(exp(-HcalHit_eta)));
	    double HcalHit_phi=theGeometry->getPosition(hbheIter->id()).phi();
	    if(HcalHit_phi<0) HcalHit_phi+=6.283185308;
	    float deltaeta=fabs(HcalHit_eta-(scIter->position().eta()));
	    float deltaphi=fabs(HcalHit_phi-(scIter->position().phi()));
	    if(deltaphi>6.283185308) deltaphi-=6.283185308;
	    if(deltaphi>3.141592654) deltaphi=6.283185308-deltaphi;   
	    float newDelta= (deltaphi*deltaphi+ deltaeta*deltaeta);
	    if(newDelta<0.04) hcalIsol+=HcalHit_pth;        
	  }      
	  } */
	
	
	// tracker based isolation studies   
	TrackerIsolation trackIsolation( &(*bestDr_Track), theTracks);
	trackIsolation.setIntRadius(0.02);
	trackIsolation.setExtRadius(0.3);
	float sumPt03 = trackIsolation.getPtTracks(true);
	trackIsolation.setExtRadius(0.4);
	float sumPt04 = trackIsolation.getPtTracks(true);
	trackIsolation.setExtRadius(0.5);
	float sumPt05 = trackIsolation.getPtTracks(true);

	// calo towers based isolation in ECAL and HCAL   
	CalotowerIsolation calotowerIsolation(scIterRef, theCaloTowers);
	calotowerIsolation.setIntRadius(0.1);
	calotowerIsolation.setExtRadius(0.3);	
	float sumHadEt03 = calotowerIsolation.getEtHcal(true);
	float sumEmEt03  = calotowerIsolation.getEtEcal(true);	
	calotowerIsolation.setExtRadius(0.4);	
	float sumHadEt04 = calotowerIsolation.getEtHcal(true);
	float sumEmEt04  = calotowerIsolation.getEtEcal(true);	
	calotowerIsolation.setExtRadius(0.5);
	float sumHadEt05 = calotowerIsolation.getEtHcal(true);
	float sumEmEt05  = calotowerIsolation.getEtEcal(true);


	// filling the tree with infos on all the superclusters
	if (bestDr_q>0) scMatchedPlusNumber++;
	if (bestDr_q<0) scMatchedMinusNumber++; 
	OutputTree->fillSuperclusters( bestDr, bestDrNotCorr, bestDr_q, ftsPV.momentum().x(), ftsPV.momentum().y(), ftsPV.momentum().z(), ftsPV.momentum().eta(), ftsPV.momentum().phi(), ftsOr.momentum().x(), ftsOr.momentum().y(), ftsOr.momentum().z(), ftsOr.momentum().eta(), ftsOr.momentum().phi(), bestDr_Track->px(), bestDr_Track->py(), bestDr_Track->pz(), bestDr_Track->eta(), bestDr_Track->phi(), scEnergy, scEt, scEta, etaScCorr, s9s25, sigmaEE, HoE, dEtaForSelected, dPhiForSelected, bestDr_eop, sumPt03, sumPt04, sumPt05, sumHadEt03, sumHadEt04, sumHadEt05, sumEmEt03, sumEmEt04, sumEmEt05); 
      }
      
    } // loop over superclusters
    
    // counting the matched superclusters
    scMatchedAllNumber = scMatchedPlusNumber+scMatchedMinusNumber;
    
    



    // --------------------------------------------------------------------
    // analysis based on gsf electrons

    int countMP = 0;
    numberOfElectrons = gsfElectrons->size();
    GsfElectronCollection::const_iterator eleIter;
    for (eleIter=gsfElectrons->begin(); eleIter != gsfElectrons->end() ; eleIter++) { 
      
      SuperClusterRef theScRef   = eleIter->get<SuperClusterRef>();
      BasicClusterRef theSeedRef = theScRef->seed();

      // looking to cluster shape variables
      const EcalRecHitCollection *rechits = 0;
      float seedEta = theSeedRef->position().eta();      
      if( fabs(seedEta) < 1.479 ) rechits = rechitsEB;
      else rechits = rechitsEE;     
      float e3x3   = EcalClusterTools::e3x3( *theSeedRef, &(*rechits), theTopology );
      float e5x5   = EcalClusterTools::e5x5( *theSeedRef, &(*rechits), theTopology );
      float s9s25  = e3x3/e5x5;
      std::vector<float> covMatrix = EcalClusterTools::covariances( *theSeedRef, rechits, theTopology, theGeometry );
      float sigmaEE = sqrt(covMatrix[0]); 
      if( fabs(seedEta) >= 1.479 ) sigmaEE = sigmaEE - 0.02*(fabs(seedEta) - 2.3);
      
      // tracker based isolation studies   
      TrackerIsolation trackIsolation( &(*eleIter->gsfTrack()), theTracks );
      trackIsolation.setIntRadius(0.02);
      trackIsolation.setExtRadius(0.3);
      float sumPt03 = trackIsolation.getPtTracks(true);
      trackIsolation.setExtRadius(0.4);
      float sumPt04 = trackIsolation.getPtTracks(true);
      trackIsolation.setExtRadius(0.5);
      float sumPt05 = trackIsolation.getPtTracks(true);
      
      // calo towers based isolation in ECAL and HCAL   
      CalotowerIsolation calotowerIsolation( &(*theScRef), theCaloTowers);
      calotowerIsolation.setIntRadius(0.1);
      calotowerIsolation.setExtRadius(0.3);	
      float sumHadEt03 = calotowerIsolation.getEtHcal(true);
      float sumEmEt03  = calotowerIsolation.getEtEcal(true);	
      calotowerIsolation.setExtRadius(0.4);	
      float sumHadEt04 = calotowerIsolation.getEtHcal(true);
      float sumEmEt04  = calotowerIsolation.getEtEcal(true);	
      calotowerIsolation.setExtRadius(0.5);
      float sumHadEt05 = calotowerIsolation.getEtHcal(true);
      float sumEmEt05  = calotowerIsolation.getEtEcal(true);

      // jurassic isolation
      const reco::GsfElectronRef eleRef(gsfElectrons,countMP);
      const isoFromDepositsMap & eIsoFromDepsTkVal   = *( (*eIsoFromDepsValueMap_)[0] );
      const isoFromDepositsMap & eIsoFromDepsEcalVal = *( (*eIsoFromDepsValueMap_)[1] );
      const isoFromDepositsMap & eIsoFromDepsHcalVal = *( (*eIsoFromDepsValueMap_)[2] );
      float  jurTrackerEle = eIsoFromDepsTkVal[eleRef];
      float  jurECALEle    = eIsoFromDepsEcalVal[eleRef];
      float  jurHCALEle    = eIsoFromDepsHcalVal[eleRef];

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
      float derre     = sqrt(deta*deta + dphi*dphi);
      float eleEoP    = eleIter->eSuperClusterOverP();
      OutputTree->fillElectrons( derre, eleCharge, elePx, elePy, elePz, eleEta, elePhi, eleEnergy, eleEt, s9s25, sigmaEE, HoE, deta, dphi, eleEoP, sumPt03, sumPt04, sumPt05, sumHadEt03, sumHadEt04, sumHadEt05, sumEmEt03, sumEmEt04, sumEmEt05, jurTrackerEle, jurECALEle, jurHCALEle);

      countMP++;
      
    } // loop over electrons

   
    
    // filling the tree: some summary numbers
    int signal = 0;
    if (isSignal_) signal = 1;
    OutputTree->fillGeneral(signal, hVtx->size(), numberMcParticle, allRecoScGt4, scMatchedPlusNumber, scMatchedMinusNumber, scMatchedAllNumber, numberOfElectrons, intHlt29, gpVertexPos.x(), gpVertexPos.y(), gpVertexPos.z());
    
    OutputTree->store();

    // deleting
    delete eIsoFromDepsValueMap_;
    
  } // segnale o vero fondo


}


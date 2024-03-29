#ifndef EMEnrichingFilterAlgo_h
#define EMEnrichingFilterAlgo_h

/** \class EMEnrichingFilter
 *
 *  EMEnrichingFilter 
 *
 * \author J Lamb, UCSB
 *
 ************************************************************/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class EMEnrichingFilterAlgo {
 public:
  EMEnrichingFilterAlgo(const edm::ParameterSet&);
  ~EMEnrichingFilterAlgo();
  
  bool filter(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  bool hasBCAncestors(reco::GenParticle gp);

 private:
  int filterPhotonElectronSeed(float clusterthreshold,
			       float seedthreshold,
			       float isoConeSize,
			       float hOverEMax,
			       float tkIsoMax,
			       float caloIsoMax,
			       bool requiretrackmatch,
			       const std::vector<reco::GenParticle> &genPars,
			       const std::vector<reco::GenParticle> &genParsCurved);



  std::vector<reco::GenParticle> applyBFieldCurv(const std::vector<reco::GenParticle> &genPars, const edm::EventSetup& iSetup);
  int filterIsoGenPar(float etmin, float conesize,const reco::GenParticleCollection &gph,
		      const reco::GenParticleCollection &gphCurved);
  float deltaRxyAtEE(const reco::GenParticle &gp1, const reco::GenParticle &gp2);

  bool isBCHadron(reco::GenParticle gp);
  bool isBCMeson(reco::GenParticle gp);
  bool isBCBaryon(reco::GenParticle gp);
    
		       
 private:

  //constants
  float FILTER_TKISOCUT_;
  float FILTER_CALOISOCUT_;
  float FILTER_ETA_MIN_;
  float FILTER_ETA_MAX_;
  float ECALBARRELMAXETA_;
  float ECALBARRELRADIUS_;
  float ECALENDCAPZ_;

  // parameters of the filter
  float isoGenParETMin_;
  float isoGenParConeSize_;
  float clusterThreshold_;
  float seedThreshold_;
  float isoConeSize_;
  float hOverEMax_;
  float tkIsoMax_;
  float caloIsoMax_;
  float eTThreshold_;    // from bctoe
  bool requireTrackMatch_;
  edm::InputTag genParSource_;
  
  // for double em object
  std::vector<reco::GenParticle> sel1seeds;
  std::vector<reco::GenParticle> sel2seeds;
  std::vector<reco::GenParticle> selBCtoEseeds;

};
#endif

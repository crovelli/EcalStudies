#ifndef JPsieeAnalyzerSignal_h
#define JPsieeAnalyzerSignal_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/FTSFromVertexToPointFactory.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include <map>
#include <cmath>
#include <string>

#include "JPsiTree.h"

class TFile;
class TH1D;
class TH2D;
class TTree;

using namespace edm;
using namespace std;


class JPsieeAnalyzerSignal : public edm::EDAnalyzer 
{
 public:
   
 explicit JPsieeAnalyzerSignal(const edm::ParameterSet& conf);
   
 virtual void beginJob(edm::EventSetup const& iSetup);
 virtual void endJob();
 virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

 private:

 // collections
 edm::InputTag electronCollection_;
 edm::InputTag superclusterCollectionEB_;
 edm::InputTag superclusterCollectionEE_;
 edm::InputTag ecalRechitsCollectionEB_; 
 edm::InputTag ecalRechitsCollectionEE_; 
 edm::InputTag hcalRechitsCollectionHBHE_;
 edm::InputTag tracksCollection_;
 edm::InputTag calotowersCollection_;
 edm::InputTag vertices_;
 edm::InputTag triggerResults_;

 bool isSignal_;  
 
 FTSFromVertexToPointFactory myftsPV;
 FTSFromVertexToPointFactory myftsOr;
 HBHERecHitMetaCollection *mhbhe_;

 math::XYZPoint xyzVertexPos;
 GlobalPoint gpVertexPos;
 GlobalPoint origin;
 
 edm::TriggerNames triggerNames_;
 int intHlt29, intHlt31;

 typedef edm::ValueMap<double> isoFromDepositsMap;
 typedef std::vector< edm::Handle<isoFromDepositsMap> > isoContainer;
 isoContainer *eIsoFromDepsValueMap_;

 // tree
 std::string fOutFileTreeName_;
 JPsiTree* OutputTree;

 TH1F *H_deltaR, *H_deltaRNotCorr;
};

#endif

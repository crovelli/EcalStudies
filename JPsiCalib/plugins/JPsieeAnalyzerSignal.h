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
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"

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
 edm::InputTag ecalRechitsCollectionEB_; 
 edm::InputTag ecalRechitsCollectionEE_; 
 edm::InputTag triggerResults_;

 bool isSignal_;  
 bool isUpsiAnalysis_;  
 
 edm::TriggerNames triggerNames_;
 int intHltJPsi, intHltUpsilon, intHltBoth;
 
 typedef edm::ValueMap<float> eleIdMap;
 typedef std::vector< edm::Handle<eleIdMap> > eleIdContainer;
 eleIdContainer *eleIdResults_;

 // tree
 std::string fOutFileTreeName_;
 JPsiTree* OutputTree;
};

#endif

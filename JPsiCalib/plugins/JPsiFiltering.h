#ifndef JPsiFiltering_h
#define JPsiFiltering_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TMath.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include <map>
#include <cmath>
#include <string>

#include "FilteringTree.h"

class TFile;
class TTree;

using namespace edm;
using namespace std;


class JPsiFiltering : public edm::EDAnalyzer 
{
 public:
  
  explicit JPsiFiltering(const edm::ParameterSet& conf);
  
  virtual void beginJob(edm::EventSetup const& iSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
 private:
  
  // collections
  edm::InputTag electronCollection_;
  edm::InputTag triggerResults_;

  typedef edm::ValueMap<float> eleIdMap;
  typedef std::vector< edm::Handle<eleIdMap> > eleIdContainer;
  eleIdContainer *eleIdResults_;

  edm::TriggerNames triggerNames_;
  
  // tree
  std::string fOutFileTreeName_;
  FilteringTree* OutputTree;
};

#endif

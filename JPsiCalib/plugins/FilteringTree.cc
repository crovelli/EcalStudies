#include "FilteringTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"


FilteringTree::FilteringTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","filter studies tree");
  
  // initialization
  myNumberOfGen=0;
  myNumberOfEle=0;
  
  // numbers
  myTree->Branch("numberOfGenerated",  &numberOfGenerated,  "numberOfGenerated/I");
  myTree->Branch("numberOfElectrons",  &numberOfElectrons,  "numberOfElectrons/I");
  myTree->Branch("passedFilter1",      &passedFilter1,      "passedFilter1/I");
  myTree->Branch("passedFilter2",      &passedFilter2,      "passedFilter2/I");
  myTree->Branch("passedFilter3",      &passedFilter3,      "passedFilter3/I");
  myTree->Branch("passedFilter4",      &passedFilter4,      "passedFilter4/I");
  myTree->Branch("passedFilter5",      &passedFilter5,      "passedFilter5/I");
  myTree->Branch("passedFilter6",      &passedFilter6,      "passedFilter6/I");
  myTree->Branch("passedFilter7",      &passedFilter7,      "passedFilter7/I");
  myTree->Branch("passedFilter8",      &passedFilter8,      "passedFilter8/I");

  // run infos
  myTree->Branch("runNumber",     &runNumber,     "runNumber/I");
  myTree->Branch("eventNumber",   &eventNumber,   "eventNumber/I");
  myTree->Branch("lumiSection",   &lumiSection,   "lumiSection/I");

  // generated electrons
  myTree->Branch("pxGen",       &pxGen,       "pxGen[numberOfGenerated]/F");
  myTree->Branch("pyGen",       &pyGen,       "pyGen[numberOfGenerated]/F");
  myTree->Branch("pzGen",       &pzGen,       "pzGen[numberOfGenerated]/F");
  myTree->Branch("eneGen",      &eneGen,      "eneGen[numberOfGenerated]/F");
  myTree->Branch("etaGen",      &etaGen,      "etaGen[numberOfGenerated]/F");
  myTree->Branch("phiGen",      &phiGen,      "phiGen[numberOfGenerated]/F");
  myTree->Branch("idGen",       &idGen,       "idGen[numberOfGenerated]/I");
  myTree->Branch("statusGen",   &statusGen,   "statusGen[numberOfGenerated]/I");
  myTree->Branch("motherIdGen", &motherIdGen, "motherIdGen[numberOfGenerated]/I");
  myTree->Branch("motherGen",   &motherGen,   "motherGen[numberOfGenerated]/I");

  // electrons
  myTree->Branch("chargeEle",               &chargeEle,              "chargeEle[numberOfElectrons]/I");
  myTree->Branch("xRecoEle",                &xRecoEle,               "xRecoEle[numberOfElectrons]/F");
  myTree->Branch("yRecoEle",                &yRecoEle,               "yRecoEle[numberOfElectrons]/F");
  myTree->Branch("zRecoEle",                &zRecoEle,               "zRecoEle[numberOfElectrons]/F");
  myTree->Branch("etaRecoEle",              &etaRecoEle,             "etaRecoEle[numberOfElectrons]/F");
  myTree->Branch("phiRecoEle",              &phiRecoEle,             "phiRecoEle[numberOfElectrons]/F");
  myTree->Branch("eneRecoEle",              &eneRecoEle,             "eneRecoEle[numberOfElectrons]/F");
  myTree->Branch("etRecoEle",               &etRecoEle,              "etRecoEle[numberOfElectrons]/F");
  myTree->Branch("HoverERecoEle",           &HoverERecoEle,          "HoverERecoEle[numberOfElectrons]/F");  
  myTree->Branch("dEtaWithTrackerRecoEle",  &dEtaWithTrackerRecoEle, "dEtaWithTrackerRecoEle[numberOfElectrons]/F");
  myTree->Branch("dPhiWithTrackerRecoEle",  &dPhiWithTrackerRecoEle, "dPhiWithTrackerRecoEle[numberOfElectrons]/F"); 
  myTree->Branch("EoverPRecoEle",  	    &EoverPRecoEle, 	     "EoverPRecoEle[numberOfElectrons]/F");  
  myTree->Branch("trkIsolRecoEle_03", 	    &trkIsolRecoEle03,	     "trkIsolRecoEle_03[numberOfElectrons]/F");
  myTree->Branch("trkIsolRecoEle_04", 	    &trkIsolRecoEle04,	     "trkIsolRecoEle_04[numberOfElectrons]/F");
  myTree->Branch("trkIsolRecoEle_05", 	    &trkIsolRecoEle05,	     "trkIsolRecoEle_05[numberOfElectrons]/F");
  myTree->Branch("jurTrkIsolEle", 	    &jurTrkIsolEle,	     "jurTrkIsolEle[numberOfElectrons]/F");
  myTree->Branch("jurEmIsolEle", 	    &jurEmIsolEle,	     "jurEmIsolEle[numberOfElectrons]/F");
  myTree->Branch("jurHadIsolEle", 	    &jurHadIsolEle,	     "jurHadIsolEle[numberOfElectrons]/F");
  myTree->Branch("eleIdLoose", 	            &eleIdLoose,	     "eleIdLoose[numberOfElectrons]/I");
  myTree->Branch("eleIdRobLoose", 	    &eleIdRobLoose,	     "eleIdRobLoose[numberOfElectrons]/I");
  myTree->Branch("eleIdRobTight", 	    &eleIdRobTight,	     "eleIdRobTight[numberOfElectrons]/I");
  myTree->Branch("eleIdTight", 	            &eleIdTight,	     "eleIdTight[numberOfElectrons]/I");
}

FilteringTree::~FilteringTree() {
  
  myFile->cd();
  myTree->Write();
  myFile->Close();
  
  delete myFile;
}

void FilteringTree::store() {

  myTree->Fill();

  myNumberOfGen = 0;
  myNumberOfEle = 0;
}

void FilteringTree::save() {
  
  myFile->cd();
  myTree->Write();
  myFile->Close();
}


void FilteringTree::fillRunInfos( int run, int event, int lumis) {
  
  runNumber   = run;
  eventNumber = event;
  lumiSection = lumis;
}

void FilteringTree::fillGeneral( int nmc, int ea, int f1, int f2, int f3, int f4, int f5, int f6, int f7, int f8 ) {
  
  numberOfGenerated = nmc;
  numberOfElectrons = ea;
  passedFilter1     = f1;
  passedFilter2     = f2;
  passedFilter3     = f3;
  passedFilter4     = f4;
  passedFilter5     = f5;
  passedFilter6     = f6;
  passedFilter7     = f7;
  passedFilter8     = f8;
}

void FilteringTree::fillGenerated( float gen_x, float gen_y, float gen_z, float gen_ene, float gen_eta, float gen_phi, int gen_id, int gen_status, int gen_mid, int gen_mother ){

  if (myNumberOfGen < NGENMAX) {
    pxGen[myNumberOfGen]       = gen_x;
    pyGen[myNumberOfGen]       = gen_y;
    pzGen[myNumberOfGen]       = gen_z;
    eneGen[myNumberOfGen]      = gen_ene;
    etaGen[myNumberOfGen]      = gen_eta;
    phiGen[myNumberOfGen]      = gen_phi;
    idGen[myNumberOfGen]       = gen_id;
    statusGen[myNumberOfGen]   = gen_status;
    motherIdGen[myNumberOfGen] = gen_mid;
    motherGen[myNumberOfGen]   = gen_mother;
    myNumberOfGen++;
  }
}

void FilteringTree::fillElectrons( int ele_q, float ele_x, float ele_y, float ele_z, float ele_eta, float ele_phi, float ele_ene, float ele_et, float ele_hoe, float ele_deta, float ele_dphi, float ele_eop, float ele_trk03, float ele_trk04, float ele_trk05, float ele_jtrk, float ele_jem, float ele_jhad, int loose, int robloose, int robtight, int tight ){       
  
  if (myNumberOfEle < NMAX) {
    chargeEle[myNumberOfEle]              = ele_q;
    xRecoEle[myNumberOfEle]               = ele_x;
    yRecoEle[myNumberOfEle]               = ele_y;
    zRecoEle[myNumberOfEle]               = ele_z;
    etaRecoEle[myNumberOfEle]             = ele_eta;
    phiRecoEle[myNumberOfEle]             = ele_phi;
    eneRecoEle[myNumberOfEle]             = ele_ene;
    etRecoEle[myNumberOfEle]              = ele_et;
    HoverERecoEle[myNumberOfEle]          = ele_hoe;
    dEtaWithTrackerRecoEle[myNumberOfEle] = ele_deta;
    dPhiWithTrackerRecoEle[myNumberOfEle] = ele_dphi;
    EoverPRecoEle[myNumberOfEle]	  = ele_eop;
    trkIsolRecoEle03[myNumberOfEle] 	  = ele_trk03;
    trkIsolRecoEle04[myNumberOfEle] 	  = ele_trk04;
    trkIsolRecoEle05[myNumberOfEle] 	  = ele_trk05;
    jurTrkIsolEle[myNumberOfEle] 	  = ele_jtrk;
    jurEmIsolEle[myNumberOfEle] 	  = ele_jem;
    jurHadIsolEle[myNumberOfEle] 	  = ele_jhad;
    eleIdLoose[myNumberOfEle] 	          = loose;
    eleIdRobLoose[myNumberOfEle]          = robloose;
    eleIdRobTight[myNumberOfEle] 	  = robtight;
    eleIdTight[myNumberOfEle] 	          = tight;
    myNumberOfEle++;
  }
}



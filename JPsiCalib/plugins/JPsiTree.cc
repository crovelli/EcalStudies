#include "JPsiTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"


JPsiTree::JPsiTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","J/Psi tree");

  // initialization
  myNumberOfGen=0;
  myNumberOfSc=0;
  myNumberOfEle=0;

  // numbers
  myTree->Branch("signal",                  &signal,                 "signal/I");
  myTree->Branch("numberOfPrimaryVtx",      &numberOfPrimaryVtx,     "numberOfPrimaryVtx/I");
  myTree->Branch("numberOfGenerated",       &numberOfGenerated,      "numberOfGenerated/I");
  myTree->Branch("numberOfScGt4",           &numberOfScGt4,          "numberOfScGt4/I");
  myTree->Branch("numberOfScMatchedPlus",   &numberOfScMatchedPlus,  "numberOfScMatchedPlus/I");
  myTree->Branch("numberOfScMatchedMinus",  &numberOfScMatchedMinus, "numberOfScMatchedMinus/I");  
  myTree->Branch("numberOfScMatchedAll",    &numberOfScMatchedAll,   "numberOfScMatchedAll/I");  
  myTree->Branch("numberOfElectrons",       &numberOfElectrons,      "numberOfElectrons/I");
  myTree->Branch("hlt29",                   &hlt29,                  "hlt29/I");  
  myTree->Branch("vertexX",                 &vertexX,                "vertexX/F");  
  myTree->Branch("vertexY",                 &vertexY,                "vertexY/F");  
  myTree->Branch("vertexZ",                 &vertexZ,                "vertexZ/F");  

  // run infos
  myTree->Branch("runNumber",     &runNumber,     "runNumber/I");
  myTree->Branch("eventNumber",   &eventNumber,   "eventNumber/I");
  myTree->Branch("lumiSection",   &lumiSection,   "lumiSection/I");

  // generated electrons
  myTree->Branch("chargeGenEle",  &chargeGenEle, "chargeEl[numberOfGenerated]/I");
  myTree->Branch("pxGenEle",      &pxGenEle,     "pxGenEle[numberOfGenerated]/F");
  myTree->Branch("pyGenEle",      &pyGenEle,     "pyGenEle[numberOfGenerated]/F");
  myTree->Branch("pzGenEle",      &pzGenEle,     "pzGenEle[numberOfGenerated]/F");
  myTree->Branch("eneGenEle",     &eneGenEle,    "eneGenEle[numberOfGenerated]/F");

  // electrons
  myTree->Branch("dRWithTrackerRecoEle",    &dRWithTrackerRecoEle,   "dRWithTrackerRecoEle[numberOfElectrons]/F");
  myTree->Branch("chargeEle",               &chargeEle,              "chargeEle[numberOfElectrons]/I");
  myTree->Branch("xRecoEle",                &xRecoEle,               "xRecoEle[numberOfElectrons]/F");
  myTree->Branch("yRecoEle",                &yRecoEle,               "yRecoEle[numberOfElectrons]/F");
  myTree->Branch("zRecoEle",                &zRecoEle,               "zRecoEle[numberOfElectrons]/F");
  myTree->Branch("etaRecoEle",              &etaRecoEle,             "etaRecoEle[numberOfElectrons]/F");
  myTree->Branch("phiRecoEle",              &phiRecoEle,             "phiRecoEle[numberOfElectrons]/F");
  myTree->Branch("eneRecoEle",              &eneRecoEle,             "eneRecoEle[numberOfElectrons]/F");
  myTree->Branch("etRecoEle",               &etRecoEle,              "etRecoEle[numberOfElectrons]/F");
  myTree->Branch("e9e25RecoEle",            &e9e25RecoEle,           "e9e25RecoEle[numberOfElectrons]/F");
  myTree->Branch("sigmaEtaEtaRecoEle",      &sigmaEtaEtaRecoEle,     "sigmaEtaEtaRecoEle[numberOfElectrons]/F");
  myTree->Branch("HoverERecoEle",           &HoverERecoEle,          "HoverERecoEle[numberOfElectrons]/F");  
  myTree->Branch("dEtaWithTrackerRecoEle",  &dEtaWithTrackerRecoEle, "dEtaWithTrackerRecoEle[numberOfElectrons]/F");
  myTree->Branch("dPhiWithTrackerRecoEle",  &dPhiWithTrackerRecoEle, "dPhiWithTrackerRecoEle[numberOfElectrons]/F"); 
  myTree->Branch("EoverPRecoEle",  	    &EoverPRecoEle, 	     "EoverPRecoEle[numberOfElectrons]/F");  
  myTree->Branch("trkIsolRecoEle_03", 	    &trkIsolRecoEle03,	     "trkIsolRecoEle_03[numberOfElectrons]/F");
  myTree->Branch("trkIsolRecoEle_04", 	    &trkIsolRecoEle04,	     "trkIsolRecoEle_04[numberOfElectrons]/F");
  myTree->Branch("trkIsolRecoEle_05", 	    &trkIsolRecoEle05,	     "trkIsolRecoEle_05[numberOfElectrons]/F");
  myTree->Branch("hadIsolRecoEle_03", 	    &hadIsolRecoEle03,	     "hadIsolRecoEle_03[numberOfElectrons]/F");
  myTree->Branch("hadIsolRecoEle_04", 	    &hadIsolRecoEle04,	     "hadIsolRecoEle_04[numberOfElectrons]/F");
  myTree->Branch("hadIsolRecoEle_05", 	    &hadIsolRecoEle05,	     "hadIsolRecoEle_05[numberOfElectrons]/F");
  myTree->Branch("emIsolRecoEle_03", 	    &emIsolRecoEle03,	     "emIsolRecoEle_03[numberOfElectrons]/F");
  myTree->Branch("emIsolRecoEle_04", 	    &emIsolRecoEle04,	     "emIsolRecoEle_04[numberOfElectrons]/F");
  myTree->Branch("emIsolRecoEle_05", 	    &emIsolRecoEle05,	     "emIsolRecoEle_05[numberOfElectrons]/F");
  myTree->Branch("jurTrkIsolEle", 	    &jurTrkIsolEle,	     "jurTrkIsolEle[numberOfElectrons]/F");
  myTree->Branch("jurEmIsolEle", 	    &jurEmIsolEle,	     "jurEmIsolEle[numberOfElectrons]/F");
  myTree->Branch("jurHadIsolEle", 	    &jurHadIsolEle,	     "jurHadIsolEle[numberOfElectrons]/F");

  // superclusters
  myTree->Branch("dRWithTrackerRecoSc",        &dRWithTrackerRecoSc,        "dRWithTrackerRecoSc[numberOfScGt4]/F");
  myTree->Branch("dRWithTrackerRecoNotCorrSc", &dRWithTrackerRecoNotCorrSc, "dRWithTrackerRecoNotCorrSc[numberOfScGt4]/F");
  myTree->Branch("chargeSc",                   &chargeSc,                   "chargeSc[numberOfScGt4]/I");
  myTree->Branch("xPVRecoSc",                  &xPVRecoSc,                  "xPVRecoSc[numberOfScGt4]/F");
  myTree->Branch("yPVRecoSc",                  &yPVRecoSc,                  "yPVRecoSc[numberOfScGt4]/F");
  myTree->Branch("zPVRecoSc",                  &zPVRecoSc,                  "zPVRecoSc[numberOfScGt4]/F");
  myTree->Branch("etaPVRecoSc",                &etaPVRecoSc,                "etaPVRecoSc[numberOfScGt4]/F");
  myTree->Branch("phiPVRecoSc",                &phiPVRecoSc,                "phiPVRecoSc[numberOfScGt4]/F");
  myTree->Branch("xOrRecoSc",                  &xOrRecoSc,                  "xOrRecoSc[numberOfScGt4]/F");
  myTree->Branch("yOrRecoSc",                  &yOrRecoSc,                  "yOrRecoSc[numberOfScGt4]/F");
  myTree->Branch("zOrRecoSc",                  &zOrRecoSc,                  "zOrRecoSc[numberOfScGt4]/F");
  myTree->Branch("etaOrRecoSc",                &etaOrRecoSc,                "etaOrRecoSc[numberOfScGt4]/F");
  myTree->Branch("phiOrRecoSc",                &phiOrRecoSc,                "phiOrRecoSc[numberOfScGt4]/F");
  myTree->Branch("pxFromTrackerSc",            &pxFromTrackerSc,            "pxFromTrackerSc[numberOfScGt4]/F");
  myTree->Branch("pyFromTrackerSc",            &pyFromTrackerSc,            "pyFromTrackerSc[numberOfScGt4]/F");
  myTree->Branch("pzFromTrackerSc",            &pzFromTrackerSc,            "pzFromTrackerSc[numberOfScGt4]/F");
  myTree->Branch("etaFromTrackerSc",           &etaFromTrackerSc,           "etaFromTrackerSc[numberOfScGt4]/F");
  myTree->Branch("phiFromTrackerSc",           &phiFromTrackerSc,           "phiFromTrackerSc[numberOfScGt4]/F");
  myTree->Branch("eneRecoSc",                  &eneRecoSc,                  "eneRecoSc[numberOfScGt4]/F");
  myTree->Branch("etRecoSc",                   &etRecoSc,                   "etRecoSc[numberOfScGt4]/F");
  myTree->Branch("etaRecoSc",                  &etaRecoSc,                  "etaRecoSc[numberOfScGt4]/F");
  myTree->Branch("etaCorrRecoSc",              &etaCorrRecoSc,              "etaCorrRecoSc[numberOfScGt4]/F");
  myTree->Branch("e9e25RecoSc",                &e9e25RecoSc,                "e9e25RecoSc[numberOfScGt4]/F");
  myTree->Branch("sigmaEtaEtaRecoSc",          &sigmaEtaEtaRecoSc,          "sigmaEtaEtaRecoSc[numberOfScGt4]/F");
  myTree->Branch("HoverERecoSc",               &HoverERecoSc,               "HoverERecoSc[numberOfScGt4]/F");  
  myTree->Branch("dEtaWithTrackerRecoSc",      &dEtaWithTrackerRecoSc,      "dEtaWithTrackerRecoSc[numberOfScGt4]/F");
  myTree->Branch("dPhiWithTrackerRecoSc",      &dPhiWithTrackerRecoSc,      "dPhiWithTrackerRecoSc[numberOfScGt4]/F");
  myTree->Branch("EoverPRecoSc",     	       &EoverPRecoSc,	            "EoverPRecoSc[numberOfScGt4]/F");
  myTree->Branch("trkIsolRecoSc_03", 	       &trkIsolRecoSc03, 	    "trkIsolRecoSc_03[numberOfScGt4]/F");
  myTree->Branch("trkIsolRecoSc_04", 	       &trkIsolRecoSc04,	    "trkIsolRecoSc_04[numberOfScGt4]/F");
  myTree->Branch("trkIsolRecoSc_05", 	       &trkIsolRecoSc05,	    "trkIsolRecoSc_05[numberOfScGt4]/F");
  myTree->Branch("hadIsolRecoSc_03",           &hadIsolRecoSc03,	    "hadIsolRecoSc_03[numberOfScGt4]/F");
  myTree->Branch("hadIsolRecoSc_04", 	       &hadIsolRecoSc04,	    "hadIsolRecoSc_04[numberOfScGt4]/F");
  myTree->Branch("hadIsolRecoSc_05", 	       &hadIsolRecoSc05,	    "hadIsolRecoSc_05[numberOfScGt4]/F");
  myTree->Branch("emIsolRecoSc_03", 	       &emIsolRecoSc03,	            "emIsolRecoSc_03[numberOfScGt4]/F");
  myTree->Branch("emIsolRecoSc_04", 	       &emIsolRecoSc04,	            "emIsolRecoSc_04[numberOfScGt4]/F");
  myTree->Branch("emIsolRecoSc_05", 	       &emIsolRecoSc05,	            "emIsolRecoSc_05[numberOfScGt4]/F");
}

JPsiTree::~JPsiTree() {
  
  myFile->cd();
  myTree->Write();
  myFile->Close();
  
  delete myFile;
}

void JPsiTree::store() {

  myTree->Fill();

  myNumberOfGen = 0;
  myNumberOfSc  = 0;
  myNumberOfEle = 0;
}

void JPsiTree::save() {

  myFile->cd();
  myTree->Write();
  myFile->Close();
}


void JPsiTree::fillRunInfos( int run, int event, int lumis) {
  
  runNumber   = run;
  eventNumber = event;
  lumiSection = lumis;
}

void JPsiTree::fillGeneral( int sgn, int nvtx, int nmc, int ngt, int mp, int mm, int ma, int ea, int tr29, float pvX, float pvY, float pvZ ) {

  signal                 = sgn;			    
  numberOfPrimaryVtx     = nvtx;
  numberOfGenerated      = nmc;
  numberOfScGt4          = ngt;
  numberOfScMatchedPlus  = mp;
  numberOfScMatchedMinus = mm;
  numberOfScMatchedAll   = ma;
  numberOfElectrons      = ea;
  hlt29                  = tr29;
  vertexX                = pvX;
  vertexY                = pvY;
  vertexZ                = pvZ;
}

void JPsiTree::fillGenerated( int gen_q, float gen_x, float gen_y, float gen_z, float gen_ene ){

  chargeGenEle[myNumberOfGen] = gen_q;
  pxGenEle[myNumberOfGen]     = gen_x;
  pyGenEle[myNumberOfGen]     = gen_y;
  pzGenEle[myNumberOfGen]     = gen_z;
  eneGenEle[myNumberOfGen]    = gen_ene;
  myNumberOfGen++;
}

void JPsiTree::fillSuperclusters(float sc_dr, float sc_drnc, int sc_q, float sc_pvx, float sc_pvy, float sc_pvz, float sc_pveta, float sc_pvphi, float sc_orx, float sc_ory, float sc_orz, float sc_oreta, float sc_orphi, float sc_trpx, float sc_trpy, float sc_trpz, float sc_treta, float sc_trphi, float sc_ene, float sc_et, float sc_eta, float sc_etac, float sc_s9s25, float sc_see, float sc_hoe, float sc_deta, float sc_dphi, float sc_eop, float sc_trk03, float sc_trk04, float sc_trk05, float sc_had03, float sc_had04, float sc_had05, float sc_em03, float sc_em04, float sc_em05) {
  
  if (myNumberOfSc < NMAX) {

    dRWithTrackerRecoSc[myNumberOfSc]        = sc_dr;
    dRWithTrackerRecoNotCorrSc[myNumberOfSc] = sc_drnc;
    chargeSc[myNumberOfSc]                   = sc_q;
    xPVRecoSc[myNumberOfSc]                  = sc_pvx;
    yPVRecoSc[myNumberOfSc]                  = sc_pvy;
    zPVRecoSc[myNumberOfSc]                  = sc_pvz;
    etaPVRecoSc[myNumberOfSc]                = sc_pveta;
    phiPVRecoSc[myNumberOfSc]                = sc_pvphi;
    xOrRecoSc[myNumberOfSc]                  = sc_orx;
    yOrRecoSc[myNumberOfSc]                  = sc_ory;
    zOrRecoSc[myNumberOfSc]                  = sc_orz;
    etaOrRecoSc[myNumberOfSc]                = sc_oreta;
    phiOrRecoSc[myNumberOfSc]                = sc_orphi;
    pxFromTrackerSc[myNumberOfSc]            = sc_trpx;
    pyFromTrackerSc[myNumberOfSc]            = sc_trpy;
    pzFromTrackerSc[myNumberOfSc]            = sc_trpz;
    etaFromTrackerSc[myNumberOfSc]           = sc_treta;
    phiFromTrackerSc[myNumberOfSc]           = sc_trphi;
    eneRecoSc[myNumberOfSc]                  = sc_ene;
    etRecoSc[myNumberOfSc]                   = sc_et;
    etaRecoSc[myNumberOfSc]                  = sc_eta;
    etaCorrRecoSc[myNumberOfSc]              = sc_etac;
    e9e25RecoSc[myNumberOfSc]                = sc_s9s25;
    sigmaEtaEtaRecoSc[myNumberOfSc]          = sc_see;
    HoverERecoSc[myNumberOfSc]               = sc_hoe;
    dEtaWithTrackerRecoSc[myNumberOfSc]      = sc_deta;
    dPhiWithTrackerRecoSc[myNumberOfSc]      = sc_dphi;
    EoverPRecoSc[myNumberOfSc]		     = sc_eop;
    trkIsolRecoSc03[myNumberOfSc]	     = sc_trk03;
    trkIsolRecoSc04[myNumberOfSc]	     = sc_trk04;
    trkIsolRecoSc05[myNumberOfSc]  	     = sc_trk05;
    hadIsolRecoSc03[myNumberOfSc]	     = sc_had03;
    hadIsolRecoSc04[myNumberOfSc]	     = sc_had04;
    hadIsolRecoSc05[myNumberOfSc]	     = sc_had05;
    emIsolRecoSc03[myNumberOfSc]	     = sc_em03;
    emIsolRecoSc04[myNumberOfSc]	     = sc_em04;
    emIsolRecoSc05[myNumberOfSc]	     = sc_em05;
    myNumberOfSc++;
  }
}


void JPsiTree::fillElectrons( float ele_dr, int ele_q, float ele_x, float ele_y, float ele_z, float ele_eta, float ele_phi, float ele_ene, float ele_et, float ele_s9s25, float ele_see, float ele_hoe, float ele_deta, float ele_dphi, float ele_eop, float ele_trk03, float ele_trk04, float ele_trk05, float ele_had03, float ele_had04, float ele_had05, float ele_em03, float ele_em04, float ele_em05, float ele_jtrk, float ele_jem, float ele_jhad){       

  if (myNumberOfEle < NMAX) {
    dRWithTrackerRecoEle[myNumberOfEle]   = ele_dr;
    chargeEle[myNumberOfEle]              = ele_q;
    xRecoEle[myNumberOfEle]               = ele_x;
    yRecoEle[myNumberOfEle]               = ele_y;
    zRecoEle[myNumberOfEle]               = ele_z;
    etaRecoEle[myNumberOfEle]             = ele_eta;
    phiRecoEle[myNumberOfEle]             = ele_phi;
    eneRecoEle[myNumberOfEle]             = ele_ene;
    etRecoEle[myNumberOfEle]              = ele_et;
    e9e25RecoEle[myNumberOfEle]           = ele_s9s25;
    sigmaEtaEtaRecoEle[myNumberOfEle]     = ele_see;
    HoverERecoEle[myNumberOfEle]          = ele_hoe;
    dEtaWithTrackerRecoEle[myNumberOfEle] = ele_deta;
    dPhiWithTrackerRecoEle[myNumberOfEle] = ele_dphi;
    EoverPRecoEle[myNumberOfEle]	  = ele_eop;
    trkIsolRecoEle03[myNumberOfEle] 	  = ele_trk03;
    trkIsolRecoEle04[myNumberOfEle] 	  = ele_trk04;
    trkIsolRecoEle05[myNumberOfEle] 	  = ele_trk05;
    hadIsolRecoEle03[myNumberOfEle] 	  = ele_had03;
    hadIsolRecoEle04[myNumberOfEle] 	  = ele_had04;
    hadIsolRecoEle05[myNumberOfEle] 	  = ele_had05;
    emIsolRecoEle03[myNumberOfEle] 	  = ele_em03;
    emIsolRecoEle04[myNumberOfEle] 	  = ele_em04;
    emIsolRecoEle05[myNumberOfEle] 	  = ele_em05;
    jurTrkIsolEle[myNumberOfEle] 	  = ele_jtrk;
    jurEmIsolEle[myNumberOfEle] 	  = ele_jem;
    jurHadIsolEle[myNumberOfEle] 	  = ele_jhad;
    myNumberOfEle++;
  }
}



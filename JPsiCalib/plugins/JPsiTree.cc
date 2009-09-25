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
  myTree = new TTree("T1","lowEnergyElectronsTree");

  // initialization
  myNumberOfGen=0;
  myNumberOfEle=0;

  // run infos
  myTree->Branch("runNumber",     &runNumber,     "runNumber/I");
  myTree->Branch("eventNumber",   &eventNumber,   "eventNumber/I");
  myTree->Branch("lumiSection",   &lumiSection,   "lumiSection/I");
  myTree->Branch("pthat",         &pthat,         "pthat/F");

  // summary numbers
  myTree->Branch("signal",               &signal,              "signal/I");
  myTree->Branch("numberOfGenerated",    &numberOfGenerated,   "numberOfGenerated/I");
  myTree->Branch("numberOfElectrons",    &numberOfElectrons,   "numberOfElectrons/I");
  myTree->Branch("hltJpsi",              &hltJpsi,             "hltJpsi/I");  
  myTree->Branch("hltUpsilon",           &hltUpsilon,          "hltUpsilon/I");  
  myTree->Branch("hltBoth",              &hltBoth,             "hltBoth/I");  

  // generated electrons
  myTree->Branch("chargeGenEle",  &chargeGenEle, "chargeEl[numberOfGenerated]/I");
  myTree->Branch("pxGenEle",      &pxGenEle,     "pxGenEle[numberOfGenerated]/F");
  myTree->Branch("pyGenEle",      &pyGenEle,     "pyGenEle[numberOfGenerated]/F");
  myTree->Branch("pzGenEle",      &pzGenEle,     "pzGenEle[numberOfGenerated]/F");
  myTree->Branch("eneGenEle",     &eneGenEle,    "eneGenEle[numberOfGenerated]/F");

  // reconstructed electrons
  myTree->Branch("e9RecoEle",                 &e9RecoEle,                 "e9RecoEle[numberOfElectrons]/F");
  myTree->Branch("e25RecoEle",                &e25RecoEle,                "e25RecoEle[numberOfElectrons]/F");
  myTree->Branch("e9e25RecoEle",              &e9e25RecoEle,              "e9e25RecoEle[numberOfElectrons]/F");
  myTree->Branch("sigmaEtaEtaRecoEle",        &sigmaEtaEtaRecoEle,        "sigmaEtaEtaRecoEle[numberOfElectrons]/F");
  myTree->Branch("sigmaIetaIetaRecoEle",      &sigmaIetaIetaRecoEle,      "sigmaIetaIetaRecoEle[numberOfElectrons]/F");
  myTree->Branch("HoverERecoEle",             &HoverERecoEle,             "HoverERecoEle[numberOfElectrons]/F");  
  myTree->Branch("eleIdLooseRecoEle",         &eleIdLooseRecoEle,         "eleIdLooseRecoEle[numberOfElectrons]/I");  
  myTree->Branch("eleIdRobLooseRecoEle",      &eleIdRobLooseRecoEle,      "eleIdRobLooseRecoEle[numberOfElectrons]/I");  
  myTree->Branch("eleIdRobTightRecoEle",      &eleIdRobTightRecoEle,      "eleIdRobTightRecoEle[numberOfElectrons]/I");  
  myTree->Branch("eleIdTightRecoEle",         &eleIdTightRecoEle,         "eleIdTightRecoEle[numberOfElectrons]/I");  
  myTree->Branch("pFlowMvaRecoEle",           &pFlowMvaRecoEle,           "pFlowMvaRecoEle[numberOfElectrons]/F");  
  myTree->Branch("EoverPRecoEle",    	      &EoverPRecoEle, 	          "EoverPRecoEle[numberOfElectrons]/F");  
  myTree->Branch("EoverPoutRecoEle",   	      &EoverPoutRecoEle,          "EoverPoutRecoEle[numberOfElectrons]/F");  
  myTree->Branch("dEtaAtVtxRecoEle",          &dEtaAtVtxRecoEle,          "dEtaAtVtxRecoEle[numberOfElectrons]/F");
  myTree->Branch("dEtaAtCaloRecoEle",         &dEtaAtCaloRecoEle,         "dEtaAtCaloRecoEle[numberOfElectrons]/F");
  myTree->Branch("dPhiAtVtxRecoEle",          &dPhiAtVtxRecoEle,          "dPhiAtVtxRecoEle[numberOfElectrons]/F");
  myTree->Branch("dPhiAtCaloRecoEle",         &dPhiAtCaloRecoEle,         "dPhiAtCaloRecoEle[numberOfElectrons]/F");
  myTree->Branch("fBremRecoEle",              &fBremRecoEle,              "fBremRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr03TkSumPtRecoEle",        &dr03TkSumPtRecoEle,        "dr03TkSumPtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr04TkSumPtRecoEle",        &dr04TkSumPtRecoEle,        "dr04TkSumPtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr03EcalSumEtRecoEle",      &dr03EcalSumEtRecoEle,      "dr03EcalSumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr04EcalSumEtRecoEle",      &dr04EcalSumEtRecoEle,      "dr04EcalSumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr03Hcal1SumEtRecoEle",     &dr03Hcal1SumEtRecoEle,     "dr03Hcal1SumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr04Hcal1SumEtRecoEle",     &dr04Hcal1SumEtRecoEle,     "dr04Hcal1SumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr03Hcal2SumEtRecoEle",     &dr03Hcal2SumEtRecoEle,     "dr03Hcal2SumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("dr04Hcal2SumEtRecoEle",     &dr04Hcal2SumEtRecoEle,     "dr04Hcal2SumEtRecoEle[numberOfElectrons]/F");
  myTree->Branch("rawSCenergyRecoEle",        &rawSCenergyRecoEle,        "rawSCenergyRecoEle[numberOfElectrons]/F");
  myTree->Branch("rawESenergyRecoEle",        &rawESenergyRecoEle,        "rawESenergyRecoEle[numberOfElectrons]/F");
  myTree->Branch("ecalEnergyRecoEle",         &ecalEnergyRecoEle,         "ecalEnergyRecoEle[numberOfElectrons]/F");
  myTree->Branch("ecalEnergyErrorRecoEle",    &ecalEnergyErrorRecoEle,    "ecalEnergyErrorRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackPxRecoEle",            &trackPxRecoEle,            "trackPxRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackPyRecoEle",            &trackPyRecoEle,            "trackPyRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackPzRecoEle",            &trackPzRecoEle,            "trackPzRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackEtaRecoEle",           &trackEtaRecoEle,           "trackEtaRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackPhiRecoEle",           &trackPhiRecoEle,           "trackPhiRecoEle[numberOfElectrons]/F");
  myTree->Branch("trackMomentumErrorRecoEle", &trackMomentumErrorRecoEle, "trackMomentumErrorRecoEle[numberOfElectrons]/F");
  myTree->Branch("isEcalEneCorrectedRecoEle", &isEcalEneCorrectedRecoEle, "isEcalEneCorrectedRecoEle[numberOfElectrons]/I");
  myTree->Branch("isPCorrectedRecoEle",       &isPCorrectedRecoEle,       "isPCorrectedRecoEle[numberOfElectrons]/I");
  myTree->Branch("isEcalDrivenRecoEle",       &isEcalDrivenRecoEle,       "isEcalDrivenRecoEle[numberOfElectrons]/I");
  myTree->Branch("isPFlowRecoEle",            &isPFlowRecoEle,            "isPFlowRecoEle[numberOfElectrons]/I");
  myTree->Branch("classRecoEle",              &classRecoEle,              "classRecoEle[numberOfElectrons]/I");
  myTree->Branch("chargeRecoEle",             &chargeRecoEle,             "chargeRecoEle[numberOfElectrons]/I");
  myTree->Branch("pxRecoEle",                 &pxRecoEle,                 "pxRecoEle[numberOfElectrons]/F");
  myTree->Branch("pyRecoEle",                 &pyRecoEle,                 "pyRecoEle[numberOfElectrons]/F");
  myTree->Branch("pzRecoEle",                 &pzRecoEle,                 "pzRecoEle[numberOfElectrons]/F");
  myTree->Branch("etaRecoEle",                &etaRecoEle,                "etaRecoEle[numberOfElectrons]/F");
  myTree->Branch("phiRecoEle",                &phiRecoEle,                "phiRecoEle[numberOfElectrons]/F");
  myTree->Branch("eneRecoEle",                &eneRecoEle,                "eneRecoEle[numberOfElectrons]/F");
  myTree->Branch("etRecoEle",                 &etRecoEle,                 "etRecoEle[numberOfElectrons]/F");
  myTree->Branch("momentumErrorRecoEle",      &momentumErrorRecoEle,      "momentumErrorRecoEle[numberOfElectrons]/F");
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
  myNumberOfEle = 0;
}

void JPsiTree::save() {

  myFile->cd();
  myTree->Write();
  myFile->Close();
}


void JPsiTree::fillRunInfos( int run, int event, int lumis, float pth) {
  
  runNumber   = run;
  eventNumber = event;
  lumiSection = lumis;
  pthat       = pth;
}

void JPsiTree::fillGeneral( int sgn, int nmc, int nreco, int trJ, int trU, int trB ){

  signal            = sgn;			    
  numberOfGenerated = nmc;
  numberOfElectrons = nreco;
  hltJpsi           = trJ;
  hltUpsilon        = trU;
  hltBoth           = trB;
}

void JPsiTree::fillGenerated( int gen_q, float gen_x, float gen_y, float gen_z, float gen_ene ){

  chargeGenEle[myNumberOfGen] = gen_q;
  pxGenEle[myNumberOfGen]     = gen_x;
  pyGenEle[myNumberOfGen]     = gen_y;
  pzGenEle[myNumberOfGen]     = gen_z;
  eneGenEle[myNumberOfGen]    = gen_ene;
  myNumberOfGen++;
}

void JPsiTree::fillElectrons( float e3x3, float e5x5, float s9s25, float sigmaEtaEta, float sigmaIetaIeta, float hcalOverEcal, int eleIdLoose, int eleIdRobLoose, int eleIdRobTight, int eleIdTight, float pfMva, float eSuperClusterOverP, float eSeedClusterOverPout, float deltaEtaSuperClusterAtVtx, float deltaEtaSeedClusterAtCalo, float deltaPhiSuperClusterAtVtx, float deltaPhiSeedClusterAtCalo, float fbrem, float dr03TkSumPt, float dr04TkSumPt, float dr03EcalRecHitSumEt, float dr04EcalRecHitSumEt, float dr03HcalDepth1TowerSumEt, float dr04HcalDepth1TowerSumEt, float dr03HcalDepth2TowerSumEt, float dr04HcalDepth2TowerSumEt, float rawSCenergy, float rawESenergy, float ecalEnergy, float ecalEnergyError, float trackPx, float trackPy, float trackPz, float trackEta, float trackPhi, float trackMomentumError, int isEcalEnergyCorrected, int isMomentumCorrected, int isEcalDriven, int isParticleFlow, int eleClass, int eleCharge, float elePx, float elePy, float elePz, float eleEta, float elePhi, float eleEnergy, float eleEt, float electronMomentumError) {

  if (myNumberOfEle < NMAX) {

    e9RecoEle[myNumberOfEle]                 = e3x3;
    e25RecoEle[myNumberOfEle]                = e5x5;
    e9e25RecoEle[myNumberOfEle]              = s9s25;
    sigmaEtaEtaRecoEle[myNumberOfEle]        = sigmaEtaEta;
    sigmaIetaIetaRecoEle[myNumberOfEle]      = sigmaIetaIeta;
    HoverERecoEle[myNumberOfEle]             = hcalOverEcal; 
    eleIdLooseRecoEle[myNumberOfEle]         = eleIdLoose;  
    eleIdRobLooseRecoEle[myNumberOfEle]      = eleIdRobLoose;
    eleIdRobTightRecoEle[myNumberOfEle]      = eleIdRobTight;
    eleIdTightRecoEle[myNumberOfEle]         = eleIdTight;
    pFlowMvaRecoEle[myNumberOfEle]           = pfMva;
    EoverPRecoEle[myNumberOfEle]             = eSuperClusterOverP;
    EoverPoutRecoEle[myNumberOfEle]          = eSeedClusterOverPout;  
    dEtaAtVtxRecoEle[myNumberOfEle]          = deltaEtaSuperClusterAtVtx;
    dEtaAtCaloRecoEle[myNumberOfEle]         = deltaEtaSeedClusterAtCalo; 
    dPhiAtVtxRecoEle[myNumberOfEle]          = deltaPhiSuperClusterAtVtx;
    dPhiAtCaloRecoEle[myNumberOfEle]         = deltaPhiSeedClusterAtCalo;
    fBremRecoEle[myNumberOfEle]              = fbrem; 
    dr03TkSumPtRecoEle[myNumberOfEle]        = dr03TkSumPt;
    dr04TkSumPtRecoEle[myNumberOfEle]        = dr04TkSumPt;
    dr03EcalSumEtRecoEle[myNumberOfEle]      = dr03EcalRecHitSumEt;
    dr04EcalSumEtRecoEle[myNumberOfEle]      = dr04EcalRecHitSumEt;
    dr03Hcal1SumEtRecoEle[myNumberOfEle]     = dr03HcalDepth1TowerSumEt;
    dr04Hcal1SumEtRecoEle[myNumberOfEle]     = dr04HcalDepth1TowerSumEt;
    dr03Hcal2SumEtRecoEle[myNumberOfEle]     = dr03HcalDepth2TowerSumEt;
    dr04Hcal2SumEtRecoEle[myNumberOfEle]     = dr04HcalDepth2TowerSumEt;
    rawSCenergyRecoEle[myNumberOfEle]        = rawSCenergy;
    rawESenergyRecoEle[myNumberOfEle]        = rawESenergy;
    ecalEnergyRecoEle[myNumberOfEle]         = ecalEnergy;
    ecalEnergyErrorRecoEle[myNumberOfEle]    = ecalEnergyError;
    trackPxRecoEle[myNumberOfEle]            = trackPx;
    trackPyRecoEle[myNumberOfEle]            = trackPy;
    trackPzRecoEle[myNumberOfEle]            = trackPz;
    trackEtaRecoEle[myNumberOfEle]           = trackEta;
    trackPhiRecoEle[myNumberOfEle]           = trackPhi;
    trackMomentumErrorRecoEle[myNumberOfEle] = trackMomentumError;
    isEcalEneCorrectedRecoEle[myNumberOfEle] = isEcalEnergyCorrected;
    isPCorrectedRecoEle[myNumberOfEle]       = isMomentumCorrected;
    isEcalDrivenRecoEle[myNumberOfEle]       = isEcalDriven;
    isPFlowRecoEle[myNumberOfEle]            = isParticleFlow;
    classRecoEle[myNumberOfEle]              = eleClass;
    chargeRecoEle[myNumberOfEle]             = eleCharge;
    pxRecoEle[myNumberOfEle]                 = elePx;
    pyRecoEle[myNumberOfEle]                 = elePy;
    pzRecoEle[myNumberOfEle]                 = elePz;
    etaRecoEle[myNumberOfEle]                = eleEta;
    phiRecoEle[myNumberOfEle]                = elePhi;
    eneRecoEle[myNumberOfEle]                = eleEnergy;
    etRecoEle[myNumberOfEle]                 = eleEt;
    momentumErrorRecoEle[myNumberOfEle]      = electronMomentumError;
    myNumberOfEle++;
  }
}



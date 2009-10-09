//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 28 11:12:17 2009 by ROOT version 5.22/00a
// from TTree T1/lowEnergyElectronsTree
// found on file: /tmp/arcidiac/signalJPsi10TeV_1.root
//////////////////////////////////////////////////////////

#ifndef JPsiTreeBase_h
#define JPsiTreeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class JPsiTreeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiSection;
   Float_t         pthat;
   Int_t           signal;
   Int_t           numberOfGenerated;
   Int_t           numberOfElectrons;
   Int_t           hltJpsi;
   Int_t           hltUpsilon;
   Int_t           hltBoth;
   Int_t           chargeGenEle[50];   //[numberOfGenerated]
   Float_t         pxGenEle[50];   //[numberOfGenerated]
   Float_t         pyGenEle[50];   //[numberOfGenerated]
   Float_t         pzGenEle[50];   //[numberOfGenerated]
   Float_t         eneGenEle[50];   //[numberOfGenerated]
   Float_t         e9RecoEle[50];   //[numberOfElectrons]
   Float_t         e25RecoEle[50];   //[numberOfElectrons]
   Float_t         e9e25RecoEle[50];   //[numberOfElectrons]
   Float_t         sigmaEtaEtaRecoEle[50];   //[numberOfElectrons]
   Float_t         sigmaIetaIetaRecoEle[50];   //[numberOfElectrons]
   Float_t         HoverERecoEle[50];   //[numberOfElectrons]
   Int_t           eleIdLooseRecoEle[50];   //[numberOfElectrons]
   Int_t           eleIdRobLooseRecoEle[50];   //[numberOfElectrons]
   Int_t           eleIdRobTightRecoEle[50];   //[numberOfElectrons]
   Int_t           eleIdTightRecoEle[50];   //[numberOfElectrons]
   Float_t         pFlowMvaRecoEle[50];   //[numberOfElectrons]
   Float_t         EoverPRecoEle[50];   //[numberOfElectrons]
   Float_t         EoverPoutRecoEle[50];   //[numberOfElectrons]
   Float_t         dEtaAtVtxRecoEle[50];   //[numberOfElectrons]
   Float_t         dEtaAtCaloRecoEle[50];   //[numberOfElectrons]
   Float_t         dPhiAtVtxRecoEle[50];   //[numberOfElectrons]
   Float_t         dPhiAtCaloRecoEle[50];   //[numberOfElectrons]
   Float_t         fBremRecoEle[50];   //[numberOfElectrons]
   Float_t         dr03TkSumPtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr04TkSumPtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr03EcalSumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr04EcalSumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr03Hcal1SumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr04Hcal1SumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr03Hcal2SumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         dr04Hcal2SumEtRecoEle[50];   //[numberOfElectrons]
   Float_t         rawSCenergyRecoEle[50];   //[numberOfElectrons]
   Float_t         rawESenergyRecoEle[50];   //[numberOfElectrons]
   Float_t         ecalEnergyRecoEle[50];   //[numberOfElectrons]
   Float_t         ecalEnergyErrorRecoEle[50];   //[numberOfElectrons]
   Float_t         trackPxRecoEle[50];   //[numberOfElectrons]
   Float_t         trackPyRecoEle[50];   //[numberOfElectrons]
   Float_t         trackPzRecoEle[50];   //[numberOfElectrons]
   Float_t         trackEtaRecoEle[50];   //[numberOfElectrons]
   Float_t         trackPhiRecoEle[50];   //[numberOfElectrons]
   Float_t         trackMomentumErrorRecoEle[50];   //[numberOfElectrons]
   Int_t           isEcalEneCorrectedRecoEle[50];   //[numberOfElectrons]
   Int_t           isPCorrectedRecoEle[50];   //[numberOfElectrons]
   Int_t           isEcalDrivenRecoEle[50];   //[numberOfElectrons]
   Int_t           isPFlowRecoEle[50];   //[numberOfElectrons]
   Int_t           classRecoEle[50];   //[numberOfElectrons]
   Int_t           chargeRecoEle[50];   //[numberOfElectrons]
   Float_t         pxRecoEle[50];   //[numberOfElectrons]
   Float_t         pyRecoEle[50];   //[numberOfElectrons]
   Float_t         pzRecoEle[50];   //[numberOfElectrons]
   Float_t         etaRecoEle[50];   //[numberOfElectrons]
   Float_t         phiRecoEle[50];   //[numberOfElectrons]
   Float_t         eneRecoEle[50];   //[numberOfElectrons]
   Float_t         etRecoEle[50];   //[numberOfElectrons]
   Float_t         momentumErrorRecoEle[50];   //[numberOfElectrons]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_signal;   //!
   TBranch        *b_numberOfGenerated;   //!
   TBranch        *b_numberOfElectrons;   //!
   TBranch        *b_hltJpsi;   //!
   TBranch        *b_hltUpsilon;   //!
   TBranch        *b_hltBoth;   //!
   TBranch        *b_chargeGenEle;   //!
   TBranch        *b_pxGenEle;   //!
   TBranch        *b_pyGenEle;   //!
   TBranch        *b_pzGenEle;   //!
   TBranch        *b_eneGenEle;   //!
   TBranch        *b_e9RecoEle;   //!
   TBranch        *b_e25RecoEle;   //!
   TBranch        *b_e9e25RecoEle;   //!
   TBranch        *b_sigmaEtaEtaRecoEle;   //!
   TBranch        *b_sigmaIetaIetaRecoEle;   //!
   TBranch        *b_HoverERecoEle;   //!
   TBranch        *b_eleIdLooseRecoEle;   //!
   TBranch        *b_eleIdRobLooseRecoEle;   //!
   TBranch        *b_eleIdRobTightRecoEle;   //!
   TBranch        *b_eleIdTightRecoEle;   //!
   TBranch        *b_pFlowMvaRecoEle;   //!
   TBranch        *b_EoverPRecoEle;   //!
   TBranch        *b_EoverPoutRecoEle;   //!
   TBranch        *b_dEtaAtVtxRecoEle;   //!
   TBranch        *b_dEtaAtCaloRecoEle;   //!
   TBranch        *b_dPhiAtVtxRecoEle;   //!
   TBranch        *b_dPhiAtCaloRecoEle;   //!
   TBranch        *b_fBremRecoEle;   //!
   TBranch        *b_dr03TkSumPtRecoEle;   //!
   TBranch        *b_dr04TkSumPtRecoEle;   //!
   TBranch        *b_dr03EcalSumEtRecoEle;   //!
   TBranch        *b_dr04EcalSumEtRecoEle;   //!
   TBranch        *b_dr03Hcal1SumEtRecoEle;   //!
   TBranch        *b_dr04Hcal1SumEtRecoEle;   //!
   TBranch        *b_dr03Hcal2SumEtRecoEle;   //!
   TBranch        *b_dr04Hcal2SumEtRecoEle;   //!
   TBranch        *b_rawSCenergyRecoEle;   //!
   TBranch        *b_rawESenergyRecoEle;   //!
   TBranch        *b_ecalEnergyRecoEle;   //!
   TBranch        *b_ecalEnergyErrorRecoEle;   //!
   TBranch        *b_trackPxRecoEle;   //!
   TBranch        *b_trackPyRecoEle;   //!
   TBranch        *b_trackPzRecoEle;   //!
   TBranch        *b_trackEtaRecoEle;   //!
   TBranch        *b_trackPhiRecoEle;   //!
   TBranch        *b_trackMomentumErrorRecoEle;   //!
   TBranch        *b_isEcalEneCorrectedRecoEle;   //!
   TBranch        *b_isPCorrectedRecoEle;   //!
   TBranch        *b_isEcalDrivenRecoEle;   //!
   TBranch        *b_isPFlowRecoEle;   //!
   TBranch        *b_classRecoEle;   //!
   TBranch        *b_chargeRecoEle;   //!
   TBranch        *b_pxRecoEle;   //!
   TBranch        *b_pyRecoEle;   //!
   TBranch        *b_pzRecoEle;   //!
   TBranch        *b_etaRecoEle;   //!
   TBranch        *b_phiRecoEle;   //!
   TBranch        *b_eneRecoEle;   //!
   TBranch        *b_etRecoEle;   //!
   TBranch        *b_momentumErrorRecoEle;   //!

   JPsiTreeBase(TTree *tree=0);
   virtual ~JPsiTreeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JPsiTreeBase_cxx
JPsiTreeBase::JPsiTreeBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/arcidiac/signalJPsi10TeV_1.root");
      if (!f) {
         f = new TFile("/tmp/arcidiac/signalJPsi10TeV_1.root");
      }
      tree = (TTree*)gDirectory->Get("T1");

   }
   Init(tree);
}

JPsiTreeBase::~JPsiTreeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JPsiTreeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JPsiTreeBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JPsiTreeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("signal", &signal, &b_signal);
   fChain->SetBranchAddress("numberOfGenerated", &numberOfGenerated, &b_numberOfGenerated);
   fChain->SetBranchAddress("numberOfElectrons", &numberOfElectrons, &b_numberOfElectrons);
   fChain->SetBranchAddress("hltJpsi", &hltJpsi, &b_hltJpsi);
   fChain->SetBranchAddress("hltUpsilon", &hltUpsilon, &b_hltUpsilon);
   fChain->SetBranchAddress("hltBoth", &hltBoth, &b_hltBoth);
   fChain->SetBranchAddress("chargeGenEle", chargeGenEle, &b_chargeGenEle);
   fChain->SetBranchAddress("pxGenEle", pxGenEle, &b_pxGenEle);
   fChain->SetBranchAddress("pyGenEle", pyGenEle, &b_pyGenEle);
   fChain->SetBranchAddress("pzGenEle", pzGenEle, &b_pzGenEle);
   fChain->SetBranchAddress("eneGenEle", eneGenEle, &b_eneGenEle);
   fChain->SetBranchAddress("e9RecoEle", e9RecoEle, &b_e9RecoEle);
   fChain->SetBranchAddress("e25RecoEle", e25RecoEle, &b_e25RecoEle);
   fChain->SetBranchAddress("e9e25RecoEle", e9e25RecoEle, &b_e9e25RecoEle);
   fChain->SetBranchAddress("sigmaEtaEtaRecoEle", sigmaEtaEtaRecoEle, &b_sigmaEtaEtaRecoEle);
   fChain->SetBranchAddress("sigmaIetaIetaRecoEle", sigmaIetaIetaRecoEle, &b_sigmaIetaIetaRecoEle);
   fChain->SetBranchAddress("HoverERecoEle", HoverERecoEle, &b_HoverERecoEle);
   fChain->SetBranchAddress("eleIdLooseRecoEle", eleIdLooseRecoEle, &b_eleIdLooseRecoEle);
   fChain->SetBranchAddress("eleIdRobLooseRecoEle", eleIdRobLooseRecoEle, &b_eleIdRobLooseRecoEle);
   fChain->SetBranchAddress("eleIdRobTightRecoEle", eleIdRobTightRecoEle, &b_eleIdRobTightRecoEle);
   fChain->SetBranchAddress("eleIdTightRecoEle", eleIdTightRecoEle, &b_eleIdTightRecoEle);
   fChain->SetBranchAddress("pFlowMvaRecoEle", pFlowMvaRecoEle, &b_pFlowMvaRecoEle);
   fChain->SetBranchAddress("EoverPRecoEle", EoverPRecoEle, &b_EoverPRecoEle);
   fChain->SetBranchAddress("EoverPoutRecoEle", EoverPoutRecoEle, &b_EoverPoutRecoEle);
   fChain->SetBranchAddress("dEtaAtVtxRecoEle", dEtaAtVtxRecoEle, &b_dEtaAtVtxRecoEle);
   fChain->SetBranchAddress("dEtaAtCaloRecoEle", dEtaAtCaloRecoEle, &b_dEtaAtCaloRecoEle);
   fChain->SetBranchAddress("dPhiAtVtxRecoEle", dPhiAtVtxRecoEle, &b_dPhiAtVtxRecoEle);
   fChain->SetBranchAddress("dPhiAtCaloRecoEle", dPhiAtCaloRecoEle, &b_dPhiAtCaloRecoEle);
   fChain->SetBranchAddress("fBremRecoEle", fBremRecoEle, &b_fBremRecoEle);
   fChain->SetBranchAddress("dr03TkSumPtRecoEle", dr03TkSumPtRecoEle, &b_dr03TkSumPtRecoEle);
   fChain->SetBranchAddress("dr04TkSumPtRecoEle", dr04TkSumPtRecoEle, &b_dr04TkSumPtRecoEle);
   fChain->SetBranchAddress("dr03EcalSumEtRecoEle", dr03EcalSumEtRecoEle, &b_dr03EcalSumEtRecoEle);
   fChain->SetBranchAddress("dr04EcalSumEtRecoEle", dr04EcalSumEtRecoEle, &b_dr04EcalSumEtRecoEle);
   fChain->SetBranchAddress("dr03Hcal1SumEtRecoEle", dr03Hcal1SumEtRecoEle, &b_dr03Hcal1SumEtRecoEle);
   fChain->SetBranchAddress("dr04Hcal1SumEtRecoEle", dr04Hcal1SumEtRecoEle, &b_dr04Hcal1SumEtRecoEle);
   fChain->SetBranchAddress("dr03Hcal2SumEtRecoEle", dr03Hcal2SumEtRecoEle, &b_dr03Hcal2SumEtRecoEle);
   fChain->SetBranchAddress("dr04Hcal2SumEtRecoEle", dr04Hcal2SumEtRecoEle, &b_dr04Hcal2SumEtRecoEle);
   fChain->SetBranchAddress("rawSCenergyRecoEle", rawSCenergyRecoEle, &b_rawSCenergyRecoEle);
   fChain->SetBranchAddress("rawESenergyRecoEle", rawESenergyRecoEle, &b_rawESenergyRecoEle);
   fChain->SetBranchAddress("ecalEnergyRecoEle", ecalEnergyRecoEle, &b_ecalEnergyRecoEle);
   fChain->SetBranchAddress("ecalEnergyErrorRecoEle", ecalEnergyErrorRecoEle, &b_ecalEnergyErrorRecoEle);
   fChain->SetBranchAddress("trackPxRecoEle", trackPxRecoEle, &b_trackPxRecoEle);
   fChain->SetBranchAddress("trackPyRecoEle", trackPyRecoEle, &b_trackPyRecoEle);
   fChain->SetBranchAddress("trackPzRecoEle", trackPzRecoEle, &b_trackPzRecoEle);
   fChain->SetBranchAddress("trackEtaRecoEle", trackEtaRecoEle, &b_trackEtaRecoEle);
   fChain->SetBranchAddress("trackPhiRecoEle", trackPhiRecoEle, &b_trackPhiRecoEle);
   fChain->SetBranchAddress("trackMomentumErrorRecoEle", trackMomentumErrorRecoEle, &b_trackMomentumErrorRecoEle);
   fChain->SetBranchAddress("isEcalEneCorrectedRecoEle", isEcalEneCorrectedRecoEle, &b_isEcalEneCorrectedRecoEle);
   fChain->SetBranchAddress("isPCorrectedRecoEle", isPCorrectedRecoEle, &b_isPCorrectedRecoEle);
   fChain->SetBranchAddress("isEcalDrivenRecoEle", isEcalDrivenRecoEle, &b_isEcalDrivenRecoEle);
   fChain->SetBranchAddress("isPFlowRecoEle", isPFlowRecoEle, &b_isPFlowRecoEle);
   fChain->SetBranchAddress("classRecoEle", classRecoEle, &b_classRecoEle);
   fChain->SetBranchAddress("chargeRecoEle", chargeRecoEle, &b_chargeRecoEle);
   fChain->SetBranchAddress("pxRecoEle", pxRecoEle, &b_pxRecoEle);
   fChain->SetBranchAddress("pyRecoEle", pyRecoEle, &b_pyRecoEle);
   fChain->SetBranchAddress("pzRecoEle", pzRecoEle, &b_pzRecoEle);
   fChain->SetBranchAddress("etaRecoEle", etaRecoEle, &b_etaRecoEle);
   fChain->SetBranchAddress("phiRecoEle", phiRecoEle, &b_phiRecoEle);
   fChain->SetBranchAddress("eneRecoEle", eneRecoEle, &b_eneRecoEle);
   fChain->SetBranchAddress("etRecoEle", etRecoEle, &b_etRecoEle);
   fChain->SetBranchAddress("momentumErrorRecoEle", momentumErrorRecoEle, &b_momentumErrorRecoEle);
   Notify();
}

Bool_t JPsiTreeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JPsiTreeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JPsiTreeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JPsiTreeBase_cxx

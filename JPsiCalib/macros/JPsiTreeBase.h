//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 16 15:22:14 2009 by ROOT version 5.18/00a
// from TTree T1/J/Psi tree
// found on file: /u1/crovelli/dataJPSI/analisi_conHLT/signalPrompt/signalPrompt_notHLTskimmed_mergedTree.root
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
   Int_t           signal;
   Int_t           numberOfPrimaryVtx;
   Int_t           numberOfGenerated;
   Int_t           numberOfScGt4;
   Int_t           numberOfScMatchedPlus;
   Int_t           numberOfScMatchedMinus;
   Int_t           numberOfScMatchedAll;
   Int_t           numberOfElectrons;
   Int_t           hlt29;
   Float_t         vertexX;
   Float_t         vertexY;
   Float_t         vertexZ;
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiSection;
   Int_t           chargeGenEle[4];   //[numberOfGenerated]
   Float_t         pxGenEle[4];   //[numberOfGenerated]
   Float_t         pyGenEle[4];   //[numberOfGenerated]
   Float_t         pzGenEle[4];   //[numberOfGenerated]
   Float_t         eneGenEle[4];   //[numberOfGenerated]
   Float_t         dRWithTrackerRecoEle[5];   //[numberOfElectrons]
   Int_t           chargeEle[5];   //[numberOfElectrons]
   Float_t         xRecoEle[5];   //[numberOfElectrons]
   Float_t         yRecoEle[5];   //[numberOfElectrons]
   Float_t         zRecoEle[5];   //[numberOfElectrons]
   Float_t         etaRecoEle[5];   //[numberOfElectrons]
   Float_t         phiRecoEle[5];   //[numberOfElectrons]
   Float_t         eneRecoEle[5];   //[numberOfElectrons]
   Float_t         etRecoEle[5];   //[numberOfElectrons]
   Float_t         e9e25RecoEle[5];   //[numberOfElectrons]
   Float_t         sigmaEtaEtaRecoEle[5];   //[numberOfElectrons]
   Float_t         HoverERecoEle[5];   //[numberOfElectrons]
   Float_t         dEtaWithTrackerRecoEle[5];   //[numberOfElectrons]
   Float_t         dPhiWithTrackerRecoEle[5];   //[numberOfElectrons]
   Float_t         EoverPRecoEle[5];   //[numberOfElectrons]
   Float_t         trkIsolRecoEle_03[5];   //[numberOfElectrons]
   Float_t         trkIsolRecoEle_04[5];   //[numberOfElectrons]
   Float_t         trkIsolRecoEle_05[5];   //[numberOfElectrons]
   Float_t         hadIsolRecoEle_03[5];   //[numberOfElectrons]
   Float_t         hadIsolRecoEle_04[5];   //[numberOfElectrons]
   Float_t         hadIsolRecoEle_05[5];   //[numberOfElectrons]
   Float_t         emIsolRecoEle_03[5];   //[numberOfElectrons]
   Float_t         emIsolRecoEle_04[5];   //[numberOfElectrons]
   Float_t         emIsolRecoEle_05[5];   //[numberOfElectrons]
   Float_t         jurTrkIsolEle[5];   //[numberOfElectrons]
   Float_t         jurEmIsolEle[5];   //[numberOfElectrons]
   Float_t         jurHadIsolEle[5];   //[numberOfElectrons]
   Float_t         dRWithTrackerRecoSc[11];   //[numberOfScGt4]
   Float_t         dRWithTrackerRecoNotCorrSc[11];   //[numberOfScGt4]
   Int_t           chargeSc[11];   //[numberOfScGt4]
   Float_t         xPVRecoSc[11];   //[numberOfScGt4]
   Float_t         yPVRecoSc[11];   //[numberOfScGt4]
   Float_t         zPVRecoSc[11];   //[numberOfScGt4]
   Float_t         etaPVRecoSc[11];   //[numberOfScGt4]
   Float_t         phiPVRecoSc[11];   //[numberOfScGt4]
   Float_t         xOrRecoSc[11];   //[numberOfScGt4]
   Float_t         yOrRecoSc[11];   //[numberOfScGt4]
   Float_t         zOrRecoSc[11];   //[numberOfScGt4]
   Float_t         etaOrRecoSc[11];   //[numberOfScGt4]
   Float_t         phiOrRecoSc[11];   //[numberOfScGt4]
   Float_t         pxFromTrackerSc[11];   //[numberOfScGt4]
   Float_t         pyFromTrackerSc[11];   //[numberOfScGt4]
   Float_t         pzFromTrackerSc[11];   //[numberOfScGt4]
   Float_t         etaFromTrackerSc[11];   //[numberOfScGt4]
   Float_t         phiFromTrackerSc[11];   //[numberOfScGt4]
   Float_t         eneRecoSc[11];   //[numberOfScGt4]
   Float_t         etRecoSc[11];   //[numberOfScGt4]
   Float_t         etaRecoSc[11];   //[numberOfScGt4]
   Float_t         etaCorrRecoSc[11];   //[numberOfScGt4]
   Float_t         e9e25RecoSc[11];   //[numberOfScGt4]
   Float_t         sigmaEtaEtaRecoSc[11];   //[numberOfScGt4]
   Float_t         HoverERecoSc[11];   //[numberOfScGt4]
   Float_t         dEtaWithTrackerRecoSc[11];   //[numberOfScGt4]
   Float_t         dPhiWithTrackerRecoSc[11];   //[numberOfScGt4]
   Float_t         EoverPRecoSc[11];   //[numberOfScGt4]
   Float_t         trkIsolRecoSc_03[11];   //[numberOfScGt4]
   Float_t         trkIsolRecoSc_04[11];   //[numberOfScGt4]
   Float_t         trkIsolRecoSc_05[11];   //[numberOfScGt4]
   Float_t         hadIsolRecoSc_03[11];   //[numberOfScGt4]
   Float_t         hadIsolRecoSc_04[11];   //[numberOfScGt4]
   Float_t         hadIsolRecoSc_05[11];   //[numberOfScGt4]
   Float_t         emIsolRecoSc_03[11];   //[numberOfScGt4]
   Float_t         emIsolRecoSc_04[11];   //[numberOfScGt4]
   Float_t         emIsolRecoSc_05[11];   //[numberOfScGt4]

   // List of branches
   TBranch        *b_signal;   //!
   TBranch        *b_numberOfPrimaryVtx;   //!
   TBranch        *b_numberOfGenerated;   //!
   TBranch        *b_numberOfScGt4;   //!
   TBranch        *b_numberOfScMatchedPlus;   //!
   TBranch        *b_numberOfScMatchedMinus;   //!
   TBranch        *b_numberOfScMatchedAll;   //!
   TBranch        *b_numberOfElectrons;   //!
   TBranch        *b_hlt29;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_chargeGenEle;   //!
   TBranch        *b_pxGenEle;   //!
   TBranch        *b_pyGenEle;   //!
   TBranch        *b_pzGenEle;   //!
   TBranch        *b_eneGenEle;   //!
   TBranch        *b_dRWithTrackerRecoEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_xRecoEle;   //!
   TBranch        *b_yRecoEle;   //!
   TBranch        *b_zRecoEle;   //!
   TBranch        *b_etaRecoEle;   //!
   TBranch        *b_phiRecoEle;   //!
   TBranch        *b_eneRecoEle;   //!
   TBranch        *b_etRecoEle;   //!
   TBranch        *b_e9e25RecoEle;   //!
   TBranch        *b_sigmaEtaEtaRecoEle;   //!
   TBranch        *b_HoverERecoEle;   //!
   TBranch        *b_dEtaWithTrackerRecoEle;   //!
   TBranch        *b_dPhiWithTrackerRecoEle;   //!
   TBranch        *b_EoverPRecoEle;   //!
   TBranch        *b_trkIsolRecoEle_03;   //!
   TBranch        *b_trkIsolRecoEle_04;   //!
   TBranch        *b_trkIsolRecoEle_05;   //!
   TBranch        *b_hadIsolRecoEle_03;   //!
   TBranch        *b_hadIsolRecoEle_04;   //!
   TBranch        *b_hadIsolRecoEle_05;   //!
   TBranch        *b_emIsolRecoEle_03;   //!
   TBranch        *b_emIsolRecoEle_04;   //!
   TBranch        *b_emIsolRecoEle_05;   //!
   TBranch        *b_jurTrkIsolEle;   //!
   TBranch        *b_jurEmIsolEle;   //!
   TBranch        *b_jurHadIsolEle;   //!
   TBranch        *b_dRWithTrackerRecoSc;   //!
   TBranch        *b_dRWithTrackerRecoNotCorrSc;   //!
   TBranch        *b_chargeSc;   //!
   TBranch        *b_xPVRecoSc;   //!
   TBranch        *b_yPVRecoSc;   //!
   TBranch        *b_zPVRecoSc;   //!
   TBranch        *b_etaPVRecoSc;   //!
   TBranch        *b_phiPVRecoSc;   //!
   TBranch        *b_xOrRecoSc;   //!
   TBranch        *b_yOrRecoSc;   //!
   TBranch        *b_zOrRecoSc;   //!
   TBranch        *b_etaOrRecoSc;   //!
   TBranch        *b_phiOrRecoSc;   //!
   TBranch        *b_pxFromTrackerSc;   //!
   TBranch        *b_pyFromTrackerSc;   //!
   TBranch        *b_pzFromTrackerSc;   //!
   TBranch        *b_etaFromTrackerSc;   //!
   TBranch        *b_phiFromTrackerSc;   //!
   TBranch        *b_eneRecoSc;   //!
   TBranch        *b_etRecoSc;   //!
   TBranch        *b_etaRecoSc;   //!
   TBranch        *b_etaCorrRecoSc;   //!
   TBranch        *b_e9e25RecoSc;   //!
   TBranch        *b_sigmaEtaEtaRecoSc;   //!
   TBranch        *b_HoverERecoSc;   //!
   TBranch        *b_dEtaWithTrackerRecoSc;   //!
   TBranch        *b_dPhiWithTrackerRecoSc;   //!
   TBranch        *b_EoverPRecoSc;   //!
   TBranch        *b_trkIsolRecoSc_03;   //!
   TBranch        *b_trkIsolRecoSc_04;   //!
   TBranch        *b_trkIsolRecoSc_05;   //!
   TBranch        *b_hadIsolRecoSc_03;   //!
   TBranch        *b_hadIsolRecoSc_04;   //!
   TBranch        *b_hadIsolRecoSc_05;   //!
   TBranch        *b_emIsolRecoSc_03;   //!
   TBranch        *b_emIsolRecoSc_04;   //!
   TBranch        *b_emIsolRecoSc_05;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u1/crovelli/dataJPSI/analisi_conHLT/signalPrompt/signalPrompt_notHLTskimmed_mergedTree.root");
      if (!f) {
         f = new TFile("/u1/crovelli/dataJPSI/analisi_conHLT/signalPrompt/signalPrompt_notHLTskimmed_mergedTree.root");
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

   fChain->SetBranchAddress("signal", &signal, &b_signal);
   fChain->SetBranchAddress("numberOfPrimaryVtx", &numberOfPrimaryVtx, &b_numberOfPrimaryVtx);
   fChain->SetBranchAddress("numberOfGenerated", &numberOfGenerated, &b_numberOfGenerated);
   fChain->SetBranchAddress("numberOfScGt4", &numberOfScGt4, &b_numberOfScGt4);
   fChain->SetBranchAddress("numberOfScMatchedPlus", &numberOfScMatchedPlus, &b_numberOfScMatchedPlus);
   fChain->SetBranchAddress("numberOfScMatchedMinus", &numberOfScMatchedMinus, &b_numberOfScMatchedMinus);
   fChain->SetBranchAddress("numberOfScMatchedAll", &numberOfScMatchedAll, &b_numberOfScMatchedAll);
   fChain->SetBranchAddress("numberOfElectrons", &numberOfElectrons, &b_numberOfElectrons);
   fChain->SetBranchAddress("hlt29", &hlt29, &b_hlt29);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("chargeGenEle", chargeGenEle, &b_chargeGenEle);
   fChain->SetBranchAddress("pxGenEle", pxGenEle, &b_pxGenEle);
   fChain->SetBranchAddress("pyGenEle", pyGenEle, &b_pyGenEle);
   fChain->SetBranchAddress("pzGenEle", pzGenEle, &b_pzGenEle);
   fChain->SetBranchAddress("eneGenEle", eneGenEle, &b_eneGenEle);
   fChain->SetBranchAddress("dRWithTrackerRecoEle", dRWithTrackerRecoEle, &b_dRWithTrackerRecoEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("xRecoEle", xRecoEle, &b_xRecoEle);
   fChain->SetBranchAddress("yRecoEle", yRecoEle, &b_yRecoEle);
   fChain->SetBranchAddress("zRecoEle", zRecoEle, &b_zRecoEle);
   fChain->SetBranchAddress("etaRecoEle", etaRecoEle, &b_etaRecoEle);
   fChain->SetBranchAddress("phiRecoEle", phiRecoEle, &b_phiRecoEle);
   fChain->SetBranchAddress("eneRecoEle", eneRecoEle, &b_eneRecoEle);
   fChain->SetBranchAddress("etRecoEle", etRecoEle, &b_etRecoEle);
   fChain->SetBranchAddress("e9e25RecoEle", e9e25RecoEle, &b_e9e25RecoEle);
   fChain->SetBranchAddress("sigmaEtaEtaRecoEle", sigmaEtaEtaRecoEle, &b_sigmaEtaEtaRecoEle);
   fChain->SetBranchAddress("HoverERecoEle", HoverERecoEle, &b_HoverERecoEle);
   fChain->SetBranchAddress("dEtaWithTrackerRecoEle", dEtaWithTrackerRecoEle, &b_dEtaWithTrackerRecoEle);
   fChain->SetBranchAddress("dPhiWithTrackerRecoEle", dPhiWithTrackerRecoEle, &b_dPhiWithTrackerRecoEle);
   fChain->SetBranchAddress("EoverPRecoEle", EoverPRecoEle, &b_EoverPRecoEle);
   fChain->SetBranchAddress("trkIsolRecoEle_03", trkIsolRecoEle_03, &b_trkIsolRecoEle_03);
   fChain->SetBranchAddress("trkIsolRecoEle_04", trkIsolRecoEle_04, &b_trkIsolRecoEle_04);
   fChain->SetBranchAddress("trkIsolRecoEle_05", trkIsolRecoEle_05, &b_trkIsolRecoEle_05);
   fChain->SetBranchAddress("hadIsolRecoEle_03", hadIsolRecoEle_03, &b_hadIsolRecoEle_03);
   fChain->SetBranchAddress("hadIsolRecoEle_04", hadIsolRecoEle_04, &b_hadIsolRecoEle_04);
   fChain->SetBranchAddress("hadIsolRecoEle_05", hadIsolRecoEle_05, &b_hadIsolRecoEle_05);
   fChain->SetBranchAddress("emIsolRecoEle_03", emIsolRecoEle_03, &b_emIsolRecoEle_03);
   fChain->SetBranchAddress("emIsolRecoEle_04", emIsolRecoEle_04, &b_emIsolRecoEle_04);
   fChain->SetBranchAddress("emIsolRecoEle_05", emIsolRecoEle_05, &b_emIsolRecoEle_05);
   fChain->SetBranchAddress("jurTrkIsolEle", jurTrkIsolEle, &b_jurTrkIsolEle);
   fChain->SetBranchAddress("jurEmIsolEle", jurEmIsolEle, &b_jurEmIsolEle);
   fChain->SetBranchAddress("jurHadIsolEle", jurHadIsolEle, &b_jurHadIsolEle);
   fChain->SetBranchAddress("dRWithTrackerRecoSc", dRWithTrackerRecoSc, &b_dRWithTrackerRecoSc);
   fChain->SetBranchAddress("dRWithTrackerRecoNotCorrSc", dRWithTrackerRecoNotCorrSc, &b_dRWithTrackerRecoNotCorrSc);
   fChain->SetBranchAddress("chargeSc", chargeSc, &b_chargeSc);
   fChain->SetBranchAddress("xPVRecoSc", xPVRecoSc, &b_xPVRecoSc);
   fChain->SetBranchAddress("yPVRecoSc", yPVRecoSc, &b_yPVRecoSc);
   fChain->SetBranchAddress("zPVRecoSc", zPVRecoSc, &b_zPVRecoSc);
   fChain->SetBranchAddress("etaPVRecoSc", etaPVRecoSc, &b_etaPVRecoSc);
   fChain->SetBranchAddress("phiPVRecoSc", phiPVRecoSc, &b_phiPVRecoSc);
   fChain->SetBranchAddress("xOrRecoSc", xOrRecoSc, &b_xOrRecoSc);
   fChain->SetBranchAddress("yOrRecoSc", yOrRecoSc, &b_yOrRecoSc);
   fChain->SetBranchAddress("zOrRecoSc", zOrRecoSc, &b_zOrRecoSc);
   fChain->SetBranchAddress("etaOrRecoSc", etaOrRecoSc, &b_etaOrRecoSc);
   fChain->SetBranchAddress("phiOrRecoSc", phiOrRecoSc, &b_phiOrRecoSc);
   fChain->SetBranchAddress("pxFromTrackerSc", pxFromTrackerSc, &b_pxFromTrackerSc);
   fChain->SetBranchAddress("pyFromTrackerSc", pyFromTrackerSc, &b_pyFromTrackerSc);
   fChain->SetBranchAddress("pzFromTrackerSc", pzFromTrackerSc, &b_pzFromTrackerSc);
   fChain->SetBranchAddress("etaFromTrackerSc", etaFromTrackerSc, &b_etaFromTrackerSc);
   fChain->SetBranchAddress("phiFromTrackerSc", phiFromTrackerSc, &b_phiFromTrackerSc);
   fChain->SetBranchAddress("eneRecoSc", eneRecoSc, &b_eneRecoSc);
   fChain->SetBranchAddress("etRecoSc", etRecoSc, &b_etRecoSc);
   fChain->SetBranchAddress("etaRecoSc", etaRecoSc, &b_etaRecoSc);
   fChain->SetBranchAddress("etaCorrRecoSc", etaCorrRecoSc, &b_etaCorrRecoSc);
   fChain->SetBranchAddress("e9e25RecoSc", e9e25RecoSc, &b_e9e25RecoSc);
   fChain->SetBranchAddress("sigmaEtaEtaRecoSc", sigmaEtaEtaRecoSc, &b_sigmaEtaEtaRecoSc);
   fChain->SetBranchAddress("HoverERecoSc", HoverERecoSc, &b_HoverERecoSc);
   fChain->SetBranchAddress("dEtaWithTrackerRecoSc", dEtaWithTrackerRecoSc, &b_dEtaWithTrackerRecoSc);
   fChain->SetBranchAddress("dPhiWithTrackerRecoSc", dPhiWithTrackerRecoSc, &b_dPhiWithTrackerRecoSc);
   fChain->SetBranchAddress("EoverPRecoSc", EoverPRecoSc, &b_EoverPRecoSc);
   fChain->SetBranchAddress("trkIsolRecoSc_03", trkIsolRecoSc_03, &b_trkIsolRecoSc_03);
   fChain->SetBranchAddress("trkIsolRecoSc_04", trkIsolRecoSc_04, &b_trkIsolRecoSc_04);
   fChain->SetBranchAddress("trkIsolRecoSc_05", trkIsolRecoSc_05, &b_trkIsolRecoSc_05);
   fChain->SetBranchAddress("hadIsolRecoSc_03", hadIsolRecoSc_03, &b_hadIsolRecoSc_03);
   fChain->SetBranchAddress("hadIsolRecoSc_04", hadIsolRecoSc_04, &b_hadIsolRecoSc_04);
   fChain->SetBranchAddress("hadIsolRecoSc_05", hadIsolRecoSc_05, &b_hadIsolRecoSc_05);
   fChain->SetBranchAddress("emIsolRecoSc_03", emIsolRecoSc_03, &b_emIsolRecoSc_03);
   fChain->SetBranchAddress("emIsolRecoSc_04", emIsolRecoSc_04, &b_emIsolRecoSc_04);
   fChain->SetBranchAddress("emIsolRecoSc_05", emIsolRecoSc_05, &b_emIsolRecoSc_05);
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

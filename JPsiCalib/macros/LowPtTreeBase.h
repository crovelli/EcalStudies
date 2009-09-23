//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 14 11:45:42 2009 by ROOT version 5.22/00a
// from TTree T1/filter studies tree
// found on file: /tmp/arcidiac/filteringMB_tree_1.root
//////////////////////////////////////////////////////////

#ifndef LowPtTreeBase_h
#define LowPtTreeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class LowPtTreeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           numberOfGenerated;
   Int_t           numberOfElectrons;
   Int_t           passedFilter1;
   Int_t           passedFilter2;
   Int_t           passedFilter3;
   Int_t           passedFilter4;
   Int_t           passedFilter5;
   Int_t           passedFilter6;
   Int_t           passedFilter7;
   Int_t           passedFilter8;
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiSection;
   Float_t         ptHat;
   Float_t         pxGen[1500];   //[numberOfGenerated]
   Float_t         pyGen[1500];   //[numberOfGenerated]
   Float_t         pzGen[1500];   //[numberOfGenerated]
   Float_t         eneGen[1500];   //[numberOfGenerated]
   Float_t         etaGen[1500];   //[numberOfGenerated]
   Float_t         phiGen[1500];   //[numberOfGenerated]
   Int_t           idGen[1500];   //[numberOfGenerated]
   Int_t           statusGen[1500];   //[numberOfGenerated]
   Int_t           motherIdGen[1500];   //[numberOfGenerated]
   Int_t           motherGen[1500];   //[numberOfGenerated]
   Int_t           chargeEle[50];   //[numberOfElectrons]
   Float_t         xRecoEle[50];   //[numberOfElectrons]
   Float_t         yRecoEle[50];   //[numberOfElectrons]
   Float_t         zRecoEle[50];   //[numberOfElectrons]
   Float_t         etaRecoEle[50];   //[numberOfElectrons]
   Float_t         phiRecoEle[50];   //[numberOfElectrons]
   Float_t         eneRecoEle[50];   //[numberOfElectrons]
   Float_t         etRecoEle[50];   //[numberOfElectrons]
   Float_t         HoverERecoEle[50];   //[numberOfElectrons]
   Float_t         dEtaWithTrackerRecoEle[50];   //[numberOfElectrons]
   Float_t         dPhiWithTrackerRecoEle[50];   //[numberOfElectrons]
   Float_t         EoverPRecoEle[50];   //[numberOfElectrons]
   Int_t           eleIdLoose[50];   //[numberOfElectrons]
   Int_t           eleIdRobLoose[50];   //[numberOfElectrons]
   Int_t           eleIdRobTight[50];   //[numberOfElectrons]
   Int_t           eleIdTight[50];   //[numberOfElectrons]
   Float_t         sumPt03[50];   //[numberOfElectrons]
   Float_t         sumPt04[50];   //[numberOfElectrons]
   Float_t         sumEtEcal03[50];   //[numberOfElectrons]
   Float_t         sumEtEcal04[50];   //[numberOfElectrons]
   Float_t         sumEtHcalD103[50];   //[numberOfElectrons]
   Float_t         sumEtHcalD104[50];   //[numberOfElectrons]
   Float_t         sumEtHcalD203[50];   //[numberOfElectrons]
   Float_t         sumEtHcalD204[50];   //[numberOfElectrons]

   // List of branches
   TBranch        *b_numberOfGenerated;   //!
   TBranch        *b_numberOfElectrons;   //!
   TBranch        *b_passedFilter1;   //!
   TBranch        *b_passedFilter2;   //!
   TBranch        *b_passedFilter3;   //!
   TBranch        *b_passedFilter4;   //!
   TBranch        *b_passedFilter5;   //!
   TBranch        *b_passedFilter6;   //!
   TBranch        *b_passedFilter7;   //!
   TBranch        *b_passedFilter8;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_pxGen;   //!
   TBranch        *b_pyGen;   //!
   TBranch        *b_pzGen;   //!
   TBranch        *b_eneGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_idGen;   //!
   TBranch        *b_statusGen;   //!
   TBranch        *b_motherIdGen;   //!
   TBranch        *b_motherGen;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_xRecoEle;   //!
   TBranch        *b_yRecoEle;   //!
   TBranch        *b_zRecoEle;   //!
   TBranch        *b_etaRecoEle;   //!
   TBranch        *b_phiRecoEle;   //!
   TBranch        *b_eneRecoEle;   //!
   TBranch        *b_etRecoEle;   //!
   TBranch        *b_HoverERecoEle;   //!
   TBranch        *b_dEtaWithTrackerRecoEle;   //!
   TBranch        *b_dPhiWithTrackerRecoEle;   //!
   TBranch        *b_EoverPRecoEle;   //!
   TBranch        *b_eleIdLoose;   //!
   TBranch        *b_eleIdRobLoose;   //!
   TBranch        *b_eleIdRobTight;   //!
   TBranch        *b_eleIdTight;   //!
   TBranch        *b_sumPt03;   //!
   TBranch        *b_sumPt04;   //!
   TBranch        *b_sumEtEcal03;   //!
   TBranch        *b_sumEtEcal04;   //!
   TBranch        *b_sumEtHcalD103;   //!
   TBranch        *b_sumEtHcalD104;   //!
   TBranch        *b_sumEtHcalD203;   //!
   TBranch        *b_sumEtHcalD204;   //!

   LowPtTreeBase(TTree *tree=0);
   virtual ~LowPtTreeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LowPtTreeBase_cxx
LowPtTreeBase::LowPtTreeBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/arcidiac/filteringMB_tree_1.root");
      if (!f) {
         f = new TFile("/tmp/arcidiac/filteringMB_tree_1.root");
      }
      tree = (TTree*)gDirectory->Get("T1");

   }
   Init(tree);
}

LowPtTreeBase::~LowPtTreeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LowPtTreeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LowPtTreeBase::LoadTree(Long64_t entry)
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

void LowPtTreeBase::Init(TTree *tree)
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

   fChain->SetBranchAddress("numberOfGenerated", &numberOfGenerated, &b_numberOfGenerated);
   fChain->SetBranchAddress("numberOfElectrons", &numberOfElectrons, &b_numberOfElectrons);
   fChain->SetBranchAddress("passedFilter1", &passedFilter1, &b_passedFilter1);
   fChain->SetBranchAddress("passedFilter2", &passedFilter2, &b_passedFilter2);
   fChain->SetBranchAddress("passedFilter3", &passedFilter3, &b_passedFilter3);
   fChain->SetBranchAddress("passedFilter4", &passedFilter4, &b_passedFilter4);
   fChain->SetBranchAddress("passedFilter5", &passedFilter5, &b_passedFilter5);
   fChain->SetBranchAddress("passedFilter6", &passedFilter6, &b_passedFilter6);
   fChain->SetBranchAddress("passedFilter7", &passedFilter7, &b_passedFilter7);
   fChain->SetBranchAddress("passedFilter8", &passedFilter8, &b_passedFilter8);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("pxGen", pxGen, &b_pxGen);
   fChain->SetBranchAddress("pyGen", pyGen, &b_pyGen);
   fChain->SetBranchAddress("pzGen", pzGen, &b_pzGen);
   fChain->SetBranchAddress("eneGen", eneGen, &b_eneGen);
   fChain->SetBranchAddress("etaGen", etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", phiGen, &b_phiGen);
   fChain->SetBranchAddress("idGen", idGen, &b_idGen);
   fChain->SetBranchAddress("statusGen", statusGen, &b_statusGen);
   fChain->SetBranchAddress("motherIdGen", motherIdGen, &b_motherIdGen);
   fChain->SetBranchAddress("motherGen", motherGen, &b_motherGen);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("xRecoEle", xRecoEle, &b_xRecoEle);
   fChain->SetBranchAddress("yRecoEle", yRecoEle, &b_yRecoEle);
   fChain->SetBranchAddress("zRecoEle", zRecoEle, &b_zRecoEle);
   fChain->SetBranchAddress("etaRecoEle", etaRecoEle, &b_etaRecoEle);
   fChain->SetBranchAddress("phiRecoEle", phiRecoEle, &b_phiRecoEle);
   fChain->SetBranchAddress("eneRecoEle", eneRecoEle, &b_eneRecoEle);
   fChain->SetBranchAddress("etRecoEle", etRecoEle, &b_etRecoEle);
   fChain->SetBranchAddress("HoverERecoEle", HoverERecoEle, &b_HoverERecoEle);
   fChain->SetBranchAddress("dEtaWithTrackerRecoEle", dEtaWithTrackerRecoEle, &b_dEtaWithTrackerRecoEle);
   fChain->SetBranchAddress("dPhiWithTrackerRecoEle", dPhiWithTrackerRecoEle, &b_dPhiWithTrackerRecoEle);
   fChain->SetBranchAddress("EoverPRecoEle", EoverPRecoEle, &b_EoverPRecoEle);
   fChain->SetBranchAddress("eleIdLoose", eleIdLoose, &b_eleIdLoose);
   fChain->SetBranchAddress("eleIdRobLoose", eleIdRobLoose, &b_eleIdRobLoose);
   fChain->SetBranchAddress("eleIdRobTight", eleIdRobTight, &b_eleIdRobTight);
   fChain->SetBranchAddress("eleIdTight", eleIdTight, &b_eleIdTight);
   fChain->SetBranchAddress("sumPt03", sumPt03, &b_sumPt03);
   fChain->SetBranchAddress("sumPt04", sumPt04, &b_sumPt04);
   fChain->SetBranchAddress("sumEtEcal03", sumEtEcal03, &b_sumEtEcal03);
   fChain->SetBranchAddress("sumEtEcal04", sumEtEcal04, &b_sumEtEcal04);
   fChain->SetBranchAddress("sumEtHcalD103", sumEtHcalD103, &b_sumEtHcalD103);
   fChain->SetBranchAddress("sumEtHcalD104", sumEtHcalD104, &b_sumEtHcalD104);
   fChain->SetBranchAddress("sumEtHcalD203", sumEtHcalD203, &b_sumEtHcalD203);
   fChain->SetBranchAddress("sumEtHcalD204", sumEtHcalD204, &b_sumEtHcalD204);
   Notify();
}

Bool_t LowPtTreeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LowPtTreeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LowPtTreeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LowPtTreeBase_cxx

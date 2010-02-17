#include <vector>
#include <iostream>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <THStack.h>

#include "myFunctions.cc"

#define NDATA 4

// vectors to be used in the following
TFile *file[NDATA];
float expected[NDATA];

// histos declaration
TH1F *ScHisto_invMassHighestEt[NDATA];
TH1F *ScHisto_invMassHighestEt_EB[NDATA];
TH1F *ScHisto_invMassHighestEt_EE[NDATA];
TH1F *ScHisto_invMassHighestEt_EBEE[NDATA];
TH1F *ScHisto_invMassHighestEt_norm1[NDATA];
TH1F *ScHisto_invMassHighestEt_EB_norm1[NDATA];
TH1F *ScHisto_invMassHighestEt_EE_norm1[NDATA];
TH1F *ScHisto_invMassHighestEt_EBEE_norm1[NDATA];
TH1F *ScHisto_invMassHighestEt_normL[NDATA];
TH1F *ScHisto_invMassHighestEt_EB_normL[NDATA];
TH1F *ScHisto_invMassHighestEt_EE_normL[NDATA];
TH1F *ScHisto_invMassHighestEt_EBEE_normL[NDATA];
TH1F *ScHisto_invMassHighestEt_normL_all;

TH1F *EleHisto_detaVtx[NDATA];
TH1F *EleHisto_dphiVtx[NDATA];
TH1F *EleHisto_EoP[NDATA];
TH1F *EleHisto_dr03ecal[NDATA];
TH1F *EleHisto_dr03tkCorr[NDATA];
TH1F *EleHisto_detaVtx_normL[NDATA];
TH1F *EleHisto_dphiVtx_normL[NDATA];
TH1F *EleHisto_EoP_normL[NDATA];
TH1F *EleHisto_dr03ecal_normL[NDATA];
TH1F *EleHisto_dr03tkCorr_normL[NDATA];


void chargeFiles() {
  file[0] = new TFile("Outfile_Jpsi.root");
  file[1] = new TFile("Outfile_bce2030.root");
  file[2] = new TFile("Outfile_bce3080.root");
  file[3] = new TFile("Outfile_bce80170.root");
  file[4] = new TFile("Outfile_EM20to30.root");
  file[5] = new TFile("Outfile_EM30to80.root");
  file[6] = new TFile("Outfile_EM80to170.root");
  file[7] = new TFile("Outfile_doubleEM6to20.root");
  file[8] = new TFile("Outfile_DYee.root");
}

void chargeFilesMB() {
  file[0] = new TFile("Outfile_Jpsi.root");
  file[1] = new TFile("Outfile_doubleEM6to20.root");
  file[2] = new TFile("Outfile_doubleEMgt20.root");
  file[3] = new TFile("Outfile_DYee.root");
}

void setExpected(int theStep) {

  float lumi = 1.;   // pb-1   !! number are already related to 10pb-1
  
  // float exp_         =           kineEff  * filterEff_prod*crossSection;
  // cout << "after HLT  : " << 10.*exp_hlt      << endl; 
  
  if (theStep==1) {     // after HLT only
    expected[0] =  *lumi;
    expected[1] =  *lumi;
    expected[2] =  *lumi;
    expected[3] =  *lumi;
    expected[4] =  *lumi;
    expected[5] =  *lumi;
    expected[6] =  *lumi;
    expected[7] =  *lumi;
    expected[8] =  *lumi;
  }
  
  if (theStep==2) {     // after HLT + reco
    expected[0] = *lumi;
    expected[1] = *lumi;
    expected[2] = *lumi;
    expected[3] = *lumi;
    expected[4] = *lumi;
    expected[5] = *lumi;
    expected[6] = *lumi;
    expected[7] = *lumi;
    expected[8] = *lumi;
  }

  if (theStep==3) {     // after HLT + reco + Et>4 
    expected[0] = *lumi;
    expected[1] = *lumi;
    expected[2] = *lumi;
    expected[3] = *lumi;
    expected[4] = *lumi;
    expected[5] = *lumi;
    expected[6] = *lumi;
    expected[7] = *lumi;
    expected[8] = *lumi;
  }

  if (theStep==4) {     // after HLT + reco + Et>4 + charge 
    expected[0] = 7808*lumi;
    expected[1] =  7251.*lumi;
    expected[2] = 46979.8*lumi;
    expected[3] =  9749.7*lumi;
    expected[4] = 41525.*lumi;
    expected[5] = 314940 *lumi; 
    expected[6] =  46827.*lumi; 
    expected[7] = 112821*lumi;
    expected[8] = 39.8*lumi;
  }

  if (theStep==5) {     // after HLT + reco + Et>4 + charge + eleID 
    expected[0] = 5048*lumi;
    expected[1] = 2234*lumi;
    expected[2] = 8590.*lumi;
    expected[3] = 864.*lumi;
    expected[4] = 4856*lumi;
    expected[5] = 25530*lumi;
    expected[6] =  2517*lumi;
    expected[7] = 38236*lumi;
    expected[8] = 29*lumi;
  }

  if (theStep==6) {     // after HLT + reco + Et>4 + charge + eleID + tk isol 
    expected[0] = 4431.76*lumi;
    expected[1] = 585.8*lumi;
    expected[2] = 835  *lumi;
    expected[3] = 29.1 *lumi;
    expected[4] = 1508.*lumi;
    expected[5] = 3077 *lumi;
    expected[6] = 155.8*lumi;
    expected[7] = 15105*lumi; 
    expected[8] = 26.6 *lumi;
  }
}

void setExpectedMB(int theStep) {

  float lumi = 1.;   // pb-1   !! number are already related to 10pb-1
  
  // float exp_         =           kineEff  * filterEff_prod*crossSection;
  // cout << "after HLT  : " << 10.*exp_hlt      << endl; 
  

  if (theStep==4) {     // after HLT + reco + Et>4 + charge 
    expected[0] = 7808*lumi;
    expected[1] = 112821*lumi;
    expected[2] = 1363340*lumi; 
    expected[3] = 39.8*lumi;
  }

  if (theStep==5) {     // after HLT + reco + Et>4 + charge + eleID 
    expected[0] = 5048*lumi;
    expected[1] = 38236*lumi;
    expected[2] = 174937 *lumi; 
    expected[3] = 29*lumi;
  }

  if (theStep==6) {     // after HLT + reco + Et>4 + charge + eleID + isol 
    expected[0] = 4431.76*lumi;
    expected[1] = 15105*lumi; 
    expected[2] = 24926 *lumi; 
    expected[3] = 26.6 *lumi;
  }
}

void readHistos(int theStep) {
  for(int ifile=0; ifile<NDATA; ifile++){
    ostringstream name;
    name << "ScHisto_invMassHighestEt" << theStep ;
    string strname = name.str();
    ScHisto_invMassHighestEt[ifile]      = (TH1F*)file[ifile]->Get(strname.c_str());
    ScHisto_invMassHighestEt_EB[ifile]   = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EB");
    ScHisto_invMassHighestEt_EE[ifile]   = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EE");
    ScHisto_invMassHighestEt_EBEE[ifile] = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EBEE");
  }
}

void normalizeTo1() { 

  char title[700];

  for(int ifile=0; ifile<NDATA; ifile++){
    sprintf(title,"ScHisto_invMassHighestEt_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_norm1[ifile]      = (TH1F*)ScHisto_invMassHighestEt[ifile]->Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EB_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EB_norm1[ifile]   = (TH1F*)ScHisto_invMassHighestEt_EB[ifile] -> Clone(title);  
    sprintf(title,"ScHisto_invMassHighestEt_EE_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EE_norm1[ifile]   = (TH1F*)ScHisto_invMassHighestEt_EE[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EBEE_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EBEE_norm1[ifile] = (TH1F*)ScHisto_invMassHighestEt_EBEE[ifile] -> Clone(title);
  }

  for(int ifile=0; ifile<NDATA; ifile++){
    if (ScHisto_invMassHighestEt_norm1[ifile]    -> Integral()>0) 
      ScHisto_invMassHighestEt_norm1[ifile]      -> Scale(1./(ScHisto_invMassHighestEt_norm1[ifile]     -> Integral()));
    if (ScHisto_invMassHighestEt_EB_norm1[ifile] -> Integral()>0) 
      ScHisto_invMassHighestEt_EB_norm1[ifile]   -> Scale(1./(ScHisto_invMassHighestEt_EB_norm1[ifile]  -> Integral()));
    if (ScHisto_invMassHighestEt_EE_norm1[ifile] -> Integral()>0) 
      ScHisto_invMassHighestEt_EE_norm1[ifile]   -> Scale(1./(ScHisto_invMassHighestEt_EE_norm1[ifile]  -> Integral()));
    if (ScHisto_invMassHighestEt_EBEE_norm1[ifile] -> Integral()>0) 
      ScHisto_invMassHighestEt_EBEE_norm1[ifile] -> Scale(1./(ScHisto_invMassHighestEt_EBEE_norm1[ifile]-> Integral()));
  }
}


void normalizeToL() { 

  char title[500];
  for(int ifile=0; ifile<NDATA; ifile++){
    sprintf(title,"ScHisto_invMassHighestEt_normL[%d]", ifile);
    ScHisto_invMassHighestEt_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EB_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EB_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EB[ifile] -> Clone(title);  
    sprintf(title,"ScHisto_invMassHighestEt_EE_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EE_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EE[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EBEE_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EBEE_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EBEE[ifile] -> Clone(title);
  }

  for(int ifile=0; ifile<NDATA; ifile++){
    if (ScHisto_invMassHighestEt_normL[ifile]   -> GetEntries()>0) 
      ScHisto_invMassHighestEt_normL[ifile]     -> Scale(expected[ifile]/(ScHisto_invMassHighestEt_normL[ifile] -> GetEntries()));
    if (ScHisto_invMassHighestEt_EB_normL[ifile]-> GetEntries()>0) 
      ScHisto_invMassHighestEt_EB_normL[ifile]  -> Scale(expected[ifile]/(ScHisto_invMassHighestEt_EB_normL[ifile] -> GetEntries()));
    if (ScHisto_invMassHighestEt_EE_normL[ifile]-> GetEntries()>0) 
      ScHisto_invMassHighestEt_EE_normL[ifile]  -> Scale(expected[ifile]/(ScHisto_invMassHighestEt_EE_normL[ifile] -> GetEntries()));
    if (ScHisto_invMassHighestEt_EBEE_normL[ifile]-> GetEntries()>0) 
      ScHisto_invMassHighestEt_EBEE_normL[ifile]  -> Scale(expected[ifile]/(ScHisto_invMassHighestEt_EBEE_normL[ifile] -> GetEntries()));
  }
}

void addStoB() {

  ScHisto_invMassHighestEt_normL_all = (TH1F*)ScHisto_invMassHighestEt_normL[0]->Clone("ScHisto_invMassHighestEt_normL_all");
  for(int ifile=1; ifile<NDATA; ifile++){
    ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[ifile]);
  }
}

void cosmeticsSvsB() {
  int color;
  for(int ifile=0; ifile<NDATA; ifile++){
    color = ifile+2;
    if (ifile > 7 ) {color = 20 + ifile;}
    
    ScHisto_invMassHighestEt_normL[ifile]      -> SetFillColor(color);
    ScHisto_invMassHighestEt_EB_normL[ifile]   -> SetLineColor(color);
    ScHisto_invMassHighestEt_EE_normL[ifile]   -> SetLineColor(color); 
    ScHisto_invMassHighestEt_EBEE_normL[ifile] -> SetLineColor(color);
    ScHisto_invMassHighestEt_norm1[ifile]      -> SetLineColor(color);
    ScHisto_invMassHighestEt_EB_norm1[ifile]   -> SetLineColor(color);
    ScHisto_invMassHighestEt_EE_norm1[ifile]   -> SetLineColor(color); 
    ScHisto_invMassHighestEt_EBEE_norm1[ifile] -> SetLineColor(color);
    
    ScHisto_invMassHighestEt_normL[ifile]      -> SetLineWidth(2);
    ScHisto_invMassHighestEt_EB_normL[ifile]   -> SetLineWidth(2);
    ScHisto_invMassHighestEt_EE_normL[ifile]   -> SetLineWidth(2); 
    ScHisto_invMassHighestEt_EBEE_normL[ifile] -> SetLineWidth(2);
    ScHisto_invMassHighestEt_norm1[ifile]      -> SetLineWidth(2);
    ScHisto_invMassHighestEt_EB_norm1[ifile]   -> SetLineWidth(2);
    ScHisto_invMassHighestEt_EE_norm1[ifile]   -> SetLineWidth(2); 
    ScHisto_invMassHighestEt_EBEE_norm1[ifile] -> SetLineWidth(2);
  }

  ScHisto_invMassHighestEt_normL_all -> SetLineColor(1);
  ScHisto_invMassHighestEt_normL_all -> SetLineWidth(2);
}

// analysis: S vs B
void drawSvsB(int theStep) {

  chargeFiles();
  setExpected(theStep);
  readHistos(theStep);
  normalizeTo1();
  normalizeToL();
  addStoB();
  cosmeticsSvsB();

  TStyle *tesiStyle = new TStyle("tesiStyle","");
  tesiStyle->SetCanvasColor(0);
  tesiStyle->SetFrameFillColor(0);
  tesiStyle->SetStatColor(0);
  tesiStyle->SetOptStat(000);
  tesiStyle->SetOptFit(0011);
  tesiStyle->SetTitleFillColor(0);
  tesiStyle->SetCanvasBorderMode(0);
  tesiStyle->SetPadBorderMode(0);
  tesiStyle->SetFrameBorderMode(0);
  tesiStyle->cd();

  TLegend leg(0.7,0.8,0.89,0.99);
  leg.AddEntry(ScHisto_invMassHighestEt_normL[0], "Jpsi prompt signal",  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[1], "bc->e, 20-30"      ,  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[2], "bc->e, 30-80"      ,  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[3], "bc->e, 80-170",       "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[4], "em.enriched,  20-30", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[5], "em.enriched,  30-80", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[6], "em.enriched, 80-170", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[7], "double em.enr  6-20", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[8], "DYee  1-10 GeV     ", "f");
  leg.SetFillColor(0);
  leg.SetBorderSize(0.4);
  leg.SetTextSize(0.03);
  leg.SetMargin(0.35);

  TCanvas c1("c1","all",1);  
  ScHisto_invMassHighestEt_normL[0]->Draw();
  for(int ifile=0; ifile<NDATA; ifile++) ScHisto_invMassHighestEt_normL[ifile] -> Draw("same");
  leg.Draw();
  c1.SaveAs("normL.eps");
  c1.SaveAs("normL.root");

  TCanvas c11("c11","all",1);  
  ScHisto_invMassHighestEt_norm1[0]->Draw();
  for(int ifile=0; ifile<NDATA; ifile++) ScHisto_invMassHighestEt_norm1[ifile] -> Draw("same");
  leg.Draw();
  c11.SaveAs("norm1.eps");
  c11.SaveAs("norm1.root");

  TCanvas c21("c21","all",1);  
  ScHisto_invMassHighestEt_normL_all->Draw();
  c21.SaveAs("insieme.eps");
  

  // summing up backgrounds and signals
  THStack staB_invMassHighestEt("staB_invMassHighestEt","staB_invMassHighestEt");
  staB_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[8]);
  for (int ii=1; ii<8; ii++) {
    staB_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  THStack staS_invMassHighestEt("staS_invMassHighestEt","staS_invMassHighestEt");
  for (int ii=0; ii<1; ii++) {
    staS_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  TCanvas c31("c31","all",1);  
  staB_invMassHighestEt.Draw();
  staS_invMassHighestEt.Draw("same");
  leg.Draw();
  c31.SaveAs("stack.eps");
  c31.SaveAs("stack.root");

  THStack staA_invMassHighestEt("staA_invMassHighestEt","staA_invMassHighestEt");
  staA_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[8]);
  for (int ii=0; ii<8; ii++) {
    staA_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  TCanvas c32("c32","all",1);  
  staA_invMassHighestEt.Draw();
  leg.Draw();
  c32.SaveAs("SandB.root");

}

// analysis: S vs B for MinBias samples
void drawSvsB_MB(int theStep) {

  chargeFilesMB();
  setExpectedMB(theStep);
  readHistos(theStep);
  normalizeTo1();
  normalizeToL();
  addStoB();
  cosmeticsSvsB();

  TStyle *tesiStyle = new TStyle("tesiStyle","");
  tesiStyle->SetCanvasColor(0);
  tesiStyle->SetFrameFillColor(0);
  tesiStyle->SetStatColor(0);
  tesiStyle->SetOptStat(000);
  tesiStyle->SetOptFit(0011);
  tesiStyle->SetTitleFillColor(0);
  tesiStyle->SetCanvasBorderMode(0);
  tesiStyle->SetPadBorderMode(0);
  tesiStyle->SetFrameBorderMode(0);
  tesiStyle->cd();

  TLegend leg(0.7,0.8,0.85,0.92);
  leg.AddEntry(ScHisto_invMassHighestEt_normL[0], "Jpsi prompt signal",  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[1], "double em.enr  6-20", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[2], "double em.enr  > 20", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[3], "DYee  1-10 GeV     ", "f");
  leg.SetFillColor(0);
  leg.SetBorderSize(0.4);
  

  // summing up backgrounds and signals
  THStack staB_invMassHighestEt("staB_invMassHighestEt","staB_invMassHighestEt");
  for (int ii=1; ii<4; ii++) {
    staB_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  THStack staS_invMassHighestEt("staS_invMassHighestEt","staS_invMassHighestEt");
  for (int ii=0; ii<1; ii++) {
    staS_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  TCanvas c31("c31","all",1);  
  staB_invMassHighestEt.Draw();
  staS_invMassHighestEt.Draw("same");
  leg.Draw();
  c31.SaveAs("stackMB.eps");
  c31.SaveAs("stackMB.root");

  THStack staA_invMassHighestEt("staA_invMassHighestEt","staA_invMassHighestEt");
  //staA_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[3]);
  for (int ii=0; ii<4; ii++) {
    staA_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  TCanvas c32("c32","all",1);  
  staA_invMassHighestEt.Draw();
  leg.Draw();
  c32.SaveAs("SandBMB.root");

}

// analysis: S vs B for MinBias samples
void drawCuts() {

  chargeFiles();
  setExpected(4);

  TH1F *EleHisto_detaVtx_normL[NDATA];
  TH1F *EleHisto_dphiVtx_normL[NDATA];
  TH1F *EleHisto_EoP_normL[NDATA];
  TH1F *EleHisto_dr03ecal_normL[NDATA];
  TH1F *EleHisto_dr03tkCorr_normL[NDATA];

  for(int ifile=0; ifile<NDATA; ifile++){
    EleHisto_detaVtx[ifile]      = (TH1F*)file[ifile]->Get("EleHisto_detaVtx");
    EleHisto_dphiVtx[ifile]      = (TH1F*)file[ifile]->Get("EleHisto_dphiVtx");
    EleHisto_EoP[ifile]          = (TH1F*)file[ifile]->Get("EleHisto_EoP");
    EleHisto_dr03ecal[ifile]     = (TH1F*)file[ifile]->Get("EleHisto_dr03ecal");
    EleHisto_dr03tkCorr[ifile]    = (TH1F*)file[ifile]->Get("EleHisto_dr03tkCorr");
  }
  char title[500];
  for(int ifile=0; ifile<NDATA; ifile++){
    sprintf(title,"EleHisto_detaVtx_normL[%d]", ifile);
    EleHisto_detaVtx_normL[ifile] = (TH1F*)EleHisto_detaVtx[ifile] -> Clone(title);
    sprintf(title,"EleHisto_dphiVtx_normL[%d]", ifile);
    EleHisto_dphiVtx_normL[ifile] = (TH1F*)EleHisto_dphiVtx[ifile] -> Clone(title);
    sprintf(title,"EleHisto_EoP_normL[%d]", ifile);
    EleHisto_EoP_normL[ifile] = (TH1F*)EleHisto_EoP[ifile] -> Clone(title);
    sprintf(title,"EleHisto_dr03ecal_normL[%d]", ifile);
    EleHisto_dr03ecal_normL[ifile] = (TH1F*)EleHisto_dr03ecal[ifile] -> Clone(title);
    sprintf(title,"EleHisto_dr03tkCorr_normL[%d]", ifile);
    EleHisto_dr03tkCorr_normL[ifile] = (TH1F*)EleHisto_dr03tkCorr[ifile] -> Clone(title);
  }

  for(int ifile=0; ifile<NDATA; ifile++){
    if (EleHisto_detaVtx_normL[ifile]   -> GetEntries()>0) 
      EleHisto_detaVtx_normL[ifile]     -> Scale(expected[ifile]/(EleHisto_detaVtx_normL[ifile] -> GetEntries()));
    if (EleHisto_dphiVtx_normL[ifile]   -> GetEntries()>0) 
      EleHisto_dphiVtx_normL[ifile]     -> Scale(expected[ifile]/(EleHisto_dphiVtx_normL[ifile] -> GetEntries()));
    if (EleHisto_EoP_normL[ifile]   -> GetEntries()>0) 
      EleHisto_EoP_normL[ifile]     -> Scale(expected[ifile]/(EleHisto_EoP_normL[ifile] -> GetEntries()));
    if (EleHisto_dr03ecal_normL[ifile]   -> GetEntries()>0) 
      EleHisto_dr03ecal_normL[ifile]     -> Scale(expected[ifile]/(EleHisto_dr03ecal_normL[ifile] -> GetEntries()));
    if (EleHisto_detaVtx_normL[ifile]   -> GetEntries()>0) 
      EleHisto_dr03tkCorr_normL[ifile]     -> Scale(expected[ifile]/(EleHisto_dr03tkCorr_normL[ifile] -> GetEntries()));
 
    color = 2;
    if (ifile > 0 ) {color = 4;}
    
    EleHisto_detaVtx_normL[ifile]     -> SetLineColor(color);
    EleHisto_dphiVtx_normL[ifile]     -> SetLineColor(color);
    EleHisto_EoP_normL[ifile]         -> SetLineColor(color);
    EleHisto_dr03ecal_normL[ifile]    -> SetLineColor(color);
    EleHisto_dr03tkCorr_normL[ifile]  -> SetLineColor(color);
    
    EleHisto_detaVtx_normL[ifile]     -> SetLineWidth(2);
    EleHisto_dphiVtx_normL[ifile]     -> SetLineWidth(2);
    EleHisto_EoP_normL[ifile]         -> SetLineWidth(2);
    EleHisto_dr03ecal_normL[ifile]    -> SetLineWidth(2);
    EleHisto_dr03tkCorr_normL[ifile]  -> SetLineWidth(2);
  }



  TStyle *tesiStyle = new TStyle("tesiStyle","");
  tesiStyle->SetCanvasColor(0);
  tesiStyle->SetFrameFillColor(0);
  tesiStyle->SetStatColor(0);
  tesiStyle->SetOptStat(000);
  tesiStyle->SetOptFit(0011);
  tesiStyle->SetTitleFillColor(0);
  tesiStyle->SetCanvasBorderMode(0);
  tesiStyle->SetPadBorderMode(0);
  tesiStyle->SetFrameBorderMode(0);
  tesiStyle->cd();

  TLegend leg(0.7,0.8,0.85,0.92);
  leg.AddEntry( EleHisto_detaVtx_normL[0], "Jpsi prompt signal",  "f");
  leg.AddEntry( EleHisto_detaVtx_normL[1], "background", "f");
  leg.SetFillColor(0);
  leg.SetBorderSize(0.4);
  

  // summing up backgrounds and signals
  TH1F  *staB_deta = (TH1F*) EleHisto_detaVtx_normL[1]->Clone( "BEleHisto_detaVtx");
  for (int ii=1; ii<NDATA; ii++) {
    staB_deta->Add(EleHisto_detaVtx_normL[ii]);
  }
  TH1F  *staS_deta = (TH1F*) EleHisto_detaVtx_normL[0]->Clone( "SEleHisto_detaVtx");
  for (int ii=0; ii<1; ii++) {
    staS_deta->Add(EleHisto_detaVtx_normL[ii]);
  }
  TCanvas c31("c31","all",1);  
  staB_deta->Draw();
  staS_deta->Draw("same");
  leg.Draw();

  c31.SaveAs("detacut.root");

}


// analisi signal: EB vs EE vs EBEB vs all
void drawRegions() {

  chargeFiles();
  readHistos(theStep);
  // cosmeticsSignal();
  normalizeTo1();

  TStyle *tesiStyle = new TStyle("tesiStyle","");
  tesiStyle->SetCanvasColor(0);
  tesiStyle->SetFrameFillColor(0);
  tesiStyle->SetStatColor(0);
  tesiStyle->SetOptStat(000);
  tesiStyle->SetOptFit(1111);
  tesiStyle->SetTitleFillColor(0);
  tesiStyle->SetCanvasBorderMode(0);
  tesiStyle->SetPadBorderMode(0);
  tesiStyle->SetFrameBorderMode(0);
  tesiStyle->cd();

  TLegend leg(0.7,0.8,0.85,0.92);
  leg.AddEntry(ScHisto_invMassHighestEt[0],      "prompt signal, all",       "l");
  leg.AddEntry(ScHisto_invMassHighestEt_EB[0],   "prompt signal, barrel",    "l");
  leg.AddEntry(ScHisto_invMassHighestEt_EE[0],   "prompt signal, endcap",    "l");
  leg.AddEntry(ScHisto_invMassHighestEt_EBEE[0], "prompt signal, 1 barrel, 1 endcap", "l");
  leg.SetFillColor(0);
  leg.SetBorderSize(0.4);
  
  TCanvas c1("c1","signal",1);  
  ScHisto_invMassHighestEt_EB[0]   -> Draw();
  ScHisto_invMassHighestEt[0]      -> Draw("same");
  ScHisto_invMassHighestEt_EE[0]   -> Draw("same");
  ScHisto_invMassHighestEt_EBEE[0] -> Draw("same");
  leg.Draw();
  c1.SaveAs("regions.eps");
  c1.SaveAs("regions.root");


  // crujif fit
  TH1F *myH[4];
  myH[0] = ScHisto_invMassHighestEt[0];
  myH[1] = ScHisto_invMassHighestEt_EB[0];
  myH[2] = ScHisto_invMassHighestEt_EE[0];
  myH[3] = ScHisto_invMassHighestEt_EBEE[0];

  for (int ii=0; ii<4; ii++){ 
    int peakBin   = myH[ii]->GetMaximumBin();
    double h_norm = myH[ii]->GetMaximum();
    double h_rms  = myH[ii]->GetRMS();
    double h_peak = myH[ii]->GetBinCenter(peakBin);

    // gaussian fit to initialize
    TF1 *gausa = new TF1 ("gausa","[0]*exp(-1*(x-[1])*(x-[1])/2/[2]/[2])",0.6);
    gausa->SetParameters(h_norm,h_peak,h_rms);
    myH[ii]->Fit("gausa","","",h_peak-3*h_rms,h_peak+3*h_rms);
    double gausNorm  = gausa->GetParameter(0);
    double gausMean  = gausa->GetParameter(1);
    double gausSigma = fabs(gausa->GetParameter(2));
    double gausChi2  = gausa->GetChisquare()/gausa->GetNDF();
    if (gausChi2>100){ gausMean = h_peak; gausSigma = h_rms; }

    TF1 *crf = new TF1("crf",cruijffFunction,0.,6.,5);
    crf->SetParNames ("Mean","Sigma","alpha1","alpha2","Norm","Constant");
    crf->SetParameter(0, gausMean);
    crf->SetParameter(1, gausSigma);
    crf->SetParameter(2, 1.); //2
    crf->SetParameter(3, 0.); //2
    crf->SetParameter(4, gausNorm);
    myH[ii]->Fit("crf","lR","",1.,5.);
  }

  TCanvas c2("c2","signal",1);  
  c2.Divide(2,2);
  c2.cd(1); myH[0]->Draw();
  c2.cd(2); myH[1]->Draw();
  c2.cd(3); myH[2]->Draw();
  c2.cd(4); myH[3]->Draw();
  c2.SaveAs("regionsFit.eps");
  c2.SaveAs("regionsFit.root");
}


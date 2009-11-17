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

// vectors to be used in the following
TFile *file[8];
float expected[8];

// histos declaration
TH1F *ScHisto_invMassHighestEt[8];
TH1F *ScHisto_invMassHighestEt_EB[8];
TH1F *ScHisto_invMassHighestEt_EE[8];
TH1F *ScHisto_invMassHighestEt_EBEE[8];
TH1F *ScHisto_invMassHighestEt_norm1[8];
TH1F *ScHisto_invMassHighestEt_EB_norm1[8];
TH1F *ScHisto_invMassHighestEt_EE_norm1[8];
TH1F *ScHisto_invMassHighestEt_EBEE_norm1[8];
TH1F *ScHisto_invMassHighestEt_normL[8];
TH1F *ScHisto_invMassHighestEt_EB_normL[8];
TH1F *ScHisto_invMassHighestEt_EE_normL[8];
TH1F *ScHisto_invMassHighestEt_EBEE_normL[8];
TH1F *ScHisto_invMassHighestEt_normL_all;

void chargeFiles() {
  file[0] = new TFile("Outfile_Jpsi.root");
  file[1] = new TFile("Outfile_bce2030.root");
  file[2] = new TFile("Outfile_bce3080.root");
  file[3] = new TFile("Outfile_bce80170.root");
  file[4] = new TFile("Outfile_EM20to30.root");
  file[5] = new TFile("Outfile_EM30to80.root");
  file[6] = new TFile("Outfile_EM80to170.root");
  file[7] = new TFile("Outfile_doubleEM.root");
}

void setExpected(int theStep) {

  float lumi = 1.;   // pb-1

  // float exp_         =           kineEff  * filterEff_prod*crossSection;
  // cout << "after HLT  : " << 10.*exp_hlt      << endl; 

  if (theStep==1) {     // after HLT only
    expected[0] =  9411.04*lumi;
    expected[1] =  160481.*lumi;
    expected[2] =  205563.*lumi;
    expected[3] =  36137.3*lumi;
    expected[4] =  733553.*lumi;
    expected[5] =  3746500.*lumi;
    expected[6] =  381960.*lumi;
    expected[7] =  4402330.*lumi;
  }

  if (theStep==2) {     // after HLT + reco
    expected[0] = 8192.73*lumi;
    expected[1] = 73696.3*lumi;
    expected[2] = 94399.1*lumi;
    expected[3] = 20059.1*lumi;
    expected[4] = 102448*lumi;
    expected[5] = 785246*lumi;
    expected[6] = 120682*lumi;
    expected[7] = 296849*lumi;
  }

  if (theStep==3) {     // after HLT + reco + Et>4 
    expected[0] = 8006.41*lumi;
    expected[1] = 61305.3*lumi;
    expected[2] = 78527.2*lumi;
    expected[3] = 16578.6*lumi;
    expected[4] = 68837.1*lumi;
    expected[5] = 557055*lumi;
    expected[6] = 85397*lumi;
    expected[7] = 229638*lumi;
  }

  if (theStep==4) {     // after HLT + reco + Et>4 + charge 
    expected[0] = 7831.54*lumi;
    expected[1] = 38579.5*lumi;
    expected[2] = 49417.2*lumi;
    expected[3] = 10395.7*lumi;
    expected[4] = 42469*lumi;
    expected[5] = 327059*lumi; 
    expected[6] = 50312.9*lumi; 
    expected[7] = 162427*lumi;
  }

  if (theStep==5) {     // after HLT + reco + Et>4 + charge + eleID 
    expected[0] =   3769.58*lumi;
    expected[1] = 1002.55*lumi;
    expected[2] = 1284.1*lumi;
    expected[3] = 58.98*lumi;
    expected[4] = 3216.6*lumi;
    expected[5] = 5860.14*lumi;
    expected[6] = 295.439*lumi;
    expected[7] = 33605.6*lumi;
  }

  if (theStep==6) {     // after HLT + reco + Et>4 + charge + eleID + isol 
    expected[0] = 2350.*lumi;
    expected[1] =  55.9*lumi;
    expected[2] =  71.6*lumi;
    expected[3] =  2.73*lumi;
    expected[4] = 547.0*lumi;
    expected[5] = 582.28*lumi;
    expected[6] = 25.97*lumi;
    expected[7] = 5600.9*lumi;
  }
}

void readHistos() {
  for(int ifile=0; ifile<8; ifile++){
    ScHisto_invMassHighestEt[ifile]      = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt");
    ScHisto_invMassHighestEt_EB[ifile]   = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EB");
    ScHisto_invMassHighestEt_EE[ifile]   = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EE");
    ScHisto_invMassHighestEt_EBEE[ifile] = (TH1F*)file[ifile]->Get("ScHisto_invMassHighestEt_EBEE");
  }
}

void normalizeTo1() { 

  char title[700];

  for(int ifile=0; ifile<8; ifile++){
    sprintf(title,"ScHisto_invMassHighestEt_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_norm1[ifile]      = (TH1F*)ScHisto_invMassHighestEt[ifile]->Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EB_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EB_norm1[ifile]   = (TH1F*)ScHisto_invMassHighestEt_EB[ifile] -> Clone(title);  
    sprintf(title,"ScHisto_invMassHighestEt_EE_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EE_norm1[ifile]   = (TH1F*)ScHisto_invMassHighestEt_EE[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EBEE_norm1[%d]", ifile);
    ScHisto_invMassHighestEt_EBEE_norm1[ifile] = (TH1F*)ScHisto_invMassHighestEt_EBEE[ifile] -> Clone(title);
  }

  for(int ifile=0; ifile<8; ifile++){
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
  for(int ifile=0; ifile<8; ifile++){
    sprintf(title,"ScHisto_invMassHighestEt_normL[%d]", ifile);
    ScHisto_invMassHighestEt_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EB_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EB_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EB[ifile] -> Clone(title);  
    sprintf(title,"ScHisto_invMassHighestEt_EE_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EE_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EE[ifile] -> Clone(title);
    sprintf(title,"ScHisto_invMassHighestEt_EBEE_normL[%d]", ifile);
    ScHisto_invMassHighestEt_EBEE_normL[ifile] = (TH1F*)ScHisto_invMassHighestEt_EBEE[ifile] -> Clone(title);
  }

  for(int ifile=0; ifile<8; ifile++){
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
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[1]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[2]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[3]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[4]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[5]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[6]);
  ScHisto_invMassHighestEt_normL_all -> Add(ScHisto_invMassHighestEt_normL[7]);
}

void cosmeticsSvsB() {
  int color;
  for(int ifile=0; ifile<8; ifile++){
    color = ifile+2;
    
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
  readHistos();
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
  leg.AddEntry(ScHisto_invMassHighestEt_normL[1], "bc->e, 20-30"      ,  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[2], "bc->e, 30-80"      ,  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[3], "bc->e, 80-170",       "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[4], "em.enriched,  20-30", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[5], "em.enriched,  30-80", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[6], "em.enriched, 80-170", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[7], "double em.enr  6-20", "f");
  leg.SetFillColor(0);
  leg.SetBorderSize(0.4);

  TCanvas c1("c1","all",1);  
  ScHisto_invMassHighestEt_normL[0]->Draw();
  for(int ifile=0; ifile<8; ifile++) ScHisto_invMassHighestEt_normL[ifile] -> Draw("same");
  leg.Draw();
  c1.SaveAs("normL.eps");
  c1.SaveAs("normL.root");

  TCanvas c11("c11","all",1);  
  ScHisto_invMassHighestEt_norm1[0]->Draw();
  for(int ifile=0; ifile<8; ifile++) ScHisto_invMassHighestEt_norm1[ifile] -> Draw("same");
  leg.Draw();
  c11.SaveAs("norm1.eps");
  c11.SaveAs("norm1.root");

  TCanvas c21("c21","all",1);  
  ScHisto_invMassHighestEt_normL_all->Draw();
  c21.SaveAs("insieme.eps");
  

  // summing up backgrounds and signals
  THStack staB_invMassHighestEt("staB_invMassHighestEt","staB_invMassHighestEt");
  staB_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[7]);
  for (int ii=1; ii<7; ii++) {
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

}


// analisi signal: EB vs EE vs EBEB vs all
void drawRegions() {

  chargeFiles();
  readHistos();
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


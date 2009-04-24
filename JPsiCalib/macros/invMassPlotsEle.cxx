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
  file[0] = new TFile("signalPrompt.root");
  file[1] = new TFile("signalNonPrompt.root");
  file[2] = new TFile("bce2030.root");
  file[3] = new TFile("bce3080.root");
  file[4] = new TFile("bce80170.root");
  file[5] = new TFile("emEnrich3080.root");
  file[6] = new TFile("emEnrich80170.root");
  file[7] = new TFile("minBias.root");
}

void setExpected(int theStep) {

  float lumi = 100.;   // pb-1

  if (theStep==1) {     // after HLT only
    expected[0] =    1644.*lumi;
    expected[1] =     286.*lumi;
    expected[2] =    6624.*lumi;
    expected[3] =   26856.*lumi;
    expected[4] =    4676.*lumi;
    expected[5] =  321480.*lumi;
    expected[6] =   36275.*lumi;
    expected[7] = 5947120.*lumi;
  }

  if (theStep==2) {     // after HLT + reco
    expected[0] =   1168.*lumi;
    expected[1] =    168.*lumi;
    expected[2] =   1259.*lumi;
    expected[3] =   7174.*lumi;
    expected[4] =   1750.*lumi;
    expected[5] =  32066.*lumi;
    expected[6] =   8430.*lumi;
    expected[7] = 471329.*lumi;
  }

  if (theStep==3) {     // after HLT + reco + Et>4 
    expected[0] =   1167.*lumi;
    expected[1] =    167.*lumi;
    expected[2] =   1245.*lumi;
    expected[3] =   7105.*lumi;
    expected[4] =   1738.*lumi;
    expected[5] =  31518.*lumi;
    expected[6] =   8287.*lumi;
    expected[7] = 460410.*lumi;
  }

  if (theStep==4) {     // after HLT + reco + Et>4 + charge 
    expected[0] =   1142.*lumi;
    expected[1] =    158.*lumi;
    expected[2] =    819.*lumi;
    expected[3] =   4389.*lumi;
    expected[4] =   1069.*lumi;
    expected[5] =  18197.*lumi; 
    expected[6] =   4782.*lumi; 
    expected[7] = 298448.*lumi;
  }

  if (theStep==5) {     // after HLT + reco + Et>4 + charge + eleID 
    expected[0] =   803.*lumi;
    expected[1] =    94.*lumi;
    expected[2] =   222.*lumi;
    expected[3] =   740.*lumi;
    expected[4] =    86.*lumi;
    expected[5] =   659.*lumi;
    expected[6] =   115.*lumi;
    expected[7] = 10919.*lumi;
  }

  if (theStep==6) {     // after HLT + reco + Et>4 + charge + eleID + isol 
    expected[0] = 710.*lumi;
    expected[1] =  50.*lumi;
    expected[2] =  80.*lumi;
    expected[3] = 129.*lumi;
    expected[4] =   5.*lumi;
    expected[5] =  76.*lumi;
    expected[6] =   9.*lumi;
    expected[7] =   0.*lumi;
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
  leg.AddEntry(ScHisto_invMassHighestEt_normL[0], "prompt signal",      "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[1], "non prompt signal",  "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[2], "bc->e, 20-30",       "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[3], "bc->e, 30-80",       "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[4], "bc->e, 80-170",      "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[5], "em.enriched, 30-80", "f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[6], "em.enriched, 80-170","f");
  leg.AddEntry(ScHisto_invMassHighestEt_normL[7], "minimum bias","f");
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
  for (int ii=2; ii<7; ii++) {
    staB_invMassHighestEt.Add(ScHisto_invMassHighestEt_normL[ii]);
  }
  THStack staS_invMassHighestEt("staS_invMassHighestEt","staS_invMassHighestEt");
  for (int ii=0; ii<2; ii++) {
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


//// optimizer 

void OptimizeSingleCut(char basename[100]="hnuMsqr2",char file[100]="NuMass", bool addHistos = false, float xmin = 0.0, float xmax = 100.0)
{
  char namegif[100];
  char namefile[100];
  char source[10];
  char str[20];
  char title[100];
  char all[110];
  
  gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);

  // signal prompt / signal non-prompt / bce 2030 / bce 3080 / bce 80170 
  Float_t lumi[5] = {25.70,33.91,10.9,8.80,47.3};

  sprintf(namefile,"Outfile_signalprompt.root");
  TFile SPfile(namefile);
  sprintf(namefile,"Outfile_signalnonprompt.root");
  TFile SNPfile(namefile);
  sprintf(namefile,"Outfile_bce2030.root");
  TFile Bkgfile1(namefile);
  sprintf(namefile,"Outfile_bce3080.root");
  TFile Bkgfile2(namefile);
  sprintf(namefile,"Outfile_bce80170.root");
  TFile Bkgfile3(namefile);
  
  ///////////// SP ////////////////
  SPfile.cd();
  
  TH1F *hS =((TH1F*)gDirectory->Get(basename))->Clone();
  hS->SetName("hSP");
  hS->Scale(10./lumi[0]);
  
  ///////////// SNP ////////////////
  SNPfile.cd();
  
  TH1F *hSNP =((TH1F*)gDirectory->Get(basename))->Clone();
  hSNP->SetName("hSNP");
  hSNP->Scale(10./lumi[1]);

  ///////////// bkg ////////////////
  Bkgfile1.cd();
  
  TH1F *hBkg =((TH1F*)gDirectory->Get(basename))->Clone();
  hBkg->SetName("hBkg");
  hBkg->Scale(10./lumi[2]);

  Bkgfile2.cd();
  
  TH1F *hBkg2 =((TH1F*)gDirectory->Get(basename))->Clone();
  hBkg2->SetName("hBkg2");
  hBkg2->Scale(10./lumi[3]);

  Bkgfile3.cd();
  
  TH1F *hBkg3 =((TH1F*)gDirectory->Get(basename))->Clone();
  hBkg3->SetName("hBkg3");
  hBkg3->Scale(10./lumi[4]); 
  
  ////////////////// Sum histos////////////////////// 
  hS->Add(hSNP);
  hBkg->Add(hBkg2);
  hBkg->Add(hBkg3);

  TCanvas *c1 = new TCanvas();
 
  if (addHistos) {
    int firstbin = 0;
    int lastbin = 100;
    for (int ii=1; ii < hS->GetNbinsX()-1; ii++) {
      if (hS->GetBinLowEdge(ii) <= xmin && hS->GetBinLowEdge(ii+1) >= xmin) firstbin = ii;
      if (hS->GetBinLowEdge(ii) <= xmax && hS->GetBinLowEdge(ii+1) >= xmax) lastbin = ii; 
    } 
    cout << "Signal events in region: " << hS->Integral(firstbin,lastbin) << endl;
    cout << "Background events in region: " << hBkg->Integral(firstbin,lastbin) << endl; 
    hS->Add(hBkg);
    hS->SetLineWidth(2);
    hS->SetLineColor(kRed);
    hBkg->SetLineWidth(2);
    hBkg->SetLineColor(kBlack); 
    hS->Draw();
    hBkg->Draw("same");
  } else {
    hS->SetLineWidth(2);
    hS->SetLineColor(kRed);
    hBkg->SetLineWidth(2);
    hBkg->SetLineColor(kBlack); 
    hS->DrawNormalized();
    hBkg->DrawNormalized("same");
  }
 
  leg = new TLegend(0.65, 0.6, 0.85, 0.8);
  // leg->AddEntry(hData,"Data","f");
  leg->AddEntry(hS,"Signal J/#psi","f");
  leg->AddEntry(hBkg,"b,c -> e","f");
  leg->Draw("same"); 

  sprintf(namegif,"%s.gif",file);
  c1->Print(namegif); 

  /* SPfile.Close();
  SNPfile.Close();
  Bkgfile1.Close();
  Bkgfile2.Close();
  Bkgfile3.Close();*/
 
  //////////////// CUT "OPTIMIZATION" ////////////////
  if (!addHistos) {
    for (int i=1; i < hS->GetNbinsX(); i++) {
      float seg = hS->Integral(1,i);
      float fon = hBkg->Integral(1,i);
      cout << "Cut: " << hS->GetBinLowEdge(i) << " \t S = " << seg << " \t B = " << fon << " \t Signif = " << ( seg+fon>0 ? seg/sqrt(seg+fon) : 0.) << endl;
      // cout << "Cut: " << hS->GetBinLowEdge(i) << " \t S = " << seg << " \t B = " << fon << " \t Purity = " << ( fon>0 ? seg/fon : 0.) << endl;
    }  
  }
}

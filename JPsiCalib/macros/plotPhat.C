void plotPhat() {
  gStyle->SetOptStat(0); 

  TFile *_file0 = TFile::Open("Outfile_qcdpt15.root");
  //TFile *_file0 = TFile::Open("Outfile_minbias.root");
  //HepHisto_ptHatAll->Draw();
  //c1->SetLogy(1);
  HepHisto_ptHat1->SetLineColor(2);
  HepHisto_ptHat1->Rebin(4);
  HepHisto_ptHat1->Draw();
  //HepHisto_ptHat1->Draw("same");
  c1->SetLogy(1);
  HepHisto_ptHat2->SetLineColor(4);
  HepHisto_ptHat2->Rebin(4);
  HepHisto_ptHat2->Draw("same");
  HepHisto_ptHat3->SetLineColor(5);
  HepHisto_ptHat3->Rebin(4);
  HepHisto_ptHat3->Draw("same");
  HepHisto_ptHat4->SetLineColor(6);
  HepHisto_ptHat4->Rebin(4);
  HepHisto_ptHat4->Draw("same");
  HepHisto_ptHat5->SetLineColor(8);
  HepHisto_ptHat5->Rebin(4);
  HepHisto_ptHat5->Draw("same");
  HepHisto_ptHat6->SetLineColor(9);
  HepHisto_ptHat6->Rebin(4);
  HepHisto_ptHat6->Draw("same");
  c1->BuildLegend();
}

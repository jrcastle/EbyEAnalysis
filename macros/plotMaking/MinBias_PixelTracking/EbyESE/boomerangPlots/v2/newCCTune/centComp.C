#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void centComp(){

  bool test = 0;

  TFile * fDefault;
  TTree * treeDefault;
  double centvalDefault;
  TH1D * hCentDefault;

  TFile * fNewCC;
  TTree * treeNewCC;
  double centvalNewCC;
  TH1D * hCentNewCC;

  TH1D * hCentRatio;

  //
  // MAIN
  //
  setTDRStyle();
  TH1D::SetDefaultSumw2();

  fDefault = new TFile("MultCent_NoCCFilter.root");
  treeDefault = (TTree*) fDefault->Get("multcentana/tree");
  treeDefault->SetBranchAddress("Cent", &centvalDefault);

  fNewCC = new TFile("MultCent_NewCCFilter_2pct.root");
  treeNewCC = (TTree*) fNewCC->Get("multcentana/tree");
  treeNewCC->SetBranchAddress("Cent", &centvalNewCC);

  hCentDefault = new TH1D("hCentDefault", "hCentDefault", 200, 0, 100);
  hCentDefault->GetXaxis()->SetTitle("Centrality %");
  hCentDefault->GetYaxis()->SetTitle("Events");
  hCentDefault->SetLineColor(1);
  hCentDefault->SetMarkerColor(1);
  hCentDefault->SetMarkerStyle(20);

  hCentNewCC = new TH1D("hCentNewCC", "hCentNewCC", 200, 0, 100);
  hCentNewCC->GetXaxis()->SetTitle("Centrality %");
  hCentNewCC->GetYaxis()->SetTitle("Events");
  hCentNewCC->SetLineColor(2);
  hCentNewCC->SetMarkerColor(2);
  hCentNewCC->SetMarkerStyle(21);


  int N;
  if(test) N = 1000000;
  //else     N = treeNewCC->GetEntries();
  else     N = treeDefault->GetEntries();

  for(int ievent = 0; ievent < N; ievent++){

    if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<100.*(ievent+1)/(double)N<<"% Completed"<<endl;

    treeDefault->GetEntry(ievent);
    hCentDefault->Fill( centvalDefault );

    if( ievent < treeNewCC->GetEntries() ){  
      treeNewCC->GetEntry(ievent);
      hCentNewCC->Fill( centvalNewCC );
    }

  }

  //-- Normalize
  hCentDefault->Scale( 1./hCentDefault->Integral() );
  hCentNewCC->Scale( 1./hCentNewCC->Integral() );

  hCentRatio = (TH1D*) hCentNewCC->Clone("hCentRatio");
  hCentRatio->GetXaxis()->SetTitle("Centrality %");
  hCentRatio->GetYaxis()->SetTitle("Ratio: NewTune/OldTune");
  hCentRatio->Divide( hCentDefault );
  hCentRatio->SetMinimum(0.8);
  hCentRatio->SetMaximum(1.2);

  //-- Draw
  TLegend * leg = new TLegend(0.1775, 0.2366, 0.4781, 0.3380);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hCentDefault, "Old Tune", "lp");
  leg->AddEntry(hCentNewCC,   "New Tune", "lp");

  TCanvas * c = new TCanvas("c", "c", 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  c->cd(1)->SetLogy();
  hCentDefault->Draw();
  hCentNewCC->Draw("same");
  leg->Draw("same");
  c->cd(2);
  hCentRatio->Draw();
  c->SaveAs("centCompare.pdf");


}

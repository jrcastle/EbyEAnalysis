#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void clusTree(){

  const int NCENT = 11; 
  const int cmin[NCENT]  = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int cmax[NCENT]  = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const double cbins[]   = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const int ccol[NCENT]  = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};
  const int Nbins[NCENT] = {300, 250, 250, 250, 175, 125, 100, 100, 100, 75, 50};
  const int mmax[NCENT]  = {6000, 5000, 4500, 3500, 2750, 2250, 1000, 500, 250, 150, 100};
  const double CCMax     = 8.0;
  TH1D * hCent = new TH1D("hCent", "hCent", NCENT, cbins);


  TFile * f;
  TTree * tree;
  double centval;
  unsigned int Mult_pt03_3_eta24;
  unsigned int Mult_pt03_3_eta10;

  TH1D * hMult24[NCENT];
  TH1D * hMult10[NCENT];

  TCanvas * cMultCent;
  TH1D * hCent2;

  TH2D * hCCvsNpix;
  TF1 * polyCut;
  TF1 * lineCut;
  TF1 * lineCutHigh;

  TLatex latex;

  //
  // MAIN
  //
  TH1D::SetDefaultSumw2();
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  f = new TFile("MultCent_NewCCFilter_2pct.root");
  //f = new TFile("MultCent_NoCCFilter.root");
  tree = (TTree*) f->Get("multcentana/tree");
  tree->SetBranchAddress("Cent",              &centval);
  tree->SetBranchAddress("Mult_pt03_3_eta24", &Mult_pt03_3_eta24);
  tree->SetBranchAddress("Mult_pt03_3_eta10", &Mult_pt03_3_eta10);

  hCent2 = new TH1D("hCent2", "hCent2", 200, 0, 100);
  hCent2->GetXaxis()->SetTitle("Centrality %");
  hCent2->GetYaxis()->SetTitle("Events");

  hCCvsNpix = (TH2D*) f->Get("multcentana/hClusVtxCompInteg");
  hCCvsNpix->GetYaxis()->SetTitle("Clus-Vtx Compatibility");
  hCCvsNpix->GetYaxis()->SetTitleSize(0.06);
  hCCvsNpix->GetYaxis()->SetLabelSize(0.05);
  hCCvsNpix->GetXaxis()->SetTitle("N_{Pixel}");
  hCCvsNpix->GetXaxis()->SetTitleSize(0.06);
  hCCvsNpix->GetXaxis()->SetLabelSize(0.05);
  hCCvsNpix->GetXaxis()->SetNdivisions(509);
  int bccmax = hCCvsNpix->GetYaxis()->FindBin(CCMax);
  hCCvsNpix->GetYaxis()->SetRange(1, bccmax);

  polyCut = new TF1("polyCut", "pol5", 1000, 60000);
  polyCut->SetParameters(3.01672, 5.63152e-05, -3.82712e-09, 1.31546e-13, -2.2803e-18, 1.57971e-23);
  polyCut->SetLineWidth(2);
  polyCut->SetLineColor(8);

  lineCut = new TF1("lineCut", "pol1", 150, 1000);
  lineCut->SetParameters(2.14116, 0.000928176);
  lineCut->SetLineWidth(2);
  lineCut->SetLineColor(8);

  lineCutHigh = new TF1("lineCutHigh", "pol1", 60000, 100000);
  lineCutHigh->SetParameters(polyCut->Eval(60000), 0.);
  lineCutHigh->SetLineWidth(2);
  lineCutHigh->SetLineColor(8); 


  //-- Initialize Event Counters
  for(int c = 0; c < NCENT; c++){
    hMult24[c] = new TH1D(Form("hMult24_c%i_%i", cmin[c], cmax[c]), Form("hMult24_c%i_%i", cmin[c], cmax[c]), Nbins[c], 0, mmax[c]);
    hMult24[c]->GetXaxis()->SetNdivisions(509);
    hMult24[c]->GetXaxis()->SetTitle("Event Multiplicity");
    hMult24[c]->GetYaxis()->SetTitle("N_{Events}");
    hMult24[c]->SetMarkerColor(1);
    hMult24[c]->SetMarkerStyle(6);

    hMult10[c] = new TH1D(Form("hMult10_c%i_%i", cmin[c], cmax[c]), Form("hMult10_c%i_%i", cmin[c], cmax[c]), Nbins[c], 0, mmax[c]);
    hMult10[c]->GetXaxis()->SetNdivisions(509);
    hMult10[c]->GetXaxis()->SetTitle("Event Multiplicity");
    hMult10[c]->GetYaxis()->SetTitle("N_{Events}");
    hMult10[c]->SetMarkerColor(1);
    hMult10[c]->SetMarkerStyle(6);
  }

  int N = tree->GetEntries();
  for(int ievent = 0; ievent < N; ievent++) {
        
    if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<100.*(ievent+1)/(double)N<<"% Completed"<<endl;
    tree->GetEntry(ievent);
    int icent = hCent->FindBin(centval)-1;
    hCent2->Fill(centval);
    hMult24[icent]->Fill(Mult_pt03_3_eta24);
    hMult10[icent]->Fill(Mult_pt03_3_eta10);

  } //-- End event loop


  //-- Make plots

  //-- Group by Centrality
  TLegend * legCent = new TLegend(0.75, 0.2, 0.95, 0.5);
  legCent->SetBorderSize(0);
  legCent->SetFillStyle(0);

  cMultCent = new TCanvas("cMultCent", "cMultCent", 500, 500);
  cMultCent->SetLogy();

  for(int c = 0; c < NCENT; c++){
    hMult24[c]->SetLineColor(ccol[c]);
    hMult24[c]->SetMaximum(20.*hMult24[0]->GetMaximum());
    legCent->AddEntry(hMult24[c], Form("Cent %i-%i%s", cmin[c], cmax[c], "%"), "lp");

    cMultCent->cd();
    if(c==0) hMult24[c]->Draw();
    else     hMult24[c]->Draw("same");
  }

  cMultCent->cd();
  legCent->Draw("same");
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  cMultCent->SaveAs( "cMultCent.png" );

  TCanvas * cCent = new TCanvas("cCent", "cCent", 500, 500);
  cCent->cd();
  cCent->SetLogy();
  hCent2->Draw();
  cCent->SaveAs("hCent.png");

  TCanvas * cCCvsNpix = new TCanvas("cCCvsNpix", "cCCvsNpix", 500, 500);
  cCCvsNpix->cd()->SetRightMargin(0.1);
  hCCvsNpix->Draw();
  polyCut->Draw("same");
  lineCut->Draw("same");
  lineCutHigh->Draw("same");

  cCCvsNpix->Update();
  TPaletteAxis* palette = (TPaletteAxis*) hCCvsNpix->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.75);
  palette->SetY1NDC(0.6);
  palette->SetX2NDC(0.8);
  palette->SetY2NDC(0.9);

  cCCvsNpix->SaveAs("cCCvsNpix_NewTune.pdf");

}

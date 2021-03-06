#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void illustrateFilterCut(){

  const int NCLUS = 8;
  string clusTrunc[NCLUS]          = {"2p0", "2p3", "2p6", "2p9", "3p2", "3p5", "3p8", "4p1"};
  const double clusTruncVal[NCLUS] = {2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1};
  const double clusterPar_         = 0.006;
  const int col[]                  = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};

  double NpixMax = 5000;
  double CCMax   = 5.0;

  TFile * f;
  TH2D * hCCvsNpix;
  TLine * l0[NCLUS];
  TLine * l1[NCLUS];
  TLine * l2[NCLUS];
  TLine * l3[NCLUS];

  //
  // MAIN
  //
  setTDRStyle();

  f = new TFile("MultCent_Pixel.root");
  hCCvsNpix = (TH2D*) f->Get("multcentana/hClusVtxCompInteg");
  int binNpixMax = hCCvsNpix->GetXaxis()->FindBin(NpixMax);
  int binCCMax   = hCCvsNpix->GetYaxis()->FindBin(CCMax);
  hCCvsNpix->GetXaxis()->SetRange(1, binNpixMax);
  hCCvsNpix->GetYaxis()->SetRange(1, binCCMax);

  TCanvas * c = new TCanvas("c", "c", 500, 500);
  c->cd();
  c->SetLogz();
  hCCvsNpix->Draw();


  for(int i = 0; i < NCLUS; i++){

    c->cd();
    l0[i] = new TLine(0.0, 0.0, 150.0, 0.0);
    l0[i]->SetLineColor( col[i] );
    l0[i]->SetLineStyle( i );
    l0[i]->SetLineWidth( 2 );
    l0[i]->Draw("same");

    l1[i] = new TLine(150.0, 0.0, 150.0, 0.9);
    l1[i]->SetLineColor( col[i] );
    l1[i]->SetLineStyle( i );
    l1[i]->SetLineWidth( 2 );
    l1[i]->Draw("same");

    l2[i] = new TLine(150., 0.9, clusTruncVal[i]/clusterPar_, clusTruncVal[i]);
    l2[i]->SetLineColor( col[i] );
    l2[i]->SetLineStyle( i );
    l2[i]->SetLineWidth( 2 );
    l2[i]->Draw("same");

    l3[i] = new TLine(clusTruncVal[i]/clusterPar_, clusTruncVal[i], NpixMax, clusTruncVal[i]);
    l3[i]->SetLineColor( col[i] );
    l3[i]->SetLineStyle( i );
    l3[i]->SetLineWidth( 2 );
    l3[i]->Draw("same");

  }

  TLegend * leg = new TLegend(0.4415, 0.1949, 0.9718, 0.4131);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  for(int i = 0; i < NCLUS; i++) leg->AddEntry( l0[i], Form("ClusterTrunc %.1f", clusTruncVal[i]), "l" );
  c->cd();
  leg->Draw("same");
  c->SaveAs("illustrateCuts.pdf");

}

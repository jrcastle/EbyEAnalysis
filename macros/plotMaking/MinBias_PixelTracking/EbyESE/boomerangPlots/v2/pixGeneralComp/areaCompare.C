#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void areaCompare(){

  const int N = 11;
  const int cmin[N] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int cmax[N] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const int ccol[N] = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};

  const int NCLUS = 7;
  string clusTrunc[NCLUS]     = {"2p0", "2p3", "2p6", "2p9", "3p2", "3p5", "3p8"};
  double clusTruncVal[NCLUS]  = {2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8};

  TFile * fMult[NCLUS];
  TH2D * hMult2D[NCLUS];
  TH1D * hMultProj[NCLUS][N];

  TCanvas * cDist[N];

  TLatex latex;

  //
  // MAIN
  //
  TH1D::SetDefaultSumw2();
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  for(int i = 0; i < NCLUS; i++){

    fMult[i] = new TFile( Form("MultCent_ClusVtx%s.root", clusTrunc[i].data()) );
    hMult2D[i] = (TH2D*) fMult[i]->Get( "multcentana/hMultCent_eta24_pt03_3_hiPixelAndGeneral" );

    for(int icent = 0; icent < N; icent++){
      hMultProj[i][icent] = (TH1D*) hMult2D[i]->ProjectionY( Form("hMultProj%s_c%i", clusTrunc[i].data(), icent), 2*cmin[icent], 2*cmax[icent] );
      hMultProj[i][icent]->Scale(1./hMultProj[i][icent]->Integral());
      hMultProj[i][icent]->SetMarkerColor(1);
      hMultProj[i][icent]->SetLineColor(ccol[i]);
    }

  }


  //-- Draw Distributions
  TLegend * leg = new TLegend(0.7, 0.2, 0.9, 0.5);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  for(int i =0; i < NCLUS; i++) leg->AddEntry(hMultProj[i][0], Form("Trunc %.1f", clusTruncVal[i]), "lp");

  for(int icent = 0; icent < N; icent++){

    cDist[icent] = new TCanvas( Form("cDist_c%i_%i", cmin[icent], cmax[icent]), Form("cDist_c%i_%i", cmin[icent], cmax[icent]), 500, 500);
    cDist[icent]->cd();
    cDist[icent]->SetLogy();
    for(int i = 0; i < NCLUS; i++){
      if( i == 0 ) hMultProj[NCLUS-1-i][icent]->Draw();
      else         hMultProj[NCLUS-1-i][icent]->Draw("same");
    }
    latex.DrawLatex(0.65, 0.88, Form("Cent %i-%i%s", cmin[icent], cmax[icent], "%"));
    leg->Draw("same");
    cDist[icent]->SaveAs( Form("cDist_c%i_%i.png", cmin[icent], cmax[icent]) );
  }


}

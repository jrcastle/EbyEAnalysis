#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TDecompSVD.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void svdUnfoldResults(){

  int norder_ = 3;
  double tkEta = 1.0;

  double cumuMin = 0.0;
  double cumuMax = 0.15;

  double g1eMin   = -1.0;
  double g1eMax   = 0.5;

  double ratioMin       = 0.9;
  double ratioMax       = 1.03;
  double ratioMinVn8Vn6 = 0.97;
  double ratioMaxVn8Vn6 = 1.002;

  TFile * fSVD;
  TH1D * finalUnfm1[NCENT];
  TH1D * finalUnf[NCENT];
  TH1D * finalUnfp1[NCENT];

  TFile * fStatSVD;
  TH1D * hVarianceOfMean_Vn2;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
  TH1D * hVarianceOfMean_Gamma1Exp;
  TH1D * hVarianceOfMean_Vn6Vn4;
  TH1D * hVarianceOfMean_Vn8Vn4;
  TH1D * hVarianceOfMean_Vn8Vn6;

  double Vn2[NCENT];
  double Vn4[NCENT];
  double Vn6[NCENT];
  double Vn8[NCENT];

  double Vn6Vn4[NCENT];
  double Vn8Vn4[NCENT];
  double Vn8Vn6[NCENT];
  double G1E[NCENT];

  double Vn2_err[NCENT];
  double Vn4_err[NCENT];
  double Vn6_err[NCENT];
  double Vn8_err[NCENT];

  double Vn6Vn4_err[NCENT];
  double Vn8Vn4_err[NCENT];
  double Vn8Vn6_err[NCENT];
  double G1E_err[NCENT];

  TGraphErrors * grVn2;
  TGraphErrors * grVn4;
  TGraphErrors * grVn6;
  TGraphErrors * grVn8;
  TGraphErrors * grVn6Vn4;
  TGraphErrors * grVn8Vn4;
  TGraphErrors * grVn8Vn6;
  TGraphErrors * grG1E;

  TFile * fOut;

  TH1D * hObs[NCENT];
  TH1D * hrefoldm1[NCENT];
  TH1D * hrefold[NCENT];
  TH1D * hrefoldp1[NCENT];
  TH2D * hResponse[NCENT];
  TH2D * hCovMatrix[NCENT];
  TH2D * hCorrMatrix[NCENT];
  TCanvas * cCovRespUnf[NCENT];

  TH1D * hKreg;
  TH1D * hdi[NCENT];
  TLine * lKreg[NCENT];
  TLine * lOne[NCENT];

  TLatex latex;
  TLatex latex2;
  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);
  latex2.SetNDC();

  //-- Get Stat Errors
  fStatSVD = new TFile( Form("../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertaintiesSVD_v%i.root", norder_, tkEta, norder_)  );
  hVarianceOfMean_Vn2       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn2");
  hVarianceOfMean_Vn4       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn4");
  hVarianceOfMean_Vn6       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn6");
  hVarianceOfMean_Vn8       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8");
  hVarianceOfMean_Gamma1Exp = (TH1D*) fStatSVD->Get("hVarianceOfMean_Gamma1Exp");
  hVarianceOfMean_Vn6Vn4    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn6Vn4");
  hVarianceOfMean_Vn8Vn4    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8Vn4");
  hVarianceOfMean_Vn8Vn6    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8Vn6");

  //-- Get unfold file
  fSVD = new TFile("../UnfoldResults/dataResp/data2_svd.root");

  //-- Get kregs
  hKreg = (TH1D*) fSVD->Get( "hKreg" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the observed distn
    hObs[icent] = (TH1D*) fSVD->Get( Form("hVnFull_c%i", icent) );
    hObs[icent]->SetLineColor(1);
    hObs[icent]->SetMarkerColor(1);
    hObs[icent]->SetMarkerStyle(24);

    //-- Get the Response matrix
    hResponse[icent] = (TH2D*) fSVD->Get( Form("hresp_c%i", icent) );

    //-- Grab the cov matrix
    hCovMatrix[icent] = (TH2D*) fSVD->Get( Form("hCovMatkreg4_c%i", icent) );

    //-- Grab the di vectors
    hdi[icent] = (TH1D*) fSVD->Get( Form("hdi_c%i", icent) );
    hdi[icent]->GetXaxis()->SetRange(1,25);
    hdi[icent]->GetXaxis()->SetTitle( "i" );
    hdi[icent]->GetYaxis()->SetTitle( "|d_{i}|" );

    int nb = hdi[icent]->GetNbinsX();
    double dimin = 0.5*hdi[icent]->GetMinimum();
    double dimax = 2.*hdi[icent]->GetMaximum();
    hdi[icent]->SetMinimum(dimin);
    hdi[icent]->SetMaximum(dimax);
    double lmin = hdi[icent]->GetBinLowEdge(1);
    double lmax = 25;
    double kreg = hKreg->GetBinContent(icent+1);
    lOne[icent] = new TLine(lmin, 1., lmax, 1.);
    lOne[icent]->SetLineColor(1);
    lOne[icent]->SetLineStyle(2);
    lOne[icent]->SetLineWidth(2);

    lKreg[icent] = new TLine(kreg, dimin, kreg, dimax);
    lKreg[icent]->SetLineColor(2);
    lKreg[icent]->SetLineWidth(2);

    //-- refold
    hrefoldm1[icent] = (TH1D*) fSVD->Get( Form("hrefoldkreg3_c%i", icent) );
    hrefoldm1[icent]->SetLineColor(2);
    hrefoldm1[icent]->SetLineWidth(2);

    hrefold[icent] = (TH1D*) fSVD->Get( Form("hrefoldkreg4_c%i", icent) );
    hrefold[icent]->SetLineColor(4);
    hrefold[icent]->SetLineWidth(2);

    hrefoldp1[icent] = (TH1D*) fSVD->Get( Form("hrefoldkreg5_c%i", icent) );
    hrefoldp1[icent]->SetLineColor(8);
    hrefoldp1[icent]->SetLineWidth(2);

    finalUnfm1[icent] = (TH1D*) fSVD->Get( Form("hrecokreg3_c%i", icent) );
    finalUnfm1[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    finalUnfm1[icent]->SetLineColor(2);
    finalUnfm1[icent]->SetMarkerColor(2);
    FixUnfold( finalUnfm1[icent] );

    finalUnf[icent] = (TH1D*) fSVD->Get( Form("hrecokreg4_c%i", icent) );
    finalUnf[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    finalUnf[icent]->SetLineColor(4);
    finalUnf[icent]->SetMarkerColor(4);
    FixUnfold( finalUnf[icent] );

    finalUnfp1[icent] = (TH1D*) fSVD->Get( Form("hrecokreg5_c%i", icent) ); 
    finalUnfp1[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    finalUnfp1[icent]->SetLineColor(8);
    finalUnfp1[icent]->SetMarkerColor(8);
    FixUnfold( finalUnfp1[icent] );

    EbyECumu cumu( finalUnf[icent] );
    double vn2 = cumu.GetCumu_vn2();
    double vn4 = cumu.GetCumu_vn4();
    double vn6 = cumu.GetCumu_vn6();
    double vn8 = cumu.GetCumu_vn8();
    double g1e = cumu.GetGamma1Exp();

    double vn6vn4;
    double vn8vn4;
    double vn8vn6;

    if(vn4 == 0 || vn6 == 0) vn6vn4 = -1;
    else                     vn6vn4 = vn6 / vn4;
    if(vn4 == 0 || vn8 == 0) vn8vn4 = -1;
    else                     vn8vn4 = vn8 / vn4;
    if(vn6 == 0 || vn8 == 0) vn8vn6 = -1;
    else                     vn8vn6 = vn8 / vn6;

    double vn2e    = sqrt( hVarianceOfMean_Vn2->GetBinContent(icent+1) );
    double vn4e    = sqrt( hVarianceOfMean_Vn4->GetBinContent(icent+1) );
    double vn6e    = sqrt( hVarianceOfMean_Vn6->GetBinContent(icent+1) );
    double vn8e    = sqrt( hVarianceOfMean_Vn8->GetBinContent(icent+1) );
    double vn6vn4e = sqrt( hVarianceOfMean_Vn6Vn4->GetBinContent(icent+1) );
    double vn8vn4e = sqrt( hVarianceOfMean_Vn8Vn4->GetBinContent(icent+1) );
    double vn8vn6e = sqrt( hVarianceOfMean_Vn8Vn6->GetBinContent(icent+1) );
    double g1ee    = sqrt( hVarianceOfMean_Gamma1Exp->GetBinContent(icent+1) );

    Vn2[icent] = vn2;
    Vn4[icent] = vn4;
    Vn6[icent] = vn6;
    Vn8[icent] = vn8;

    Vn6Vn4[icent] = vn6vn4;
    Vn8Vn4[icent] = vn8vn4;
    Vn8Vn6[icent] = vn8vn6;
    G1E[icent]    = g1e;

    Vn2_err[icent] = vn2e;
    Vn4_err[icent] = vn4e;
    Vn6_err[icent] = vn6e;
    Vn8_err[icent] = vn8e;

    Vn6Vn4_err[icent] = vn6vn4e;
    Vn8Vn4_err[icent] = vn8vn4e;
    Vn8Vn6_err[icent] = vn8vn6e;
    G1E_err[icent]    = g1ee;

  }

  //- Tgraph time
  grVn2    = new TGraphErrors(NCENT, centBinCenter, Vn2,    CERR, Vn2_err);
  grVn4    = new TGraphErrors(NCENT, centBinCenter, Vn4,    CERR, Vn4_err);
  grVn6    = new TGraphErrors(NCENT, centBinCenter, Vn6,    CERR, Vn6_err);
  grVn8    = new TGraphErrors(NCENT, centBinCenter, Vn8,    CERR, Vn8_err);
  grVn6Vn4 = new TGraphErrors(NCENT, centBinCenter, Vn6Vn4, CERR, Vn6Vn4_err);
  grVn8Vn4 = new TGraphErrors(NCENT, centBinCenter, Vn8Vn4, CERR, Vn8Vn4_err);
  grVn8Vn6 = new TGraphErrors(NCENT, centBinCenter, Vn8Vn6, CERR, Vn8Vn6_err);
  grG1E    = new TGraphErrors(NCENT, centBinCenter, G1E,    CERR, G1E_err);

  formatGraph(grVn2,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{2}", norder_),                      9,         24, "grVn2");
  formatGraph(grVn4,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{4}", norder_),                      kSpring+4, 24, "grVn4");
  formatGraph(grVn6,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{6}", norder_),                      6,         24, "grVn6");
  formatGraph(grVn8,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{8}", norder_),                      kOrange+7, 24, "grVn8");
  formatGraph(grVn6Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_), 4,         24, "grVn6Vn4");
  formatGraph(grVn8Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_), kGreen+2,  24, "grVn8Vn4");
  formatGraph(grVn8Vn6, "Centrality %", ratioMinVn8Vn6, ratioMaxVn8Vn6, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_), kViolet-1, 24, "grVn8Vn6");
  formatGraph(grG1E,    "Centrality %", g1eMin,         g1eMax,         "#gamma_{1}^{Exp}",                              2,         24, "grG1E");

  //Save to a file
  fOut = new TFile("SVDPhysics.root", "recreate");
  grVn2->Write();
  grVn4->Write();
  grVn6->Write();
  grVn8->Write();
  grVn6Vn4->Write();
  grVn8Vn4->Write();
  grVn8Vn6->Write();
  grG1E->Write();
  for(int icent = 0; icent < NCENT; icent++) finalUnf[icent]->Write();

  //-- Draw the resp and cov response matrices with the final unfolded dstn
  TPaletteAxis * palette[NCENT];
  TPaletteAxis * paletteCov[NCENT];

  TLegend * legUnf = new TLegend(0.5708, 0.6684, 0.9701, 0.8666);
  legUnf->SetBorderSize(0);
  legUnf->SetFillStyle(0);
  legUnf->AddEntry(hObs[0],       "Observed", "lp");
  legUnf->AddEntry(finalUnfm1[0], "Unfolded, kreg-1 ","lp");
  legUnf->AddEntry(finalUnf[0],   "Unfolded, kreg",   "lp");
  legUnf->AddEntry(finalUnfp1[0], "Unfolded, kreg+1", "lp");

  TLegend * legRef = new TLegend(0.5708, 0.6684, 0.9701, 0.8666);
  legRef->SetBorderSize(0);
  legRef->SetFillStyle(0);
  legRef->AddEntry(hObs[0],      "Observed", "lp");
  legRef->AddEntry(hrefoldm1[0], "Refolded, kreg-1", "l");
  legRef->AddEntry(hrefold[0],   "Refolded, kreg",   "l");
  legRef->AddEntry(hrefoldp1[0], "Refolded, kreg+1", "l");

  for(int icent = 0; icent < NCENT; icent++){
    cCovRespUnf[icent] = new TCanvas(Form("cCovRespUnf_c%i", icent), Form("cCovRespUnf_c%i", icent), 1500, 1000);
    cCovRespUnf[icent]->Divide(3,2);

    cCovRespUnf[icent]->cd(1);
    cCovRespUnf[icent]->cd(1)->SetLogz();
    hResponse[icent]->GetZaxis()->SetNdivisions(508);
    //hResponse[icent]->Scale(1./hResponse[icent]->Integral());
    hResponse[icent]->Draw("colz");
    cCovRespUnf[icent]->Update();
    palette[icent] = (TPaletteAxis*) hResponse[icent]->GetListOfFunctions()->FindObject("palette");
    palette[icent]->SetX1NDC(0.75);
    palette[icent]->SetY1NDC(0.2);
    palette[icent]->SetX2NDC(0.8);
    palette[icent]->SetY2NDC(0.5);
    latex.DrawLatex(0.6, 0.89, "Response Matrix");
    latex.DrawLatex(0.18, 0.89, "A");

    cCovRespUnf[icent]->cd(2);
    cCovRespUnf[icent]->cd(2)->SetLogy();
    hObs[icent]->Draw();
    finalUnfm1[icent]->Draw("same");
    finalUnf[icent]->Draw("same");
    finalUnfp1[icent]->Draw("same");
    legUnf->Draw("same");
    latex.DrawLatex(0.53, 0.89, "Unfolded Distribution");
    latex.DrawLatex(0.18, 0.89, "B");

    //-- Convert cov matrix to corr matrix
    hCorrMatrix[icent] = new TH2D( Form("hCorrMatrix_c%i", icent), Form("hCorrMatrix_c%i", icent), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_] );
    hCorrMatrix[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hCorrMatrix[icent]->GetYaxis()->SetTitle( Form("v_{%i}", norder_) );
    hCorrMatrix[icent]->SetOption("colz");

    TMatrixD MCov = H2M( hCovMatrix[icent] );
    TMatrixD S(NBins, NBins);
    for(int i = 0; i<NBins; i++){
      S(i,i) = sqrt( MCov(i,i) );
    }
    TDecompSVD svd(S);
    TMatrixD Sinv = svd.Invert();
    TMatrixD MCov2 = MCov;
    MCov2 *= Sinv;
    TMatrixD MCorr = Sinv;
    MCorr *= MCov2;
    M2H(MCorr, hCorrMatrix[icent]);

    cCovRespUnf[icent]->cd(3);
    hCorrMatrix[icent]->GetXaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->GetYaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->GetZaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->Draw("colz");
    cCovRespUnf[icent]->Update();
    paletteCov[icent] = (TPaletteAxis*) hCorrMatrix[icent]->GetListOfFunctions()->FindObject("palette");
    paletteCov[icent]->SetX1NDC(0.75);
    paletteCov[icent]->SetY1NDC(0.2);
    paletteCov[icent]->SetX2NDC(0.8);
    paletteCov[icent]->SetY2NDC(0.5);
    latex.DrawLatex(0.58, 0.89, "Correlation Matrix");
    latex.DrawLatex(0.18, 0.89, "C");

    cCovRespUnf[icent]->cd(4);
    cCovRespUnf[icent]->cd(4)->SetLogy();
    hdi[icent]->Draw();
    lOne[icent]->Draw("same");
    lKreg[icent]->Draw("same");
    latex.DrawLatex(0.65, 0.89,  "Regularization");
    latex.DrawLatex(0.18, 0.89, "D");

    cCovRespUnf[icent]->cd(5);
    cCovRespUnf[icent]->cd(5)->SetLogy();
    hObs[icent]->Draw();
    hrefoldm1[icent]->Draw("same");
    hrefold[icent]->Draw("same");
    hrefoldp1[icent]->Draw("same");
    legRef->Draw();
    latex.DrawLatex(0.75, 0.89,  "Refolding");
    latex.DrawLatex(0.18, 0.89, "E");

    cCovRespUnf[icent]->cd(6);
    latex2.DrawLatex(0.4, 0.6,  Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));

    cCovRespUnf[icent]->SaveAs( Form("../plots/unfolding/cSVDCovRespUnf_cent%i_%i.pdf", cent_min[icent], cent_max[icent]) );
  }


}

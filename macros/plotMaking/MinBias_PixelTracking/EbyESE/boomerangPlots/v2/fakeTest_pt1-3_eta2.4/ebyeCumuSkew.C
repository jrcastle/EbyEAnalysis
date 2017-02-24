#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void ebyeCumuSkew(){

  TH1D::SetDefaultSumw2();

  int centbin  = 4;
  int norder_  = 2;
  double tkEta = 2.4;
  double ptMin = 1.0;
  double ptMax = 3.0;

  bool dosys     = 0;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double momentMin = 0.00;
  double momentMax = 0.29;
  double chi2Min   = 0.1;
  double chi2Max   = 1000000.;

  TLatex latex;

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Unfolding output
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- Chi Squares for iteration cutoffs
  double chi2NDF_Refold[NCENT][NITER];
  int iterCutoff[NCENT];

  //-- Gamma1Exp
  double gamma1Exp[NCENT];
  double gamma1Expe[NCENT];
  TH1D * grGamma1Exp;

  //-- vn{6} / vn{4} ratio
  double ratio_vn6_vn4[NCENT];
  double ratio_vn6_vn4e[NCENT];
  TH1D * grvn6vn4Ratio;

  //-- vn{8} / vn{4} ratio
  double ratio_vn8_vn4[NCENT];
  double ratio_vn8_vn4e[NCENT];
  TH1D * grvn8vn4Ratio;

  //-- vn{8} / vn{6} ratio
  double ratio_vn8_vn6[NCENT];
  double ratio_vn8_vn6e[NCENT];
  TH1D * grvn8vn6Ratio;

  //-- Statistical uncertainties
  TFile * fStat;
  TH1D * hVarianceOfMean_Gamma1Exp;
  TH1D * hVarianceOfMean_Vn6Vn4;
  TH1D * hVarianceOfMean_Vn8Vn4;
  TH1D * hVarianceOfMean_Vn8Vn6;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get the Analyzer output file
  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );

  //-- Get the unfolding output file
  bool R1 = dataResp && studTResp && gaussResp;
  bool R2 = dataResp && studTResp;
  bool R3 = dataResp && gaussResp;
  bool R4 = studTResp && gaussResp;

  if( R1 || R2 || R3 || R4){
    std::cout<<"WARNING! More than one response function defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  fUnfold = 0;
  if( gaussResp ){
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/gaussResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/gaussResp/data%i_dosys.root", norder_) );
  }
  if( studTResp ){
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/studTResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/studTResp/data%i_dosys.root", norder_) );
  }
  if( dataResp ){
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i_dosys.root", norder_) );
  }

  //-- Get the statistical errors
  fStat = new TFile( Form("../../statErrorHandle/v%i/fakeTest_pt03-3_eta1.0/StatisticalUncertainties_v%i.root ", norder_, norder_) );
  hVarianceOfMean_Gamma1Exp = (TH1D*) fStat->Get("hVarianceOfMean_Gamma1Exp");
  hVarianceOfMean_Vn6Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn6Vn4");
  hVarianceOfMean_Vn8Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn4");
  hVarianceOfMean_Vn8Vn6    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn6");

  //-- Initialize cumu graphs
  grGamma1Exp = new TH1D("grGamma1Exp", "grGamma1Exp", NCENT, centbinsDefault);
  grGamma1Exp->SetLineColor(2);
  grGamma1Exp->SetMarkerColor(2);
  grGamma1Exp->SetMarkerStyle(20);
  grGamma1Exp->GetXaxis()->SetTitle( "Centrality %");
  grGamma1Exp->GetYaxis()->SetTitle( "#gamma_{1}^{exp}");
  grGamma1Exp->SetMinimum(-1.0);
  grGamma1Exp->SetMaximum(0.5);

  grvn6vn4Ratio = new TH1D("grvn6vn4Ratio","grvn6vn4Ratio", NCENT, centbinsDefault);
  grvn6vn4Ratio->SetLineColor(4);
  grvn6vn4Ratio->SetMarkerColor(4);
  grvn6vn4Ratio->SetMarkerStyle(21);
  grvn6vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
  grvn6vn4Ratio->SetMinimum(0.95);
  grvn6vn4Ratio->SetMaximum(1.03);

  grvn8vn4Ratio = new TH1D("grvn8vn4Ratio","grvn8vn4Ratio", NCENT, centbinsDefault);
  grvn8vn4Ratio->SetLineColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerStyle(34);
  grvn8vn4Ratio->SetMarkerSize(1.2);
  grvn8vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
  grvn8vn4Ratio->SetMinimum(0.95);
  grvn8vn4Ratio->SetMaximum(1.03);

  grvn8vn6Ratio = new TH1D("grvn8vn6Ratio","grvn8vn6Ratio", NCENT, centbinsDefault);
  grvn8vn6Ratio->SetLineColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerStyle(33);
  grvn8vn6Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_) );
  grvn8vn6Ratio->SetMinimum(0.99);
  grvn8vn6Ratio->SetMaximum(1.002);


  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfold[icent][i] = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold[icent][i]->SetLineColor(col[i]);
      hUnfold[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefold[icent][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold[icent][i]->SetLineWidth(2);
      hRefold[icent][i]->SetLineColor(col[i]);
      hRefold[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      chi2NDF_Refold[icent][i] = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF"); 

      if( chi2NDF_Refold[icent][i] < 1.2 ){
	iterCutoff[icent] = i;
	break;
      }
      if( i == NITER - 1 ) iterCutoff[icent] = i;

    } //-- End unfold iteration loop

  } //-- End cent loop

  //-- Calculate physics results here, that way warnings aren't buried in the Chi2Test warnings
  for(int icent = 0; icent < NCENT; icent++){

    std::cout << Form("============ Cent Bin %i ============",icent) << std::endl;

    int iter = iterCutoff[icent];
    FixUnfold( hUnfold[icent][iter] ) ;
    EbyECumu cumu(hUnfold[icent][iter]);

    //-- Gamma1Exp
    gamma1Exp[icent]  = cumu.GetGamma1Exp();
    gamma1Expe[icent] = sqrt( hVarianceOfMean_Gamma1Exp->GetBinContent(icent+1) );
    grGamma1Exp->SetBinContent(icent+1, gamma1Exp[icent]);
    grGamma1Exp->SetBinError(icent+1, gamma1Expe[icent]);

    //-- vn{*} / vn{4} ratio
    double vn8 = cumu.GetCumu_vn8();
    double vn6 = cumu.GetCumu_vn6();
    double vn4 = cumu.GetCumu_vn4();
    double vn2 = cumu.GetCumu_vn2();

    //-- Cumu ratio quality check vn{6} / vn{4}
    if( vn4 == 0|| vn6 == 0){
      ratio_vn6_vn4[icent]  = 0.; 
      ratio_vn6_vn4e[icent] = 1000;
    }
    else{
      ratio_vn6_vn4[icent]  = vn6 / vn4;
      ratio_vn6_vn4e[icent] = sqrt( hVarianceOfMean_Vn6Vn4->GetBinContent(icent+1) );
    }
 
    //-- Cumu ratio quality check vn{8} / vn{4}
    if( vn4 == 0|| vn8 == 0){
      ratio_vn8_vn4[icent]  = 0.;
      ratio_vn8_vn4e[icent] = 1000;
    }
    else{
      ratio_vn8_vn4[icent]  = vn8 / vn4;
      ratio_vn8_vn4e[icent] = sqrt( hVarianceOfMean_Vn8Vn4->GetBinContent(icent+1) );
    }

    //-- Cumu ratio quality check vn{8} / vn{6}
    if( vn6 == 0|| vn8 == 0){
      ratio_vn8_vn6[icent]  = 0.;
      ratio_vn8_vn6e[icent] = 1000;
    }
    else{
      ratio_vn8_vn6[icent]  = vn8 / vn6;
      ratio_vn8_vn6e[icent] = sqrt( hVarianceOfMean_Vn8Vn6->GetBinContent(icent+1) );
    }

    //-- Fill
    grvn6vn4Ratio->SetBinContent(icent+1, ratio_vn6_vn4[icent]);
    grvn6vn4Ratio->SetBinError(icent+1, ratio_vn6_vn4e[icent]);

    grvn8vn4Ratio->SetBinContent(icent+1, ratio_vn8_vn4[icent]);
    grvn8vn4Ratio->SetBinError(icent+1, ratio_vn8_vn4e[icent]);

    grvn8vn6Ratio->SetBinContent(icent+1, ratio_vn8_vn6[icent]);
    grvn8vn6Ratio->SetBinError(icent+1, ratio_vn8_vn6e[icent]);

  }

  //-- DRAW!
  TLine * line0 = new TLine(centbinsDefault[0], 0.0, centbinsDefault[NCENT], 0.0);
  line0->SetLineColor(1);
  line0->SetLineStyle(2);
  line0->SetLineWidth(2);

  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  grGamma1Exp->Draw();
  line0->Draw("same");
  latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.21, Form("%.2f < p_{T} < %.2f GeV/c", ptMin, ptMax));
  cGamma1Exp->SaveAs("cGamma1Exp.pdf");

  TLine * line = new TLine(centbinsDefault[0], 1.0, centbinsDefault[NCENT], 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * cvn6vn4Ratio = new TCanvas("cvn6vn4Ratio", "cvn6vn4Ratio", 500, 500);
  cvn6vn4Ratio->cd();
  grvn6vn4Ratio->Draw();
  line->Draw("same");
  latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.21, Form("%.2f < p_{T} < %.2f GeV/c", ptMin, ptMax));
  cvn6vn4Ratio->SaveAs("cvn6vn4Ratio.pdf");

  TCanvas * cvn8vn4Ratio = new TCanvas("cvn8vn4Ratio", "cvn8vn4Ratio", 500, 500);
  cvn8vn4Ratio->cd();
  grvn8vn4Ratio->Draw();
  line->Draw("same");
  latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.21, Form("%.2f < p_{T} < %.2f GeV/c", ptMin, ptMax));
  cvn8vn4Ratio->SaveAs("cvn8vn4Ratio.pdf");

  TCanvas * cvn8vn6Ratio = new TCanvas("cvn8vn6Ratio", "cvn8vn6Ratio", 500, 500);
  cvn8vn6Ratio->cd();
  grvn8vn6Ratio->Draw();
  line->Draw("same");
  latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.21, Form("%.2f < p_{T} < %.2f GeV/c", ptMin, ptMax));
  cvn8vn6Ratio->SaveAs("cvn8vn6Ratio.pdf");



}

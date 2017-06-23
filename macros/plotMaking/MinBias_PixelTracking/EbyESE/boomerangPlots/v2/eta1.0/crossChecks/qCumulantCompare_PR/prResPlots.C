#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/SysTablesEbyESE.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;
using namespace sysebyese;

void prResPlots(){

  TH1D::SetDefaultSumw2();

  int centbin     = 4;
  int norder_     = 2;
  double tkEta    = 1.0;

  bool dosys     = 0;
  bool compATLAS = 1;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double cumuMin = 0.0;
  double cumuMax = 0.25;

  double gamm1expMin = -1.0;
  double gamm1expMax = 0.5;
  double ratioMin    = 0.95;
  double ratioMax    = 1.03;
  double ratioMinVn8Vn6 = 0.99;
  double ratioMaxVn8Vn6 = 1.002;

  double titleOffset = 1.4;

  TLatex latex;

  //-- Save results to file
  TFile * fOut;

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Unfolding output
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- Chi Squares for iteration cutoffs
  int iterCutoff[NCENT];

  //--vn{2}
  double vn2Raw[NCENT];
  double vn2RawState[NCENT];
  TGraphErrors * grVn2Raw;

  //--vn{4}
  double vn4Raw[NCENT];
  double vn4RawState[NCENT];
  TGraphErrors * grVn4Raw;

  //--vn{6}
  double vn6Raw[NCENT];
  double vn6RawState[NCENT];
  TGraphErrors * grVn6Raw;

  //--vn{8}
  double vn8Raw[NCENT];
  double vn8RawState[NCENT];
  TGraphErrors * grVn8Raw;

  //-- Gamma1Exp
  double gamma1Exp[NCENT];
  double gamma1ExpState[NCENT];
  TGraphErrors * grGamma1Exp;

  //-- vn{6} / vn{4} ratio
  double ratio_vn6_vn4[NCENT];
  double ratio_vn6_vn4State[NCENT];
  TGraphErrors * grvn6vn4Ratio;

  TGraphErrors * grvn6vn4Ratio_Npart;
  TGraphAsymmErrors * grvn6vn4Ratio_ATLASNpart;

  //-- vn{8} / vn{4} ratio
  double ratio_vn8_vn4[NCENT];
  double ratio_vn8_vn4State[NCENT];
  TGraphErrors * grvn8vn4Ratio;

  TGraphErrors * grvn8vn4Ratio_Npart;
  TGraphAsymmErrors * grvn8vn4Ratio_ATLASNpart;

  //-- vn{8} / vn{6} ratio
  double ratio_vn8_vn6[NCENT];
  double ratio_vn8_vn6State[NCENT];
  TGraphErrors * grvn8vn6Ratio;

  //-- Statistical Errors
  TFile * fStat;
  TH1D * hVarianceOfMean_Vn2;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
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
  fAna = new TFile( "../../../promptRECOComp_pt1-3_eta1.0/AnalyzerResults/CastleEbyE.root" );

  //-- Get the unfolding output file
  fUnfold = new TFile( Form("../../../promptRECOComp_pt1-3_eta1.0/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Statistical Errors
  fStat = new TFile( Form("../../../../statErrorHandle/v%i/promptRECOComp_pt1-3_eta1.0/StatisticalUncertainties_v%i.root ", norder_,  norder_) );
  hVarianceOfMean_Vn2       = (TH1D*) fStat->Get("hVarianceOfMean_Vn2");
  hVarianceOfMean_Vn4       = (TH1D*) fStat->Get("hVarianceOfMean_Vn4");
  hVarianceOfMean_Vn6       = (TH1D*) fStat->Get("hVarianceOfMean_Vn6");
  hVarianceOfMean_Vn8       = (TH1D*) fStat->Get("hVarianceOfMean_Vn8");
  hVarianceOfMean_Gamma1Exp = (TH1D*) fStat->Get("hVarianceOfMean_Gamma1Exp");
  hVarianceOfMean_Vn6Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn6Vn4");
  hVarianceOfMean_Vn8Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn4");
  hVarianceOfMean_Vn8Vn6    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn6");

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
      double chi2NDF_Refold = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF"); 

      if( chi2NDF_Refold < 1.2 ){
	iterCutoff[icent] = i;
	break;
      }
      if( i == NITER - 1 ) iterCutoff[icent] = i;

    } //-- End unfold iteration loop

  } //-- End cent loop

  //-- Calculate physics results here, that way warnings aren't buried in the Chi2Test warnings
  for(int icent = 0; icent < NCENT; icent++){

    std::cout << Form("============ Cent Bin %i ============",icent) << std::endl;

    int ii = iterCutoff[icent];
    std::cout<<ii<<std::endl;
    FixUnfold( hUnfold[icent][ii] );
    EbyECumu cumu(hUnfold[icent][ii]);
    double vn2  = cumu.GetCumu_vn2();
    double vn4  = cumu.GetCumu_vn4();
    double vn6  = cumu.GetCumu_vn6();
    double vn8  = cumu.GetCumu_vn8();
    double gamma1exp = cumu.GetGamma1Exp();
    double vn6vn4;
    double vn8vn4;
    if( vn4 == 0 ){
      vn6vn4 = -1;
      vn8vn4 = -1;
    }
    else{
      vn6vn4 = vn6 / vn4;
      vn8vn4 = vn8 / vn4;
    }
    double vn8vn6;
    if( vn6 == 0 ) vn8vn6 = -1;
    else           vn8vn6 = vn8 / vn6;

    if( vn2 == 0 ) vn2 = -1;
    if( vn4 == 0 ) vn4 = -1;
    if( vn6 == 0 ) vn6 = -1;
    if( vn8 == 0 ) vn8 = -1;

    double vn2e       = sqrt( hVarianceOfMean_Vn2->GetBinContent(icent+1) );
    double vn4e       = sqrt( hVarianceOfMean_Vn4->GetBinContent(icent+1) );
    double vn6e       = sqrt( hVarianceOfMean_Vn6->GetBinContent(icent+1) );
    double vn8e       = sqrt( hVarianceOfMean_Vn8->GetBinContent(icent+1) );
    double gamma1expe = sqrt( hVarianceOfMean_Gamma1Exp->GetBinContent(icent+1) );
    double vn6vn4e    = sqrt( hVarianceOfMean_Vn6Vn4->GetBinContent(icent+1) );
    double vn8vn4e    = sqrt( hVarianceOfMean_Vn8Vn4->GetBinContent(icent+1) );
    double vn8vn6e    = sqrt( hVarianceOfMean_Vn8Vn6->GetBinContent(icent+1) );

    //-- vn{2}
    vn2Raw[icent]      = vn2;
    vn2RawState[icent] = vn2e;

    //-- vn{4}
    vn4Raw[icent]      = vn4;
    vn4RawState[icent] = vn4e;

    //-- vn{6}
    vn6Raw[icent]      = vn6;
    vn6RawState[icent] = vn6e;

    //-- vn{8}
    vn8Raw[icent]      = vn8;
    vn8RawState[icent] = vn8e;

    //-- Gamma1Exp
    gamma1Exp[icent]      = gamma1exp;
    gamma1ExpState[icent] = gamma1expe;

    //-- vn{6} / vn{4} ratio
    ratio_vn6_vn4[icent]      = vn6vn4;
    ratio_vn6_vn4State[icent] = vn6vn4e;

    //-- vn{8} / vn{4} ratio
    ratio_vn8_vn4[icent]      = vn8vn4;
    ratio_vn8_vn4State[icent] = vn8vn4e;

    //-- vn{8} / vn{6} ratio
    ratio_vn8_vn6[icent]      = vn8vn6;
    ratio_vn8_vn6State[icent] = vn8vn6e;

  }

  double nullCentErr[NCENT];
  for(int i = 0; i < NCENT; i++) nullCentErr[i] = 0.;

  //-- Initialize cumu graphs

  //-- vn{2}
  grVn2Raw = new TGraphErrors(NCENT, centBinCenter, vn2Raw, nullCentErr, vn2RawState);
  grVn2Raw->SetLineColor(9);
  grVn2Raw->SetMarkerColor(9);
  grVn2Raw->SetMarkerStyle(21);
  grVn2Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn2Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn2Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{4}
  grVn4Raw = new TGraphErrors(NCENT, centBinCenter, vn4Raw, nullCentErr, vn4RawState);
  grVn4Raw->SetLineColor(kSpring+4);
  grVn4Raw->SetMarkerColor(kSpring+4);
  grVn4Raw->SetMarkerStyle(20);
  grVn4Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn4Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn4Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{6}
  grVn6Raw = new TGraphErrors(NCENT, centBinCenter, vn6Raw, nullCentErr, vn6RawState);
  grVn6Raw->SetLineColor(6);
  grVn6Raw->SetMarkerColor(6);
  grVn6Raw->SetMarkerStyle(25);
  grVn6Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn6Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn6Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{8}
  grVn8Raw = new TGraphErrors(NCENT, centBinCenter, vn8Raw, nullCentErr, vn8RawState);
  grVn8Raw->SetLineColor(kOrange+7);
  grVn8Raw->SetMarkerColor(kOrange+7);
  grVn8Raw->SetMarkerStyle(28);
  grVn8Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn8Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn8Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- Gamma1Exp
  grGamma1Exp = new TGraphErrors(NCENT, centBinCenter, gamma1Exp, nullCentErr, gamma1ExpState) ;
  grGamma1Exp->SetLineColor(2);
  grGamma1Exp->SetMarkerColor(2);
  grGamma1Exp->SetMarkerStyle(20);
  grGamma1Exp->GetXaxis()->SetTitle( "Centrality %");
  grGamma1Exp->GetYaxis()->SetTitle( "#gamma_{1}^{exp}");
  grGamma1Exp->GetYaxis()->SetRangeUser(gamm1expMin, gamm1expMax);

  //-- vn{6} / vn{4}
  grvn6vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn6_vn4, nullCentErr, ratio_vn6_vn4State);
  grvn6vn4Ratio->SetLineColor(4);
  grvn6vn4Ratio->SetMarkerColor(4);
  grvn6vn4Ratio->SetMarkerStyle(21);
  grvn6vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
  grvn6vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn6vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);

  //-- vn{8} / vn{4}  
  grvn8vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn4, nullCentErr, ratio_vn8_vn4State);
  grvn8vn4Ratio->SetLineColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerStyle(34);
  grvn8vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
  grvn8vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);

  //-- vn{8} / vn{6}
  grvn8vn6Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn6, nullCentErr, ratio_vn8_vn6State);
  grvn8vn6Ratio->SetLineColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerStyle(33);
  grvn8vn6Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_) );
  grvn8vn6Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn6Ratio->GetYaxis()->SetTitleOffset(titleOffset+0.25);

  //-- DRAW!
  TLegend * legCumu = new TLegend(0.1915, 0.7182, 0.4355, 0.9322);
  legCumu->SetBorderSize(0);
  legCumu->SetFillStyle(0);
  legCumu->AddEntry(grVn2Raw, "k = 1", "lp");
  legCumu->AddEntry(grVn4Raw, "k = 2", "lp");
  legCumu->AddEntry(grVn6Raw, "k = 3", "lp");
  legCumu->AddEntry(grVn8Raw, "k = 4", "lp");

  TCanvas * cCumuRaw = new TCanvas("cCumuRaw", "cCumuRaw", 500, 500);
  cCumuRaw->cd();
  grVn2Raw->Draw("ap");
  grVn4Raw->Draw("psame");
  grVn6Raw->Draw("psame");
  grVn8Raw->Draw("psame");
  legCumu->Draw("same");
  latex.DrawLatex(0.53, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.53, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.53, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.53, 0.21, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

  TLegend * legg1e = new TLegend(0.3880, 0.7638, 0.9418, 0.9021);
  legInit( legg1e );
  legg1e->AddEntry(grGamma1Exp,       "CMS",        "lp");

  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  grGamma1Exp->Draw("ap");
  legg1e->Draw("same");
  latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.21, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));


  TLine * line2 = new TLine(centbinsDefault[0], 1.0, grvn8vn4Ratio->GetXaxis()->GetXmax(), 1.0);
  line2->SetLineColor(1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);

  TLegend * leg64 = new TLegend(0.3880, 0.7638, 0.9418, 0.9021);
  legInit( leg64 );
  leg64->AddEntry(grvn6vn4Ratio,       "CMS",        "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  
  cCumuRatio->cd(1);
  cCumuRatio->cd(1)->SetLeftMargin(0.2);
  grvn6vn4Ratio->Draw("ap");
  leg64->Draw("same");
  //latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  //latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  //latex.DrawLatex(0.2, 0.21, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

  cCumuRatio->cd(2);
  cCumuRatio->cd(2)->SetLeftMargin(0.2);
  grvn8vn4Ratio->Draw("ap");
  line2->Draw("same");
  latex.DrawLatex(0.25, 0.39, "CMS #it{Preliminary}");
  latex.DrawLatex(0.25, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.25, 0.27, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.25, 0.21, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

  cCumuRatio->cd(3);
  cCumuRatio->cd(3)->SetLeftMargin(0.2);
  grvn8vn6Ratio->Draw("ap");
  line2->Draw("same");
  //latex.DrawLatex(0.2, 0.39, "CMS #it{Preliminary}");
  //latex.DrawLatex(0.2, 0.33, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  //latex.DrawLatex(0.2, 0.27, Form("|#eta| < %.1f", tkEta));
  //latex.DrawLatex(0.2, 0.21, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  cCumuRatio->Update();

  //-- Save results to file
  fOut = new TFile("PhysicsResults.root", "recreate");
  fOut->cd();
  grVn2Raw->Write("grVn2Raw");
  grVn4Raw->Write("grVn4Raw");
  grVn6Raw->Write("grVn6Raw");
  grVn8Raw->Write("grVn8Raw");
  grvn6vn4Ratio->Write("grvn6vn4Ratio");
  grvn8vn4Ratio->Write("grvn8vn4Ratio");
  grvn8vn6Ratio->Write("grvn8vn6Ratio");
  grGamma1Exp->Write("grGamma1Exp");


}

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

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void sysFakes(){

  const int norder_ = 2;

  double vn6vn4Min = 0.95;
  double vn6vn4Max = 1.03;
  double vn8vn4Min = 0.95;
  double vn8vn4Max = 1.03;
  double vn8vn6Min = 0.99;
  double vn8vn6Max = 1.002;
  double g1eMin    = -1.0;
  double g1eMax    = 0.5;

  double rMin = 0.95;
  double rMax = 1.05;
  double g1rMin = 0.;
  double g1rMax = 5;

  TLatex latex;

  //-- 0.3 < pT < 3.0, |eta| < 2.4 (Default)
  TFile * fAnaPt03_3Eta2p4;
  TFile * fUnfPt03_3Eta2p4;
  TH1D * hObsPt03_3Eta2p4[NCENT];
  TH1D * hUnfoldPt03_3Eta2p4[NCENT][NITER];
  TH1D * hRefoldPt03_3Eta2p4[NCENT][NITER];
  bool stopFoundPt03_3Eta2p4;
  int iterCutPt03_3Eta2p4[NCENT];

  TFile * fStatPt03_3Eta2p4;
  TH1D * hVarianceOfMean_Vn6Vn4_Pt03_3Eta2p4;
  TH1D * hVarianceOfMean_Vn8Vn4_Pt03_3Eta2p4;
  TH1D * hVarianceOfMean_Vn8Vn6_Pt03_3Eta2p4;
  TH1D * hVarianceOfMean_Gamma1Exp_Pt03_3Eta2p4;

  double vn6vn4Pt03_3Eta2p4[NCENT];
  double vn8vn4Pt03_3Eta2p4[NCENT];
  double vn8vn6Pt03_3Eta2p4[NCENT];
  double g1ePt03_3Eta2p4[NCENT];

  double vn6vn4Pt03_3Eta2p4_err[NCENT];
  double vn8vn4Pt03_3Eta2p4_err[NCENT];
  double vn8vn6Pt03_3Eta2p4_err[NCENT];
  double g1ePt03_3Eta2p4_err[NCENT];

  TGraphErrors * grVn6Vn4Pt03_3Eta2p4;
  TGraphErrors * grVn8Vn4Pt03_3Eta2p4;
  TGraphErrors * grVn8Vn6Pt03_3Eta2p4;
  TGraphErrors * grG1EPt03_3Eta2p4;

  //-- 1.0 < pT < 3.0, |eta| < 2.4
  TFile * fAnaPt1_3Eta2p4;
  TFile * fUnfPt1_3Eta2p4;
  TH1D * hObsPt1_3Eta2p4[NCENT];
  TH1D * hUnfoldPt1_3Eta2p4[NCENT][NITER];
  TH1D * hRefoldPt1_3Eta2p4[NCENT][NITER];
  bool stopFoundPt1_3Eta2p4;
  int iterCutPt1_3Eta2p4[NCENT];

  TFile * fStatPt1_3Eta2p4;
  TH1D * hVarianceOfMean_Vn6Vn4_Pt1_3Eta2p4;
  TH1D * hVarianceOfMean_Vn8Vn4_Pt1_3Eta2p4;
  TH1D * hVarianceOfMean_Vn8Vn6_Pt1_3Eta2p4;
  TH1D * hVarianceOfMean_Gamma1Exp_Pt1_3Eta2p4;

  double vn6vn4Pt1_3Eta2p4[NCENT];
  double vn8vn4Pt1_3Eta2p4[NCENT];
  double vn8vn6Pt1_3Eta2p4[NCENT];
  double g1ePt1_3Eta2p4[NCENT];

  double vn6vn4Pt1_3Eta2p4_err[NCENT];
  double vn8vn4Pt1_3Eta2p4_err[NCENT];
  double vn8vn6Pt1_3Eta2p4_err[NCENT];
  double g1ePt1_3Eta2p4_err[NCENT];

  double vn6vn4Pt1_3Eta2p4_RatioToDefault[NCENT];
  double vn8vn4Pt1_3Eta2p4_RatioToDefault[NCENT];
  double vn8vn6Pt1_3Eta2p4_RatioToDefault[NCENT];
  double g1ePt1_3Eta2p4_RatioToDefault[NCENT];

  double vn6vn4Pt1_3Eta2p4_RatioToDefault_err[NCENT];
  double vn8vn4Pt1_3Eta2p4_RatioToDefault_err[NCENT];
  double vn8vn6Pt1_3Eta2p4_RatioToDefault_err[NCENT];
  double g1ePt1_3Eta2p4_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn4Pt1_3Eta2p4;
  TGraphErrors * grVn8Vn4Pt1_3Eta2p4;
  TGraphErrors * grVn8Vn6Pt1_3Eta2p4;
  TGraphErrors * grG1EPt1_3Eta2p4;

  TGraphErrors * grVn6Vn4Pt1_3Eta2p4_RatioToDefault;
  TGraphErrors * grVn8Vn4Pt1_3Eta2p4_RatioToDefault;
  TGraphErrors * grVn8Vn6Pt1_3Eta2p4_RatioToDefault;
  TGraphErrors * grG1EPt1_3Eta2p4_RatioToDefault;

  //-- 0.3 < pT < 3.0, |eta| < 1.0
  TFile * fAnaPt03_3Eta1p0;
  TFile * fUnfPt03_3Eta1p0;
  TH1D * hObsPt03_3Eta1p0[NCENT];
  TH1D * hUnfoldPt03_3Eta1p0[NCENT][NITER];
  TH1D * hRefoldPt03_3Eta1p0[NCENT][NITER];
  bool stopFoundPt03_3Eta1p0;
  int iterCutPt03_3Eta1p0[NCENT];

  TFile * fStatPt03_3Eta1p0;
  TH1D * hVarianceOfMean_Vn6Vn4_Pt03_3Eta1p0;
  TH1D * hVarianceOfMean_Vn8Vn4_Pt03_3Eta1p0;
  TH1D * hVarianceOfMean_Vn8Vn6_Pt03_3Eta1p0;
  TH1D * hVarianceOfMean_Gamma1Exp_Pt03_3Eta1p0;

  double vn6vn4Pt03_3Eta1p0[NCENT];
  double vn8vn4Pt03_3Eta1p0[NCENT];
  double vn8vn6Pt03_3Eta1p0[NCENT];
  double g1ePt03_3Eta1p0[NCENT];

  double vn6vn4Pt03_3Eta1p0_err[NCENT];
  double vn8vn4Pt03_3Eta1p0_err[NCENT];
  double vn8vn6Pt03_3Eta1p0_err[NCENT];
  double g1ePt03_3Eta1p0_err[NCENT];

  double vn6vn4Pt03_3Eta1p0_RatioToDefault[NCENT];
  double vn8vn4Pt03_3Eta1p0_RatioToDefault[NCENT];
  double vn8vn6Pt03_3Eta1p0_RatioToDefault[NCENT];
  double g1ePt03_3Eta1p0_RatioToDefault[NCENT];

  double vn6vn4Pt03_3Eta1p0_RatioToDefault_err[NCENT];
  double vn8vn4Pt03_3Eta1p0_RatioToDefault_err[NCENT];
  double vn8vn6Pt03_3Eta1p0_RatioToDefault_err[NCENT];
  double g1ePt03_3Eta1p0_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn4Pt03_3Eta1p0;
  TGraphErrors * grVn8Vn4Pt03_3Eta1p0;
  TGraphErrors * grVn8Vn6Pt03_3Eta1p0;
  TGraphErrors * grG1EPt03_3Eta1p0;

  TGraphErrors * grVn6Vn4Pt03_3Eta1p0_RatioToDefault;
  TGraphErrors * grVn8Vn4Pt03_3Eta1p0_RatioToDefault;
  TGraphErrors * grVn8Vn6Pt03_3Eta1p0_RatioToDefault;
  TGraphErrors * grG1EPt03_3Eta1p0_RatioToDefault;

  //-- Prompt RECO (1.0 < pT < 3.0, |eta| < 2.4)
  TFile * fAnaPromptRECO;
  TFile * fUnfPromptRECO;
  TH1D * hObsPromptRECO[NCENT];
  TH1D * hUnfoldPromptRECO[NCENT][NITER];
  TH1D * hRefoldPromptRECO[NCENT][NITER];
  bool stopFoundPromptRECO;
  int iterCutPromptRECO[NCENT];

  TFile * fStatPromptRECO;
  TH1D * hVarianceOfMean_Vn6Vn4_PromptRECO;
  TH1D * hVarianceOfMean_Vn8Vn4_PromptRECO;
  TH1D * hVarianceOfMean_Vn8Vn6_PromptRECO;
  TH1D * hVarianceOfMean_Gamma1Exp_PromptRECO;

  double vn6vn4PromptRECO[NCENT];
  double vn8vn4PromptRECO[NCENT];
  double vn8vn6PromptRECO[NCENT];
  double g1ePromptRECO[NCENT];

  double vn6vn4PromptRECO_err[NCENT];
  double vn8vn4PromptRECO_err[NCENT];
  double vn8vn6PromptRECO_err[NCENT];
  double g1ePromptRECO_err[NCENT];

  double vn6vn4PromptRECO_RatioToDefault[NCENT];
  double vn8vn4PromptRECO_RatioToDefault[NCENT];
  double vn8vn6PromptRECO_RatioToDefault[NCENT];
  double g1ePromptRECO_RatioToDefault[NCENT];

  double vn6vn4PromptRECO_RatioToDefault_err[NCENT];
  double vn8vn4PromptRECO_RatioToDefault_err[NCENT];
  double vn8vn6PromptRECO_RatioToDefault_err[NCENT];
  double g1ePromptRECO_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn4PromptRECO;
  TGraphErrors * grVn8Vn4PromptRECO;
  TGraphErrors * grVn8Vn6PromptRECO;
  TGraphErrors * grG1EPromptRECO;

  TGraphErrors * grVn6Vn4PromptRECO_RatioToDefault;
  TGraphErrors * grVn8Vn4PromptRECO_RatioToDefault;
  TGraphErrors * grVn8Vn6PromptRECO_RatioToDefault;
  TGraphErrors * grG1EPromptRECO_RatioToDefault;

  //-- qCumu 1.0 < pT < 3.0 |eta| < 2.4 
  TFile * fQuan_pt1_3_eta2p4;
  TGraphErrors * grQuanVn6Vn4pt1_3_eta2p4;
  TGraphErrors * grQuanVn8Vn4pt1_3_eta2p4;
  TGraphErrors * grQuanVn8Vn6pt1_3_eta2p4;

  //-- qCumu 0.3 < pT < 3.0 |eta| < 2.4
  TFile * fQuan_pt03_3_eta2p4;
  TGraphErrors * grQuanVn6Vn4pt03_3_eta2p4;
  TGraphErrors * grQuanVn8Vn4pt03_3_eta2p4;
  TGraphErrors * grQuanVn8Vn6pt03_3_eta2p4;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get Analyzer files
  fAnaPt03_3Eta2p4 = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fAnaPt1_3Eta2p4  = new TFile( "../../../fakeTest_pt1-3_eta2.4/AnalyzerResults/CastleEbyE.root" );
  fAnaPt03_3Eta1p0 = new TFile( "../../../fakeTest_pt03-3_eta1.0/AnalyzerResults/CastleEbyE.root" );
  fAnaPromptRECO   = new TFile( "../../../promptRECOComp_pt1-3_eta2.4/AnalyzerResults/CastleEbyE.root" );

  //-- Get Unfold files
  fUnfPt03_3Eta2p4 = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnfPt1_3Eta2p4  = new TFile( Form("../../../fakeTest_pt1-3_eta2.4/UnfoldResults/dataResp/data%i.root", norder_)  );
  fUnfPt03_3Eta1p0 = new TFile( Form("../../../fakeTest_pt03-3_eta1.0/UnfoldResults/dataResp/data%i.root", norder_)  );
  fUnfPromptRECO   = new TFile( Form("../../../promptRECOComp_pt1-3_eta2.4/UnfoldResults/dataResp/data%i.root", norder_)  );

  //-- Get Statistical Uncertainties
  fStatPt03_3Eta2p4 = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_Pt03_3Eta2p4    = (TH1D*) fStatPt03_3Eta2p4->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Pt03_3Eta2p4    = (TH1D*) fStatPt03_3Eta2p4->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Pt03_3Eta2p4    = (TH1D*) fStatPt03_3Eta2p4->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_Pt03_3Eta2p4 = (TH1D*) fStatPt03_3Eta2p4->Get( "hVarianceOfMean_Gamma1Exp" );

  fStatPt1_3Eta2p4  = new TFile( Form("../../../../statErrorHandle/v%i/fakeTest_pt1-3_eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_Pt1_3Eta2p4    = (TH1D*) fStatPt1_3Eta2p4->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Pt1_3Eta2p4    = (TH1D*) fStatPt1_3Eta2p4->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Pt1_3Eta2p4    = (TH1D*) fStatPt1_3Eta2p4->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_Pt1_3Eta2p4 = (TH1D*) fStatPt1_3Eta2p4->Get( "hVarianceOfMean_Gamma1Exp" );

  fStatPt03_3Eta1p0 = new TFile( Form("../../../../statErrorHandle/v%i/fakeTest_pt03-3_eta1.0/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_Pt03_3Eta1p0    = (TH1D*) fStatPt03_3Eta1p0->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Pt03_3Eta1p0    = (TH1D*) fStatPt03_3Eta1p0->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Pt03_3Eta1p0    = (TH1D*) fStatPt03_3Eta1p0->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_Pt03_3Eta1p0 = (TH1D*) fStatPt03_3Eta1p0->Get( "hVarianceOfMean_Gamma1Exp" );

  fStatPromptRECO = new TFile( Form("../../../../statErrorHandle/v%i/promptRECOComp_pt1-3_eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_PromptRECO    = (TH1D*) fStatPromptRECO->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_PromptRECO    = (TH1D*) fStatPromptRECO->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_PromptRECO    = (TH1D*) fStatPromptRECO->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_PromptRECO = (TH1D*) fStatPromptRECO->Get( "hVarianceOfMean_Gamma1Exp" );


  //-- Start Grabbing Histos
  for(int icent = 0; icent < NCENT; icent++){

    //-- Vn Observed 
    hObsPt03_3Eta2p4[icent] = (TH1D*) fAnaPt03_3Eta2p4->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsPt03_3Eta2p4[icent]->SetName( Form("hVnFullPt03_3Eta2p4_c%i", icent) );
    hObsPt1_3Eta2p4[icent]  = (TH1D*) fAnaPt1_3Eta2p4->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsPt1_3Eta2p4[icent]->SetName( Form("hVnFullPt1_3Eta2p4_c%i", icent) );
    hObsPt03_3Eta1p0[icent] = (TH1D*) fAnaPt03_3Eta1p0->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsPt03_3Eta1p0[icent]->SetName( Form("hVnFullPt03_3Eta1p0_c%i", icent) );
    hObsPromptRECO[icent] = (TH1D*) fAnaPromptRECO->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsPromptRECO[icent]->SetName( Form("hVnFullPromptRECO_c%i", icent) );

    //-- Reset stopFound boolean
    stopFoundPt03_3Eta2p4 = 0;
    stopFoundPt1_3Eta2p4  = 0;
    stopFoundPt03_3Eta1p0 = 0;
    stopFoundPromptRECO = 0;
    //-- Grab un/refolded distns
    for(int i = 0; i < NITER; i++){

      //-- 0.3 < pT < 3.0, |eta| < 2.4 (Default)
      hUnfoldPt03_3Eta2p4[icent][i] = (TH1D*) fUnfPt03_3Eta2p4->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldPt03_3Eta2p4[icent][i]->SetName( Form("hrecoPt03_3Eta2p4%i_c%i", iter[i], icent) );
      hRefoldPt03_3Eta2p4[icent][i] = (TH1D*) fUnfPt03_3Eta2p4->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldPt03_3Eta2p4[icent][i]->SetName( Form("hrefolPt03_3Eta2p4%i_c%i", iter[i], icent) );

      double chi2Pt03_3Eta2p4 = hRefoldPt03_3Eta2p4[icent][i]->Chi2Test(hObsPt03_3Eta2p4[icent], "CHI2/NDF");
      if( chi2Pt03_3Eta2p4 <= 1.2 && !stopFoundPt03_3Eta2p4 ){
	stopFoundPt03_3Eta2p4 = 1;
	iterCutPt03_3Eta2p4[icent] = i;
      }
      if( i == NITER-1 && !stopFoundPt03_3Eta2p4 ){
	stopFoundPt03_3Eta2p4 = 1;
        iterCutPt03_3Eta2p4[icent] = i;
      }

      //-- 1.0 < pT < 3.0, |eta| < 2.4
      hUnfoldPt1_3Eta2p4[icent][i] = (TH1D*) fUnfPt1_3Eta2p4->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldPt1_3Eta2p4[icent][i]->SetName( Form("hrecoPt1_3Eta2p4%i_c%i", iter[i], icent) );
      hRefoldPt1_3Eta2p4[icent][i] = (TH1D*) fUnfPt1_3Eta2p4->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldPt1_3Eta2p4[icent][i]->SetName( Form("hrefolPt1_3Eta2p4%i_c%i", iter[i], icent) );

      double chi2Pt1_3Eta2p4 = hRefoldPt1_3Eta2p4[icent][i]->Chi2Test(hObsPt1_3Eta2p4[icent], "CHI2/NDF");
      if( chi2Pt1_3Eta2p4 <= 1.2 && !stopFoundPt1_3Eta2p4 ){
        stopFoundPt1_3Eta2p4 = 1;
        iterCutPt1_3Eta2p4[icent] = i;
      } 
      if( i == NITER-1 && !stopFoundPt1_3Eta2p4 ){
        stopFoundPt1_3Eta2p4 = 1;
        iterCutPt1_3Eta2p4[icent] = i;
      }

      //-- 0.3 < pT < 3.0, |eta| < 1.0
      hUnfoldPt03_3Eta1p0[icent][i] = (TH1D*) fUnfPt03_3Eta1p0->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldPt03_3Eta1p0[icent][i]->SetName( Form("hrecoPt03_3Eta1p0%i_c%i", iter[i], icent) );
      hRefoldPt03_3Eta1p0[icent][i] = (TH1D*) fUnfPt03_3Eta1p0->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldPt03_3Eta1p0[icent][i]->SetName( Form("hrefolPt03_3Eta1p0%i_c%i", iter[i], icent) );

      double chi2Pt03_3Eta1p0 = hRefoldPt03_3Eta1p0[icent][i]->Chi2Test(hObsPt03_3Eta1p0[icent], "CHI2/NDF");
      if( chi2Pt03_3Eta1p0 <= 1.2 && !stopFoundPt03_3Eta1p0 ){
        stopFoundPt03_3Eta1p0 = 1;
        iterCutPt03_3Eta1p0[icent] = i;
      } 
      if( i == NITER-1 && !stopFoundPt03_3Eta1p0 ){
        stopFoundPt03_3Eta1p0 = 1;
        iterCutPt03_3Eta1p0[icent] = i;
      } 

      //-- Prompt RECO
      hUnfoldPromptRECO[icent][i] = (TH1D*) fUnfPromptRECO->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldPromptRECO[icent][i]->SetName( Form("hrecoPromptRECO%i_c%i", iter[i], icent) );
      hRefoldPromptRECO[icent][i] = (TH1D*) fUnfPromptRECO->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldPromptRECO[icent][i]->SetName( Form("hrefolPromptRECO%i_c%i", iter[i], icent) );

      double chi2PromptRECO = hRefoldPromptRECO[icent][i]->Chi2Test(hObsPromptRECO[icent], "CHI2/NDF");
      if( chi2PromptRECO <= 1.2 && !stopFoundPromptRECO ){
        stopFoundPromptRECO = 1;
        iterCutPromptRECO[icent] = i;
      }
      if( i == NITER-1 && !stopFoundPromptRECO ){
        stopFoundPromptRECO = 1;
        iterCutPromptRECO[icent] = i;
      }

    } //-- End iter loop

  } //-- End cent loop


  //-- Now that we're done fetching things, let's build some arrays!
  for(int icent = 0; icent < NCENT; icent++){

    //-- 0.3 < pT < 3.0, |eta| < 2.4 (Default)
    int iPt03_3Eta2p4 = iterCutPt03_3Eta2p4[icent];
    FixUnfold(hUnfoldPt03_3Eta2p4[icent][iPt03_3Eta2p4]);
    EbyECumu cumuPt03_3Eta2p4(hUnfoldPt03_3Eta2p4[icent][iPt03_3Eta2p4]);
    double vn4Pt03_3Eta2p4  = cumuPt03_3Eta2p4.GetCumu_vn4();
    double vn6Pt03_3Eta2p4  = cumuPt03_3Eta2p4.GetCumu_vn6();
    double vn8Pt03_3Eta2p4  = cumuPt03_3Eta2p4.GetCumu_vn8();
    double g1exPt03_3Eta2p4 = cumuPt03_3Eta2p4.GetGamma1Exp();

    if(vn4Pt03_3Eta2p4 == 0 || vn6Pt03_3Eta2p4 == 0) vn6vn4Pt03_3Eta2p4[icent]= -1000.;
    else                                             vn6vn4Pt03_3Eta2p4[icent] = vn6Pt03_3Eta2p4 / vn4Pt03_3Eta2p4;
    if(vn4Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn4Pt03_3Eta2p4[icent]= -1000.;
    else                                             vn8vn4Pt03_3Eta2p4[icent] = vn8Pt03_3Eta2p4 / vn4Pt03_3Eta2p4;
    if(vn6Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn6Pt03_3Eta2p4[icent]= -1000.;
    else                                             vn8vn6Pt03_3Eta2p4[icent] = vn8Pt03_3Eta2p4 / vn6Pt03_3Eta2p4;
    g1ePt03_3Eta2p4[icent] = g1exPt03_3Eta2p4;

    double vn6vn4Pt03_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn6Vn4_Pt03_3Eta2p4->GetBinContent(icent+1) );
    double vn8vn4Pt03_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn8Vn4_Pt03_3Eta2p4->GetBinContent(icent+1) );
    double vn8vn6Pt03_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn8Vn6_Pt03_3Eta2p4->GetBinContent(icent+1) );
    double g1ePt03_3Eta2p4_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_Pt03_3Eta2p4->GetBinContent(icent+1) );

    vn6vn4Pt03_3Eta2p4_err[icent] = vn6vn4Pt03_3Eta2p4_staterr;
    vn8vn4Pt03_3Eta2p4_err[icent] = vn8vn4Pt03_3Eta2p4_staterr;
    vn8vn6Pt03_3Eta2p4_err[icent] = vn8vn6Pt03_3Eta2p4_staterr;
    g1ePt03_3Eta2p4_err[icent]    = g1ePt03_3Eta2p4_staterr;

    //-- 1.0 < pT < 3.0, |eta| < 2.4
    //-- Raw
    int iPt1_3Eta2p4 = iterCutPt1_3Eta2p4[icent];
    FixUnfold(hUnfoldPt1_3Eta2p4[icent][iPt1_3Eta2p4]);
    EbyECumu cumuPt1_3Eta2p4(hUnfoldPt1_3Eta2p4[icent][iPt1_3Eta2p4]);
    double vn4Pt1_3Eta2p4  = cumuPt1_3Eta2p4.GetCumu_vn4();
    double vn6Pt1_3Eta2p4  = cumuPt1_3Eta2p4.GetCumu_vn6();
    double vn8Pt1_3Eta2p4  = cumuPt1_3Eta2p4.GetCumu_vn8();
    double g1exPt1_3Eta2p4 = cumuPt1_3Eta2p4.GetGamma1Exp();

    if(vn4Pt1_3Eta2p4 == 0 || vn6Pt1_3Eta2p4 == 0) vn6vn4Pt1_3Eta2p4[icent] = -1000.;
    else                                           vn6vn4Pt1_3Eta2p4[icent] = vn6Pt1_3Eta2p4 / vn4Pt1_3Eta2p4;
    if(vn4Pt1_3Eta2p4 == 0 || vn8Pt1_3Eta2p4 == 0) vn8vn4Pt1_3Eta2p4[icent] = -1000.;
    else                                           vn8vn4Pt1_3Eta2p4[icent] = vn8Pt1_3Eta2p4 / vn4Pt1_3Eta2p4;
    if(vn6Pt1_3Eta2p4 == 0 || vn8Pt1_3Eta2p4 == 0) vn8vn6Pt1_3Eta2p4[icent] = -1000.;
    else                                           vn8vn6Pt1_3Eta2p4[icent] = vn8Pt1_3Eta2p4 / vn6Pt1_3Eta2p4;
    g1ePt1_3Eta2p4[icent] = g1exPt1_3Eta2p4;

    double vn6vn4Pt1_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn6Vn4_Pt1_3Eta2p4->GetBinContent(icent+1) );
    double vn8vn4Pt1_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn8Vn4_Pt1_3Eta2p4->GetBinContent(icent+1) );
    double vn8vn6Pt1_3Eta2p4_staterr = sqrt( hVarianceOfMean_Vn8Vn6_Pt1_3Eta2p4->GetBinContent(icent+1) );
    double g1ePt1_3Eta2p4_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_Pt1_3Eta2p4->GetBinContent(icent+1) );

    vn6vn4Pt1_3Eta2p4_err[icent] = vn6vn4Pt1_3Eta2p4_staterr;
    vn8vn4Pt1_3Eta2p4_err[icent] = vn8vn4Pt1_3Eta2p4_staterr;
    vn8vn6Pt1_3Eta2p4_err[icent] = vn8vn6Pt1_3Eta2p4_staterr;
    g1ePt1_3Eta2p4_err[icent]    = g1ePt1_3Eta2p4_staterr;

    //-- Ratio to default
    if(vn4Pt1_3Eta2p4 == 0 || vn6Pt1_3Eta2p4 == 0 || vn4Pt03_3Eta2p4 == 0 || vn6Pt03_3Eta2p4 == 0) vn6vn4Pt1_3Eta2p4_RatioToDefault[icent] = -1000.;
    else                                                                                           vn6vn4Pt1_3Eta2p4_RatioToDefault[icent] = vn6vn4Pt1_3Eta2p4[icent] / vn6vn4Pt03_3Eta2p4[icent];
    if(vn4Pt1_3Eta2p4 == 0 || vn8Pt1_3Eta2p4 == 0 || vn4Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn4Pt1_3Eta2p4_RatioToDefault[icent] = -1000.;
    else                                                                                           vn8vn4Pt1_3Eta2p4_RatioToDefault[icent] = vn8vn4Pt1_3Eta2p4[icent] / vn8vn4Pt03_3Eta2p4[icent];
    if(vn6Pt1_3Eta2p4 == 0 || vn8Pt1_3Eta2p4 == 0 || vn6Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn6Pt1_3Eta2p4_RatioToDefault[icent] = -1000.;
    else                                                                                           vn8vn6Pt1_3Eta2p4_RatioToDefault[icent] = vn8vn6Pt1_3Eta2p4[icent] / vn8vn6Pt03_3Eta2p4[icent];
    if( g1exPt1_3Eta2p4 == -10000. || g1exPt03_3Eta2p4 == -10000.) g1ePt1_3Eta2p4_RatioToDefault[icent] = -1000.;
    else                                                           g1ePt1_3Eta2p4_RatioToDefault[icent] = g1ePt1_3Eta2p4[icent] / g1ePt03_3Eta2p4[icent];

    vn6vn4Pt1_3Eta2p4_RatioToDefault_err[icent] = sqrt( pow(vn6vn4Pt1_3Eta2p4_err[icent]/vn6vn4Pt03_3Eta2p4[icent],2) + pow(vn6vn4Pt1_3Eta2p4[icent]*vn6vn4Pt03_3Eta2p4_err[icent]/vn6vn4Pt03_3Eta2p4[icent]/vn6vn4Pt03_3Eta2p4[icent],2) );
    vn8vn4Pt1_3Eta2p4_RatioToDefault_err[icent] = sqrt( pow(vn8vn4Pt1_3Eta2p4_err[icent]/vn8vn4Pt03_3Eta2p4[icent],2) + pow(vn8vn4Pt1_3Eta2p4[icent]*vn8vn4Pt03_3Eta2p4_err[icent]/vn8vn4Pt03_3Eta2p4[icent]/vn8vn4Pt03_3Eta2p4[icent],2) );
    vn8vn6Pt1_3Eta2p4_RatioToDefault_err[icent] = sqrt( pow(vn8vn6Pt1_3Eta2p4_err[icent]/vn8vn6Pt03_3Eta2p4[icent],2) + pow(vn8vn6Pt1_3Eta2p4[icent]*vn8vn6Pt03_3Eta2p4_err[icent]/vn8vn6Pt03_3Eta2p4[icent]/vn8vn6Pt03_3Eta2p4[icent],2) );
    g1ePt1_3Eta2p4_RatioToDefault_err[icent]    = sqrt( pow(g1ePt1_3Eta2p4_err[icent]/g1ePt03_3Eta2p4[icent],2) + pow(g1ePt1_3Eta2p4[icent]*g1ePt03_3Eta2p4_err[icent]/g1ePt03_3Eta2p4[icent]/g1ePt03_3Eta2p4[icent],2) );

    //-- 0.3 < pT < 3.0, |eta| < 1.0
    //-- Raw
    int iPt03_3Eta1p0 = iterCutPt03_3Eta1p0[icent];
    FixUnfold(hUnfoldPt03_3Eta1p0[icent][iPt03_3Eta1p0]);
    EbyECumu cumuPt03_3Eta1p0(hUnfoldPt03_3Eta1p0[icent][iPt03_3Eta1p0]);
    double vn4Pt03_3Eta1p0  = cumuPt03_3Eta1p0.GetCumu_vn4();
    double vn6Pt03_3Eta1p0  = cumuPt03_3Eta1p0.GetCumu_vn6();
    double vn8Pt03_3Eta1p0  = cumuPt03_3Eta1p0.GetCumu_vn8();
    double g1exPt03_3Eta1p0 = cumuPt03_3Eta1p0.GetGamma1Exp();

    if(vn4Pt03_3Eta1p0 == 0 || vn6Pt03_3Eta1p0 == 0) vn6vn4Pt03_3Eta1p0[icent] = -1000.;
    else                                             vn6vn4Pt03_3Eta1p0[icent] = vn6Pt03_3Eta1p0 / vn4Pt03_3Eta1p0;
    if(vn4Pt03_3Eta1p0 == 0 || vn8Pt03_3Eta1p0 == 0) vn8vn4Pt03_3Eta1p0[icent] = -1000.;
    else                                             vn8vn4Pt03_3Eta1p0[icent] = vn8Pt03_3Eta1p0 / vn4Pt03_3Eta1p0;
    if(vn6Pt03_3Eta1p0 == 0 || vn8Pt03_3Eta1p0 == 0) vn8vn6Pt03_3Eta1p0[icent] = -1000.;
    else                                             vn8vn6Pt03_3Eta1p0[icent] = vn8Pt03_3Eta1p0 / vn6Pt03_3Eta1p0;
    g1ePt03_3Eta1p0[icent] = g1exPt03_3Eta1p0;

    double vn6vn4Pt03_3Eta1p0_staterr = sqrt( hVarianceOfMean_Vn6Vn4_Pt03_3Eta1p0->GetBinContent(icent+1) );
    double vn8vn4Pt03_3Eta1p0_staterr = sqrt( hVarianceOfMean_Vn8Vn4_Pt03_3Eta1p0->GetBinContent(icent+1) );
    double vn8vn6Pt03_3Eta1p0_staterr = sqrt( hVarianceOfMean_Vn8Vn6_Pt03_3Eta1p0->GetBinContent(icent+1) );
    double g1ePt03_3Eta1p0_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_Pt03_3Eta1p0->GetBinContent(icent+1) );

    vn6vn4Pt03_3Eta1p0_err[icent] = vn6vn4Pt03_3Eta1p0_staterr;
    vn8vn4Pt03_3Eta1p0_err[icent] = vn8vn4Pt03_3Eta1p0_staterr;
    vn8vn6Pt03_3Eta1p0_err[icent] = vn8vn6Pt03_3Eta1p0_staterr;
    g1ePt03_3Eta1p0_err[icent]    = g1ePt03_3Eta1p0_staterr;

    //-- Ratio to default
    if(vn4Pt03_3Eta1p0 == 0 || vn6Pt03_3Eta1p0 == 0 || vn4Pt03_3Eta2p4 == 0 || vn6Pt03_3Eta2p4 == 0) vn6vn4Pt03_3Eta1p0_RatioToDefault[icent] = -1000.;
    else                                                                                             vn6vn4Pt03_3Eta1p0_RatioToDefault[icent] = vn6vn4Pt03_3Eta1p0[icent] / vn6vn4Pt03_3Eta2p4[icent];
    if(vn4Pt03_3Eta1p0 == 0 || vn8Pt03_3Eta1p0 == 0 || vn4Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn4Pt03_3Eta1p0_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn4Pt03_3Eta1p0_RatioToDefault[icent] = vn8vn4Pt03_3Eta1p0[icent] / vn8vn4Pt03_3Eta2p4[icent];
    if(vn6Pt03_3Eta1p0 == 0 || vn8Pt03_3Eta1p0 == 0 || vn6Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn6Pt03_3Eta1p0_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn6Pt03_3Eta1p0_RatioToDefault[icent] = vn8vn6Pt03_3Eta1p0[icent] / vn8vn6Pt03_3Eta2p4[icent];
    if(g1exPt03_3Eta1p0== -10000. || g1exPt03_3Eta2p4 == -10000.) g1ePt03_3Eta1p0_RatioToDefault[icent] = -1000.;
    else                                                          g1ePt03_3Eta1p0_RatioToDefault[icent]= g1ePt03_3Eta1p0[icent]/ g1ePt03_3Eta2p4[icent];

    vn6vn4Pt03_3Eta1p0_RatioToDefault_err[icent] = sqrt( pow(vn6vn4Pt03_3Eta1p0_err[icent]/vn6vn4Pt03_3Eta2p4[icent],2) + pow(vn6vn4Pt03_3Eta1p0[icent]*vn6vn4Pt03_3Eta2p4_err[icent]/vn6vn4Pt03_3Eta2p4[icent]/vn6vn4Pt03_3Eta2p4[icent],2) );
    vn8vn4Pt03_3Eta1p0_RatioToDefault_err[icent] = sqrt( pow(vn8vn4Pt03_3Eta1p0_err[icent]/vn8vn4Pt03_3Eta2p4[icent],2) + pow(vn8vn4Pt03_3Eta1p0[icent]*vn8vn4Pt03_3Eta2p4_err[icent]/vn8vn4Pt03_3Eta2p4[icent]/vn8vn4Pt03_3Eta2p4[icent],2) );
    vn8vn6Pt03_3Eta1p0_RatioToDefault_err[icent] = sqrt( pow(vn8vn6Pt03_3Eta1p0_err[icent]/vn8vn6Pt03_3Eta2p4[icent],2) + pow(vn8vn6Pt03_3Eta1p0[icent]*vn8vn6Pt03_3Eta2p4_err[icent]/vn8vn6Pt03_3Eta2p4[icent]/vn8vn6Pt03_3Eta2p4[icent],2) );
    g1ePt03_3Eta1p0_RatioToDefault_err[icent]    = sqrt( pow(g1ePt03_3Eta1p0_err[icent]/g1ePt03_3Eta2p4[icent],2) + pow(g1ePt03_3Eta1p0[icent]*g1ePt03_3Eta2p4_err[icent]/g1ePt03_3Eta2p4[icent]/g1ePt03_3Eta2p4[icent],2) );

    //-- Prompt RECO
    //-- Raw
    int iPromptRECO = iterCutPromptRECO[icent];
    FixUnfold(hUnfoldPromptRECO[icent][iPromptRECO]);
    EbyECumu cumuPromptRECO(hUnfoldPromptRECO[icent][iPromptRECO]);
    double vn4PromptRECO  = cumuPromptRECO.GetCumu_vn4();
    double vn6PromptRECO  = cumuPromptRECO.GetCumu_vn6();
    double vn8PromptRECO  = cumuPromptRECO.GetCumu_vn8();
    double g1exPromptRECO = cumuPromptRECO.GetGamma1Exp();

    if(vn4PromptRECO == 0 || vn6PromptRECO == 0) vn6vn4PromptRECO[icent] = -1000.;
    else                                             vn6vn4PromptRECO[icent] = vn6PromptRECO / vn4PromptRECO;
    if(vn4PromptRECO == 0 || vn8PromptRECO == 0) vn8vn4PromptRECO[icent] = -1000.;
    else                                             vn8vn4PromptRECO[icent] = vn8PromptRECO / vn4PromptRECO;
    if(vn6PromptRECO == 0 || vn8PromptRECO == 0) vn8vn6PromptRECO[icent] = -1000.;
    else                                             vn8vn6PromptRECO[icent] = vn8PromptRECO / vn6PromptRECO;
    g1ePromptRECO[icent] = g1exPromptRECO;

    double vn6vn4PromptRECO_staterr = sqrt( hVarianceOfMean_Vn6Vn4_PromptRECO->GetBinContent(icent+1) );
    double vn8vn4PromptRECO_staterr = sqrt( hVarianceOfMean_Vn8Vn4_PromptRECO->GetBinContent(icent+1) );
    double vn8vn6PromptRECO_staterr = sqrt( hVarianceOfMean_Vn8Vn6_PromptRECO->GetBinContent(icent+1) );
    double g1ePromptRECO_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_PromptRECO->GetBinContent(icent+1) );

    vn6vn4PromptRECO_err[icent] = vn6vn4PromptRECO_staterr;
    vn8vn4PromptRECO_err[icent] = vn8vn4PromptRECO_staterr;
    vn8vn6PromptRECO_err[icent] = vn8vn6PromptRECO_staterr;
    g1ePromptRECO_err[icent]    = g1ePromptRECO_staterr;

    //-- Ratio to default
    if(vn4PromptRECO == 0 || vn6PromptRECO == 0 || vn4Pt03_3Eta2p4 == 0 || vn6Pt03_3Eta2p4 == 0) vn6vn4PromptRECO_RatioToDefault[icent] = -1000.;
    else                                                                                             vn6vn4PromptRECO_RatioToDefault[icent] = vn6vn4PromptRECO[icent] / vn6vn4Pt03_3Eta2p4[icent];
    if(vn4PromptRECO == 0 || vn8PromptRECO == 0 || vn4Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn4PromptRECO_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn4PromptRECO_RatioToDefault[icent] = vn8vn4PromptRECO[icent] / vn8vn4Pt03_3Eta2p4[icent];
    if(vn6PromptRECO == 0 || vn8PromptRECO == 0 || vn6Pt03_3Eta2p4 == 0 || vn8Pt03_3Eta2p4 == 0) vn8vn6PromptRECO_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn6PromptRECO_RatioToDefault[icent] = vn8vn6PromptRECO[icent] / vn8vn6Pt03_3Eta2p4[icent];
    if(g1exPromptRECO== -10000. || g1exPt03_3Eta2p4 == -10000.) g1ePromptRECO_RatioToDefault[icent] = -1000.;
    else                                                          g1ePromptRECO_RatioToDefault[icent]= g1ePromptRECO[icent]/ g1ePt03_3Eta2p4[icent];

    vn6vn4PromptRECO_RatioToDefault_err[icent] = sqrt( pow(vn6vn4PromptRECO_err[icent]/vn6vn4Pt03_3Eta2p4[icent],2) + pow(vn6vn4PromptRECO[icent]*vn6vn4Pt03_3Eta2p4_err[icent]/vn6vn4Pt03_3Eta2p4[icent]/vn6vn4Pt03_3Eta2p4[icent],2) );
    vn8vn4PromptRECO_RatioToDefault_err[icent] = sqrt( pow(vn8vn4PromptRECO_err[icent]/vn8vn4Pt03_3Eta2p4[icent],2) + pow(vn8vn4PromptRECO[icent]*vn8vn4Pt03_3Eta2p4_err[icent]/vn8vn4Pt03_3Eta2p4[icent]/vn8vn4Pt03_3Eta2p4[icent],2) );
    vn8vn6PromptRECO_RatioToDefault_err[icent] = sqrt( pow(vn8vn6PromptRECO_err[icent]/vn8vn6Pt03_3Eta2p4[icent],2) + pow(vn8vn6PromptRECO[icent]*vn8vn6Pt03_3Eta2p4_err[icent]/vn8vn6Pt03_3Eta2p4[icent]/vn8vn6Pt03_3Eta2p4[icent],2) );
    g1ePromptRECO_RatioToDefault_err[icent]    = sqrt( pow(g1ePromptRECO_err[icent]/g1ePt03_3Eta2p4[icent],2) + pow(g1ePromptRECO[icent]*g1ePt03_3Eta2p4_err[icent]/g1ePt03_3Eta2p4[icent]/g1ePt03_3Eta2p4[icent],2) );

  } //-- End cent loop

  //-- TGraph time!
  double c_err[NCENT];
  for(int i = 0; i < NCENT; i++) c_err[i] = 0.;

  //-- 0.3 < pT < 3.0, |eta| < 2.4 (Default)
  grVn6Vn4Pt03_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn6vn4Pt03_3Eta2p4, c_err, vn6vn4Pt03_3Eta2p4_err);
  grVn8Vn4Pt03_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn8vn4Pt03_3Eta2p4, c_err, vn8vn4Pt03_3Eta2p4_err);
  grVn8Vn6Pt03_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn8vn6Pt03_3Eta2p4, c_err, vn8vn6Pt03_3Eta2p4_err);
  grG1EPt03_3Eta2p4    = new TGraphErrors(NCENT, centBinCenter, g1ePt03_3Eta2p4,    c_err, g1ePt03_3Eta2p4_err);

  formatGraph(grVn6Vn4Pt03_3Eta2p4, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         21, "grVn6Vn4Pt03_3Eta2p4");
  formatGraph(grVn8Vn4Pt03_3Eta2p4, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  34, "grVn8Vn4Pt03_3Eta2p4");
  formatGraph(grVn8Vn6Pt03_3Eta2p4, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 33, "grVn8Vn6Pt03_3Eta2p4");
  formatGraph(grG1EPt03_3Eta2p4,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         20, "grG1EPt03_3Eta2p4");

  //-- 1.0 < pT < 3.0, |eta| < 2.4
  grVn6Vn4Pt1_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn6vn4Pt1_3Eta2p4, c_err, vn6vn4Pt1_3Eta2p4_err);
  grVn8Vn4Pt1_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn8vn4Pt1_3Eta2p4, c_err, vn8vn4Pt1_3Eta2p4_err);
  grVn8Vn6Pt1_3Eta2p4 = new TGraphErrors(NCENT, centBinCenter, vn8vn6Pt1_3Eta2p4, c_err, vn8vn6Pt1_3Eta2p4_err);
  grG1EPt1_3Eta2p4    = new TGraphErrors(NCENT, centBinCenter, g1ePt1_3Eta2p4,    c_err, g1ePt1_3Eta2p4_err);

  formatGraph(grVn6Vn4Pt1_3Eta2p4, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         24, "grVn6Vn4Pt1_3Eta2p4");
  formatGraph(grVn8Vn4Pt1_3Eta2p4, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  25, "grVn8Vn4Pt1_3Eta2p4");
  formatGraph(grVn8Vn6Pt1_3Eta2p4, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 28, "grVn8Vn6Pt1_3Eta2p4");
  formatGraph(grG1EPt1_3Eta2p4,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         26, "grG1EPt1_3Eta2p4");

  //-- Ratio to Default
  grVn6Vn4Pt1_3Eta2p4_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn4Pt1_3Eta2p4_RatioToDefault, c_err, vn6vn4Pt1_3Eta2p4_RatioToDefault_err);
  grVn8Vn4Pt1_3Eta2p4_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn4Pt1_3Eta2p4_RatioToDefault, c_err, vn8vn4Pt1_3Eta2p4_RatioToDefault_err);
  grVn8Vn6Pt1_3Eta2p4_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn6Pt1_3Eta2p4_RatioToDefault, c_err, vn8vn6Pt1_3Eta2p4_RatioToDefault_err);
  grG1EPt1_3Eta2p4_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1ePt1_3Eta2p4_RatioToDefault,    c_err, g1ePt1_3Eta2p4_RatioToDefault_err);

  formatGraph(grVn6Vn4Pt1_3Eta2p4_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         24, "grVn6Vn4Pt1_3Eta2p4_RatioToDefault");
  formatGraph(grVn8Vn4Pt1_3Eta2p4_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  25, "grVn8Vn4Pt1_3Eta2p4_RatioToDefault");
  formatGraph(grVn8Vn6Pt1_3Eta2p4_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 28, "grVn8Vn6Pt1_3Eta2p4_RatioToDefault");
  formatGraph(grG1EPt1_3Eta2p4_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         26, "grG1EPt1_3Eta2p4_RatioToDefault");

  //-- 0.3 < pT < 3.0, |eta| < 1.0
  grVn6Vn4Pt03_3Eta1p0 = new TGraphErrors(NCENT, centBinCenter, vn6vn4Pt03_3Eta1p0, c_err, vn6vn4Pt03_3Eta1p0_err);
  grVn8Vn4Pt03_3Eta1p0 = new TGraphErrors(NCENT, centBinCenter, vn8vn4Pt03_3Eta1p0, c_err, vn8vn4Pt03_3Eta1p0_err);
  grVn8Vn6Pt03_3Eta1p0 = new TGraphErrors(NCENT, centBinCenter, vn8vn6Pt03_3Eta1p0, c_err, vn8vn6Pt03_3Eta1p0_err);
  grG1EPt03_3Eta1p0    = new TGraphErrors(NCENT, centBinCenter, g1ePt03_3Eta1p0,    c_err, g1ePt03_3Eta1p0_err);

  formatGraph(grVn6Vn4Pt03_3Eta1p0, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         32, "grVn6Vn4Pt03_3Eta1p0");
  formatGraph(grVn8Vn4Pt03_3Eta1p0, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  27, "grVn8Vn4Pt03_3Eta1p0");
  formatGraph(grVn8Vn6Pt03_3Eta1p0, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 30, "grVn8Vn6Pt03_3Eta1p0");
  formatGraph(grG1EPt03_3Eta1p0,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         31, "grG1EPt03_3Eta1p0");

  //-- Ratio to Default
  grVn6Vn4Pt03_3Eta1p0_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn4Pt03_3Eta1p0_RatioToDefault, c_err, vn6vn4Pt03_3Eta1p0_RatioToDefault_err);
  grVn8Vn4Pt03_3Eta1p0_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn4Pt03_3Eta1p0_RatioToDefault, c_err, vn8vn4Pt03_3Eta1p0_RatioToDefault_err);
  grVn8Vn6Pt03_3Eta1p0_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn6Pt03_3Eta1p0_RatioToDefault, c_err, vn8vn6Pt03_3Eta1p0_RatioToDefault_err);
  grG1EPt03_3Eta1p0_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1ePt03_3Eta1p0_RatioToDefault,    c_err, g1ePt03_3Eta1p0_RatioToDefault_err);

  formatGraph(grVn6Vn4Pt03_3Eta1p0_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         32, "grVn6Vn4Pt03_3Eta1p0_RatioToDefault");
  formatGraph(grVn8Vn4Pt03_3Eta1p0_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  27, "grVn8Vn4Pt03_3Eta1p0_RatioToDefault");
  formatGraph(grVn8Vn6Pt03_3Eta1p0_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 30, "grVn8Vn6Pt03_3Eta1p0_RatioToDefault");
  formatGraph(grG1EPt03_3Eta1p0_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         31, "grG1EPt03_3Eta1p0_RatioToDefault");

  //-- Prompt RECO
  grVn6Vn4PromptRECO = new TGraphErrors(NCENT, centBinCenter, vn6vn4PromptRECO, c_err, vn6vn4PromptRECO_err);
  grVn8Vn4PromptRECO = new TGraphErrors(NCENT, centBinCenter, vn8vn4PromptRECO, c_err, vn8vn4PromptRECO_err);
  grVn8Vn6PromptRECO = new TGraphErrors(NCENT, centBinCenter, vn8vn6PromptRECO, c_err, vn8vn6PromptRECO_err);
  grG1EPromptRECO    = new TGraphErrors(NCENT, centBinCenter, g1ePromptRECO,    c_err, g1ePromptRECO_err);

  formatGraph(grVn6Vn4PromptRECO, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 1, 24, "grVn6Vn4PromptRECO");
  formatGraph(grVn8Vn4PromptRECO, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", 1, 25, "grVn8Vn4PromptRECO");
  formatGraph(grVn8Vn6PromptRECO, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", 1, 28, "grVn8Vn6PromptRECO");
  formatGraph(grG1EPromptRECO,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  1, 26, "grG1EPromptRECO");

  //-- Ratio to Default
  grVn6Vn4PromptRECO_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn4PromptRECO_RatioToDefault, c_err, vn6vn4PromptRECO_RatioToDefault_err);
  grVn8Vn4PromptRECO_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn4PromptRECO_RatioToDefault, c_err, vn8vn4PromptRECO_RatioToDefault_err);
  grVn8Vn6PromptRECO_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn6PromptRECO_RatioToDefault, c_err, vn8vn6PromptRECO_RatioToDefault_err);
  grG1EPromptRECO_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1ePromptRECO_RatioToDefault,    c_err, g1ePromptRECO_RatioToDefault_err);

  formatGraph(grVn6Vn4PromptRECO_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 1, 24, "grVn6Vn4PromptRECO_RatioToDefault");
  formatGraph(grVn8Vn4PromptRECO_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", 1, 25, "grVn8Vn4PromptRECO_RatioToDefault");
  formatGraph(grVn8Vn6PromptRECO_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", 1, 28, "grVn8Vn6PromptRECO_RatioToDefault");
  formatGraph(grG1EPromptRECO_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  1, 26, "grG1EPromptRECO_RatioToDefault");

  //-- qCumu 1.0 < pT < 3.0 |eta| < 2.4
  fQuan_pt1_3_eta2p4 = new TFile("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/boomerangPlots/v2/eta1.0/crossChecks/qCumulantCompare/fQuan_pt1_3_eta2p4.root");
  grQuanVn6Vn4pt1_3_eta2p4 = (TGraphErrors*) fQuan_pt1_3_eta2p4->Get("grQuanVn6Vn4");
  grQuanVn8Vn4pt1_3_eta2p4 = (TGraphErrors*) fQuan_pt1_3_eta2p4->Get("grQuanVn8Vn4");
  grQuanVn8Vn6pt1_3_eta2p4 = (TGraphErrors*) fQuan_pt1_3_eta2p4->Get("grQuanVn8Vn6");

  grQuanVn6Vn4pt1_3_eta2p4->SetLineColor(2);
  grQuanVn6Vn4pt1_3_eta2p4->SetMarkerColor(2);
  grQuanVn6Vn4pt1_3_eta2p4->SetMarkerStyle(24);

  grQuanVn8Vn4pt1_3_eta2p4->SetLineColor(2);
  grQuanVn8Vn4pt1_3_eta2p4->SetMarkerColor(2);
  grQuanVn8Vn4pt1_3_eta2p4->SetMarkerStyle(25);

  grQuanVn8Vn6pt1_3_eta2p4->SetLineColor(2);
  grQuanVn8Vn6pt1_3_eta2p4->SetMarkerColor(2);
  grQuanVn8Vn6pt1_3_eta2p4->SetMarkerStyle(28);

  //-- qCumu 0.3 < pT < 3.0 |eta| < 2.4
  fQuan_pt03_3_eta2p4 = new TFile("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/boomerangPlots/v2/eta2.4/crossChecks/qCumulantCompare/fQuan_pt03_3_eta2p4.root");
  grQuanVn6Vn4pt03_3_eta2p4 = (TGraphErrors*) fQuan_pt03_3_eta2p4->Get("grQuanVn6Vn4");
  grQuanVn8Vn4pt03_3_eta2p4 = (TGraphErrors*) fQuan_pt03_3_eta2p4->Get("grQuanVn8Vn4");
  grQuanVn8Vn6pt03_3_eta2p4 = (TGraphErrors*) fQuan_pt03_3_eta2p4->Get("grQuanVn8Vn6");

  grQuanVn6Vn4pt03_3_eta2p4->SetLineColor(2);
  grQuanVn6Vn4pt03_3_eta2p4->SetMarkerColor(2);
  grQuanVn6Vn4pt03_3_eta2p4->SetMarkerStyle(32);

  grQuanVn8Vn4pt03_3_eta2p4->SetLineColor(2);
  grQuanVn8Vn4pt03_3_eta2p4->SetMarkerColor(2);
  grQuanVn8Vn4pt03_3_eta2p4->SetMarkerStyle(27);

  grQuanVn8Vn6pt03_3_eta2p4->SetLineColor(2);
  grQuanVn8Vn6pt03_3_eta2p4->SetMarkerColor(2);
  grQuanVn8Vn6pt03_3_eta2p4->SetMarkerStyle(30);


  //-- DRAW

  //-- cumuRatio
  TLegend * legvn6vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Pt03_3Eta2p4,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn6vn4->AddEntry(grVn6Vn4Pt03_3Eta1p0,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 1.0",              "lp");
  legvn6vn4->AddEntry(grQuanVn6Vn4pt03_3_eta2p4, "qCumu: 0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");
  legvn6vn4->AddEntry(grVn6Vn4Pt1_3Eta2p4,       "1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn6vn4->AddEntry(grVn6Vn4PromptRECO,        "Prompt RECO: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4", "lp");
  legvn6vn4->AddEntry(grQuanVn6Vn4pt1_3_eta2p4,  "qCumu: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");

  TLegend * legvn8vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Pt03_3Eta2p4,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn8vn4->AddEntry(grVn8Vn4Pt03_3Eta1p0,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 1.0",              "lp");
  legvn8vn4->AddEntry(grQuanVn8Vn4pt03_3_eta2p4, "qCumu: 0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");
  legvn8vn4->AddEntry(grVn8Vn4Pt1_3Eta2p4,       "1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn8vn4->AddEntry(grVn8Vn4PromptRECO,        "Prompt RECO: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4", "lp");
  legvn8vn4->AddEntry(grQuanVn8Vn4pt1_3_eta2p4,  "qCumu: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");

  TLegend * legvn8vn6 = new TLegend(0.4520, 0.1973, 0.94, 0.3703);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Pt03_3Eta2p4,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn8vn6->AddEntry(grVn8Vn6Pt03_3Eta1p0,      "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 1.0",              "lp");
  legvn8vn6->AddEntry(grQuanVn8Vn6pt03_3_eta2p4, "qCumu: 0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");
  legvn8vn6->AddEntry(grVn8Vn6Pt1_3Eta2p4,       "1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legvn8vn6->AddEntry(grVn8Vn6PromptRECO,        "Prompt RECO: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4", "lp");
  legvn8vn6->AddEntry(grQuanVn8Vn6pt1_3_eta2p4,  "qCumu: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",       "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  cCumuRatio->cd(1);
  grVn6Vn4Pt03_3Eta2p4->Draw("ap");
  grVn6Vn4Pt1_3Eta2p4->Draw("psame");
  grVn6Vn4Pt03_3Eta1p0->Draw("psame");
  grVn6Vn4PromptRECO->Draw("psame");
  grQuanVn6Vn4pt03_3_eta2p4->Draw("psame");
  grQuanVn6Vn4pt1_3_eta2p4->Draw("psame");
  legvn6vn4->Draw("same");

  cCumuRatio->cd(2);
  grVn8Vn4Pt03_3Eta2p4->Draw("ap");
  grVn8Vn4Pt1_3Eta2p4->Draw("psame");
  grVn8Vn4Pt03_3Eta1p0->Draw("psame");
  grVn8Vn4PromptRECO->Draw("psame");
  grQuanVn8Vn4pt03_3_eta2p4->Draw("psame");
  grQuanVn8Vn4pt1_3_eta2p4->Draw("psame");
  legvn8vn4->Draw("same");
  cCumuRatio->cd(3);

  grVn8Vn6Pt03_3Eta2p4->Draw("ap");
  grVn8Vn6Pt1_3Eta2p4->Draw("psame");
  grVn8Vn6Pt03_3Eta1p0->Draw("psame");
  grVn8Vn6PromptRECO->Draw("psame");
  grQuanVn8Vn6pt03_3_eta2p4->Draw("psame");
  grQuanVn8Vn6pt1_3_eta2p4->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- cumuRatio To Default
  TCanvas * cCumuRatio_RatToDef = new TCanvas("cCumuRatio_RatToDef", "cCumuRatio_RatToDef", 1500, 500);
  cCumuRatio_RatToDef->Divide(3,1);
  cCumuRatio_RatToDef->cd(1);
  grVn6Vn4Pt1_3Eta2p4_RatioToDefault->Draw("ap");
  grVn6Vn4Pt03_3Eta1p0_RatioToDefault->Draw("psame");
  grVn6Vn4PromptRECO_RatioToDefault->Draw("psame");
  cCumuRatio_RatToDef->cd(2);
  grVn8Vn4Pt1_3Eta2p4_RatioToDefault->Draw("ap");
  grVn8Vn4Pt03_3Eta1p0_RatioToDefault->Draw("psame");
  grVn6Vn4PromptRECO_RatioToDefault->Draw("psame");
  cCumuRatio_RatToDef->cd(3);
  grVn8Vn6Pt1_3Eta2p4_RatioToDefault->Draw("ap");
  grVn8Vn6Pt03_3Eta1p0_RatioToDefault->Draw("psame");
  grVn6Vn4PromptRECO_RatioToDefault->Draw("psame");
  cCumuRatio_RatToDef->SaveAs("cCumuRatio_RatToDef.pdf");

  //-- g1e
  TLegend * legG1E = new TLegend(0.4538, 0.7421, 0.9529, 0.9172);
  legG1E->SetBorderSize(0);
  legG1E->SetFillStyle(0);
  legG1E->AddEntry(grG1EPt03_3Eta2p4, "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legG1E->AddEntry(grG1EPt03_3Eta1p0, "0.3 < p_{T} < 3.0 GeV/c, |#eta| < 1.0",              "lp");
  legG1E->AddEntry(grG1EPt1_3Eta2p4,  "1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4",              "lp");
  legG1E->AddEntry(grG1EPromptRECO,   "Prompt RECO: 1.0 < p_{T} < 3.0 GeV/c, |#eta| < 2.4", "lp");

  TCanvas* cG1e = new TCanvas("cG1e", "cG1e", 1000, 500);
  cG1e->Divide(2,1);
  cG1e->cd(1);
  grG1EPt03_3Eta2p4->Draw("ap");
  grG1EPt1_3Eta2p4->Draw("psame");
  grG1EPt03_3Eta1p0->Draw("psame");
  grG1EPromptRECO->Draw("psame");
  legG1E->Draw("same");
  cG1e->cd(2);
  grG1EPt1_3Eta2p4_RatioToDefault->Draw("ap");
  grG1EPt03_3Eta1p0_RatioToDefault->Draw("psame");
  grG1EPromptRECO_RatioToDefault->Draw("psame");
  cG1e->SaveAs("cG1e.pdf");








}

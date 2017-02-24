#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TError.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void sysTkEff(){

  int norder_ = 2;

  double ratioMin = 0.8;
  double ratioMax = 1.2;

  double ratioG1Min = 0.0;
  double ratioG1Max = 1.0;

  double ratioCumuRatioMin = 0.95;
  double ratioCumuRatioMax = 1.05;

  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Standard Unfolding
  TFile * fUnf;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- NoTeff Unfolding
  TFile * fUnfNoTeff;
  TH1D * hUnfoldNoTeff[NCENT][NITER];
  TH1D * hRefoldNoTeff[NCENT][NITER];

  //-- Iter Cutoffs
  bool iterStopFound[NCENT];
  int iterCutoff[NCENT];

  bool iterStopFoundNoTeff[NCENT];
  int iterCutoffNoTeff[NCENT];

  TLatex latex;

  //-- Systematics
  double vn2NoTeff_RatioToNominal[NCENT];
  double vn4NoTeff_RatioToNominal[NCENT];
  double vn6NoTeff_RatioToNominal[NCENT];
  double vn8NoTeff_RatioToNominal[NCENT];
  double gamma1expNoTeff_RatioToNominal[NCENT];
  double vn6vn4NoTeff_RatioToNominal[NCENT];
  double vn8vn4NoTeff_RatioToNominal[NCENT];
  double vn8vn6NoTeff_RatioToNominal[NCENT];

  double vn2NoTeff_RatioToNominal_staterr[NCENT];
  double vn4NoTeff_RatioToNominal_staterr[NCENT];
  double vn6NoTeff_RatioToNominal_staterr[NCENT];
  double vn8NoTeff_RatioToNominal_staterr[NCENT];
  double gamma1expNoTeff_RatioToNominal_staterr[NCENT];
  double vn6vn4NoTeff_RatioToNominal_staterr[NCENT];
  double vn8vn4NoTeff_RatioToNominal_staterr[NCENT];
  double vn8vn6NoTeff_RatioToNominal_staterr[NCENT];

  double vn2NoTeff_PctDiffToNominal[NCENT];
  double vn4NoTeff_PctDiffToNominal[NCENT];
  double vn6NoTeff_PctDiffToNominal[NCENT];
  double vn8NoTeff_PctDiffToNominal[NCENT];
  double gamma1expNoTeff_PctDiffToNominal[NCENT];
  double vn6vn4NoTeff_PctDiffToNominal[NCENT];
  double vn8vn4NoTeff_PctDiffToNominal[NCENT];
  double vn8vn6NoTeff_PctDiffToNominal[NCENT];

  //-- Systematic Performance Plots
  TGraphErrors * grVn2NoTeff_RatioToNominal;
  TGraphErrors * grVn4NoTeff_RatioToNominal;
  TGraphErrors * grVn6NoTeff_RatioToNominal;
  TGraphErrors * grVn8NoTeff_RatioToNominal;
  TGraphErrors * grGamma1ExpNoTeff_RatioToNominal;
  TGraphErrors * grVn6Vn4NoTeff_RatioToNominal;
  TGraphErrors * grVn8Vn4NoTeff_RatioToNominal;
  TGraphErrors * grVn8Vn6NoTeff_RatioToNominal;

  //-- Stat Errors
  TFile * fStat_Nominal;
  TH1D * hVarianceOfMean_Vn2_Nominal;
  TH1D * hVarianceOfMean_Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn6_Nominal;
  TH1D * hVarianceOfMean_Vn8_Nominal;
  TH1D * hVarianceOfMean_Gamma1Exp_Nominal;
  TH1D * hVarianceOfMean_Vn6Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn8Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn8Vn6_Nominal;

  TFile * fStat_NoTeff;
  TH1D * hVarianceOfMean_Vn2_NoTeff;
  TH1D * hVarianceOfMean_Vn4_NoTeff;
  TH1D * hVarianceOfMean_Vn6_NoTeff;
  TH1D * hVarianceOfMean_Vn8_NoTeff;
  TH1D * hVarianceOfMean_Gamma1Exp_NoTeff;
  TH1D * hVarianceOfMean_Vn6Vn4_NoTeff;
  TH1D * hVarianceOfMean_Vn8Vn4_NoTeff;
  TH1D * hVarianceOfMean_Vn8Vn6_NoTeff;

  //-- RelErrors
  TFile * fRelErr;
  TH1D * relErrVn2;
  TH1D * relErrVn4;
  TH1D * relErrVn6;
  TH1D * relErrVn8;
  TH1D * relErrGamma1Exp;
  TH1D * relErrVn6Vn4;
  TH1D * relErrVn8Vn4;
  TH1D * relErrVn8Vn6;

  //
  //-- MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  gErrorIgnoreLevel = kWarning;

  fAna      = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fUnf      = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnfNoTeff = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorTeff.root", "recreate");
  fRelErr->cd();

  relErrVn2        = new TH1D("relErrVn2",       "relErrVn2",       NCENT, centbinsDefault);
  relErrVn4        = new TH1D("relErrVn4",       "relErrVn4",       NCENT, centbinsDefault);
  relErrVn6        = new TH1D("relErrVn6",       "relErrVn6",       NCENT, centbinsDefault);
  relErrVn8        = new TH1D("relErrVn8",       "relErrVn8",       NCENT, centbinsDefault);
  relErrGamma1Exp  = new TH1D("relErrGamma1Exp", "relErrGamma1Exp", NCENT, centbinsDefault);
  relErrVn6Vn4     = new TH1D("relErrVn6Vn4",    "relErrVn6Vn4",    NCENT, centbinsDefault);
  relErrVn8Vn4     = new TH1D("relErrVn8Vn4",    "relErrVn8Vn4",    NCENT, centbinsDefault);
  relErrVn8Vn6     = new TH1D("relErrVn8Vn6",    "relErrVn8Vn6",    NCENT, centbinsDefault);

  //-- Stat errors
  fStat_Nominal = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn2_Nominal       = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn2" );
  hVarianceOfMean_Vn4_Nominal       = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn4" );
  hVarianceOfMean_Vn6_Nominal       = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn6" );
  hVarianceOfMean_Vn8_Nominal       = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn8" );
  hVarianceOfMean_Gamma1Exp_Nominal = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Gamma1Exp" );
  hVarianceOfMean_Vn6Vn4_Nominal    = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Nominal    = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Nominal    = (TH1D*) fStat_Nominal->Get( "hVarianceOfMean_Vn8Vn6" );

  fStat_NoTeff = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/tkEff/StatUncertTkEff_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn2_NoTeff       = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn2_TkEff" );
  hVarianceOfMean_Vn4_NoTeff       = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn4_TkEff" );
  hVarianceOfMean_Vn6_NoTeff       = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn6_TkEff" );
  hVarianceOfMean_Vn8_NoTeff       = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn8_TkEff" );
  hVarianceOfMean_Gamma1Exp_NoTeff = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Gamma1Exp_TkEff" );
  hVarianceOfMean_Vn6Vn4_NoTeff    = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn6Vn4_TkEff" );
  hVarianceOfMean_Vn8Vn4_NoTeff    = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn8Vn4_TkEff" );
  hVarianceOfMean_Vn8Vn6_NoTeff    = (TH1D*) fStat_NoTeff->Get( "hVarianceOfMean_Vn8Vn6_TkEff" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    //-- Initialize iteration cutoff values/booleans
    iterStopFound[icent] = 0;
    iterCutoff[icent]    = 0;

    iterStopFoundNoTeff[icent] = 0;
    iterCutoffNoTeff[icent]    = 0;

    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfold[icent][i] = (TH1D*) fUnf->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold[icent][i]->SetLineColor(col[i]);
      hUnfold[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefold[icent][i] = (TH1D*) fUnf->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold[icent][i]->SetLineWidth(2);
      hRefold[icent][i]->SetLineColor(col[i]);
      hRefold[icent][i]->SetMarkerColor(col[i]);

      //-- Get the NoTeff unfolded histograms
      hUnfoldNoTeff[icent][i] = (TH1D*) fUnfNoTeff->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldNoTeff[icent][i]->SetLineColor(col[i]);
      hUnfoldNoTeff[icent][i]->SetMarkerColor(col[i]);

      //-- Get the NoTeff refolded histograms
      hRefoldNoTeff[icent][i] = (TH1D*) fUnfNoTeff->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldNoTeff[icent][i]->SetLineColor(col[i]);
      hRefoldNoTeff[icent][i]->SetMarkerColor(col[i]);


      //-- Chi squares
      double chi2NDF_Refold       = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      double chi2NDF_Refold_NoTeff = hRefoldNoTeff[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");


      //-- Normal unfolding chi2 check
      if( chi2NDF_Refold < 1.2 && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }

      //-- NoTeff unfolding chi2 check
      if( chi2NDF_Refold_NoTeff < 1.2 && !iterStopFoundNoTeff[icent] ){
        iterStopFoundNoTeff[icent] = 1;
        iterCutoffNoTeff[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFoundNoTeff[icent] ){
	iterStopFoundNoTeff[icent] = 1;
        iterCutoffNoTeff[icent]    = i;
      }

    } //-- End iter loop

  } //-- End cent loop

  //-- Figure out the response matrix uncertainty component
  for(int icent = 0; icent < NCENT; icent++){

    std::cout << Form("=================== Cent %i ===================", icent) << std::endl;

    int i      = iterCutoff[icent];
    int iNoTeff = iterCutoffNoTeff[icent];

    FixUnfold( hUnfold[icent][i] );
    FixUnfold( hUnfoldNoTeff[icent][iNoTeff] );

    EbyECumu cumu( hUnfold[icent][i] );
    EbyECumu cumuNoTeff( hUnfoldNoTeff[icent][iNoTeff] );

    double vn2        = cumu.GetCumu_vn2();
    double vn4        = cumu.GetCumu_vn4();
    double vn6        = cumu.GetCumu_vn6();
    double vn8        = cumu.GetCumu_vn8();
    double vn2NoTeff   = cumuNoTeff.GetCumu_vn2();
    double vn4NoTeff   = cumuNoTeff.GetCumu_vn4();
    double vn6NoTeff   = cumuNoTeff.GetCumu_vn6();
    double vn8NoTeff   = cumuNoTeff.GetCumu_vn8();

    std::cout << "vn2 = " << vn2 << "\tvn2NoTeff = " << vn2NoTeff << "\tratio = " << vn2NoTeff / vn2 << std::endl;

    double vn2_staterr         = sqrt( hVarianceOfMean_Vn2_Nominal->GetBinContent(icent+1) );
    double vn4_staterr         = sqrt( hVarianceOfMean_Vn4_Nominal->GetBinContent(icent+1) );
    double vn6_staterr         = sqrt( hVarianceOfMean_Vn6_Nominal->GetBinContent(icent+1) );
    double vn8_staterr         = sqrt( hVarianceOfMean_Vn8_Nominal->GetBinContent(icent+1) );
    double vn2NoTeff_staterr   = sqrt( hVarianceOfMean_Vn2_NoTeff->GetBinContent(icent+1) );
    double vn4NoTeff_staterr   = sqrt( hVarianceOfMean_Vn4_NoTeff->GetBinContent(icent+1) );
    double vn6NoTeff_staterr   = sqrt( hVarianceOfMean_Vn6_NoTeff->GetBinContent(icent+1) );
    double vn8NoTeff_staterr   = sqrt( hVarianceOfMean_Vn8_NoTeff->GetBinContent(icent+1) );

    double gamma1exp      = cumu.GetGamma1Exp();
    double gamma1expNoTeff = cumuNoTeff.GetGamma1Exp();

    double gamma1exp_staterr       = sqrt( hVarianceOfMean_Gamma1Exp_Nominal->GetBinContent(icent+1) );
    double gamma1expNoTeff_staterr = sqrt( hVarianceOfMean_Gamma1Exp_NoTeff->GetBinContent(icent+1) );

    double vn6vn4;
    double vn6vn4NoTeff;
    if( vn6 == 0 || vn4 == 0 ) vn6vn4 = 0;
    else                       vn6vn4 = vn6 / vn4;
    if( vn6NoTeff == 0 || vn4NoTeff == 0) vn6vn4NoTeff = 0;
    else                                  vn6vn4NoTeff = vn6NoTeff / vn4NoTeff;

    double vn8vn4;
    double vn8vn4NoTeff;
    if( vn8 == 0 || vn4 == 0 ) vn8vn4 = 0;
    else                       vn8vn4 = vn8 / vn4;
    if( vn8NoTeff == 0 || vn4NoTeff == 0) vn8vn4NoTeff = 0;
    else                                  vn8vn4NoTeff = vn8NoTeff / vn4NoTeff;

    double vn8vn6;
    double vn8vn6NoTeff;
    if( vn8 == 0 || vn6 == 0 ) vn8vn6 = 0;
    else                       vn8vn6 = vn8 / vn6;
    if( vn8NoTeff == 0 || vn6NoTeff == 0) vn8vn6NoTeff = 0;
    else                                  vn8vn6NoTeff = vn8NoTeff / vn6NoTeff;

    double vn6vn4_staterr       = sqrt( hVarianceOfMean_Vn6Vn4_Nominal->GetBinContent(icent+1) );
    double vn6vn4NoTeff_staterr = sqrt( hVarianceOfMean_Vn6Vn4_NoTeff->GetBinContent(icent+1) );
    double vn8vn4_staterr       = sqrt( hVarianceOfMean_Vn8Vn4_Nominal->GetBinContent(icent+1) );
    double vn8vn4NoTeff_staterr = sqrt( hVarianceOfMean_Vn8Vn4_NoTeff->GetBinContent(icent+1) );
    double vn8vn6_staterr       = sqrt( hVarianceOfMean_Vn8Vn6_Nominal->GetBinContent(icent+1) );
    double vn8vn6NoTeff_staterr = sqrt( hVarianceOfMean_Vn8Vn6_NoTeff->GetBinContent(icent+1) );

    //-- Calculate ratios
    if( vn2 == 0 || vn2NoTeff == 0 ) vn2NoTeff_RatioToNominal[icent] = 0;
    else                             vn2NoTeff_RatioToNominal[icent] = vn2NoTeff / vn2;

    if( vn4 == 0 || vn4NoTeff == 0 ) vn4NoTeff_RatioToNominal[icent] = 0;
    else                             vn4NoTeff_RatioToNominal[icent] = vn4NoTeff / vn4;

    if( vn6 == 0 || vn6NoTeff == 0 ) vn6NoTeff_RatioToNominal[icent] = 0;
    else                             vn6NoTeff_RatioToNominal[icent] = vn6NoTeff / vn6;

    if( vn8 == 0 || vn8NoTeff == 0 ) vn8NoTeff_RatioToNominal[icent] = 0;
    else                             vn8NoTeff_RatioToNominal[icent] = vn8NoTeff / vn8;

    if( gamma1exp == 0 || gamma1expNoTeff == 0 ) gamma1expNoTeff_RatioToNominal[icent] = 0;
    else                                         gamma1expNoTeff_RatioToNominal[icent] = fabs(1. - gamma1expNoTeff / gamma1exp );

    if( vn6vn4 == 0 || vn6vn4NoTeff == 0 ) vn6vn4NoTeff_RatioToNominal[icent] = 0;
    else                                   vn6vn4NoTeff_RatioToNominal[icent] = vn6vn4NoTeff / vn6vn4;

    if( vn8vn4 == 0 || vn8vn4NoTeff == 0 ) vn8vn4NoTeff_RatioToNominal[icent] = 0;
    else                                   vn8vn4NoTeff_RatioToNominal[icent] = vn8vn4NoTeff / vn8vn4;

    if( vn8vn6 == 0 || vn8vn6NoTeff == 0 ) vn8vn6NoTeff_RatioToNominal[icent] = 0;
    else                                   vn8vn6NoTeff_RatioToNominal[icent] = vn8vn6NoTeff / vn8vn6;

    //-- Ratio errors
    vn2NoTeff_RatioToNominal_staterr[icent]       = sqrt( pow( vn2NoTeff_staterr/vn2,2) + pow(vn2NoTeff*vn2_staterr/vn2/vn2,2) );
    vn4NoTeff_RatioToNominal_staterr[icent]       = sqrt( pow( vn4NoTeff_staterr/vn4,2) + pow(vn4NoTeff*vn4_staterr/vn4/vn4,2) );
    vn6NoTeff_RatioToNominal_staterr[icent]       = sqrt( pow( vn6NoTeff_staterr/vn6,2) + pow(vn6NoTeff*vn6_staterr/vn6/vn6,2) );
    vn8NoTeff_RatioToNominal_staterr[icent]       = sqrt( pow( vn8NoTeff_staterr/vn8,2) + pow(vn8NoTeff*vn8_staterr/vn8/vn8,2) );
    gamma1expNoTeff_RatioToNominal_staterr[icent] = sqrt( pow( gamma1expNoTeff_staterr/gamma1exp,2) + pow(gamma1expNoTeff*gamma1exp_staterr/gamma1exp/gamma1exp,2) );
    vn6vn4NoTeff_RatioToNominal_staterr[icent]    = sqrt( pow( vn6vn4NoTeff_staterr/vn6vn4,2) + pow(vn6vn4NoTeff*vn6vn4_staterr/vn6vn4/vn6vn4,2) );
    vn8vn4NoTeff_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn4NoTeff_staterr/vn8vn4,2) + pow(vn8vn4NoTeff*vn8vn4_staterr/vn8vn4/vn8vn4,2) );
    vn8vn6NoTeff_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn6NoTeff_staterr/vn8vn6,2) + pow(vn8vn6NoTeff*vn8vn6_staterr/vn8vn6/vn8vn6,2) );

    //-- Calculate pct difference relative to nominal
    if( vn2 == 0 || vn2NoTeff == 0 || compatibleWithOne(vn2NoTeff_RatioToNominal[icent], vn2NoTeff_RatioToNominal_staterr[icent]) ) vn2NoTeff_PctDiffToNominal[icent] = 0;
    else                            vn2NoTeff_PctDiffToNominal[icent] = fabs( vn2NoTeff - vn2 ) / fabs( vn2 );

    if( vn4 == 0 || vn4NoTeff == 0 || compatibleWithOne(vn4NoTeff_RatioToNominal[icent], vn4NoTeff_RatioToNominal_staterr[icent])) vn4NoTeff_PctDiffToNominal[icent] = 0;
    else                            vn4NoTeff_PctDiffToNominal[icent] = fabs( vn4NoTeff - vn4 ) / fabs( vn4 );

    if( vn6 == 0 || vn6NoTeff == 0 || compatibleWithOne(vn6NoTeff_RatioToNominal[icent], vn6NoTeff_RatioToNominal_staterr[icent])) vn6NoTeff_PctDiffToNominal[icent] = 0;
    else                            vn6NoTeff_PctDiffToNominal[icent] = fabs( vn6NoTeff - vn6 ) / fabs( vn6 );

    if( vn8 == 0 || vn8NoTeff == 0 || compatibleWithOne(vn8NoTeff_RatioToNominal[icent], vn8NoTeff_RatioToNominal_staterr[icent])) vn8NoTeff_PctDiffToNominal[icent] = 0;
    else                            vn8NoTeff_PctDiffToNominal[icent] = fabs( vn8NoTeff - vn8 ) / fabs( vn8 );

    if( gamma1exp == 0 || gamma1expNoTeff == 0 || compatibleWithOne(gamma1expNoTeff_RatioToNominal[icent], gamma1expNoTeff_RatioToNominal_staterr[icent]) ) gamma1expNoTeff_PctDiffToNominal[icent] = 0;
    else                                        gamma1expNoTeff_PctDiffToNominal[icent] = fabs( gamma1expNoTeff - gamma1exp ) / fabs( gamma1exp );

    if( vn6vn4 == 0 || vn6vn4NoTeff == 0 || compatibleWithOne(vn6vn4NoTeff_RatioToNominal[icent], vn6vn4NoTeff_RatioToNominal_staterr[icent]) ) vn6vn4NoTeff_PctDiffToNominal[icent] = 0;
    else                                  vn6vn4NoTeff_PctDiffToNominal[icent] = fabs( vn6vn4NoTeff - vn6vn4 ) / fabs( vn6vn4 );

    if( vn8vn4 == 0 || vn8vn4NoTeff == 0 || compatibleWithOne(vn8vn4NoTeff_RatioToNominal[icent], vn8vn4NoTeff_RatioToNominal_staterr[icent]) ) vn8vn4NoTeff_PctDiffToNominal[icent] = 0;
    else                                  vn8vn4NoTeff_PctDiffToNominal[icent] = fabs( vn8vn4NoTeff - vn8vn4 ) / fabs( vn8vn4 );

    if( vn8vn6 == 0 || vn8vn6NoTeff == 0 || compatibleWithOne(vn8vn6NoTeff_RatioToNominal[icent], vn8vn6NoTeff_RatioToNominal_staterr[icent]) ) vn8vn6NoTeff_PctDiffToNominal[icent] = 0;
    else                                  vn8vn6NoTeff_PctDiffToNominal[icent] = fabs( vn8vn6NoTeff - vn8vn6 ) / fabs( vn8vn6 );

  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT];
  for(int i = 0; i < NCENT; i++) cErr[i] = 0;
  grVn2NoTeff_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn2NoTeff_RatioToNominal, cErr, vn2NoTeff_RatioToNominal_staterr);
  formatGraph(grVn2NoTeff_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} [NoTkEff] / [Nominal]", norder_), 1, 24, "grVn2NoTeff_RatioToNominal");
  grVn4NoTeff_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn4NoTeff_RatioToNominal, cErr, vn4NoTeff_RatioToNominal_staterr);
  formatGraph(grVn4NoTeff_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} [NoTkEff] / [Nominal]", norder_), kSpring+4, 25, "grVn4NoTeff_RatioToNominal");
  grVn6NoTeff_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn6NoTeff_RatioToNominal, cErr, vn6NoTeff_RatioToNominal_staterr);
  formatGraph(grVn6NoTeff_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} [NoTkEff] / [Nominal]", norder_), 6, 28, "grVn6NoTeff_RatioToNominal");
  grVn8NoTeff_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn8NoTeff_RatioToNominal, cErr, vn8NoTeff_RatioToNominal_staterr);
  formatGraph(grVn8NoTeff_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} [NoTkEff] / [Nominal]", norder_), kOrange+7, 27, "grVn8NoTeff_RatioToNominal");
  grGamma1ExpNoTeff_RatioToNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expNoTeff_RatioToNominal, cErr, gamma1expNoTeff_RatioToNominal_staterr);
  formatGraph(grGamma1ExpNoTeff_RatioToNominal, "Centrality %", ratioG1Min, ratioG1Max, "|1-#gamma_{1}^{exp} [NoTkEff] / [Nominal]|", 2, 20, "grGamma1ExpNoTeff_RatioToNominal");
  grVn6Vn4NoTeff_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4NoTeff_RatioToNominal, cErr, vn6vn4NoTeff_RatioToNominal_staterr);
  formatGraph(grVn6Vn4NoTeff_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{6}/v_{%i}{4} [NoTkEff] / [Nominal]", norder_, norder_), 4, 21, "grVn6Vn4NoTeff_RatioToNominal");
  grVn8Vn4NoTeff_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4NoTeff_RatioToNominal, cErr, vn8vn4NoTeff_RatioToNominal_staterr);
  formatGraph(grVn8Vn4NoTeff_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{4} [NoTkEff] / [Nominal]", norder_, norder_), kGreen+2, 34, "grVn8Vn4NoTeff_RatioToNominal");
  grVn8Vn6NoTeff_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6NoTeff_RatioToNominal, cErr, vn8vn6NoTeff_RatioToNominal_staterr);
  formatGraph(grVn8Vn6NoTeff_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{6} [NoTkEff] / [Nominal]", norder_, norder_), kViolet-1, 33, "grVn8Vn6NoTeff_RatioToNominal");


  TLine * line = new TLine(grVn2NoTeff_RatioToNominal->GetXaxis()->GetXmin(), 1.0, grVn2NoTeff_RatioToNominal->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);


  TCanvas * cCumuSys = new TCanvas("cCumuSys", "cCumuSys", 1000, 1000);
  cCumuSys->Divide(2,2);
  cCumuSys->cd(1);
  grVn2NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(2);
  grVn4NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(3);
  grVn6NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(4);
  grVn8NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->SaveAs("../../plots/systematicStudies/cSysTkEff_CumuCent.pdf");

  TCanvas * cGamma1ExpSys = new TCanvas("cGamma1ExpSys", "cGamma1ExpSys", 500, 500);
  cGamma1ExpSys->cd();
  grGamma1ExpNoTeff_RatioToNominal->Draw("ap");
  //line->Draw("same");
  cGamma1ExpSys->SaveAs("../../plots/systematicStudies/cSysTkEff_Gamma1ExpCent.pdf");

  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  grVn6Vn4NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6NoTeff_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/cSysTkEff_CumuRatioCent.pdf");

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){

    std::cout<<"---------------"<<std::endl;

    relErrVn2->SetBinContent(icent+1, vn2NoTeff_PctDiffToNominal[icent]);
    relErrVn4->SetBinContent(icent+1, vn4NoTeff_PctDiffToNominal[icent]);
    relErrVn6->SetBinContent(icent+1, vn6NoTeff_PctDiffToNominal[icent]);
    relErrVn8->SetBinContent(icent+1, vn6NoTeff_PctDiffToNominal[icent]);
    relErrGamma1Exp->SetBinContent(icent+1, gamma1expNoTeff_PctDiffToNominal[icent]);
    relErrVn6Vn4->SetBinContent(icent+1, vn6vn4NoTeff_PctDiffToNominal[icent]);
    relErrVn8Vn4->SetBinContent(icent+1, vn8vn4NoTeff_PctDiffToNominal[icent]);
    relErrVn8Vn6->SetBinContent(icent+1, vn8vn6NoTeff_PctDiffToNominal[icent]);

  }

  fRelErr->Write();

}

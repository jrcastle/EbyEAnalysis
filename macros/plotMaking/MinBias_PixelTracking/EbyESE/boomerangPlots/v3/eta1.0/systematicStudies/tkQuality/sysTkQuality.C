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

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void sysTkQuality(){

  int norder_ = 3;
  double tkEta = 1.0;

  double vn2Min = 0.00;
  double vn2Max = 0.40;
  double vn4Min = 0.00;
  double vn4Max = 0.40;
  double vn6Min = 0.00;
  double vn6Max = 0.40;
  double vn8Min = 0.00;
  double vn8Max = 0.40;

  double gamma1expMin = -1.0;
  double gamma1expMax = 0.5;
  double vn6vn4Min    = 0.95;
  double vn6vn4Max    = 1.02;
  double vn8vn4Min    = 0.95;
  double vn8vn4Max    = 1.02;
  double vn8vn6Min    = 0.95;
  double vn8vn6Max    = 1.02;

  double ratioMin = 0.98;
  double ratioMax = 1.02;
  double gamma1ExpRatioMin = 0.;
  double gamma1ExpRatioMax = 0.2;

  //-- Loose Cuts
  TFile * fAnaLoose;
  TFile * fUnfLoose;
  TH1D * hObsLoose[NCENT];
  TH1D * hUnfoldLoose[NCENT][NITER];
  TH1D * hRefoldLoose[NCENT][NITER];

  TFile * fStatLoose;
  TH1D * hVarianceOfMean_Vn2_Loose;
  TH1D * hVarianceOfMean_Vn4_Loose;
  TH1D * hVarianceOfMean_Vn6_Loose;
  TH1D * hVarianceOfMean_Vn8_Loose;
  TH1D * hVarianceOfMean_Gamma1Exp_Loose;
  TH1D * hVarianceOfMean_Vn6Vn4_Loose;
  TH1D * hVarianceOfMean_Vn8Vn4_Loose;
  TH1D * hVarianceOfMean_Vn8Vn6_Loose;

  TGraphErrors * grVn2Loose;
  TGraphErrors * grVn4Loose;
  TGraphErrors * grVn6Loose;
  TGraphErrors * grVn8Loose;
  TGraphErrors * grGamma1ExpLoose;
  TGraphErrors * grVn6Vn4Loose;
  TGraphErrors * grVn8Vn4Loose;
  TGraphErrors * grVn8Vn6Loose;

  double vn2Loose[NCENT];
  double vn4Loose[NCENT];
  double vn6Loose[NCENT];
  double vn8Loose[NCENT];
  double gamma1expLoose[NCENT];
  double vn6vn4Loose[NCENT];
  double vn8vn4Loose[NCENT];
  double vn8vn6Loose[NCENT];

  double vn2Loose_staterr[NCENT];
  double vn4Loose_staterr[NCENT];
  double vn6Loose_staterr[NCENT];
  double vn8Loose_staterr[NCENT];
  double gamma1expLoose_staterr[NCENT];
  double vn6vn4Loose_staterr[NCENT];
  double vn8vn4Loose_staterr[NCENT];
  double vn8vn6Loose_staterr[NCENT];

  TGraphErrors * grVn2Loose_RatiotoNominal;
  TGraphErrors * grVn4Loose_RatiotoNominal;
  TGraphErrors * grVn6Loose_RatiotoNominal;
  TGraphErrors * grVn8Loose_RatiotoNominal;
  TGraphErrors * grGamma1ExpLoose_RatiotoNominal;
  TGraphErrors * grVn6Vn4Loose_RatiotoNominal;
  TGraphErrors * grVn8Vn4Loose_RatiotoNominal;
  TGraphErrors * grVn8Vn6Loose_RatiotoNominal;

  double vn2Loose_RatiotoNominal[NCENT];
  double vn4Loose_RatiotoNominal[NCENT];
  double vn6Loose_RatiotoNominal[NCENT];
  double vn8Loose_RatiotoNominal[NCENT];
  double gamma1expLoose_RatiotoNominal[NCENT];
  double vn6vn4Loose_RatiotoNominal[NCENT];
  double vn8vn4Loose_RatiotoNominal[NCENT];
  double vn8vn6Loose_RatiotoNominal[NCENT];

  double vn2Loose_RatiotoNominal_staterr[NCENT];
  double vn4Loose_RatiotoNominal_staterr[NCENT];
  double vn6Loose_RatiotoNominal_staterr[NCENT];
  double vn8Loose_RatiotoNominal_staterr[NCENT];
  double gamma1expLoose_RatiotoNominal_staterr[NCENT];
  double vn6vn4Loose_RatiotoNominal_staterr[NCENT];
  double vn8vn4Loose_RatiotoNominal_staterr[NCENT];
  double vn8vn6Loose_RatiotoNominal_staterr[NCENT];

  int iterCutoffLoose[NCENT];

  //-- Tight Cuts
  TFile * fAnaTight;
  TFile * fUnfTight;
  TH1D * hObsTight[NCENT];
  TH1D * hUnfoldTight[NCENT][NITER];
  TH1D * hRefoldTight[NCENT][NITER];

  TFile * fStatTight;
  TH1D * hVarianceOfMean_Vn2_Tight;
  TH1D * hVarianceOfMean_Vn4_Tight;
  TH1D * hVarianceOfMean_Vn6_Tight;
  TH1D * hVarianceOfMean_Vn8_Tight;
  TH1D * hVarianceOfMean_Gamma1Exp_Tight;
  TH1D * hVarianceOfMean_Vn6Vn4_Tight;
  TH1D * hVarianceOfMean_Vn8Vn4_Tight;
  TH1D * hVarianceOfMean_Vn8Vn6_Tight;

  TGraphErrors * grVn2Tight;
  TGraphErrors * grVn4Tight;
  TGraphErrors * grVn6Tight;
  TGraphErrors * grVn8Tight;
  TGraphErrors * grGamma1ExpTight;
  TGraphErrors * grVn6Vn4Tight;
  TGraphErrors * grVn8Vn4Tight;
  TGraphErrors * grVn8Vn6Tight;

  double vn2Tight[NCENT];
  double vn4Tight[NCENT];
  double vn6Tight[NCENT];
  double vn8Tight[NCENT];
  double gamma1expTight[NCENT];
  double vn6vn4Tight[NCENT];
  double vn8vn4Tight[NCENT];
  double vn8vn6Tight[NCENT];

  double vn2Tight_staterr[NCENT];
  double vn4Tight_staterr[NCENT];
  double vn6Tight_staterr[NCENT];
  double vn8Tight_staterr[NCENT];
  double gamma1expTight_staterr[NCENT];
  double vn6vn4Tight_staterr[NCENT];
  double vn8vn4Tight_staterr[NCENT];
  double vn8vn6Tight_staterr[NCENT];

  TGraphErrors * grVn2Tight_RatiotoNominal;
  TGraphErrors * grVn4Tight_RatiotoNominal;
  TGraphErrors * grVn6Tight_RatiotoNominal;
  TGraphErrors * grVn8Tight_RatiotoNominal;
  TGraphErrors * grGamma1ExpTight_RatiotoNominal;
  TGraphErrors * grVn6Vn4Tight_RatiotoNominal;
  TGraphErrors * grVn8Vn4Tight_RatiotoNominal;
  TGraphErrors * grVn8Vn6Tight_RatiotoNominal;

  double vn2Tight_RatiotoNominal[NCENT];
  double vn4Tight_RatiotoNominal[NCENT];
  double vn6Tight_RatiotoNominal[NCENT];
  double vn8Tight_RatiotoNominal[NCENT];
  double gamma1expTight_RatiotoNominal[NCENT];
  double vn6vn4Tight_RatiotoNominal[NCENT];
  double vn8vn4Tight_RatiotoNominal[NCENT];
  double vn8vn6Tight_RatiotoNominal[NCENT];

  double vn2Tight_RatiotoNominal_staterr[NCENT];
  double vn4Tight_RatiotoNominal_staterr[NCENT];
  double vn6Tight_RatiotoNominal_staterr[NCENT];
  double vn8Tight_RatiotoNominal_staterr[NCENT];
  double gamma1expTight_RatiotoNominal_staterr[NCENT];
  double vn6vn4Tight_RatiotoNominal_staterr[NCENT];
  double vn8vn4Tight_RatiotoNominal_staterr[NCENT];
  double vn8vn6Tight_RatiotoNominal_staterr[NCENT];

  int iterCutoffTight[NCENT];

  //-- Nominal
  TFile * fAnaNominal;
  TFile * fUnfNominal;
  TH1D * hObsNominal[NCENT];
  TH1D * hUnfoldNominal[NCENT][NITER];
  TH1D * hRefoldNominal[NCENT][NITER];

  TFile * fStatNominal;
  TH1D * hVarianceOfMean_Vn2_Nominal;
  TH1D * hVarianceOfMean_Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn6_Nominal;
  TH1D * hVarianceOfMean_Vn8_Nominal;
  TH1D * hVarianceOfMean_Gamma1Exp_Nominal;
  TH1D * hVarianceOfMean_Vn6Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn8Vn4_Nominal;
  TH1D * hVarianceOfMean_Vn8Vn6_Nominal;

  TGraphErrors * grVn2Nominal;
  TGraphErrors * grVn4Nominal;
  TGraphErrors * grVn6Nominal;
  TGraphErrors * grVn8Nominal;
  TGraphErrors * grGamma1ExpNominal;
  TGraphErrors * grVn6Vn4Nominal;
  TGraphErrors * grVn8Vn4Nominal;
  TGraphErrors * grVn8Vn6Nominal;

  double vn2Nominal[NCENT];
  double vn4Nominal[NCENT];
  double vn6Nominal[NCENT];
  double vn8Nominal[NCENT];
  double gamma1expNominal[NCENT];
  double vn6vn4Nominal[NCENT];
  double vn8vn4Nominal[NCENT];
  double vn8vn6Nominal[NCENT];

  double vn2Nominal_staterr[NCENT];
  double vn4Nominal_staterr[NCENT];
  double vn6Nominal_staterr[NCENT];
  double vn8Nominal_staterr[NCENT];
  double gamma1expNominal_staterr[NCENT];
  double vn6vn4Nominal_staterr[NCENT];
  double vn8vn4Nominal_staterr[NCENT];
  double vn8vn6Nominal_staterr[NCENT];

  int iterCutoffNominal[NCENT];

  TLatex latex;

  //
  //-- MAIN
  // 

  setTDRStyle();
  latex.SetNDC();

  //-- Grab files
  fAnaLoose    = new TFile( "loose/AnalyzerResults/CastleEbyE.root" );
  fUnfLoose    = new TFile( Form("loose/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaTight = new TFile( "tight/AnalyzerResults/CastleEbyE.root" );
  fUnfTight = new TFile( Form("tight/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaNominal   = new TFile( "../../AnalyzerResults/CastleEbyE.root"  );
  fUnfNominal   = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Stat Errors
  fStatLoose = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/tkQuality/loose/StatUncertTkQualityLoose_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Loose       = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn2_TkQualityLoose" );
  hVarianceOfMean_Vn4_Loose       = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn4_TkQualityLoose" );
  hVarianceOfMean_Vn6_Loose       = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn6_TkQualityLoose" );
  hVarianceOfMean_Vn8_Loose       = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn8_TkQualityLoose" );
  hVarianceOfMean_Gamma1Exp_Loose = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Gamma1Exp_TkQualityLoose" );
  hVarianceOfMean_Vn6Vn4_Loose    = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn6Vn4_TkQualityLoose" );
  hVarianceOfMean_Vn8Vn4_Loose    = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn8Vn4_TkQualityLoose" );
  hVarianceOfMean_Vn8Vn6_Loose    = (TH1D*) fStatLoose->Get( "hVarianceOfMean_Vn8Vn6_TkQualityLoose" );

  fStatTight = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/tkQuality/tight/StatUncertTkQualityTight_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Tight       = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn2_TkQualityTight" );
  hVarianceOfMean_Vn4_Tight       = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn4_TkQualityTight" );
  hVarianceOfMean_Vn6_Tight       = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn6_TkQualityTight" );
  hVarianceOfMean_Vn8_Tight       = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn8_TkQualityTight" );
  hVarianceOfMean_Gamma1Exp_Tight = (TH1D*) fStatTight->Get( "hVarianceOfMean_Gamma1Exp_TkQualityTight" );
  hVarianceOfMean_Vn6Vn4_Tight    = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn6Vn4_TkQualityTight" );
  hVarianceOfMean_Vn8Vn4_Tight    = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn8Vn4_TkQualityTight" );
  hVarianceOfMean_Vn8Vn6_Tight    = (TH1D*) fStatTight->Get( "hVarianceOfMean_Vn8Vn6_TkQualityTight" );

  fStatNominal = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertainties_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Nominal       = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn2" );
  hVarianceOfMean_Vn4_Nominal       = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn4" );
  hVarianceOfMean_Vn6_Nominal       = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn6" );
  hVarianceOfMean_Vn8_Nominal       = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn8" );
  hVarianceOfMean_Gamma1Exp_Nominal = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Gamma1Exp" );
  hVarianceOfMean_Vn6Vn4_Nominal    = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Nominal    = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Nominal    = (TH1D*) fStatNominal->Get( "hVarianceOfMean_Vn8Vn6" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- -------------------- Loose --------------------
    hObsLoose[icent] = (TH1D*) fAnaLoose->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsLoose[icent]->SetName( Form("hVnFull_c%i_loose", icent) );

    //-- Iter loop
    iterCutoffLoose[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldLoose[icent][i] = (TH1D*) fUnfLoose->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldLoose[icent][i]->SetName( Form("hreco%i_c%i_loose", iter[i], icent) );
      hUnfoldLoose[icent][i]->SetLineColor(col[i]);
      hUnfoldLoose[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldLoose[icent][i] = (TH1D*) fUnfLoose->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldLoose[icent][i]->SetName( Form("hrefold%i_c%i_loose", iter[i], icent) );
      hRefoldLoose[icent][i]->SetLineWidth(2);
      hRefoldLoose[icent][i]->SetLineColor(col[i]);
      hRefoldLoose[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldLoose[icent][i]->Chi2Test(hObsLoose[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
	iterCutoffLoose[icent] = i;

	FixUnfold( hUnfoldLoose[icent][i] );
	EbyECumu cumu(hUnfoldLoose[icent][i]);
	vn2Loose[icent]  = cumu.GetCumu_vn2();
	vn4Loose[icent]  = cumu.GetCumu_vn4();
	vn6Loose[icent]  = cumu.GetCumu_vn6();
	vn8Loose[icent]  = cumu.GetCumu_vn8();
	gamma1expLoose[icent] = cumu.GetGamma1Exp();

	if( vn4Loose[icent] == 0 ){
	  vn6vn4Loose[icent] = 0;
	  vn8vn4Loose[icent] = 0;
	}
	else{
	  vn6vn4Loose[icent] = vn6Loose[icent] / vn4Loose[icent];
	  vn8vn4Loose[icent] = vn8Loose[icent] / vn4Loose[icent];
	}
	if( vn6Loose[icent] == 0 ) vn8vn6Loose[icent] = 0;
	else                      vn8vn6Loose[icent] = vn8Loose[icent] / vn6Loose[icent];

	vn2Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Loose->GetBinContent(icent+1) );
	vn4Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Loose->GetBinContent(icent+1) );
	vn6Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Loose->GetBinContent(icent+1) );
	vn8Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Loose->GetBinContent(icent+1) );
	gamma1expLoose_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Loose->GetBinContent(icent+1) );
	vn6vn4Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Loose->GetBinContent(icent+1) );
	vn8vn4Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Loose->GetBinContent(icent+1) );
	vn8vn6Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Loose->GetBinContent(icent+1) );

	break;
      }
      if( i == NITER - 1 ){
	iterCutoffLoose[icent] = i;

	FixUnfold( hUnfoldLoose[icent][i] );
        EbyECumu cumu(hUnfoldLoose[icent][i]);
        vn2Loose[icent]  = cumu.GetCumu_vn2();
        vn4Loose[icent]  = cumu.GetCumu_vn4();
        vn6Loose[icent]  = cumu.GetCumu_vn6();
	vn8Loose[icent]  = cumu.GetCumu_vn8();
        gamma1expLoose[icent] = cumu.GetGamma1Exp();

        if( vn4Loose[icent] == 0 ){
          vn6vn4Loose[icent] = 0;
          vn8vn4Loose[icent] = 0;
	}
        else{
          vn6vn4Loose[icent] = vn6Loose[icent] / vn4Loose[icent];
          vn8vn4Loose[icent] = vn8Loose[icent] / vn4Loose[icent];
	}
        if( vn6Loose[icent] == 0 ) vn8vn6Loose[icent] = 0;
        else                      vn8vn6Loose[icent] = vn8Loose[icent] / vn6Loose[icent];

	vn2Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Loose->GetBinContent(icent+1) );
        vn4Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Loose->GetBinContent(icent+1) );
        vn6Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Loose->GetBinContent(icent+1) );
        vn8Loose_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Loose->GetBinContent(icent+1) );
        gamma1expLoose_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Loose->GetBinContent(icent+1) );
        vn6vn4Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Loose->GetBinContent(icent+1) );
        vn8vn4Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Loose->GetBinContent(icent+1) );
        vn8vn6Loose_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Loose->GetBinContent(icent+1) );

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Tight --------------------
    hObsTight[icent] = (TH1D*) fAnaTight->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsTight[icent]->SetName( Form("hVnFull_c%i_tight", icent) );

    //-- Iter loop
    iterCutoffTight[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldTight[icent][i] = (TH1D*) fUnfTight->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldTight[icent][i]->SetName( Form("hreco%i_c%i_tight", iter[i], icent) );
      hUnfoldTight[icent][i]->SetLineColor(col[i]);
      hUnfoldTight[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldTight[icent][i] = (TH1D*) fUnfTight->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldTight[icent][i]->SetName( Form("hrefold%i_c%i_tight", iter[i], icent) );
      hRefoldTight[icent][i]->SetLineWidth(2);
      hRefoldTight[icent][i]->SetLineColor(col[i]);
      hRefoldTight[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldTight[icent][i]->Chi2Test(hObsTight[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
        iterCutoffTight[icent] = i;

	FixUnfold( hUnfoldTight[icent][i] );
        EbyECumu cumu(hUnfoldTight[icent][i]);
        vn2Tight[icent]  = cumu.GetCumu_vn2();
        vn4Tight[icent]  = cumu.GetCumu_vn4();
        vn6Tight[icent]  = cumu.GetCumu_vn6();
	vn8Tight[icent]  = cumu.GetCumu_vn8();
        gamma1expTight[icent] = cumu.GetGamma1Exp();

        if( vn4Tight[icent] == 0 ){
          vn6vn4Tight[icent] = 0;
          vn8vn4Tight[icent] = 0;
	}
        else{
          vn6vn4Tight[icent] = vn6Tight[icent] / vn4Tight[icent];
          vn8vn4Tight[icent] = vn8Tight[icent] / vn4Tight[icent];
	}
        if( vn6Tight[icent] == 0 ) vn8vn6Tight[icent] = 0;
        else                      vn8vn6Tight[icent] = vn8Tight[icent] / vn6Tight[icent];

	vn2Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Tight->GetBinContent(icent+1) );
        vn4Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Tight->GetBinContent(icent+1) );
        vn6Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Tight->GetBinContent(icent+1) );
        vn8Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Tight->GetBinContent(icent+1) );
        gamma1expTight_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Tight->GetBinContent(icent+1) );
        vn6vn4Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Tight->GetBinContent(icent+1) );
        vn8vn4Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Tight->GetBinContent(icent+1) );
        vn8vn6Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Tight->GetBinContent(icent+1) );

        break;
      }
      if( i == NITER - 1 ){
        iterCutoffTight[icent] = i;

	FixUnfold( hUnfoldTight[icent][i] );
        EbyECumu cumu(hUnfoldTight[icent][i]);
        vn2Tight[icent]  = cumu.GetCumu_vn2();
        vn4Tight[icent]  = cumu.GetCumu_vn4();
        vn6Tight[icent]  = cumu.GetCumu_vn6();
        vn8Tight[icent]  = cumu.GetCumu_vn8();
        gamma1expTight[icent] = cumu.GetGamma1Exp();

        if( vn4Tight[icent] == 0 ){
          vn6vn4Tight[icent] = 0;
          vn8vn4Tight[icent] = 0;
        }
        else{
          vn6vn4Tight[icent] = vn6Tight[icent] / vn4Tight[icent];
          vn8vn4Tight[icent] = vn8Tight[icent] / vn4Tight[icent];
        }
        if( vn6Tight[icent] == 0 ) vn8vn6Tight[icent] = 0;
        else                      vn8vn6Tight[icent] = vn8Tight[icent] / vn6Tight[icent];

	vn2Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Tight->GetBinContent(icent+1) );
        vn4Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Tight->GetBinContent(icent+1) );
        vn6Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Tight->GetBinContent(icent+1) );
        vn8Tight_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Tight->GetBinContent(icent+1) );
        gamma1expTight_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Tight->GetBinContent(icent+1) );
        vn6vn4Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Tight->GetBinContent(icent+1) );
        vn8vn4Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Tight->GetBinContent(icent+1) );
        vn8vn6Tight_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Tight->GetBinContent(icent+1) );

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Nominal --------------------
    hObsNominal[icent] = (TH1D*) fAnaNominal->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsNominal[icent]->SetName( Form("hVnFull_c%i_nominal", icent) );

    //-- Iter loop
    iterCutoffNominal[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldNominal[icent][i] = (TH1D*) fUnfNominal->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldNominal[icent][i]->SetName( Form("hreco%i_c%i_nominal", iter[i], icent) );
      hUnfoldNominal[icent][i]->SetLineColor(col[i]);
      hUnfoldNominal[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldNominal[icent][i] = (TH1D*) fUnfNominal->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldNominal[icent][i]->SetName( Form("hrefold%i_c%i_nominal", iter[i], icent) );
      hRefoldNominal[icent][i]->SetLineWidth(2);
      hRefoldNominal[icent][i]->SetLineColor(col[i]);
      hRefoldNominal[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldNominal[icent][i]->Chi2Test(hObsNominal[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
        iterCutoffNominal[icent] = i;

	FixUnfold( hUnfoldNominal[icent][i] );
        EbyECumu cumu(hUnfoldNominal[icent][i]);
        vn2Nominal[icent]  = cumu.GetCumu_vn2();
        vn4Nominal[icent]  = cumu.GetCumu_vn4();
        vn6Nominal[icent]  = cumu.GetCumu_vn6();
	vn8Nominal[icent]  = cumu.GetCumu_vn8();
        gamma1expNominal[icent] = cumu.GetGamma1Exp();

        if( vn4Nominal[icent] == 0 ){
          vn6vn4Nominal[icent] = 0;
          vn8vn4Nominal[icent] = 0;
	}
        else{
          vn6vn4Nominal[icent] = vn6Nominal[icent] / vn4Nominal[icent];
          vn8vn4Nominal[icent] = vn8Nominal[icent] / vn4Nominal[icent];
	}
        if( vn6Nominal[icent] == 0 ) vn8vn6Nominal[icent] = 0;
        else                vn8vn6Nominal[icent] = vn8Nominal[icent] / vn6Nominal[icent];

	vn2Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Nominal->GetBinContent(icent+1) );
        vn4Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Nominal->GetBinContent(icent+1) );
        vn6Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Nominal->GetBinContent(icent+1) );
        vn8Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Nominal->GetBinContent(icent+1) );
        gamma1expNominal_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Nominal->GetBinContent(icent+1) );
        vn6vn4Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Nominal->GetBinContent(icent+1) );
        vn8vn4Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Nominal->GetBinContent(icent+1) );
        vn8vn6Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Nominal->GetBinContent(icent+1) );

        break;
      }
      if( i == NITER - 1 ){
        iterCutoffNominal[icent] = i;

	FixUnfold( hUnfoldNominal[icent][i] );
        EbyECumu cumu(hUnfoldNominal[icent][i]);
        vn2Nominal[icent]  = cumu.GetCumu_vn2();
        vn4Nominal[icent]  = cumu.GetCumu_vn4();
        vn6Nominal[icent]  = cumu.GetCumu_vn6();
        vn8Nominal[icent]  = cumu.GetCumu_vn8();
        gamma1expNominal[icent] = cumu.GetGamma1Exp();

        if( vn4Nominal[icent] == 0 ){
          vn6vn4Nominal[icent] = 0;
          vn8vn4Nominal[icent] = 0;
        }
        else{
          vn6vn4Nominal[icent] = vn6Nominal[icent] / vn4Nominal[icent];
          vn8vn4Nominal[icent] = vn8Nominal[icent] / vn4Nominal[icent];
        }
        if( vn6Nominal[icent] == 0 ) vn8vn6Nominal[icent] = 0;
        else                      vn8vn6Nominal[icent] = vn8Nominal[icent] / vn6Nominal[icent];

	vn2Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Nominal->GetBinContent(icent+1) );
        vn4Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Nominal->GetBinContent(icent+1) );
        vn6Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Nominal->GetBinContent(icent+1) );
        vn8Nominal_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Nominal->GetBinContent(icent+1) );
        gamma1expNominal_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Nominal->GetBinContent(icent+1) );
        vn6vn4Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Nominal->GetBinContent(icent+1) );
        vn8vn4Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Nominal->GetBinContent(icent+1) );
        vn8vn6Nominal_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Nominal->GetBinContent(icent+1) );

        break;
      }
 
    } //-- End iter loop

  } //-- End Cent loop

  //-- Calculate ratios to the nominal case
  for(int icent = 0; icent < NCENT; icent++){
    
    double vn2L        = vn2Loose[icent];
    double vn4L        = vn4Loose[icent];
    double vn6L        = vn6Loose[icent];
    double vn8L        = vn8Loose[icent];
    double gamma1expL  = gamma1expLoose[icent];
    double vn6vn4L     = vn6vn4Loose[icent];
    double vn8vn4L     = vn8vn4Loose[icent];
    double vn8vn6L     = vn8vn6Loose[icent];

    double vn2L_e       = vn2Loose_staterr[icent];
    double vn4L_e       = vn4Loose_staterr[icent];
    double vn6L_e       = vn6Loose_staterr[icent];
    double vn8L_e       = vn8Loose_staterr[icent];
    double gamma1expL_e = gamma1expLoose_staterr[icent];
    double vn6vn4L_e    = vn6vn4Loose_staterr[icent];
    double vn8vn4L_e    = vn8vn4Loose_staterr[icent];
    double vn8vn6L_e    = vn8vn6Loose_staterr[icent];

    double vn2T        = vn2Tight[icent];
    double vn4T        = vn4Tight[icent];
    double vn6T        = vn6Tight[icent];
    double vn8T        = vn8Tight[icent];
    double gamma1expT  = gamma1expTight[icent];
    double vn6vn4T     = vn6vn4Tight[icent];
    double vn8vn4T     = vn8vn4Tight[icent];
    double vn8vn6T     = vn8vn6Tight[icent];

    double vn2T_e       = vn2Tight_staterr[icent];
    double vn4T_e       = vn4Tight_staterr[icent];
    double vn6T_e       = vn6Tight_staterr[icent];
    double vn8T_e       = vn8Tight_staterr[icent];
    double gamma1expT_e = gamma1expTight_staterr[icent];
    double vn6vn4T_e    = vn6vn4Tight_staterr[icent];
    double vn8vn4T_e    = vn8vn4Tight_staterr[icent];
    double vn8vn6T_e    = vn8vn6Tight_staterr[icent];

    double vn2N        = vn2Nominal[icent];
    double vn4N        = vn4Nominal[icent];
    double vn6N        = vn6Nominal[icent];
    double vn8N        = vn8Nominal[icent];
    double gamma1expN  = gamma1expNominal[icent];
    double vn6vn4N     = vn6vn4Nominal[icent];
    double vn8vn4N     = vn8vn4Nominal[icent];
    double vn8vn6N     = vn8vn6Nominal[icent];

    double vn2N_e       = vn2Nominal_staterr[icent];
    double vn4N_e       = vn4Nominal_staterr[icent];
    double vn6N_e       = vn6Nominal_staterr[icent];
    double vn8N_e       = vn8Nominal_staterr[icent];
    double gamma1expN_e = gamma1expNominal_staterr[icent];
    double vn6vn4N_e    = vn6vn4Nominal_staterr[icent];
    double vn8vn4N_e    = vn8vn4Nominal_staterr[icent];
    double vn8vn6N_e    = vn8vn6Nominal_staterr[icent];

    double r2L         = vn2L / vn2N;
    double r4L         = vn4L / vn4N;
    double r6L         = vn6L / vn6N;
    double r8L         = vn8L / vn8N;
    double rgamma1expL = gamma1expL / gamma1expN;
    double r64L        = vn6vn4L / vn6vn4N;
    double r84L        = vn8vn4L / vn8vn4N;
    double r86L        = vn8vn6L / vn8vn6N;

    double r2L_e         = sqrt( pow(vn2L_e/vn2N, 2) + pow(vn2L*vn2N_e/vn2N/vn2N, 2) );
    double r4L_e         = sqrt( pow(vn4L_e/vn4N, 2) + pow(vn4L*vn4N_e/vn4N/vn4N, 2) );
    double r6L_e         = sqrt( pow(vn6L_e/vn6N, 2) + pow(vn6L*vn6N_e/vn6N/vn6N, 2) );
    double r8L_e         = sqrt( pow(vn8L_e/vn8N, 2) + pow(vn8L*vn8N_e/vn8N/vn8N, 2) );
    double rgamma1expL_e = sqrt( pow(gamma1expL_e/gamma1expN, 2) + pow(gamma1expL*gamma1expN_e/gamma1expN/gamma1expN, 2) );
    double r64L_e        = sqrt( pow(vn6vn4L_e/vn6vn4N, 2) + pow(vn6vn4L*vn6vn4N_e/vn6vn4N/vn6vn4N, 2) );
    double r84L_e        = sqrt( pow(vn8vn4L_e/vn8vn4N, 2) + pow(vn8vn4L*vn8vn4N_e/vn8vn4N/vn8vn4N, 2) );
    double r86L_e        = sqrt( pow(vn8vn6L_e/vn8vn6N, 2) + pow(vn8vn6L*vn8vn6N_e/vn8vn6N/vn8vn6N, 2) );

    double r2T         = vn2T / vn2N;
    double r4T         = vn4T / vn4N;
    double r6T         = vn6T / vn6N;
    double r8T         = vn8T / vn8N;
    double rgamma1expT = gamma1expT / gamma1expN;
    double r64T        = vn6vn4T / vn6vn4N;
    double r84T        = vn8vn4T / vn8vn4N;
    double r86T        = vn8vn6T / vn8vn6N;

    double r2T_e         = sqrt( pow(vn2T_e/vn2N, 2) + pow(vn2T*vn2N_e/vn2N/vn2N, 2) );
    double r4T_e         = sqrt( pow(vn4T_e/vn4N, 2) + pow(vn4T*vn4N_e/vn4N/vn4N, 2) );
    double r6T_e         = sqrt( pow(vn6T_e/vn6N, 2) + pow(vn6T*vn6N_e/vn6N/vn6N, 2) );
    double r8T_e         = sqrt( pow(vn8T_e/vn8N, 2) + pow(vn8T*vn8N_e/vn8N/vn8N, 2) );
    double rgamma1expT_e = sqrt( pow(gamma1expT_e/gamma1expN, 2) + pow(gamma1expT*gamma1expN_e/gamma1expN/gamma1expN, 2) );
    double r64T_e        = sqrt( pow(vn6vn4T_e/vn6vn4N, 2) + pow(vn6vn4T*vn6vn4N_e/vn6vn4N/vn6vn4N, 2) );
    double r84T_e        = sqrt( pow(vn8vn4T_e/vn8vn4N, 2) + pow(vn8vn4T*vn8vn4N_e/vn8vn4N/vn8vn4N, 2) );
    double r86T_e        = sqrt( pow(vn8vn6T_e/vn8vn6N, 2) + pow(vn8vn6T*vn8vn6N_e/vn8vn6N/vn8vn6N, 2) );

    vn2Loose_RatiotoNominal[icent]       = r2L;
    vn4Loose_RatiotoNominal[icent]       = r4L;
    vn6Loose_RatiotoNominal[icent]       = r6L;
    vn8Loose_RatiotoNominal[icent]       = r8L;
    gamma1expLoose_RatiotoNominal[icent] = fabs(1.-rgamma1expL);
    vn6vn4Loose_RatiotoNominal[icent]    = r64L;
    vn8vn4Loose_RatiotoNominal[icent]    = r84L;
    vn8vn6Loose_RatiotoNominal[icent]    = r86L;

    vn2Loose_RatiotoNominal_staterr[icent]       = r2L_e;
    vn4Loose_RatiotoNominal_staterr[icent]       = r4L_e;
    vn6Loose_RatiotoNominal_staterr[icent]       = r6L_e;
    vn8Loose_RatiotoNominal_staterr[icent]       = r8L_e;
    gamma1expLoose_RatiotoNominal_staterr[icent] = rgamma1expL_e;
    vn6vn4Loose_RatiotoNominal_staterr[icent]    = r64L_e;
    vn8vn4Loose_RatiotoNominal_staterr[icent]    = r84L_e;
    vn8vn6Loose_RatiotoNominal_staterr[icent]    = r86L_e;

    vn2Tight_RatiotoNominal[icent]       = r2T;
    vn4Tight_RatiotoNominal[icent]       = r4T;
    vn6Tight_RatiotoNominal[icent]       = r6T;
    vn8Tight_RatiotoNominal[icent]       = r8T;
    gamma1expTight_RatiotoNominal[icent] = fabs(1.-rgamma1expT);
    vn6vn4Tight_RatiotoNominal[icent]    = r64T;
    vn8vn4Tight_RatiotoNominal[icent]    = r84T;
    vn8vn6Tight_RatiotoNominal[icent]    = r86T;

    vn2Tight_RatiotoNominal_staterr[icent]       = r2T_e;
    vn4Tight_RatiotoNominal_staterr[icent]       = r4T_e;
    vn6Tight_RatiotoNominal_staterr[icent]       = r6T_e;
    vn8Tight_RatiotoNominal_staterr[icent]       = r8T_e;
    gamma1expTight_RatiotoNominal_staterr[icent] = rgamma1expT_e;
    vn6vn4Tight_RatiotoNominal_staterr[icent]    = r64T_e;
    vn8vn4Tight_RatiotoNominal_staterr[icent]    = r84T_e;
    vn8vn6Tight_RatiotoNominal_staterr[icent]    = r86T_e;

  }

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT]; 
  for(int icent = 0; icent < NCENT; icent++) cErr[icent] = 0;

  //-- Ratios 
  grVn2Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Loose_RatiotoNominal,       cErr, vn2Loose_RatiotoNominal_staterr);
  grVn4Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Loose_RatiotoNominal,       cErr, vn4Loose_RatiotoNominal_staterr);
  grVn6Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Loose_RatiotoNominal,       cErr, vn6Loose_RatiotoNominal_staterr);
  grVn8Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Loose_RatiotoNominal,       cErr, vn8Loose_RatiotoNominal_staterr);
  grGamma1ExpLoose_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expLoose_RatiotoNominal, cErr, gamma1expLoose_RatiotoNominal_staterr);
  grVn6Vn4Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Loose_RatiotoNominal,    cErr, vn6vn4Loose_RatiotoNominal_staterr);
  grVn8Vn4Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Loose_RatiotoNominal,    cErr, vn6vn4Loose_RatiotoNominal_staterr);
  grVn8Vn6Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Loose_RatiotoNominal,    cErr, vn8vn6Loose_RatiotoNominal_staterr);

  grVn2Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Tight_RatiotoNominal,       cErr, vn2Tight_RatiotoNominal_staterr);
  grVn4Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Tight_RatiotoNominal,       cErr, vn4Tight_RatiotoNominal_staterr);
  grVn6Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Tight_RatiotoNominal,       cErr, vn6Tight_RatiotoNominal_staterr);
  grVn8Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Tight_RatiotoNominal,       cErr, vn8Tight_RatiotoNominal_staterr);
  grGamma1ExpTight_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expTight_RatiotoNominal, cErr, gamma1expTight_RatiotoNominal_staterr);
  grVn6Vn4Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Tight_RatiotoNominal,    cErr, vn6vn4Tight_RatiotoNominal_staterr);
  grVn8Vn4Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Tight_RatiotoNominal,    cErr, vn8vn4Tight_RatiotoNominal_staterr);
  grVn8Vn6Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Tight_RatiotoNominal,    cErr, vn8vn6Tight_RatiotoNominal_staterr);

  //-- Cumu Ratio Plots
  TLine * line = new TLine(cent_min[0], 1.0, grVn2Loose_RatiotoNominal->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * cSysCumu = new TCanvas("cSysCumu", "cSysCumu", 1000, 1000);
  cSysCumu->Divide(2,2);

  cSysCumu->cd(1);
  grVn2Loose_RatiotoNominal->Draw("ap");
  grVn2Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(2);
  grVn4Loose_RatiotoNominal->Draw("ap");
  grVn4Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(3);
  grVn6Loose_RatiotoNominal->Draw("ap");
  grVn6Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(4);
  grVn8Loose_RatiotoNominal->Draw("ap");
  grVn8Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");

  formatGraph(grVn2Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), 1, 24, "grVn2Loose_RatiotoNominal");
  formatGraph(grVn2Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), kCyan+2, 24, "grVn2Tight_RatiotoNominal");

  formatGraph(grVn4Loose_RatiotoNominal, "Centrality %", 0.90, 1.10, Form("v_{%i}{4} Ratio", norder_), kSpring+4, 25, "grVn4Loose_RatiotoNominal");
  formatGraph(grVn4Tight_RatiotoNominal, "Centrality %", 0.90, 1.10, Form("v_{%i}{4} Ratio", norder_), kViolet-1, 25, "grVn4Tight_RatiotoNominal");

  formatGraph(grVn6Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 6, 28, "grVn6Loose_RatiotoNominal");
  formatGraph(grVn6Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 4, 28, "grVn6Tight_RatiotoNominal");

  formatGraph(grVn8Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kOrange+7, 27, "grVn8Loose_RatiotoNominal");
  formatGraph(grVn8Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kGray+2, 27, "grVn8Tight_RatiotoNominal");

  TLegend * leg31 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg31->SetFillStyle(0);
  leg31->SetBorderSize(0);
  leg31->AddEntry(grVn2Loose_RatiotoNominal, "loose/nominal", "lp");
  leg31->AddEntry(grVn2Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg32 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg32->SetFillStyle(0);
  leg32->SetBorderSize(0);
  leg32->AddEntry(grVn4Loose_RatiotoNominal, "loose/nominal", "lp");
  leg32->AddEntry(grVn4Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg33 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg33->SetFillStyle(0);
  leg33->SetBorderSize(0);
  leg33->AddEntry(grVn6Loose_RatiotoNominal, "loose/nominal", "lp");
  leg33->AddEntry(grVn6Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg34 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg34->SetFillStyle(0);
  leg34->SetBorderSize(0);
  leg34->AddEntry(grVn8Loose_RatiotoNominal, "loose/nominal", "lp");
  leg34->AddEntry(grVn8Tight_RatiotoNominal, "tight/nominal", "lp");

  cSysCumu->cd(1);
  leg31->Draw("same");
  cSysCumu->cd(2);
  leg32->Draw("same");
  cSysCumu->cd(3);
  leg33->Draw("same");
  cSysCumu->cd(4);
  leg34->Draw("same");

  cSysCumu->Update();
  cSysCumu->SaveAs("../../plots/systematicStudies/SysTkQuality_CumuCent.pdf");

  //-- Gamma1Exp Ratio Plot
  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  grGamma1ExpLoose_RatiotoNominal->Draw("ap");
  grGamma1ExpTight_RatiotoNominal->Draw("psame");
  //line->Draw("same");

  formatGraph(grGamma1ExpLoose_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "|1.-#gamma_{1}^{exp} Ratio|", 2, 28, "grGamma1ExpLoose_RatiotoNominal");
  formatGraph(grGamma1ExpTight_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "|1.-#gamma_{1}^{exp} Ratio|", 1, 28, "grGamma1ExpTight_RatiotoNominal");

  TLegend * leg4 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->AddEntry(grGamma1ExpLoose_RatiotoNominal, "loose/nominal", "lp");
  leg4->AddEntry(grGamma1ExpTight_RatiotoNominal, "tight/nominal", "lp");
  leg4->Draw("same");

  cGamma1Exp->Update();
  cGamma1Exp->SaveAs("../../plots/systematicStudies/SysTkQualityGamma1ExpCent.pdf");

  //-- Cumu Ratio Plot
  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  cCumuRatioSys->cd(1)->SetLeftMargin(0.2);
  grVn6Vn4Loose_RatiotoNominal->GetYaxis()->SetTitleOffset(1.6);
  grVn6Vn4Loose_RatiotoNominal->Draw("ap");
  grVn6Vn4Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4Loose_RatiotoNominal->Draw("ap");
  grVn8Vn4Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6Loose_RatiotoNominal->Draw("ap");
  grVn8Vn6Tight_RatiotoNominal->Draw("psame");
  line->Draw("same");

  formatGraph(grVn6Vn4Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 4, 21, "grVn6Vn4Loose_RatiotoNominal");
  formatGraph(grVn6Vn4Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 2, 21, "grVn6Vn4Tight_RatiotoNominal");

  formatGraph(grVn8Vn4Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), kGreen+2, 34, "grVn8Vn4Loose_RatiotoNominal");
  formatGraph(grVn8Vn4Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), 6, 34, "grVn8Vn4Tight_RatiotoNominal");

  formatGraph(grVn8Vn6Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kViolet-1, 33, "grVn8Vn6Loose_RatiotoNominal");
  formatGraph(grVn8Vn6Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kOrange+7, 33, "grVn8Vn6Tight_RatiotoNominal");

  TLegend * leg51 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg51->SetFillStyle(0);
  leg51->SetBorderSize(0);
  leg51->AddEntry(grVn6Vn4Loose_RatiotoNominal, "loose/nominal", "lp");
  leg51->AddEntry(grVn6Vn4Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg52 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg52->SetFillStyle(0);
  leg52->SetBorderSize(0);
  leg52->AddEntry(grVn8Vn4Loose_RatiotoNominal, "loose/nominal", "lp");
  leg52->AddEntry(grVn8Vn4Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg53 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg53->SetFillStyle(0);
  leg53->SetBorderSize(0);
  leg53->AddEntry(grVn8Vn6Loose_RatiotoNominal, "loose/nominal", "lp");
  leg53->AddEntry(grVn8Vn6Tight_RatiotoNominal, "tight/nominal", "lp");

  cCumuRatioSys->cd(1);
  leg51->Draw("same");
  cCumuRatioSys->cd(2);
  leg52->Draw("same");
  cCumuRatioSys->cd(3);
  leg53->Draw("same");

  cCumuRatioSys->Update();
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/SysTkQuality_CumuRatioCent.pdf");

  //---------------------------------------------------------------------------------------------------- 
  //-- Save plots for smoothing
  TFile * fSave = new TFile("SysTkQuality.root", "recreate");
  fSave->cd();

  grVn2Loose_RatiotoNominal->Write("grVn2Loose_RatiotoNominal");
  grVn2Tight_RatiotoNominal->Write("grVn2Tight_RatiotoNominal");
  grVn4Loose_RatiotoNominal->Write("grVn4Loose_RatiotoNominal");
  grVn4Tight_RatiotoNominal->Write("grVn4Tight_RatiotoNominal");
  grVn6Loose_RatiotoNominal->Write("grVn6Loose_RatiotoNominal");
  grVn6Tight_RatiotoNominal->Write("grVn6Tight_RatiotoNominal");
  grVn8Loose_RatiotoNominal->Write("grVn8Loose_RatiotoNominal");
  grVn8Tight_RatiotoNominal->Write("grVn8Tight_RatiotoNominal");

  grVn6Vn4Loose_RatiotoNominal->Write("grVn6Vn4Loose_RatiotoNominal");
  grVn6Vn4Tight_RatiotoNominal->Write("grVn6Vn4Tight_RatiotoNominal");
  grVn8Vn4Loose_RatiotoNominal->Write("grVn8Vn4Loose_RatiotoNominal");
  grVn8Vn4Tight_RatiotoNominal->Write("grVn8Vn4Tight_RatiotoNominal");
  grVn8Vn6Loose_RatiotoNominal->Write("grVn8Vn6Loose_RatiotoNominal");
  grVn8Vn6Tight_RatiotoNominal->Write("grVn8Vn6Tight_RatiotoNominal");

  grGamma1ExpLoose_RatiotoNominal->Write("grGamma1ExpLoose_RatiotoNominal");
  grGamma1ExpTight_RatiotoNominal->Write("grGamma1ExpTight_RatiotoNominal");

  //-- Save the unfolded distns for when the cutoff is chi2=2.
  for(int icent = 0; icent < NCENT; icent++){
    int i = iterCutoffLoose[icent];
    std::cout<<i<<std::endl;
    hUnfoldLoose[icent][i]->SetLineColor(1);
    hUnfoldLoose[icent][i]->SetMarkerColor(1);
    hUnfoldLoose[icent][i]->Write( Form("hFinalUnfold_LooseSysTkQ_c%i", icent) );

    i = iterCutoffTight[icent];
    std::cout<<i<<std::endl;
    hUnfoldTight[icent][i]->SetLineColor(1);
    hUnfoldTight[icent][i]->SetMarkerColor(1);
    hUnfoldTight[icent][i]->Write( Form("hFinalUnfold_TightSysTkQ_c%i", icent) );
  }

}

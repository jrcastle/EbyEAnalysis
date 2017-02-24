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

void sysPrior(){

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
  double gamma1ExpRatioMin = 0.9;
  double gamma1ExpRatioMax = 1.1;

  TFile * fPriorNom;
  TFile * fPriorThin;
  TFile * fPriorWide;
  TH1D * hPriorNom[NCENT];
  TH1D * hPriorThin[NCENT];
  TH1D * hPriorWide[NCENT];

  //-- Thin Cuts
  TFile * fUnfThin;
  TH1D * hUnfoldThin[NCENT][NITER];
  TH1D * hRefoldThin[NCENT][NITER];

  TFile * fStatThin;
  TH1D * hVarianceOfMean_Vn2_Thin;
  TH1D * hVarianceOfMean_Vn4_Thin;
  TH1D * hVarianceOfMean_Vn6_Thin;
  TH1D * hVarianceOfMean_Vn8_Thin;
  TH1D * hVarianceOfMean_Gamma1Exp_Thin;
  TH1D * hVarianceOfMean_Vn6Vn4_Thin;
  TH1D * hVarianceOfMean_Vn8Vn4_Thin;
  TH1D * hVarianceOfMean_Vn8Vn6_Thin;

  TGraphErrors * grVn2Thin;
  TGraphErrors * grVn4Thin;
  TGraphErrors * grVn6Thin;
  TGraphErrors * grVn8Thin;
  TGraphErrors * grGamma1ExpThin;
  TGraphErrors * grVn6Vn4Thin;
  TGraphErrors * grVn8Vn4Thin;
  TGraphErrors * grVn8Vn6Thin;

  double vn2Thin[NCENT];
  double vn4Thin[NCENT];
  double vn6Thin[NCENT];
  double vn8Thin[NCENT];
  double gamma1expThin[NCENT];
  double vn6vn4Thin[NCENT];
  double vn8vn4Thin[NCENT];
  double vn8vn6Thin[NCENT];

  double vn2Thin_staterr[NCENT];
  double vn4Thin_staterr[NCENT];
  double vn6Thin_staterr[NCENT];
  double vn8Thin_staterr[NCENT];
  double gamma1expThin_staterr[NCENT];
  double vn6vn4Thin_staterr[NCENT];
  double vn8vn4Thin_staterr[NCENT];
  double vn8vn6Thin_staterr[NCENT];

  TGraphErrors * grVn2Thin_RatiotoNominal;
  TGraphErrors * grVn4Thin_RatiotoNominal;
  TGraphErrors * grVn6Thin_RatiotoNominal;
  TGraphErrors * grVn8Thin_RatiotoNominal;
  TGraphErrors * grGamma1ExpThin_RatiotoNominal;
  TGraphErrors * grVn6Vn4Thin_RatiotoNominal;
  TGraphErrors * grVn8Vn4Thin_RatiotoNominal;
  TGraphErrors * grVn8Vn6Thin_RatiotoNominal;

  double vn2Thin_RatiotoNominal[NCENT];
  double vn4Thin_RatiotoNominal[NCENT];
  double vn6Thin_RatiotoNominal[NCENT];
  double vn8Thin_RatiotoNominal[NCENT];
  double gamma1expThin_RatiotoNominal[NCENT];
  double vn6vn4Thin_RatiotoNominal[NCENT];
  double vn8vn4Thin_RatiotoNominal[NCENT];
  double vn8vn6Thin_RatiotoNominal[NCENT];

  double vn2Thin_RatiotoNominal_staterr[NCENT];
  double vn4Thin_RatiotoNominal_staterr[NCENT];
  double vn6Thin_RatiotoNominal_staterr[NCENT];
  double vn8Thin_RatiotoNominal_staterr[NCENT];
  double gamma1expThin_RatiotoNominal_staterr[NCENT];
  double vn6vn4Thin_RatiotoNominal_staterr[NCENT];
  double vn8vn4Thin_RatiotoNominal_staterr[NCENT];
  double vn8vn6Thin_RatiotoNominal_staterr[NCENT];

  int iterCutoffThin[NCENT];

  //-- Wide Cuts
  TFile * fUnfWide;
  TH1D * hUnfoldWide[NCENT][NITER];
  TH1D * hRefoldWide[NCENT][NITER];

  TFile * fStatWide;
  TH1D * hVarianceOfMean_Vn2_Wide;
  TH1D * hVarianceOfMean_Vn4_Wide;
  TH1D * hVarianceOfMean_Vn6_Wide;
  TH1D * hVarianceOfMean_Vn8_Wide;
  TH1D * hVarianceOfMean_Gamma1Exp_Wide;
  TH1D * hVarianceOfMean_Vn6Vn4_Wide;
  TH1D * hVarianceOfMean_Vn8Vn4_Wide;
  TH1D * hVarianceOfMean_Vn8Vn6_Wide;

  TGraphErrors * grVn2Wide;
  TGraphErrors * grVn4Wide;
  TGraphErrors * grVn6Wide;
  TGraphErrors * grVn8Wide;
  TGraphErrors * grGamma1ExpWide;
  TGraphErrors * grVn6Vn4Wide;
  TGraphErrors * grVn8Vn4Wide;
  TGraphErrors * grVn8Vn6Wide;

  double vn2Wide[NCENT];
  double vn4Wide[NCENT];
  double vn6Wide[NCENT];
  double vn8Wide[NCENT];
  double gamma1expWide[NCENT];
  double vn6vn4Wide[NCENT];
  double vn8vn4Wide[NCENT];
  double vn8vn6Wide[NCENT];

  double vn2Wide_staterr[NCENT];
  double vn4Wide_staterr[NCENT];
  double vn6Wide_staterr[NCENT];
  double vn8Wide_staterr[NCENT];
  double gamma1expWide_staterr[NCENT];
  double vn6vn4Wide_staterr[NCENT];
  double vn8vn4Wide_staterr[NCENT];
  double vn8vn6Wide_staterr[NCENT];

  TGraphErrors * grVn2Wide_RatiotoNominal;
  TGraphErrors * grVn4Wide_RatiotoNominal;
  TGraphErrors * grVn6Wide_RatiotoNominal;
  TGraphErrors * grVn8Wide_RatiotoNominal;
  TGraphErrors * grGamma1ExpWide_RatiotoNominal;
  TGraphErrors * grVn6Vn4Wide_RatiotoNominal;
  TGraphErrors * grVn8Vn4Wide_RatiotoNominal;
  TGraphErrors * grVn8Vn6Wide_RatiotoNominal;

  double vn2Wide_RatiotoNominal[NCENT];
  double vn4Wide_RatiotoNominal[NCENT];
  double vn6Wide_RatiotoNominal[NCENT];
  double vn8Wide_RatiotoNominal[NCENT];
  double gamma1expWide_RatiotoNominal[NCENT];
  double vn6vn4Wide_RatiotoNominal[NCENT];
  double vn8vn4Wide_RatiotoNominal[NCENT];
  double vn8vn6Wide_RatiotoNominal[NCENT];

  double vn2Wide_RatiotoNominal_staterr[NCENT];
  double vn4Wide_RatiotoNominal_staterr[NCENT];
  double vn6Wide_RatiotoNominal_staterr[NCENT];
  double vn8Wide_RatiotoNominal_staterr[NCENT];
  double gamma1expWide_RatiotoNominal_staterr[NCENT];
  double vn6vn4Wide_RatiotoNominal_staterr[NCENT];
  double vn8vn4Wide_RatiotoNominal_staterr[NCENT];
  double vn8vn6Wide_RatiotoNominal_staterr[NCENT];

  int iterCutoffWide[NCENT];

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
  fPriorNom  = new TFile("../../AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root");
  fPriorThin = new TFile("newPriors/dataDrivenResponseAndPriors_Thinner.root");
  fPriorWide = new TFile("newPriors/dataDrivenResponseAndPriors_Wider.root");

  fUnfThin = new TFile( Form("../../UnfoldResults/dataResp/data%i_Thinner.root", norder_) );
  fUnfWide = new TFile( Form("../../UnfoldResults/dataResp/data%i_Wider.root", norder_) );

  fAnaNominal   = new TFile( "../../AnalyzerResults/CastleEbyE.root"  );
  fUnfNominal   = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Stat Errors
  fStatThin = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/prior/StatUncertPriorThinner_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Thin       = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn2_PriorThinner" );
  hVarianceOfMean_Vn4_Thin       = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn4_PriorThinner" );
  hVarianceOfMean_Vn6_Thin       = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn6_PriorThinner" );
  hVarianceOfMean_Vn8_Thin       = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn8_PriorThinner" );
  hVarianceOfMean_Gamma1Exp_Thin = (TH1D*) fStatThin->Get( "hVarianceOfMean_Gamma1Exp_PriorThinner" );
  hVarianceOfMean_Vn6Vn4_Thin    = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn6Vn4_PriorThinner" );
  hVarianceOfMean_Vn8Vn4_Thin    = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn8Vn4_PriorThinner" );
  hVarianceOfMean_Vn8Vn6_Thin    = (TH1D*) fStatThin->Get( "hVarianceOfMean_Vn8Vn6_PriorThinner" );

  fStatWide = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/prior/StatUncertPriorWider_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Wide       = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn2_PriorWider" );
  hVarianceOfMean_Vn4_Wide       = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn4_PriorWider" );
  hVarianceOfMean_Vn6_Wide       = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn6_PriorWider" );
  hVarianceOfMean_Vn8_Wide       = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn8_PriorWider" );
  hVarianceOfMean_Gamma1Exp_Wide = (TH1D*) fStatWide->Get( "hVarianceOfMean_Gamma1Exp_PriorWider" );
  hVarianceOfMean_Vn6Vn4_Wide    = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn6Vn4_PriorWider" );
  hVarianceOfMean_Vn8Vn4_Wide    = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn8Vn4_PriorWider" );
  hVarianceOfMean_Vn8Vn6_Wide    = (TH1D*) fStatWide->Get( "hVarianceOfMean_Vn8Vn6_PriorWider" );

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

    hObsNominal[icent] = (TH1D*) fAnaNominal->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsNominal[icent]->SetName( Form("hVnFull_c%i_nominal", icent) );

    hPriorNom[icent] = (TH1D*) fPriorNom->Get( Form("hPrior_c%i", icent) );
    hPriorNom[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hPriorNom[icent]->GetYaxis()->SetTitle( "Events" );
    hPriorNom[icent]->SetLineColor(1);
    hPriorNom[icent]->SetMarkerColor(1);
    hPriorNom[icent]->SetMarkerStyle(20);

    hPriorThin[icent] = (TH1D*) fPriorThin->Get( Form("hPrior_c%i", icent) );
    hPriorThin[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hPriorThin[icent]->GetYaxis()->SetTitle( "Events" );
    hPriorThin[icent]->SetLineColor(2);
    hPriorThin[icent]->SetMarkerColor(2);
    hPriorThin[icent]->SetMarkerStyle(21);

    hPriorWide[icent] = (TH1D*) fPriorWide->Get( Form("hPrior_c%i", icent) );
    hPriorWide[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hPriorWide[icent]->GetYaxis()->SetTitle( "Events" );
    hPriorWide[icent]->SetLineColor(4);
    hPriorWide[icent]->SetMarkerColor(4);
    hPriorWide[icent]->SetMarkerStyle(22);

    //-- -------------------- Thin --------------------

    //-- Iter loop
    iterCutoffThin[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldThin[icent][i] = (TH1D*) fUnfThin->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldThin[icent][i]->SetName( Form("hreco%i_c%i_thin", iter[i], icent) );
      hUnfoldThin[icent][i]->SetLineColor(col[i]);
      hUnfoldThin[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldThin[icent][i] = (TH1D*) fUnfThin->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldThin[icent][i]->SetName( Form("hrefold%i_c%i_thin", iter[i], icent) );
      hRefoldThin[icent][i]->SetLineWidth(2);
      hRefoldThin[icent][i]->SetLineColor(col[i]);
      hRefoldThin[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldThin[icent][i]->Chi2Test(hObsNominal[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
	iterCutoffThin[icent] = i;

	FixUnfold( hUnfoldThin[icent][i] );
	EbyECumu cumu(hUnfoldThin[icent][i]);
	vn2Thin[icent]  = cumu.GetCumu_vn2();
	vn4Thin[icent]  = cumu.GetCumu_vn4();
	vn6Thin[icent]  = cumu.GetCumu_vn6();
	vn8Thin[icent]  = cumu.GetCumu_vn8();
	gamma1expThin[icent] = cumu.GetGamma1Exp();

	if( vn4Thin[icent] == 0 ){
	  vn6vn4Thin[icent] = 0;
	  vn8vn4Thin[icent] = 0;
	}
	else{
	  vn6vn4Thin[icent] = vn6Thin[icent] / vn4Thin[icent];
	  vn8vn4Thin[icent] = vn8Thin[icent] / vn4Thin[icent];
	}
	if( vn6Thin[icent] == 0 ) vn8vn6Thin[icent] = 0;
	else                      vn8vn6Thin[icent] = vn8Thin[icent] / vn6Thin[icent];

	vn2Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Thin->GetBinContent(icent+1) );
	vn4Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Thin->GetBinContent(icent+1) );
	vn6Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Thin->GetBinContent(icent+1) );
	vn8Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Thin->GetBinContent(icent+1) );
	gamma1expThin_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Thin->GetBinContent(icent+1) );
	vn6vn4Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Thin->GetBinContent(icent+1) );
	vn8vn4Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Thin->GetBinContent(icent+1) );
	vn8vn6Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Thin->GetBinContent(icent+1) );

	break;
      }
      if( i == NITER - 1 ){
	iterCutoffThin[icent] = i;

	FixUnfold( hUnfoldThin[icent][i] );
        EbyECumu cumu(hUnfoldThin[icent][i]);
        vn2Thin[icent]  = cumu.GetCumu_vn2();
        vn4Thin[icent]  = cumu.GetCumu_vn4();
        vn6Thin[icent]  = cumu.GetCumu_vn6();
	vn8Thin[icent]  = cumu.GetCumu_vn8();
        gamma1expThin[icent] = cumu.GetGamma1Exp();

        if( vn4Thin[icent] == 0 ){
          vn6vn4Thin[icent] = 0;
          vn8vn4Thin[icent] = 0;
	}
        else{
          vn6vn4Thin[icent] = vn6Thin[icent] / vn4Thin[icent];
          vn8vn4Thin[icent] = vn8Thin[icent] / vn4Thin[icent];
	}
        if( vn6Thin[icent] == 0 ) vn8vn6Thin[icent] = 0;
        else                      vn8vn6Thin[icent] = vn8Thin[icent] / vn6Thin[icent];

	vn2Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Thin->GetBinContent(icent+1) );
        vn4Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Thin->GetBinContent(icent+1) );
        vn6Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Thin->GetBinContent(icent+1) );
        vn8Thin_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Thin->GetBinContent(icent+1) );
        gamma1expThin_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Thin->GetBinContent(icent+1) );
        vn6vn4Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Thin->GetBinContent(icent+1) );
        vn8vn4Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Thin->GetBinContent(icent+1) );
        vn8vn6Thin_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Thin->GetBinContent(icent+1) );

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Wide --------------------

    //-- Iter loop
    iterCutoffWide[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldWide[icent][i] = (TH1D*) fUnfWide->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldWide[icent][i]->SetName( Form("hreco%i_c%i_wide", iter[i], icent) );
      hUnfoldWide[icent][i]->SetLineColor(col[i]);
      hUnfoldWide[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldWide[icent][i] = (TH1D*) fUnfWide->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldWide[icent][i]->SetName( Form("hrefold%i_c%i_wide", iter[i], icent) );
      hRefoldWide[icent][i]->SetLineWidth(2);
      hRefoldWide[icent][i]->SetLineColor(col[i]);
      hRefoldWide[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldWide[icent][i]->Chi2Test(hObsNominal[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
        iterCutoffWide[icent] = i;

	FixUnfold( hUnfoldWide[icent][i] );
        EbyECumu cumu(hUnfoldWide[icent][i]);
        vn2Wide[icent]  = cumu.GetCumu_vn2();
        vn4Wide[icent]  = cumu.GetCumu_vn4();
        vn6Wide[icent]  = cumu.GetCumu_vn6();
	vn8Wide[icent]  = cumu.GetCumu_vn8();
        gamma1expWide[icent] = cumu.GetGamma1Exp();

        if( vn4Wide[icent] == 0 ){
          vn6vn4Wide[icent] = 0;
          vn8vn4Wide[icent] = 0;
	}
        else{
          vn6vn4Wide[icent] = vn6Wide[icent] / vn4Wide[icent];
          vn8vn4Wide[icent] = vn8Wide[icent] / vn4Wide[icent];
	}
        if( vn6Wide[icent] == 0 ) vn8vn6Wide[icent] = 0;
        else                      vn8vn6Wide[icent] = vn8Wide[icent] / vn6Wide[icent];

	vn2Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Wide->GetBinContent(icent+1) );
        vn4Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Wide->GetBinContent(icent+1) );
        vn6Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Wide->GetBinContent(icent+1) );
        vn8Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Wide->GetBinContent(icent+1) );
        gamma1expWide_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Wide->GetBinContent(icent+1) );
        vn6vn4Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Wide->GetBinContent(icent+1) );
        vn8vn4Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Wide->GetBinContent(icent+1) );
        vn8vn6Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Wide->GetBinContent(icent+1) );

        break;
      }
      if( i == NITER - 1 ){
        iterCutoffWide[icent] = i;

	FixUnfold( hUnfoldWide[icent][i] );
        EbyECumu cumu(hUnfoldWide[icent][i]);
        vn2Wide[icent]  = cumu.GetCumu_vn2();
        vn4Wide[icent]  = cumu.GetCumu_vn4();
        vn6Wide[icent]  = cumu.GetCumu_vn6();
        vn8Wide[icent]  = cumu.GetCumu_vn8();
        gamma1expWide[icent] = cumu.GetGamma1Exp();

        if( vn4Wide[icent] == 0 ){
          vn6vn4Wide[icent] = 0;
          vn8vn4Wide[icent] = 0;
        }
        else{
          vn6vn4Wide[icent] = vn6Wide[icent] / vn4Wide[icent];
          vn8vn4Wide[icent] = vn8Wide[icent] / vn4Wide[icent];
        }
        if( vn6Wide[icent] == 0 ) vn8vn6Wide[icent] = 0;
        else                      vn8vn6Wide[icent] = vn8Wide[icent] / vn6Wide[icent];

	vn2Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn2_Wide->GetBinContent(icent+1) );
        vn4Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn4_Wide->GetBinContent(icent+1) );
        vn6Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn6_Wide->GetBinContent(icent+1) );
        vn8Wide_staterr[icent]       = sqrt( hVarianceOfMean_Vn8_Wide->GetBinContent(icent+1) );
        gamma1expWide_staterr[icent] = sqrt( hVarianceOfMean_Gamma1Exp_Wide->GetBinContent(icent+1) );
        vn6vn4Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn6Vn4_Wide->GetBinContent(icent+1) );
        vn8vn4Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn4_Wide->GetBinContent(icent+1) );
        vn8vn6Wide_staterr[icent]    = sqrt( hVarianceOfMean_Vn8Vn6_Wide->GetBinContent(icent+1) );

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Nominal --------------------

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
    
    double vn2T        = vn2Thin[icent];
    double vn4T        = vn4Thin[icent];
    double vn6T        = vn6Thin[icent];
    double vn8T        = vn8Thin[icent];
    double gamma1expT  = gamma1expThin[icent];
    double vn6vn4T     = vn6vn4Thin[icent];
    double vn8vn4T     = vn8vn4Thin[icent];
    double vn8vn6T     = vn8vn6Thin[icent];

    double vn2T_e       = vn2Thin_staterr[icent];
    double vn4T_e       = vn4Thin_staterr[icent];
    double vn6T_e       = vn6Thin_staterr[icent];
    double vn8T_e       = vn8Thin_staterr[icent];
    double gamma1expT_e = gamma1expThin_staterr[icent];
    double vn6vn4T_e    = vn6vn4Thin_staterr[icent];
    double vn8vn4T_e    = vn8vn4Thin_staterr[icent];
    double vn8vn6T_e    = vn8vn6Thin_staterr[icent];

    double vn2W        = vn2Wide[icent];
    double vn4W        = vn4Wide[icent];
    double vn6W        = vn6Wide[icent];
    double vn8W        = vn8Wide[icent];
    double gamma1expW  = gamma1expWide[icent];
    double vn6vn4W     = vn6vn4Wide[icent];
    double vn8vn4W     = vn8vn4Wide[icent];
    double vn8vn6W     = vn8vn6Wide[icent];

    double vn2W_e       = vn2Wide_staterr[icent];
    double vn4W_e       = vn4Wide_staterr[icent];
    double vn6W_e       = vn6Wide_staterr[icent];
    double vn8W_e       = vn8Wide_staterr[icent];
    double gamma1expW_e = gamma1expWide_staterr[icent];
    double vn6vn4W_e    = vn6vn4Wide_staterr[icent];
    double vn8vn4W_e    = vn8vn4Wide_staterr[icent];
    double vn8vn6W_e    = vn8vn6Wide_staterr[icent];

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

    double r2W         = vn2W / vn2N;
    double r4W         = vn4W / vn4N;
    double r6W         = vn6W / vn6N;
    double r8W         = vn8W / vn8N;
    double rgamma1expW = gamma1expW / gamma1expN;
    double r64W        = vn6vn4W / vn6vn4N;
    double r84W        = vn8vn4W / vn8vn4N;
    double r86W        = vn8vn6W / vn8vn6N;

    double r2W_e         = sqrt( pow(vn2W_e/vn2N, 2) + pow(vn2W*vn2N_e/vn2N/vn2N, 2) );
    double r4W_e         = sqrt( pow(vn4W_e/vn4N, 2) + pow(vn4W*vn4N_e/vn4N/vn4N, 2) );
    double r6W_e         = sqrt( pow(vn6W_e/vn6N, 2) + pow(vn6W*vn6N_e/vn6N/vn6N, 2) );
    double r8W_e         = sqrt( pow(vn8W_e/vn8N, 2) + pow(vn8W*vn8N_e/vn8N/vn8N, 2) );
    double rgamma1expW_e = sqrt( pow(gamma1expW_e/gamma1expN, 2) + pow(gamma1expW*gamma1expN_e/gamma1expN/gamma1expN, 2) );
    double r64W_e        = sqrt( pow(vn6vn4W_e/vn6vn4N, 2) + pow(vn6vn4W*vn6vn4N_e/vn6vn4N/vn6vn4N, 2) );
    double r84W_e        = sqrt( pow(vn8vn4W_e/vn8vn4N, 2) + pow(vn8vn4W*vn8vn4N_e/vn8vn4N/vn8vn4N, 2) );
    double r86W_e        = sqrt( pow(vn8vn6W_e/vn8vn6N, 2) + pow(vn8vn6W*vn8vn6N_e/vn8vn6N/vn8vn6N, 2) );

    vn2Thin_RatiotoNominal[icent]       = r2T;
    vn4Thin_RatiotoNominal[icent]       = r4T;
    vn6Thin_RatiotoNominal[icent]       = r6T;
    vn8Thin_RatiotoNominal[icent]       = r8T;
    gamma1expThin_RatiotoNominal[icent] = rgamma1expT;
    vn6vn4Thin_RatiotoNominal[icent]    = r64T;
    vn8vn4Thin_RatiotoNominal[icent]    = r84T;
    vn8vn6Thin_RatiotoNominal[icent]    = r86T;

    vn2Thin_RatiotoNominal_staterr[icent]       = r2T_e;
    vn4Thin_RatiotoNominal_staterr[icent]       = r4T_e;
    vn6Thin_RatiotoNominal_staterr[icent]       = r6T_e;
    vn8Thin_RatiotoNominal_staterr[icent]       = r8T_e;
    gamma1expThin_RatiotoNominal_staterr[icent] = rgamma1expT_e;
    vn6vn4Thin_RatiotoNominal_staterr[icent]    = r64T_e;
    vn8vn4Thin_RatiotoNominal_staterr[icent]    = r84T_e;
    vn8vn6Thin_RatiotoNominal_staterr[icent]    = r86T_e;

    vn2Wide_RatiotoNominal[icent]       = r2W;
    vn4Wide_RatiotoNominal[icent]       = r4W;
    vn6Wide_RatiotoNominal[icent]       = r6W;
    vn8Wide_RatiotoNominal[icent]       = r8W;
    gamma1expWide_RatiotoNominal[icent] = rgamma1expW;
    vn6vn4Wide_RatiotoNominal[icent]    = r64W;
    vn8vn4Wide_RatiotoNominal[icent]    = r84W;
    vn8vn6Wide_RatiotoNominal[icent]    = r86W;

    vn2Wide_RatiotoNominal_staterr[icent]       = r2W_e;
    vn4Wide_RatiotoNominal_staterr[icent]       = r4W_e;
    vn6Wide_RatiotoNominal_staterr[icent]       = r6W_e;
    vn8Wide_RatiotoNominal_staterr[icent]       = r8W_e;
    gamma1expWide_RatiotoNominal_staterr[icent] = rgamma1expW_e;
    vn6vn4Wide_RatiotoNominal_staterr[icent]    = r64W_e;
    vn8vn4Wide_RatiotoNominal_staterr[icent]    = r84W_e;
    vn8vn6Wide_RatiotoNominal_staterr[icent]    = r86W_e;

  }

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT]; 
  for(int icent = 0; icent < NCENT; icent++) cErr[icent] = 0;

  //-- Ratios 
  grVn2Thin_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Thin_RatiotoNominal,       cErr, vn2Thin_RatiotoNominal_staterr);
  grVn4Thin_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Thin_RatiotoNominal,       cErr, vn4Thin_RatiotoNominal_staterr);
  grVn6Thin_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Thin_RatiotoNominal,       cErr, vn6Thin_RatiotoNominal_staterr);
  grVn8Thin_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Thin_RatiotoNominal,       cErr, vn8Thin_RatiotoNominal_staterr);
  grGamma1ExpThin_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expThin_RatiotoNominal, cErr, gamma1expThin_RatiotoNominal_staterr);
  grVn6Vn4Thin_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Thin_RatiotoNominal,    cErr, vn6vn4Thin_RatiotoNominal_staterr);
  grVn8Vn4Thin_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Thin_RatiotoNominal,    cErr, vn6vn4Thin_RatiotoNominal_staterr);
  grVn8Vn6Thin_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Thin_RatiotoNominal,    cErr, vn8vn6Thin_RatiotoNominal_staterr);

  grVn2Wide_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Wide_RatiotoNominal,       cErr, vn2Wide_RatiotoNominal_staterr);
  grVn4Wide_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Wide_RatiotoNominal,       cErr, vn4Wide_RatiotoNominal_staterr);
  grVn6Wide_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Wide_RatiotoNominal,       cErr, vn6Wide_RatiotoNominal_staterr);
  grVn8Wide_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Wide_RatiotoNominal,       cErr, vn8Wide_RatiotoNominal_staterr);
  grGamma1ExpWide_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expWide_RatiotoNominal, cErr, gamma1expWide_RatiotoNominal_staterr);
  grVn6Vn4Wide_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Wide_RatiotoNominal,    cErr, vn6vn4Wide_RatiotoNominal_staterr);
  grVn8Vn4Wide_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Wide_RatiotoNominal,    cErr, vn8vn4Wide_RatiotoNominal_staterr);
  grVn8Vn6Wide_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Wide_RatiotoNominal,    cErr, vn8vn6Wide_RatiotoNominal_staterr);

  //-- Cumu Ratio Plots
  TLine * line = new TLine(cent_min[0], 1.0, grVn2Thin_RatiotoNominal->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * cSysCumu = new TCanvas("cSysCumu", "cSysCumu", 1000, 1000);
  cSysCumu->Divide(2,2);

  cSysCumu->cd(1);
  grVn2Thin_RatiotoNominal->Draw("ap");
  grVn2Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(2);
  grVn4Thin_RatiotoNominal->Draw("ap");
  grVn4Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(3);
  grVn6Thin_RatiotoNominal->Draw("ap");
  grVn6Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(4);
  grVn8Thin_RatiotoNominal->Draw("ap");
  grVn8Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  formatGraph(grVn2Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), 1, 24, "grVn2Thin_RatiotoNominal");
  formatGraph(grVn2Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), kCyan+2, 24, "grVn2Wide_RatiotoNominal");

  formatGraph(grVn4Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kSpring+4, 25, "grVn4Thin_RatiotoNominal");
  formatGraph(grVn4Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kViolet-1, 25, "grVn4Wide_RatiotoNominal");

  formatGraph(grVn6Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 6, 28, "grVn6Thin_RatiotoNominal");
  formatGraph(grVn6Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 4, 28, "grVn6Wide_RatiotoNominal");

  formatGraph(grVn8Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kOrange+7, 27, "grVn8Thin_RatiotoNominal");
  formatGraph(grVn8Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kGray+2, 27, "grVn8Wide_RatiotoNominal");

  TLegend * leg31 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg31->SetFillStyle(0);
  leg31->SetBorderSize(0);
  leg31->AddEntry(grVn2Thin_RatiotoNominal, "thin/nominal", "lp");
  leg31->AddEntry(grVn2Wide_RatiotoNominal, "wide/nominal", "lp");

  TLegend * leg32 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg32->SetFillStyle(0);
  leg32->SetBorderSize(0);
  leg32->AddEntry(grVn4Thin_RatiotoNominal, "thin/nominal", "lp");
  leg32->AddEntry(grVn4Wide_RatiotoNominal, "wide/nominal", "lp");

  TLegend * leg33 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg33->SetFillStyle(0);
  leg33->SetBorderSize(0);
  leg33->AddEntry(grVn6Thin_RatiotoNominal, "thin/nominal", "lp");
  leg33->AddEntry(grVn6Wide_RatiotoNominal, "wide/nominal", "lp");

  TLegend * leg34 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg34->SetFillStyle(0);
  leg34->SetBorderSize(0);
  leg34->AddEntry(grVn8Thin_RatiotoNominal, "thin/nominal", "lp");
  leg34->AddEntry(grVn8Wide_RatiotoNominal, "wide/nominal", "lp");

  cSysCumu->cd(1);
  leg31->Draw("same");
  cSysCumu->cd(2);
  leg32->Draw("same");
  cSysCumu->cd(3);
  leg33->Draw("same");
  cSysCumu->cd(4);
  leg34->Draw("same");

  cSysCumu->Update();
  cSysCumu->SaveAs("../../plots/systematicStudies/SysPrior_CumuCent.pdf");

  //-- Gamma1Exp Ratio Plot
  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  grGamma1ExpThin_RatiotoNominal->Draw("ap");
  grGamma1ExpWide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  formatGraph(grGamma1ExpThin_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 2, 28, "grGamma1ExpThin_RatiotoNominal");
  formatGraph(grGamma1ExpWide_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 1, 28, "grGamma1ExpWide_RatiotoNominal");

  TLegend * leg4 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->AddEntry(grGamma1ExpThin_RatiotoNominal, "thin/nominal", "lp");
  leg4->AddEntry(grGamma1ExpWide_RatiotoNominal, "wide/nominal", "lp");
  leg4->Draw("same");

  cGamma1Exp->Update();
  cGamma1Exp->SaveAs("../../plots/systematicStudies/SysPriorGamma1ExpCent.pdf");

  //-- Cumu Ratio Plot
  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  cCumuRatioSys->cd(1)->SetLeftMargin(0.2);
  grVn6Vn4Thin_RatiotoNominal->GetYaxis()->SetTitleOffset(1.6);
  grVn6Vn4Thin_RatiotoNominal->Draw("ap");
  grVn6Vn4Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4Thin_RatiotoNominal->Draw("ap");
  grVn8Vn4Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6Thin_RatiotoNominal->Draw("ap");
  grVn8Vn6Wide_RatiotoNominal->Draw("psame");
  line->Draw("same");

  formatGraph(grVn6Vn4Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 4, 21, "grVn6Vn4Thin_RatiotoNominal");
  formatGraph(grVn6Vn4Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 2, 21, "grVn6Vn4Wide_RatiotoNominal");

  formatGraph(grVn8Vn4Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), kGreen+2, 34, "grVn8Vn4Thin_RatiotoNominal");
  formatGraph(grVn8Vn4Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), 6, 34, "grVn8Vn4Wide_RatiotoNominal");

  formatGraph(grVn8Vn6Thin_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kViolet-1, 33, "grVn8Vn6Thin_RatiotoNominal");
  formatGraph(grVn8Vn6Wide_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kOrange+7, 33, "grVn8Vn6Wide_RatiotoNominal");

  TLegend * leg51 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg51->SetFillStyle(0);
  leg51->SetBorderSize(0);
  leg51->AddEntry(grVn6Vn4Thin_RatiotoNominal, "thin/nominal", "lp");
  leg51->AddEntry(grVn6Vn4Wide_RatiotoNominal, "wide/nominal", "lp");

  TLegend * leg52 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg52->SetFillStyle(0);
  leg52->SetBorderSize(0);
  leg52->AddEntry(grVn8Vn4Thin_RatiotoNominal, "thin/nominal", "lp");
  leg52->AddEntry(grVn8Vn4Wide_RatiotoNominal, "wide/nominal", "lp");

  TLegend * leg53 = new TLegend(0.24, 0.20, 0.76, 0.45);
  leg53->SetFillStyle(0);
  leg53->SetBorderSize(0);
  leg53->AddEntry(grVn8Vn6Thin_RatiotoNominal, "thin/nominal", "lp");
  leg53->AddEntry(grVn8Vn6Wide_RatiotoNominal, "wide/nominal", "lp");

  cCumuRatioSys->cd(1);
  leg51->Draw("same");
  cCumuRatioSys->cd(2);
  leg52->Draw("same");
  cCumuRatioSys->cd(3);
  leg53->Draw("same");

  cCumuRatioSys->Update();
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/SysPrior_CumuRatioCent.pdf");

  //-- Draw Priors
  TCanvas * cPrior[NCENT];
  TLegend * legPrior = new TLegend(0.60, 0.72, 0.99, 0.92);
  legInit( legPrior );
  legPrior->AddEntry(hPriorNom[0],  "Default Prior", "lp");
  legPrior->AddEntry(hPriorThin[0], "Thin Prior", "lp");
  legPrior->AddEntry(hPriorWide[0], "Wide Prior", "lp");

  for(int icent = 0; icent < NCENT; icent++){
    if( icent != 5 ) continue;
    cPrior[icent] = new TCanvas( Form("cPrior_c%i", icent), Form("cPrior_c%i", icent), 500, 500);
    cPrior[icent]->cd();
    cPrior[icent]->SetLogy();
    hPriorThin[icent]->Draw();
    hPriorWide[icent]->Draw("same");
    hPriorNom[icent]->Draw("same");
    legPrior->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%") );
    cPrior[icent]->SaveAs( Form("../../plots/systematicStudies/SysPriorComp_c%i.pdf", icent) );
  }

  //---------------------------------------------------------------------------------------------------- 
  //-- Save plots for smoothing
  TFile * fSave = new TFile("SysPrior.root", "recreate");
  fSave->cd();

  grVn2Thin_RatiotoNominal->Write("grVn2Thin_RatiotoNominal");
  grVn2Wide_RatiotoNominal->Write("grVn2Wide_RatiotoNominal");
  grVn4Thin_RatiotoNominal->Write("grVn4Thin_RatiotoNominal");
  grVn4Wide_RatiotoNominal->Write("grVn4Wide_RatiotoNominal");
  grVn6Thin_RatiotoNominal->Write("grVn6Thin_RatiotoNominal");
  grVn6Wide_RatiotoNominal->Write("grVn6Wide_RatiotoNominal");
  grVn8Thin_RatiotoNominal->Write("grVn8Thin_RatiotoNominal");
  grVn8Wide_RatiotoNominal->Write("grVn8Wide_RatiotoNominal");

  grVn6Vn4Thin_RatiotoNominal->Write("grVn6Vn4Thin_RatiotoNominal");
  grVn6Vn4Wide_RatiotoNominal->Write("grVn6Vn4Wide_RatiotoNominal");
  grVn8Vn4Thin_RatiotoNominal->Write("grVn8Vn4Thin_RatiotoNominal");
  grVn8Vn4Wide_RatiotoNominal->Write("grVn8Vn4Wide_RatiotoNominal");
  grVn8Vn6Thin_RatiotoNominal->Write("grVn8Vn6Thin_RatiotoNominal");
  grVn8Vn6Wide_RatiotoNominal->Write("grVn8Vn6Wide_RatiotoNominal");

  grGamma1ExpThin_RatiotoNominal->Write("grGamma1ExpThin_RatiotoNominal");
  grGamma1ExpWide_RatiotoNominal->Write("grGamma1ExpWide_RatiotoNominal");

  //-- Save the unfolded distns for when the cutoff is chi2=2.
  for(int icent = 0; icent < NCENT; icent++){
    int i = iterCutoffThin[icent];
    std::cout<<i<<std::endl;
    hUnfoldThin[icent][i]->SetLineColor(1);
    hUnfoldThin[icent][i]->SetMarkerColor(1);
    hUnfoldThin[icent][i]->Write( Form("hFinalUnfold_ThinSysPrior_c%i", icent) );

    i = iterCutoffWide[icent];
    std::cout<<i<<std::endl;
    hUnfoldWide[icent][i]->SetLineColor(1);
    hUnfoldWide[icent][i]->SetMarkerColor(1);
    hUnfoldWide[icent][i]->Write( Form("hFinalUnfold_WideSysPrior_c%i", icent) );
  }

}

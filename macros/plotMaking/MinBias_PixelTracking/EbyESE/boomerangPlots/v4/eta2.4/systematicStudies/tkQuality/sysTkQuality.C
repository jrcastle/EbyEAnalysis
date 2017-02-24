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

  int norder_ = 4;

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
  double gamma1ExpRatioMin = 0.8;
  double gamma1ExpRatioMax = 1.2;

  //-- Vertex < 3.0
  TFile * fAnaLoose;
  TFile * fUnfLoose;
  TH1D * hObsLoose[NCENT];
  TH1D * hUnfoldLoose[NCENT][NITER];
  TH1D * hRefoldLoose[NCENT][NITER];

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

  bool iterCutoffLoose[NCENT];

  //-- 3.0 < Vertex < 15.0
  TFile * fAnaTight;
  TFile * fUnfTight;
  TH1D * hObsTight[NCENT];
  TH1D * hUnfoldTight[NCENT][NITER];
  TH1D * hRefoldTight[NCENT][NITER];

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

  bool iterCutoffTight[NCENT];

  //-- Vertex < 15.0
  TFile * fAnaNominal;
  TFile * fUnfNominal;
  TH1D * hObsNominal[NCENT];
  TH1D * hUnfoldNominal[NCENT][NITER];
  TH1D * hRefoldNominal[NCENT][NITER];

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

  bool iterCutoffNominal[NCENT];

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

  fRelErr   = new TFile("relErrorTkQuality.root", "recreate");
  fRelErr->cd();
  relErrVn2        = new TH1D("relErrVn2",       "relErrVn2",       NCENT, centbinsDefault);
  relErrVn4        = new TH1D("relErrVn4",       "relErrVn4",       NCENT, centbinsDefault);
  relErrVn6        = new TH1D("relErrVn6",       "relErrVn6",       NCENT, centbinsDefault);
  relErrVn8        = new TH1D("relErrVn8",       "relErrVn8",       NCENT, centbinsDefault);
  relErrGamma1Exp  = new TH1D("relErrGamma1Exp", "relErrGamma1Exp", NCENT, centbinsDefault);
  relErrVn6Vn4     = new TH1D("relErrVn6Vn4",    "relErrVn6Vn4",    NCENT, centbinsDefault);
  relErrVn8Vn4     = new TH1D("relErrVn8Vn4",    "relErrVn8Vn4",    NCENT, centbinsDefault);
  relErrVn8Vn6     = new TH1D("relErrVn8Vn6",    "relErrVn8Vn6",    NCENT, centbinsDefault);

  for(int icent = 0; icent < NCENT; icent++){

    //-- -------------------- Vertex < 3.0 --------------------
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

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- 3.0 < Vertex < 15.0 --------------------
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

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Vertex < 15.0 --------------------
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

    double vn2T        = vn2Tight[icent];
    double vn4T        = vn4Tight[icent];
    double vn6T        = vn6Tight[icent];
    double vn8T        = vn8Tight[icent];
    double gamma1expT  = gamma1expTight[icent];
    double vn6vn4T     = vn6vn4Tight[icent];
    double vn8vn4T     = vn8vn4Tight[icent];
    double vn8vn6T     = vn8vn6Tight[icent];

    double vn2N        = vn2Nominal[icent];
    double vn4N        = vn4Nominal[icent];
    double vn6N        = vn6Nominal[icent];
    double vn8N        = vn8Nominal[icent];
    double gamma1expN  = gamma1expNominal[icent];
    double vn6vn4N     = vn6vn4Nominal[icent];
    double vn8vn4N     = vn8vn4Nominal[icent];
    double vn8vn6N     = vn8vn6Nominal[icent];

    double r2L         = vn2L / vn2N;
    double r4L         = vn4L / vn4N;
    double r6L         = vn6L / vn6N;
    double r8L         = vn8L / vn8N;
    double rgamma1expL = gamma1expL / gamma1expN;
    double r64L        = vn6vn4L / vn6vn4N;
    double r84L        = vn8vn4L / vn8vn4N;
    double r86L        = vn8vn6L / vn8vn6N;

    double r2T         = vn2T / vn2N;
    double r4T         = vn4T / vn4N;
    double r6T         = vn6T / vn6N;
    double r8T         = vn8T / vn8N;
    double rgamma1expT = gamma1expT / gamma1expN;
    double r64T        = vn6vn4T / vn6vn4N;
    double r84T        = vn8vn4T / vn8vn4N;
    double r86T        = vn8vn6T / vn8vn6N;

    vn2Loose_RatiotoNominal[icent]       = r2L;
    vn4Loose_RatiotoNominal[icent]       = r4L;
    vn6Loose_RatiotoNominal[icent]       = r6L;
    vn8Loose_RatiotoNominal[icent]       = r8L;
    gamma1expLoose_RatiotoNominal[icent] = rgamma1expL;
    vn6vn4Loose_RatiotoNominal[icent]    = r64L;
    vn8vn4Loose_RatiotoNominal[icent]    = r84L;
    vn8vn6Loose_RatiotoNominal[icent]    = r86L;

    vn2Tight_RatiotoNominal[icent]       = r2T;
    vn4Tight_RatiotoNominal[icent]       = r4T;
    vn6Tight_RatiotoNominal[icent]       = r6T;
    vn8Tight_RatiotoNominal[icent]       = r8T;
    gamma1expTight_RatiotoNominal[icent] = rgamma1expT;
    vn6vn4Tight_RatiotoNominal[icent]    = r64T;
    vn8vn4Tight_RatiotoNominal[icent]    = r84T;
    vn8vn6Tight_RatiotoNominal[icent]    = r86T;

  }

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT]; 
  for(int icent = 0; icent < NCENT; icent++) cErr[icent] = 0;

  //-- Ratios 
  grVn2Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Loose_RatiotoNominal,       cErr, cErr);
  grVn4Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Loose_RatiotoNominal,       cErr, cErr);
  grVn6Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Loose_RatiotoNominal,       cErr, cErr);
  grVn8Loose_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Loose_RatiotoNominal,       cErr, cErr);
  grGamma1ExpLoose_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expLoose_RatiotoNominal, cErr, cErr);
  grVn6Vn4Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Loose_RatiotoNominal,    cErr, cErr);
  grVn8Vn4Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Loose_RatiotoNominal,    cErr, cErr);
  grVn8Vn6Loose_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Loose_RatiotoNominal,    cErr, cErr);

  grVn2Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn2Tight_RatiotoNominal,       cErr, cErr);
  grVn4Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn4Tight_RatiotoNominal,       cErr, cErr);
  grVn6Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn6Tight_RatiotoNominal,       cErr, cErr);
  grVn8Tight_RatiotoNominal       = new TGraphErrors(NCENT, centBinCenter, vn8Tight_RatiotoNominal,       cErr, cErr);
  grGamma1ExpTight_RatiotoNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expTight_RatiotoNominal, cErr, cErr);
  grVn6Vn4Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Tight_RatiotoNominal,    cErr, cErr);
  grVn8Vn4Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Tight_RatiotoNominal,    cErr, cErr);
  grVn8Vn6Tight_RatiotoNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Tight_RatiotoNominal,    cErr, cErr);

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

  formatGraph(grVn4Loose_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kSpring+4, 25, "grVn4Loose_RatiotoNominal");
  formatGraph(grVn4Tight_RatiotoNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kViolet-1, 25, "grVn4Tight_RatiotoNominal");

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
  line->Draw("same");

  formatGraph(grGamma1ExpLoose_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 2, 28, "grGamma1ExpLoose_RatiotoNominal");
  formatGraph(grGamma1ExpTight_RatiotoNominal, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 1, 28, "grGamma1ExpTight_RatiotoNominal");

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

  TLegend * leg51 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
  leg51->SetFillStyle(0);
  leg51->SetBorderSize(0);
  leg51->AddEntry(grVn6Vn4Loose_RatiotoNominal, "loose/nominal", "lp");
  leg51->AddEntry(grVn6Vn4Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg52 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
  leg52->SetFillStyle(0);
  leg52->SetBorderSize(0);
  leg52->AddEntry(grVn8Vn4Loose_RatiotoNominal, "loose/nominal", "lp");
  leg52->AddEntry(grVn8Vn4Tight_RatiotoNominal, "tight/nominal", "lp");

  TLegend * leg53 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
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

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){

    double pctVn2_3       = fabs(vn2Loose[icent] - vn2Nominal[icent]) / fabs(vn2Nominal[icent]);
    double pctVn4_3       = fabs(vn4Loose[icent] - vn4Nominal[icent]) / fabs(vn4Nominal[icent]);
    double pctVn6_3       = fabs(vn6Loose[icent] - vn6Nominal[icent]) / fabs(vn6Nominal[icent]);
    double pctVn8_3       = fabs(vn8Loose[icent] - vn8Nominal[icent]) / fabs(vn8Nominal[icent]);
    double pctGamma1Exp_3 = fabs(gamma1expLoose[icent] - gamma1expNominal[icent]) / fabs(gamma1expNominal[icent]);
    double pctVn6Vn4_3    = fabs(vn6vn4Loose[icent] - vn6vn4Nominal[icent]) / fabs(vn6vn4Nominal[icent]);
    double pctVn8Vn4_3    = fabs(vn8vn4Loose[icent] - vn8vn4Nominal[icent]) / fabs(vn8vn4Nominal[icent]);
    double pctVn8Vn6_3    = fabs(vn8vn6Loose[icent] - vn8vn6Nominal[icent]) / fabs(vn8vn6Nominal[icent]);

    double pctVn2_3_15       = fabs(vn2Tight[icent] - vn2Nominal[icent]) / fabs(vn2Nominal[icent]);
    double pctVn4_3_15       = fabs(vn4Tight[icent] - vn4Nominal[icent]) / fabs(vn4Nominal[icent]);
    double pctVn6_3_15       = fabs(vn6Tight[icent] - vn6Nominal[icent]) / fabs(vn6Nominal[icent]);
    double pctVn8_3_15       = fabs(vn8Tight[icent] - vn8Nominal[icent]) / fabs(vn8Nominal[icent]);
    double pctGamma1Exp_3_15 = fabs(gamma1expTight[icent] - gamma1expNominal[icent]) / fabs(gamma1expNominal[icent]);
    double pctVn6Vn4_3_15    = fabs(vn6vn4Tight[icent] - vn6vn4Nominal[icent]) / fabs(vn6vn4Nominal[icent]);
    double pctVn8Vn4_3_15    = fabs(vn8vn4Tight[icent] - vn8vn4Nominal[icent]) / fabs(vn8vn4Nominal[icent]);
    double pctVn8Vn6_3_15    = fabs(vn8vn6Tight[icent] - vn8vn6Nominal[icent]) / fabs(vn8vn6Nominal[icent]);

    double reVn2       = max(pctVn2_3, pctVn2_3_15);
    double reVn4       = max(pctVn4_3, pctVn4_3_15);
    double reVn6       = max(pctVn6_3, pctVn6_3_15);
    double reVn8       = max(pctVn8_3, pctVn8_3_15);
    double reGamma1Exp = max(pctGamma1Exp_3, pctGamma1Exp_3_15);
    double reVn6Vn4    = max(pctVn6Vn4_3, pctVn6Vn4_3_15);
    double reVn8Vn4    = max(pctVn8Vn4_3, pctVn8Vn4_3_15);
    double reVn8Vn6    = max(pctVn8Vn6_3, pctVn8Vn6_3_15);

    if( TMath::IsNaN(reVn2) || std::isinf(reVn2) )             reVn2 = 0;
    if( TMath::IsNaN(reVn4) || std::isinf(reVn4) )             reVn4 = 0;
    if( TMath::IsNaN(reVn6) || std::isinf(reVn6) )             reVn6 = 0;
    if( TMath::IsNaN(reVn8) || std::isinf(reVn8) )             reVn8 = 0;
    if( TMath::IsNaN(reGamma1Exp) || std::isinf(reGamma1Exp) ) reGamma1Exp = 0;
    if( TMath::IsNaN(reVn6Vn4) || std::isinf(reVn6Vn4) )       reVn6Vn4 = 0;
    if( TMath::IsNaN(reVn8Vn4) || std::isinf(reVn8Vn4) )       reVn8Vn4 = 0;
    if( TMath::IsNaN(reVn8Vn6) || std::isinf(reVn8Vn6) )       reVn8Vn6 = 0;

    relErrVn2->SetBinContent(icent+1, reVn2);
    relErrVn4->SetBinContent(icent+1, reVn4);
    relErrVn6->SetBinContent(icent+1, reVn6);
    relErrVn8->SetBinContent(icent+1, reVn8);
    relErrGamma1Exp->SetBinContent(icent+1, reGamma1Exp);
    relErrVn6Vn4->SetBinContent(icent+1, reVn6Vn4);
    relErrVn8Vn4->SetBinContent(icent+1, reVn8Vn4);
    relErrVn8Vn6->SetBinContent(icent+1, reVn8Vn6);

  }

  fRelErr->Write();

}

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

void sysVtxCut(){

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
  TFile * fAnaVtx3;
  TFile * fUnfVtx3;
  TH1D * hObsVtx3[NCENT];
  TH1D * hUnfoldVtx3[NCENT][NITER];
  TH1D * hRefoldVtx3[NCENT][NITER];

  TGraphErrors * grVn2Vtx3;
  TGraphErrors * grVn4Vtx3;
  TGraphErrors * grVn6Vtx3;
  TGraphErrors * grVn8Vtx3;
  TGraphErrors * grGamma1ExpVtx3;
  TGraphErrors * grVn6Vn4Vtx3;
  TGraphErrors * grVn8Vn4Vtx3;
  TGraphErrors * grVn8Vn6Vtx3;

  double vn2Vtx3[NCENT];
  double vn4Vtx3[NCENT];
  double vn6Vtx3[NCENT];
  double vn8Vtx3[NCENT];
  double gamma1expVtx3[NCENT];
  double vn6vn4Vtx3[NCENT];
  double vn8vn4Vtx3[NCENT];
  double vn8vn6Vtx3[NCENT];

  TGraphErrors * grVn2Vtx3_RatiotoVtx15;
  TGraphErrors * grVn4Vtx3_RatiotoVtx15;
  TGraphErrors * grVn6Vtx3_RatiotoVtx15;
  TGraphErrors * grVn8Vtx3_RatiotoVtx15;
  TGraphErrors * grGamma1ExpVtx3_RatiotoVtx15;
  TGraphErrors * grVn6Vn4Vtx3_RatiotoVtx15;
  TGraphErrors * grVn8Vn4Vtx3_RatiotoVtx15;
  TGraphErrors * grVn8Vn6Vtx3_RatiotoVtx15;

  double vn2Vtx3_RatiotoVtx15[NCENT];
  double vn4Vtx3_RatiotoVtx15[NCENT];
  double vn6Vtx3_RatiotoVtx15[NCENT];
  double vn8Vtx3_RatiotoVtx15[NCENT];
  double gamma1expVtx3_RatiotoVtx15[NCENT];
  double vn6vn4Vtx3_RatiotoVtx15[NCENT];
  double vn8vn4Vtx3_RatiotoVtx15[NCENT];
  double vn8vn6Vtx3_RatiotoVtx15[NCENT];

  bool iterCutoffVtx3[NCENT];

  //-- 3.0 < Vertex < 15.0
  TFile * fAnaVtx3_15;
  TFile * fUnfVtx3_15;
  TH1D * hObsVtx3_15[NCENT];
  TH1D * hUnfoldVtx3_15[NCENT][NITER];
  TH1D * hRefoldVtx3_15[NCENT][NITER];

  TGraphErrors * grVn2Vtx3_15;
  TGraphErrors * grVn4Vtx3_15;
  TGraphErrors * grVn6Vtx3_15;
  TGraphErrors * grVn8Vtx3_15;
  TGraphErrors * grGamma1ExpVtx3_15;
  TGraphErrors * grVn6Vn4Vtx3_15;
  TGraphErrors * grVn8Vn4Vtx3_15;
  TGraphErrors * grVn8Vn6Vtx3_15;

  double vn2Vtx3_15[NCENT];
  double vn4Vtx3_15[NCENT];
  double vn6Vtx3_15[NCENT];
  double vn8Vtx3_15[NCENT];
  double gamma1expVtx3_15[NCENT];
  double vn6vn4Vtx3_15[NCENT];
  double vn8vn4Vtx3_15[NCENT];
  double vn8vn6Vtx3_15[NCENT];

  TGraphErrors * grVn2Vtx3_15_RatiotoVtx15;
  TGraphErrors * grVn4Vtx3_15_RatiotoVtx15;
  TGraphErrors * grVn6Vtx3_15_RatiotoVtx15;
  TGraphErrors * grVn8Vtx3_15_RatiotoVtx15;
  TGraphErrors * grGamma1ExpVtx3_15_RatiotoVtx15;
  TGraphErrors * grVn6Vn4Vtx3_15_RatiotoVtx15;
  TGraphErrors * grVn8Vn4Vtx3_15_RatiotoVtx15;
  TGraphErrors * grVn8Vn6Vtx3_15_RatiotoVtx15;

  double vn2Vtx3_15_RatiotoVtx15[NCENT];
  double vn4Vtx3_15_RatiotoVtx15[NCENT];
  double vn6Vtx3_15_RatiotoVtx15[NCENT];
  double vn8Vtx3_15_RatiotoVtx15[NCENT];
  double gamma1expVtx3_15_RatiotoVtx15[NCENT];
  double vn6vn4Vtx3_15_RatiotoVtx15[NCENT];
  double vn8vn4Vtx3_15_RatiotoVtx15[NCENT];
  double vn8vn6Vtx3_15_RatiotoVtx15[NCENT];

  bool iterCutoffVtx3_15[NCENT];

  //-- Vertex < 15.0
  TFile * fAnaVtx15;
  TFile * fUnfVtx15;
  TH1D * hObsVtx15[NCENT];
  TH1D * hUnfoldVtx15[NCENT][NITER];
  TH1D * hRefoldVtx15[NCENT][NITER];

  TGraphErrors * grVn2Vtx15;
  TGraphErrors * grVn4Vtx15;
  TGraphErrors * grVn6Vtx15;
  TGraphErrors * grVn8Vtx15;
  TGraphErrors * grGamma1ExpVtx15;
  TGraphErrors * grVn6Vn4Vtx15;
  TGraphErrors * grVn8Vn4Vtx15;
  TGraphErrors * grVn8Vn6Vtx15;

  double vn2Vtx15[NCENT];
  double vn4Vtx15[NCENT];
  double vn6Vtx15[NCENT];
  double vn8Vtx15[NCENT];
  double gamma1expVtx15[NCENT];
  double vn6vn4Vtx15[NCENT];
  double vn8vn4Vtx15[NCENT];
  double vn8vn6Vtx15[NCENT];

  bool iterCutoffVtx15[NCENT];

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
  fAnaVtx3    = new TFile( "vtx_leq_3/AnalyzerResults/CastleEbyE.root" );
  fUnfVtx3    = new TFile( Form("vtx_leq_3/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaVtx3_15 = new TFile( "vtx3_15/AnalyzerResults/CastleEbyE.root" );
  fUnfVtx3_15 = new TFile( Form("vtx3_15/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaVtx15   = new TFile( "../../AnalyzerResults/CastleEbyE.root"  );
  fUnfVtx15   = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorVtx.root", "recreate");
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
    hObsVtx3[icent] = (TH1D*) fAnaVtx3->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsVtx3[icent]->SetName( Form("hVnFull_c%i_vtx3", icent) );

    //-- Iter loop
    iterCutoffVtx3[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldVtx3[icent][i] = (TH1D*) fUnfVtx3->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldVtx3[icent][i]->SetName( Form("hreco%i_c%i_vtx3", iter[i], icent) );
      hUnfoldVtx3[icent][i]->SetLineColor(col[i]);
      hUnfoldVtx3[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldVtx3[icent][i] = (TH1D*) fUnfVtx3->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldVtx3[icent][i]->SetName( Form("hrefold%i_c%i_vtx3", iter[i], icent) );
      hRefoldVtx3[icent][i]->SetLineWidth(2);
      hRefoldVtx3[icent][i]->SetLineColor(col[i]);
      hRefoldVtx3[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldVtx3[icent][i]->Chi2Test(hObsVtx3[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
	iterCutoffVtx3[icent] = i;
	
	EbyECumu cumu(hUnfoldVtx3[icent][i]);
	vn2Vtx3[icent]  = cumu.GetCumu_vn2();
	vn4Vtx3[icent]  = cumu.GetCumu_vn4();
	vn6Vtx3[icent]  = cumu.GetCumu_vn6();
	vn8Vtx3[icent]  = cumu.GetCumu_vn8();
	gamma1expVtx3[icent] = cumu.GetGamma1Exp();

	if( vn4Vtx3[icent] == 0 ){
	  vn6vn4Vtx3[icent] = 0;
	  vn8vn4Vtx3[icent] = 0;
	}
	else{
	  vn6vn4Vtx3[icent] = vn6Vtx3[icent] / vn4Vtx3[icent];
	  vn8vn4Vtx3[icent] = vn8Vtx3[icent] / vn4Vtx3[icent];
	}
	if( vn6Vtx3[icent] == 0 ) vn8vn6Vtx3[icent] = 0;
	else                      vn8vn6Vtx3[icent] = vn8Vtx3[icent] / vn6Vtx3[icent];

	break;
      }
      if( i == NITER - 1 ){
	iterCutoffVtx3[icent] = i;

        EbyECumu cumu(hUnfoldVtx3[icent][i]);
        vn2Vtx3[icent]  = cumu.GetCumu_vn2();
        vn4Vtx3[icent]  = cumu.GetCumu_vn4();
        vn6Vtx3[icent]  = cumu.GetCumu_vn6();
	vn8Vtx3[icent]  = cumu.GetCumu_vn8();
        gamma1expVtx3[icent] = cumu.GetGamma1Exp();

        if( vn4Vtx3[icent] == 0 ){
          vn6vn4Vtx3[icent] = 0;
          vn8vn4Vtx3[icent] = 0;
	}
        else{
          vn6vn4Vtx3[icent] = vn6Vtx3[icent] / vn4Vtx3[icent];
          vn8vn4Vtx3[icent] = vn8Vtx3[icent] / vn4Vtx3[icent];
	}
        if( vn6Vtx3[icent] == 0 ) vn8vn6Vtx3[icent] = 0;
        else                      vn8vn6Vtx3[icent] = vn8Vtx3[icent] / vn6Vtx3[icent];

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- 3.0 < Vertex < 15.0 --------------------
    hObsVtx3_15[icent] = (TH1D*) fAnaVtx3_15->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsVtx3_15[icent]->SetName( Form("hVnFull_c%i_vtx3", icent) );

    //-- Iter loop
    iterCutoffVtx3_15[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldVtx3_15[icent][i] = (TH1D*) fUnfVtx3_15->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldVtx3_15[icent][i]->SetName( Form("hreco%i_c%i_vtx3", iter[i], icent) );
      hUnfoldVtx3_15[icent][i]->SetLineColor(col[i]);
      hUnfoldVtx3_15[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldVtx3_15[icent][i] = (TH1D*) fUnfVtx3_15->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldVtx3_15[icent][i]->SetName( Form("hrefold%i_c%i_vtx3", iter[i], icent) );
      hRefoldVtx3_15[icent][i]->SetLineWidth(2);
      hRefoldVtx3_15[icent][i]->SetLineColor(col[i]);
      hRefoldVtx3_15[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldVtx3_15[icent][i]->Chi2Test(hObsVtx3_15[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
        iterCutoffVtx3_15[icent] = i;

        EbyECumu cumu(hUnfoldVtx3_15[icent][i]);
        vn2Vtx3_15[icent]  = cumu.GetCumu_vn2();
        vn4Vtx3_15[icent]  = cumu.GetCumu_vn4();
        vn6Vtx3_15[icent]  = cumu.GetCumu_vn6();
	vn8Vtx3_15[icent]  = cumu.GetCumu_vn8();
        gamma1expVtx3_15[icent] = cumu.GetGamma1Exp();

        if( vn4Vtx3_15[icent] == 0 ){
          vn6vn4Vtx3_15[icent] = 0;
          vn8vn4Vtx3_15[icent] = 0;
	}
        else{
          vn6vn4Vtx3_15[icent] = vn6Vtx3_15[icent] / vn4Vtx3_15[icent];
          vn8vn4Vtx3_15[icent] = vn8Vtx3_15[icent] / vn4Vtx3_15[icent];
	}
        if( vn6Vtx3_15[icent] == 0 ) vn8vn6Vtx3_15[icent] = 0;
        else                      vn8vn6Vtx3_15[icent] = vn8Vtx3_15[icent] / vn6Vtx3_15[icent];

        break;
      }
      if( i == NITER - 1 ){
        iterCutoffVtx3_15[icent] = i;

        EbyECumu cumu(hUnfoldVtx3_15[icent][i]);
        vn2Vtx3_15[icent]  = cumu.GetCumu_vn2();
        vn4Vtx3_15[icent]  = cumu.GetCumu_vn4();
        vn6Vtx3_15[icent]  = cumu.GetCumu_vn6();
        vn8Vtx3_15[icent]  = cumu.GetCumu_vn8();
        gamma1expVtx3_15[icent] = cumu.GetGamma1Exp();

        if( vn4Vtx3_15[icent] == 0 ){
          vn6vn4Vtx3_15[icent] = 0;
          vn8vn4Vtx3_15[icent] = 0;
        }
        else{
          vn6vn4Vtx3_15[icent] = vn6Vtx3_15[icent] / vn4Vtx3_15[icent];
          vn8vn4Vtx3_15[icent] = vn8Vtx3_15[icent] / vn4Vtx3_15[icent];
        }
        if( vn6Vtx3_15[icent] == 0 ) vn8vn6Vtx3_15[icent] = 0;
        else                      vn8vn6Vtx3_15[icent] = vn8Vtx3_15[icent] / vn6Vtx3_15[icent];

        break;
      }
 
    } //-- End iter loop

    //-- -------------------- Vertex < 15.0 --------------------
    hObsVtx15[icent] = (TH1D*) fAnaVtx15->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsVtx15[icent]->SetName( Form("hVnFull_c%i_vtx3", icent) );

    //-- Iter loop
    iterCutoffVtx15[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the unfolded histograms
      hUnfoldVtx15[icent][i] = (TH1D*) fUnfVtx15->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldVtx15[icent][i]->SetName( Form("hreco%i_c%i_vtx3", iter[i], icent) );
      hUnfoldVtx15[icent][i]->SetLineColor(col[i]);
      hUnfoldVtx15[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefoldVtx15[icent][i] = (TH1D*) fUnfVtx15->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldVtx15[icent][i]->SetName( Form("hrefold%i_c%i_vtx3", iter[i], icent) );
      hRefoldVtx15[icent][i]->SetLineWidth(2);
      hRefoldVtx15[icent][i]->SetLineColor(col[i]);
      hRefoldVtx15[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold = hRefoldVtx15[icent][i]->Chi2Test(hObsVtx15[icent], "CHI2/NDF");

      if( chi2NDF_Refold < 1.2 ){
        iterCutoffVtx15[icent] = i;

        EbyECumu cumu(hUnfoldVtx15[icent][i]);
        vn2Vtx15[icent]  = cumu.GetCumu_vn2();
        vn4Vtx15[icent]  = cumu.GetCumu_vn4();
        vn6Vtx15[icent]  = cumu.GetCumu_vn6();
	vn8Vtx15[icent]  = cumu.GetCumu_vn8();
        gamma1expVtx15[icent] = cumu.GetGamma1Exp();

        if( vn4Vtx15[icent] == 0 ){
          vn6vn4Vtx15[icent] = 0;
          vn8vn4Vtx15[icent] = 0;
	}
        else{
          vn6vn4Vtx15[icent] = vn6Vtx15[icent] / vn4Vtx15[icent];
          vn8vn4Vtx15[icent] = vn8Vtx15[icent] / vn4Vtx15[icent];
	}
        if( vn6Vtx15[icent] == 0 ) vn8vn6Vtx15[icent] = 0;
        else                vn8vn6Vtx15[icent] = vn8Vtx15[icent] / vn6Vtx15[icent];

        break;
      }
      if( i == NITER - 1 ){
        iterCutoffVtx15[icent] = i;

        EbyECumu cumu(hUnfoldVtx15[icent][i]);
        vn2Vtx15[icent]  = cumu.GetCumu_vn2();
        vn4Vtx15[icent]  = cumu.GetCumu_vn4();
        vn6Vtx15[icent]  = cumu.GetCumu_vn6();
        vn8Vtx15[icent]  = cumu.GetCumu_vn8();
        gamma1expVtx15[icent] = cumu.GetGamma1Exp();

        if( vn4Vtx15[icent] == 0 ){
          vn6vn4Vtx15[icent] = 0;
          vn8vn4Vtx15[icent] = 0;
        }
        else{
          vn6vn4Vtx15[icent] = vn6Vtx15[icent] / vn4Vtx15[icent];
          vn8vn4Vtx15[icent] = vn8Vtx15[icent] / vn4Vtx15[icent];
        }
        if( vn6Vtx15[icent] == 0 ) vn8vn6Vtx15[icent] = 0;
        else                      vn8vn6Vtx15[icent] = vn8Vtx15[icent] / vn6Vtx15[icent];

        break;
      }
 
    } //-- End iter loop

  } //-- End Cent loop

  //-- Calculate ratios to the vtx < 15 case
  for(int icent = 0; icent < NCENT; icent++){
    
    double vn23        = vn2Vtx3[icent];
    double vn43        = vn4Vtx3[icent];
    double vn63        = vn6Vtx3[icent];
    double vn83        = vn8Vtx3[icent];
    double gamma1exp3  = gamma1expVtx3[icent];
    double vn6vn43     = vn6vn4Vtx3[icent];
    double vn8vn43     = vn8vn4Vtx3[icent];
    double vn8vn63     = vn8vn6Vtx3[icent];

    double vn23_15        = vn2Vtx3_15[icent];
    double vn43_15        = vn4Vtx3_15[icent];
    double vn63_15        = vn6Vtx3_15[icent];
    double vn83_15        = vn8Vtx3_15[icent];
    double gamma1exp3_15  = gamma1expVtx3_15[icent];
    double vn6vn43_15     = vn6vn4Vtx3_15[icent];
    double vn8vn43_15     = vn8vn4Vtx3_15[icent];
    double vn8vn63_15     = vn8vn6Vtx3_15[icent];

    double vn215        = vn2Vtx15[icent];
    double vn415        = vn4Vtx15[icent];
    double vn615        = vn6Vtx15[icent];
    double vn815        = vn8Vtx15[icent];
    double gamma1exp15  = gamma1expVtx15[icent];
    double vn6vn415     = vn6vn4Vtx15[icent];
    double vn8vn415     = vn8vn4Vtx15[icent];
    double vn8vn615     = vn8vn6Vtx15[icent];

    double r23         = vn23 / vn215;
    double r43         = vn43 / vn415;
    double r63         = vn63 / vn615;
    double r83         = vn83 / vn815;
    double rgamma1exp3 = gamma1exp3 / gamma1exp15;
    double r643        = vn6vn43 / vn6vn415;
    double r843        = vn8vn43 / vn8vn415;
    double r863        = vn8vn63 / vn8vn615;

    double r23_15         = vn23_15 / vn215;
    double r43_15         = vn43_15 / vn415;
    double r63_15         = vn63_15 / vn615;
    double r83_15         = vn83_15 / vn815;
    double rgamma1exp3_15 = gamma1exp3_15 / gamma1exp15;
    double r643_15        = vn6vn43_15 / vn6vn415;
    double r843_15        = vn8vn43_15 / vn8vn415;
    double r863_15        = vn8vn63_15 / vn8vn615;

    vn2Vtx3_RatiotoVtx15[icent]       = r23;
    vn4Vtx3_RatiotoVtx15[icent]       = r43;
    vn6Vtx3_RatiotoVtx15[icent]       = r63;
    vn8Vtx3_RatiotoVtx15[icent]       = r83;
    gamma1expVtx3_RatiotoVtx15[icent] = rgamma1exp3;
    vn6vn4Vtx3_RatiotoVtx15[icent]    = r643;
    vn8vn4Vtx3_RatiotoVtx15[icent]    = r843;
    vn8vn6Vtx3_RatiotoVtx15[icent]    = r863;

    vn2Vtx3_15_RatiotoVtx15[icent]       = r23_15;
    vn4Vtx3_15_RatiotoVtx15[icent]       = r43_15;
    vn6Vtx3_15_RatiotoVtx15[icent]       = r63_15;
    vn8Vtx3_15_RatiotoVtx15[icent]       = r83_15;
    gamma1expVtx3_15_RatiotoVtx15[icent] = rgamma1exp3_15;
    vn6vn4Vtx3_15_RatiotoVtx15[icent]    = r643_15;
    vn8vn4Vtx3_15_RatiotoVtx15[icent]    = r843_15;
    vn8vn6Vtx3_15_RatiotoVtx15[icent]    = r863_15;

  }

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT]; 
  for(int icent = 0; icent < NCENT; icent++) cErr[icent] = 0;

  //-- Ratios 
  grVn2Vtx3_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn2Vtx3_RatiotoVtx15,       cErr, cErr);
  grVn4Vtx3_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn4Vtx3_RatiotoVtx15,       cErr, cErr);
  grVn6Vtx3_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn6Vtx3_RatiotoVtx15,       cErr, cErr);
  grVn8Vtx3_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn8Vtx3_RatiotoVtx15,       cErr, cErr);
  grGamma1ExpVtx3_RatiotoVtx15 = new TGraphErrors(NCENT, centBinCenter, gamma1expVtx3_RatiotoVtx15, cErr, cErr);
  grVn6Vn4Vtx3_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Vtx3_RatiotoVtx15,    cErr, cErr);
  grVn8Vn4Vtx3_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Vtx3_RatiotoVtx15,    cErr, cErr);
  grVn8Vn6Vtx3_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Vtx3_RatiotoVtx15,    cErr, cErr);

  grVn2Vtx3_15_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn2Vtx3_15_RatiotoVtx15,       cErr, cErr);
  grVn4Vtx3_15_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn4Vtx3_15_RatiotoVtx15,       cErr, cErr);
  grVn6Vtx3_15_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn6Vtx3_15_RatiotoVtx15,       cErr, cErr);
  grVn8Vtx3_15_RatiotoVtx15       = new TGraphErrors(NCENT, centBinCenter, vn8Vtx3_15_RatiotoVtx15,       cErr, cErr);
  grGamma1ExpVtx3_15_RatiotoVtx15 = new TGraphErrors(NCENT, centBinCenter, gamma1expVtx3_15_RatiotoVtx15, cErr, cErr);
  grVn6Vn4Vtx3_15_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Vtx3_15_RatiotoVtx15,    cErr, cErr);
  grVn8Vn4Vtx3_15_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Vtx3_15_RatiotoVtx15,    cErr, cErr);
  grVn8Vn6Vtx3_15_RatiotoVtx15    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Vtx3_15_RatiotoVtx15,    cErr, cErr);

  //-- Cumu Ratio Plots
  TLine * line = new TLine(cent_min[0], 1.0, grVn2Vtx3_RatiotoVtx15->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * cSysCumu = new TCanvas("cSysCumu", "cSysCumu", 1000, 1000);
  cSysCumu->Divide(2,2);

  cSysCumu->cd(1);
  grVn2Vtx3_RatiotoVtx15->Draw("ap");
  grVn2Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(2);
  grVn4Vtx3_RatiotoVtx15->Draw("ap");
  grVn4Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(3);
  grVn6Vtx3_RatiotoVtx15->Draw("ap");
  grVn6Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  cSysCumu->cd(4);
  grVn8Vtx3_RatiotoVtx15->Draw("ap");
  grVn8Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  formatGraph(grVn2Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), 1, 24, "grVn2Vtx3_RatiotoVtx15");
  formatGraph(grVn2Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} Ratio", norder_), kCyan+2, 24, "grVn2Vtx3_15_RatiotoVtx15");

  formatGraph(grVn4Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kSpring+4, 25, "grVn4Vtx3_RatiotoVtx15");
  formatGraph(grVn4Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} Ratio", norder_), kViolet-1, 25, "grVn4Vtx3_15_RatiotoVtx15");

  formatGraph(grVn6Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 6, 28, "grVn6Vtx3_RatiotoVtx15");
  formatGraph(grVn6Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} Ratio", norder_), 4, 28, "grVn6Vtx3_15_RatiotoVtx15");

  formatGraph(grVn8Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kOrange+7, 27, "grVn8Vtx3_RatiotoVtx15");
  formatGraph(grVn8Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} Ratio", norder_), kGray+2, 27, "grVn8Vtx3_15_RatiotoVtx15");

  TLegend * leg31 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg31->SetFillStyle(0);
  leg31->SetBorderSize(0);
  leg31->AddEntry(grVn2Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg31->AddEntry(grVn2Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  TLegend * leg32 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg32->SetFillStyle(0);
  leg32->SetBorderSize(0);
  leg32->AddEntry(grVn4Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg32->AddEntry(grVn4Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  TLegend * leg33 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg33->SetFillStyle(0);
  leg33->SetBorderSize(0);
  leg33->AddEntry(grVn6Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg33->AddEntry(grVn6Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  TLegend * leg34 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg34->SetFillStyle(0);
  leg34->SetBorderSize(0);
  leg34->AddEntry(grVn8Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg34->AddEntry(grVn8Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  cSysCumu->cd(1);
  leg31->Draw("same");
  cSysCumu->cd(2);
  leg32->Draw("same");
  cSysCumu->cd(3);
  leg33->Draw("same");
  cSysCumu->cd(4);
  leg34->Draw("same");

  cSysCumu->Update();
  cSysCumu->SaveAs("../../plots/systematicStudies/SysVtx_CumuCent.pdf");

  //-- Gamma1Exp Ratio Plot
  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  grGamma1ExpVtx3_RatiotoVtx15->Draw("ap");
  grGamma1ExpVtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  formatGraph(grGamma1ExpVtx3_RatiotoVtx15,    "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 2, 28, "grGamma1ExpVtx3_RatiotoVtx15");
  formatGraph(grGamma1ExpVtx3_15_RatiotoVtx15, "Centrality %", gamma1ExpRatioMin, gamma1ExpRatioMax, "#gamma_{1}^{exp} Ratio", 1, 28, "grGamma1ExpVtx3_15_RatiotoVtx15");

  TLegend * leg4 = new TLegend(0.7000, 0.7034, 0.9456, 0.9449);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->AddEntry(grGamma1ExpVtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg4->AddEntry(grGamma1ExpVtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");
  leg4->Draw("same");

  cGamma1Exp->Update();
  cGamma1Exp->SaveAs("../../plots/systematicStudies/SysVtxGamma1ExpCent.pdf");

  //-- Cumu Ratio Plot
  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  grVn6Vn4Vtx3_RatiotoVtx15->Draw("ap");
  grVn6Vn4Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4Vtx3_RatiotoVtx15->Draw("ap");
  grVn8Vn4Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6Vtx3_RatiotoVtx15->Draw("ap");
  grVn8Vn6Vtx3_15_RatiotoVtx15->Draw("psame");
  line->Draw("same");

  formatGraph(grVn6Vn4Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 4, 21, "grVn6Vn4Vtx3_RatiotoVtx15");
  formatGraph(grVn6Vn4Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 2, 21, "grVn6Vn4Vtx3_15_RatiotoVtx15");

  formatGraph(grVn8Vn4Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), kGreen+2, 34, "grVn8Vn4Vtx3_RatiotoVtx15");
  formatGraph(grVn8Vn4Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), 6, 34, "grVn8Vn4Vtx3_15_RatiotoVtx15");

  formatGraph(grVn8Vn6Vtx3_RatiotoVtx15,    "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kViolet-1, 33, "grVn8Vn6Vtx3_RatiotoVtx15");
  formatGraph(grVn8Vn6Vtx3_15_RatiotoVtx15, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kOrange+7, 33, "grVn8Vn6Vtx3_15_RatiotoVtx15");

  TLegend * leg51 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
  leg51->SetFillStyle(0);
  leg51->SetBorderSize(0);
  leg51->AddEntry(grVn6Vn4Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg51->AddEntry(grVn6Vn4Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  TLegend * leg52 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
  leg52->SetFillStyle(0);
  leg52->SetBorderSize(0);
  leg52->AddEntry(grVn8Vn4Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg52->AddEntry(grVn8Vn4Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  TLegend * leg53 = new TLegend(0.1855, 0.2013, 0.4315, 0.4449);
  leg53->SetFillStyle(0);
  leg53->SetBorderSize(0);
  leg53->AddEntry(grVn8Vn6Vtx3_RatiotoVtx15,    "#frac{|v_{z}| < 3.0}{|v_{z}| < 15.0}", "lp");
  leg53->AddEntry(grVn8Vn6Vtx3_15_RatiotoVtx15, "#frac{3.0 < |v_{z}| < 15.0}{|v_{z}| < 15.0}", "lp");

  cCumuRatioSys->cd(1);
  leg51->Draw("same");
  cCumuRatioSys->cd(2);
  leg52->Draw("same");
  cCumuRatioSys->cd(3);
  leg53->Draw("same");

  cCumuRatioSys->Update();
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/SysVtx_CumuRatioCent.pdf");

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){

    double pctVn2_3       = fabs(vn2Vtx3[icent] - vn2Vtx15[icent]) / fabs(vn2Vtx15[icent]);
    double pctVn4_3       = fabs(vn4Vtx3[icent] - vn4Vtx15[icent]) / fabs(vn4Vtx15[icent]);
    double pctVn6_3       = fabs(vn6Vtx3[icent] - vn6Vtx15[icent]) / fabs(vn6Vtx15[icent]);
    double pctVn8_3       = fabs(vn8Vtx3[icent] - vn8Vtx15[icent]) / fabs(vn8Vtx15[icent]);
    double pctGamma1Exp_3 = fabs(gamma1expVtx3[icent] - gamma1expVtx15[icent]) / fabs(gamma1expVtx15[icent]);
    double pctVn6Vn4_3    = fabs(vn6vn4Vtx3[icent] - vn6vn4Vtx15[icent]) / fabs(vn6vn4Vtx15[icent]);
    double pctVn8Vn4_3    = fabs(vn8vn4Vtx3[icent] - vn8vn4Vtx15[icent]) / fabs(vn8vn4Vtx15[icent]);
    double pctVn8Vn6_3    = fabs(vn8vn6Vtx3[icent] - vn8vn6Vtx15[icent]) / fabs(vn8vn6Vtx15[icent]);

    double pctVn2_3_15       = fabs(vn2Vtx3_15[icent] - vn2Vtx15[icent]) / fabs(vn2Vtx15[icent]);
    double pctVn4_3_15       = fabs(vn4Vtx3_15[icent] - vn4Vtx15[icent]) / fabs(vn4Vtx15[icent]);
    double pctVn6_3_15       = fabs(vn6Vtx3_15[icent] - vn6Vtx15[icent]) / fabs(vn6Vtx15[icent]);
    double pctVn8_3_15       = fabs(vn8Vtx3_15[icent] - vn8Vtx15[icent]) / fabs(vn8Vtx15[icent]);
    double pctGamma1Exp_3_15 = fabs(gamma1expVtx3_15[icent] - gamma1expVtx15[icent]) / fabs(gamma1expVtx15[icent]);
    double pctVn6Vn4_3_15    = fabs(vn6vn4Vtx3_15[icent] - vn6vn4Vtx15[icent]) / fabs(vn6vn4Vtx15[icent]);
    double pctVn8Vn4_3_15    = fabs(vn8vn4Vtx3_15[icent] - vn8vn4Vtx15[icent]) / fabs(vn8vn4Vtx15[icent]);
    double pctVn8Vn6_3_15    = fabs(vn8vn6Vtx3_15[icent] - vn8vn6Vtx15[icent]) / fabs(vn8vn6Vtx15[icent]);

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

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
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/SysTablesEbyESE.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;
using namespace sysebyese;

void sysResultsPlots(){

  TH1D::SetDefaultSumw2();

  int centbin     = 4;
  int norder_     = 3;
  double tkEta    = 1.0;

  bool dosys     = 0;
  bool smoothSys = 1;
  bool compATLAS = 1;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double cumuMin = 0.0;
  double cumuMax = 0.045;

  double gamm1expMin = -1.0;
  double gamm1expMax = 0.5;
  double ratioMin    = 0.95;
  double ratioMax    = 1.01999;
  double ratioMinVn8Vn6 = 0.994;
  double ratioMaxVn8Vn6 = 1.003;

  double sysPctMin = 0.00001;
  double sysPctMax = 50.;

  double titleOffset = 1.4;

  TLatex latex;
  TLatex latex2;

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
  double vn2RawSyse[NCENT];
  TGraphErrors * grVn2Raw;
  TGraphErrors * grVn2RawSys;
  TGraphErrors * grVn2Theory;

  //--vn{4}
  double vn4Raw[NCENT];
  double vn4RawState[NCENT];
  double vn4RawSyse[NCENT];
  TGraphErrors * grVn4Raw;
  TGraphErrors * grVn4RawSys;
  TGraphErrors * grVn4Theory;

  //--vn{6}
  double vn6Raw[NCENT];
  double vn6RawState[NCENT];
  double vn6RawSyse[NCENT];
  TGraphErrors * grVn6Raw;
  TGraphErrors * grVn6RawSys;
  TGraphErrors * grVn6Theory;

  //--vn{8}
  double vn8Raw[NCENT];
  double vn8RawState[NCENT];
  double vn8RawSyse[NCENT];
  TGraphErrors * grVn8Raw;
  TGraphErrors * grVn8RawSys;
  TGraphErrors * grVn8Theory;

  //-- Gamma1Exp
  double gamma1Exp[NCENT];
  double gamma1ExpState[NCENT];
  double gamma1ExpSyse[NCENT];
  TGraphErrors * grGamma1Exp;
  TGraphErrors * grGamma1ExpSys;
  TGraphErrors * grGamma1ExpTheory;

  //-- vn{6} / vn{4} ratio
  double ratio_vn6_vn4[NCENT];
  double ratio_vn6_vn4State[NCENT];
  double ratio_vn6_vn4Syse[NCENT];
  TGraphErrors * grvn6vn4Ratio;
  TGraphErrors * grvn6vn4RatioSys;
  TGraphErrors * grvn6vn4RatioTheory;

  TGraphErrors * grvn6vn4Ratio_Npart;
  TGraphErrors * grvn6vn4RatioSys_Npart;
  TGraphErrors * grvn6vn4Ratio_ATLASNpart;
  TGraphErrors * grvn6vn4RatioSys_ATLASNpart;

  //-- vn{8} / vn{4} ratio
  double ratio_vn8_vn4[NCENT];
  double ratio_vn8_vn4State[NCENT];
  double ratio_vn8_vn4Syse[NCENT];
  TGraphErrors * grvn8vn4Ratio;
  TGraphErrors * grvn8vn4RatioSys;

  TGraphErrors * grvn8vn4Ratio_Npart;
  TGraphErrors * grvn8vn4RatioSys_Npart;
  TGraphErrors * grvn8vn4Ratio_ATLASNpart;
  TGraphErrors * grvn8vn4RatioSys_ATLASNpart;

  //-- vn{8} / vn{6} ratio
  double ratio_vn8_vn6[NCENT];
  double ratio_vn8_vn6State[NCENT];
  double ratio_vn8_vn6Syse[NCENT];
  TGraphErrors * grvn8vn6Ratio;
  TGraphErrors * grvn8vn6RatioSys;

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

  TFile * fSmoothSys;
  TH1D * SmoothSysTotVn2;
  TH1D * SmoothSysTotVn4;
  TH1D * SmoothSysTotVn6;
  TH1D * SmoothSysTotVn8;
  TH1D * SmoothSysTotVn6Vn4;
  TH1D * SmoothSysTotVn8Vn4;
  TH1D * SmoothSysTotVn8Vn6;
  TH1D * SmoothSysTotG1E;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(20);
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(23);

  //-- Get the Analyzer output file
  fAna = new TFile( "../AnalyzerResults/CastleEbyE.root" );

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
    if( !dosys ) fUnfold = new TFile( Form("../UnfoldResults/gaussResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../UnfoldResults/gaussResp/data%i_dosys.root", norder_) );
  }
  if( studTResp ){
    if( !dosys ) fUnfold = new TFile( Form("../UnfoldResults/studTResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../UnfoldResults/studTResp/data%i_dosys.root", norder_) );
  }
  if( dataResp ){
    if( !dosys ) fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i_dosys.root", norder_) );
  }

  //-- Statistical Errors
  fStat = new TFile( Form("../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertainties_v%i.root ", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2       = (TH1D*) fStat->Get("hVarianceOfMean_Vn2");
  hVarianceOfMean_Vn4       = (TH1D*) fStat->Get("hVarianceOfMean_Vn4");
  hVarianceOfMean_Vn6       = (TH1D*) fStat->Get("hVarianceOfMean_Vn6");
  hVarianceOfMean_Vn8       = (TH1D*) fStat->Get("hVarianceOfMean_Vn8");
  hVarianceOfMean_Gamma1Exp = (TH1D*) fStat->Get("hVarianceOfMean_Gamma1Exp");
  hVarianceOfMean_Vn6Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn6Vn4");
  hVarianceOfMean_Vn8Vn4    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn4");
  hVarianceOfMean_Vn8Vn6    = (TH1D*) fStat->Get("hVarianceOfMean_Vn8Vn6");

  //-- Smooth Systematics
  fSmoothSys = new TFile("SmoothSysTot.root");
  SmoothSysTotVn2    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn2");
  SmoothSysTotVn4    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn4");
  SmoothSysTotVn6    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn6");
  SmoothSysTotVn8    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8");
  SmoothSysTotVn6Vn4 = (TH1D*) fSmoothSys->Get("SmoothSysTotVn6Vn4");
  SmoothSysTotVn8Vn4 = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8Vn4");
  SmoothSysTotVn8Vn6 = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8Vn6");
  SmoothSysTotG1E    = (TH1D*) fSmoothSys->Get("SmoothSysTotG1E");

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

    FixUnfold( hUnfold[icent][ii] );
    EbyECumu cumu(hUnfold[icent][ii]);
    double vn2  = cumu.GetCumu_vn2();
    double vn4  = cumu.GetCumu_vn4();
    double vn6  = cumu.GetCumu_vn6();
    double vn8  = cumu.GetCumu_vn8();
    double gamma1exp = cumu.GetGamma1Exp();
    double vn6vn4;
    double vn8vn4;

    std::cout << "Vn{2} = " << Form("%.2f" ,vn2) << "\tVn{4} = " << Form("%.2f" ,vn4) << "\tVn{6} = " << Form("%.2f" ,vn6) << "\tVn{8} = " << Form("%.2f" ,vn8) << std::endl;
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

    //-- Sys Error
    double centval = hCentBins.GetBinCenter(icent+1); 
    int isys = hSysBins.FindBin(centval)-1;

    double sysRelErrTotVn2;
    double sysRelErrTotVn4;
    double sysRelErrTotVn6;
    double sysRelErrTotVn8;
    double sysRelErrTotGamma1Exp;
    double sysRelErrTotVn6Vn4;
    double sysRelErrTotVn8Vn4;
    double sysRelErrTotVn8Vn6;

    if( smoothSys ){
      sysRelErrTotVn2       = SmoothSysTotVn2->GetBinContent(icent+1);
      //sysRelErrTotVn4       = SmoothSysTotVn4->GetBinContent(icent+1);
      sysRelErrTotVn4       = sysTotalV3_pt0p3_3_eta1p0_Vn4[isys];
      sysRelErrTotVn6       = SmoothSysTotVn6->GetBinContent(icent+1);
      sysRelErrTotVn8       = SmoothSysTotVn8->GetBinContent(icent+1);
      sysRelErrTotGamma1Exp = SmoothSysTotG1E->GetBinContent(icent+1);
      sysRelErrTotVn6Vn4    = SmoothSysTotVn6Vn4->GetBinContent(icent+1);
      sysRelErrTotVn8Vn4    = SmoothSysTotVn8Vn4->GetBinContent(icent+1);
      sysRelErrTotVn8Vn6    = SmoothSysTotVn8Vn6->GetBinContent(icent+1);
    }
    else{
      sysRelErrTotVn2       = sysTotalV2_pt0p3_3_eta1p0_Vn2[isys];
      sysRelErrTotVn4       = sysTotalV2_pt0p3_3_eta1p0_Vn4[isys];
      sysRelErrTotVn6       = sysTotalV2_pt0p3_3_eta1p0_Vn6[isys];
      sysRelErrTotVn8       = sysTotalV2_pt0p3_3_eta1p0_Vn8[isys];
      sysRelErrTotGamma1Exp = sysTotalV2_pt0p3_3_eta1p0_Gamma1Exp[isys];
      sysRelErrTotVn6Vn4    = sysTotalV2_pt0p3_3_eta1p0_Vn6Vn4[isys];
      sysRelErrTotVn8Vn4    = sysTotalV2_pt0p3_3_eta1p0_Vn8Vn4[isys];
      sysRelErrTotVn8Vn6    = sysTotalV2_pt0p3_3_eta1p0_Vn8Vn6[isys];
    }

    double sysErrVn2       = sysRelErrTotVn2 * vn2;
    double sysErrVn4       = sysRelErrTotVn4 * vn4;
    double sysErrVn6       = sysRelErrTotVn6 * vn6;
    double sysErrVn8       = sysRelErrTotVn8 * vn8;
    double sysErrGamma1Exp = sysRelErrTotGamma1Exp * gamma1exp;
    double sysErrVn6Vn4    = sysRelErrTotVn6Vn4 * vn6vn4;
    double sysErrVn8Vn4    = sysRelErrTotVn8Vn4 * vn8vn4;
    double sysErrVn8Vn6    = sysRelErrTotVn8Vn6 * vn8vn6;

    //-- vn{2}
    vn2Raw[icent]      = vn2;
    vn2RawState[icent] = vn2e;
    vn2RawSyse[icent]  = sysErrVn2;

    //-- vn{4}
    vn4Raw[icent]      = vn4;
    vn4RawState[icent] = vn4e;
    vn4RawSyse[icent]  = sysErrVn4;

    //-- vn{6}
    vn6Raw[icent]      = vn6;
    vn6RawState[icent] = vn6e;
    vn6RawSyse[icent]  = sysErrVn6;

    //-- vn{8}
    vn8Raw[icent]      = vn8;
    vn8RawState[icent] = vn8e;
    vn8RawSyse[icent]  = sysErrVn8;

    //-- Gamma1Exp
    gamma1Exp[icent]      = gamma1exp;
    gamma1ExpState[icent] = gamma1expe;
    gamma1ExpSyse[icent]  = sysErrGamma1Exp;

    //-- vn{6} / vn{4} ratio
    ratio_vn6_vn4[icent]      = vn6vn4;
    ratio_vn6_vn4State[icent] = vn6vn4e;
    ratio_vn6_vn4Syse[icent]  = sysErrVn6Vn4;

    //-- vn{8} / vn{4} ratio
    ratio_vn8_vn4[icent]      = vn8vn4;
    ratio_vn8_vn4State[icent] = vn8vn4e;
    ratio_vn8_vn4Syse[icent]  = sysErrVn8Vn4;

    //-- vn{8} / vn{6} ratio
    ratio_vn8_vn6[icent]      = vn8vn6;
    ratio_vn8_vn6State[icent] = vn8vn6e;
    ratio_vn8_vn6Syse[icent]  = sysErrVn8Vn6; 

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
  grVn2Raw->GetYaxis()->SetTitleOffset(1.45);
  grVn2Raw->GetYaxis()->SetDecimals(2);
  grVn2Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grVn2RawSys = new TGraphErrors(NCENT, centBinCenter, vn2Raw, centBinErr, vn2RawSyse);
  grVn2RawSys->SetLineColor(9);
  grVn2RawSys->SetMarkerColor(9);
  grVn2RawSys->SetMarkerStyle(21);
  grVn2RawSys->SetFillColor(17);
  grVn2RawSys->GetXaxis()->SetTitle( "Centrality %");
  grVn2RawSys->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn2RawSys->GetYaxis()->SetTitleOffset(1.45);
  grVn2RawSys->GetYaxis()->SetDecimals(2);
  grVn2RawSys->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{4}
  grVn4Raw = new TGraphErrors(NCENT, centBinCenter, vn4Raw, nullCentErr, vn4RawState);
  grVn4Raw->SetLineColor(kSpring+4);
  grVn4Raw->SetMarkerColor(kSpring+4);
  grVn4Raw->SetMarkerStyle(20);
  grVn4Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn4Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn4Raw->GetYaxis()->SetTitleOffset(1.45);
  grVn4Raw->GetYaxis()->SetDecimals(2);
  grVn4Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grVn4RawSys = new TGraphErrors(NCENT, centBinCenter, vn4Raw, centBinErr, vn4RawSyse);
  grVn4RawSys->SetLineColor(kSpring+4);
  grVn4RawSys->SetMarkerColor(kSpring+4);
  grVn4RawSys->SetMarkerStyle(20);
  grVn4RawSys->SetFillColor(17);
  grVn4RawSys->GetXaxis()->SetTitle( "Centrality %");
  grVn4RawSys->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn4RawSys->GetYaxis()->SetTitleOffset(1.45);
  grVn4RawSys->GetYaxis()->SetDecimals(2);
  grVn4RawSys->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{6}
  grVn6Raw = new TGraphErrors(NCENT, centBinCenter, vn6Raw, nullCentErr, vn6RawState);
  grVn6Raw->SetLineColor(6);
  grVn6Raw->SetMarkerColor(6);
  grVn6Raw->SetMarkerStyle(25);
  grVn6Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn6Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn6Raw->GetYaxis()->SetTitleOffset(1.45);
  grVn6Raw->GetYaxis()->SetDecimals(2);
  grVn6Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grVn6RawSys = new TGraphErrors(NCENT, centBinCenter, vn6Raw, centBinErr, vn6RawSyse);
  grVn6RawSys->SetLineColor(6);
  grVn6RawSys->SetMarkerColor(6);
  grVn6RawSys->SetMarkerStyle(25);
  grVn6RawSys->SetFillColor(17);
  grVn6RawSys->GetXaxis()->SetTitle( "Centrality %");
  grVn6RawSys->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn6RawSys->GetYaxis()->SetTitleOffset(1.45);
  grVn6RawSys->GetYaxis()->SetDecimals(2);
  grVn6RawSys->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- vn{8}
  grVn8Raw = new TGraphErrors(NCENT, centBinCenter, vn8Raw, nullCentErr, vn8RawState);
  grVn8Raw->SetLineColor(kOrange+7);
  grVn8Raw->SetMarkerColor(kOrange+7);
  grVn8Raw->SetMarkerStyle(28);
  grVn8Raw->GetXaxis()->SetTitle( "Centrality %");
  grVn8Raw->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn8Raw->GetYaxis()->SetTitleOffset(1.45);
  grVn8Raw->GetYaxis()->SetDecimals(2);
  grVn8Raw->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grVn8RawSys = new TGraphErrors(NCENT, centBinCenter, vn8Raw, centBinErr, vn8RawSyse);
  grVn8RawSys->SetLineColor(8);
  grVn8RawSys->SetMarkerColor(kOrange+7);
  grVn8RawSys->SetMarkerStyle(28);
  grVn8RawSys->SetFillColor(17);
  grVn8RawSys->GetXaxis()->SetTitle( "Centrality %");
  grVn8RawSys->GetYaxis()->SetTitle( Form("v_{%i}{2k}", norder_) );
  grVn8RawSys->GetYaxis()->SetTitleOffset(1.45);
  grVn8RawSys->GetYaxis()->SetDecimals(2);
  grVn8RawSys->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  //-- Gamma1Exp
  grGamma1Exp = new TGraphErrors(NCENT, centBinCenter, gamma1Exp, nullCentErr, gamma1ExpState) ;
  grGamma1Exp->SetLineColor(2);
  grGamma1Exp->SetMarkerColor(2);
  grGamma1Exp->SetMarkerStyle(20);
  grGamma1Exp->GetXaxis()->SetTitle( "Centrality %");
  grGamma1Exp->GetYaxis()->SetTitle( "#gamma_{1}^{exp}");
  grGamma1Exp->GetYaxis()->SetDecimals(1);
  grGamma1Exp->GetYaxis()->SetRangeUser(gamm1expMin, gamm1expMax);

  grGamma1ExpSys = new TGraphErrors(NCENT, centBinCenter, gamma1Exp, centBinErr, gamma1ExpSyse);
  grGamma1ExpSys->SetLineColor(17);
  grGamma1ExpSys->SetMarkerColor(17);
  grGamma1ExpSys->SetMarkerStyle(20);
  grGamma1ExpSys->SetFillColor(17);
  grGamma1ExpSys->GetYaxis()->SetDecimals(1);
  grGamma1ExpSys->GetXaxis()->SetTitle( "Centrality %");
  grGamma1ExpSys->GetYaxis()->SetTitle( "#gamma_{1}^{exp}");
  grGamma1ExpSys->GetYaxis()->SetRangeUser(gamm1expMin, gamm1expMax);

  //-- vn{6} / vn{4}
  grvn6vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn6_vn4, nullCentErr, ratio_vn6_vn4State);
  grvn6vn4Ratio->SetLineColor(4);
  grvn6vn4Ratio->SetMarkerColor(4);
  grvn6vn4Ratio->SetMarkerStyle(21);
  grvn6vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
  grvn6vn4Ratio->GetYaxis()->SetNdivisions(507);
  grvn6vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn6vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);
  grvn6vn4Ratio->GetYaxis()->SetDecimals(2);

  grvn6vn4RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn6_vn4, centBinErr, ratio_vn6_vn4Syse);
  grvn6vn4RatioSys->SetLineColor(17);
  grvn6vn4RatioSys->SetMarkerColor(17);
  grvn6vn4RatioSys->SetMarkerStyle(21);
  grvn6vn4RatioSys->SetFillColor(17);
  grvn6vn4RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4RatioSys->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
  grvn6vn4RatioSys->GetYaxis()->SetNdivisions(507);
  grvn6vn4RatioSys->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn6vn4RatioSys->GetYaxis()->SetTitleOffset(titleOffset);
  grvn6vn4RatioSys->GetYaxis()->SetDecimals(2);

  //-- vn{8} / vn{4}  
  grvn8vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn4, nullCentErr, ratio_vn8_vn4State);
  grvn8vn4Ratio->SetLineColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerStyle(34);
  grvn8vn4Ratio->SetMarkerSize(1.5);
  grvn8vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
  grvn8vn4Ratio->GetYaxis()->SetNdivisions(507);
  grvn8vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);
  grvn8vn4Ratio->GetYaxis()->SetDecimals(2);

  grvn8vn4RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn4, centBinErr, ratio_vn8_vn4Syse);
  grvn8vn4RatioSys->SetLineColor(17);
  grvn8vn4RatioSys->SetMarkerColor(17);
  grvn8vn4RatioSys->SetMarkerStyle(21);
  grvn8vn4RatioSys->SetFillColor(17);
  grvn8vn4RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4RatioSys->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
  grvn8vn4RatioSys->GetYaxis()->SetNdivisions(507);
  grvn8vn4RatioSys->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn4RatioSys->GetYaxis()->SetTitleOffset(titleOffset);
  grvn8vn4RatioSys->GetYaxis()->SetDecimals(2);

  //-- vn{8} / vn{6}
  grvn8vn6Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn6, nullCentErr, ratio_vn8_vn6State);
  grvn8vn6Ratio->SetLineColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerStyle(33);
  grvn8vn6Ratio->SetMarkerSize(1.5);
  grvn8vn6Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6Ratio->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_) );
  grvn8vn6Ratio->GetYaxis()->SetNdivisions(507);
  grvn8vn6Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn6Ratio->GetYaxis()->SetTitleOffset(titleOffset+0.25);
  grvn8vn6Ratio->GetYaxis()->SetDecimals(3);

  grvn8vn6RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn6, centBinErr, ratio_vn8_vn6Syse);
  grvn8vn6RatioSys->SetLineColor(17);
  grvn8vn6RatioSys->SetMarkerColor(17);
  grvn8vn6RatioSys->SetMarkerStyle(21);
  grvn8vn6RatioSys->SetFillColor(17);
  grvn8vn6RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6RatioSys->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_) );
  grvn8vn6RatioSys->GetYaxis()->SetNdivisions(507);
  grvn8vn6RatioSys->GetYaxis()->SetRangeUser(ratioMinVn8Vn6,ratioMaxVn8Vn6);
  grvn8vn6RatioSys->GetYaxis()->SetTitleOffset(titleOffset+0.25);
  grvn8vn6RatioSys->GetYaxis()->SetDecimals(3);

  //-- DRAW!

  TLegend * legCumu = new TLegend(0.38, 0.21, 0.59, 0.33);
  legCumu->SetBorderSize(0);
  legCumu->SetFillStyle(0);
  legCumu->AddEntry(grVn2Raw, "k = 1", "lp");
  legCumu->AddEntry(grVn4Raw, "k = 2", "lp");
  //legCumu->AddEntry(grVn6Raw, "k = 3", "lp");
  //legCumu->AddEntry(grVn8Raw, "k = 4", "lp");

  //-- PAS Fig. 2
  TCanvas * cCumuRaw = new TCanvas("cCumuRaw", "cCumuRaw", 500, 500);
  cCumuRaw->cd();
  cCumuRaw->SetTopMargin(0.1);
  cCumuRaw->SetLeftMargin(0.2);
  cCumuRaw->SetRightMargin(0.1);
  grVn4RawSys->Draw("apE2");
  grVn4Raw->Draw("psame");
  grVn2RawSys->Draw("psameE2");
  grVn2Raw->Draw("psame");
  //grVn6RawSys->Draw("psameE2");
  //grVn6Raw->Draw("psame");
  //grVn8RawSys->Draw("psameE2");
  //grVn8Raw->Draw("psame");
  legCumu->Draw("same");
  legCumu->SetTextFont(43); 
  legCumu->SetTextSize(23);
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.24, 0.76, Form("|#eta| < %.1f", tkEta));
  //latex2.DrawLatex(0.57, 0.54, "#bf{v_{2}{2k}}");
  cCumuRaw->SaveAs("../plots/skew/SysCumuRaw.pdf");

  TLegend * legg1e = new TLegend(0.22, 0.22, 0.54, 0.375);
  legInit( legg1e );
  legg1e->AddEntry(grGamma1Exp,       "#gamma_{1}^{exp}", "lp");

  //-- PAS Fig. 4
  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  cGamma1Exp->SetTopMargin(0.1);
  cGamma1Exp->SetLeftMargin(0.2);
  cGamma1Exp->SetRightMargin(0.1);
  grGamma1ExpSys->Draw("apE2");
  grGamma1ExpSys->Draw("pE2same");
  grGamma1Exp->Draw("psame");
  legg1e->Draw("same");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.72, 0.76, Form("|#eta| < %.1f", tkEta));

  legg1e->SetTextFont(43);
  legg1e->SetTextSize(26);
  cGamma1Exp->Update();

  cGamma1Exp->SaveAs("../plots/skew/SysGamma1Exp.pdf");

  TLine * line = new TLine(centbinsDefault[0], 1.0, grvn6vn4RatioSys->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLine * line2 = new TLine(centbinsDefault[0], 1.0, grvn8vn4Ratio->GetXaxis()->GetXmax(), 1.0);
  line2->SetLineColor(1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);

  TLegend * leg64 = new TLegend(0.22, 0.21, 0.54, 0.335);
  legInit( leg64 );
  leg64->AddEntry(grvn6vn4Ratio,       "v_{2}{6} / v_{2}{4}", "lp");

  TLegend * leg84 = new TLegend(0.22, 0.21, 0.54, 0.26);
  legInit( leg84 );
  leg84->AddEntry(grvn8vn4Ratio, "v_{2}{8} / v_{2}{4}", "lp");

  TLegend * leg86 = new TLegend(0.22, 0.21, 0.54, 0.26);
  legInit( leg86 );
  leg86->AddEntry(grvn8vn6Ratio, "v_{2}{8} / v_{2}{6}", "lp");

  //--   PAS Fig. 3
  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  
  cCumuRatio->cd(1);
  cCumuRatio->cd(1)->SetLeftMargin(0.19);
  cCumuRatio->cd(1)->SetRightMargin(0.07);
  cCumuRatio->cd(1)->SetTopMargin(0.1); 
  grvn6vn4RatioSys->Draw("apE2");
  grvn6vn4RatioSys->Draw("pE2same");
  grvn6vn4Ratio->Draw("psame");
  line->Draw("same");
  leg64->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatio->cd(2);
  cCumuRatio->cd(2)->SetLeftMargin(0.19);
  cCumuRatio->cd(2)->SetRightMargin(0.07);
  cCumuRatio->cd(2)->SetTopMargin(0.1);
  grvn8vn4RatioSys->Draw("apE2");
  grvn8vn4Ratio->Draw("psame");
  line2->Draw("same");
  leg84->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatio->cd(3);
  cCumuRatio->cd(3)->SetLeftMargin(0.19);
  cCumuRatio->cd(3)->SetRightMargin(0.07);
  cCumuRatio->cd(3)->SetTopMargin(0.1);
  grvn8vn6RatioSys->Draw("apE2");
  grvn8vn6Ratio->Draw("psame");
  line2->Draw("same");
  leg86->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));

  leg64->SetTextFont(43);
  leg64->SetTextSize(26);

  leg84->SetTextFont(43);
  leg84->SetTextSize(26);

  leg86->SetTextFont(43);
  leg86->SetTextSize(26);

  cCumuRatio->Update();
  cCumuRatio->SaveAs("../plots/skew/SysCumuRatio.pdf");

  //-- ATLAS Comparison
  if( compATLAS ){

    ratioMin = 0.9;
    ratioMax = 1.03999;

    double npartSysWidth[NCENT+1];
    for(int icent = 0; icent < NCENT+1; icent++) npartSysWidth[icent] = 5;

    //-- CMS 5.02 TeV, 0.3 < pT < 3.0 GeV, |eta| < 2.4
    //-- vn{6} / vn{4}
    grvn6vn4Ratio_Npart = (TGraphErrors*) grvn6vn4Ratio->Clone("grvn6vn4Ratio_Npart");
    grvn6vn4Ratio_Npart->SetLineColor(2);
    grvn6vn4Ratio_Npart->SetMarkerColor(2);
    grvn6vn4Ratio_Npart->SetMarkerStyle(20);

    grvn6vn4RatioSys_Npart = (TGraphErrors*) grvn6vn4RatioSys->Clone("grvn6vn4RatioSys_Npart");
    grvn6vn4RatioSys_Npart->SetLineColor(16);
    grvn6vn4RatioSys_Npart->SetMarkerColor(16);
    grvn6vn4RatioSys_Npart->SetMarkerStyle(21);
    grvn6vn4RatioSys_Npart->SetFillColorAlpha(46, 0.6);

    //-- vn{8} / vn{4}
    grvn8vn4Ratio_Npart = (TGraphErrors*) grvn8vn4Ratio->Clone("grvn8vn4Ratio_Npart"); 
    grvn8vn4Ratio_Npart->SetLineColor(2);
    grvn8vn4Ratio_Npart->SetMarkerColor(2);
    grvn8vn4Ratio_Npart->SetMarkerStyle(21);
    grvn8vn4Ratio_Npart->SetMarkerSize(1.0);

    grvn8vn4RatioSys_Npart = (TGraphErrors*) grvn8vn4RatioSys->Clone("grvn8vn4RatioSys_Npart");
    grvn8vn4RatioSys_Npart->SetLineColor(16);
    grvn8vn4RatioSys_Npart->SetMarkerColor(16);
    grvn8vn4RatioSys_Npart->SetMarkerStyle(21);
    grvn8vn4RatioSys_Npart->SetFillColorAlpha(46, 0.6);
    grvn8vn4RatioSys_Npart->SetMarkerSize(1.0);

    //-- ATLAS 2.76 TeV, 0.5 < pT < 20 GeV, |eta| < 2.5
    //-- vn{6} / vn{4}
    Double_t ATLASVn6Vn4_xval[16] = {
      2.5,      3.5,      4.5,      7.5,      12.5,
      17.5,     22.5,     27.5,     32.5,     37.5,
      42.5,     47.5,     52.5,     57.5,     62.5,
      67.5};
    Double_t ATLASVn6Vn4_yval[16] = {
      0.83,        0.9584,      0.971,       0.997,       0.9967,
      0.9975,      0.9973,      0.9966,      0.9953,      0.9939,
      0.9918,      0.9881,      0.9836,      0.9802,      0.9815,
      0.985};
    Double_t ATLASVn6Vn4_xerr[16] = {
      0,      0,      0,      0,      0,
      0,      0,      0,      0,      0,
      0,      0,      0,      0,      0,
      0};
    Double_t ATLASVn6Vn4_ystatsys[16] = {
      0.2445,        0.05164,       0.02158,       0.007417,      0.00464,
      0.001148,      0.0002847,     0.001802,      0.004732,      0.002289,
      0.004418,      0.003787,      0.004197,      0.01227,       0.01283,
      0.003727};
    int ATLASVn6Vn4_numpoints = 16;


    grvn6vn4Ratio_ATLASNpart = new TGraphErrors(ATLASVn6Vn4_numpoints, ATLASVn6Vn4_xval, ATLASVn6Vn4_yval, ATLASVn6Vn4_xerr, ATLASVn6Vn4_ystatsys);
    grvn6vn4Ratio_ATLASNpart->SetLineColor(1);
    grvn6vn4Ratio_ATLASNpart->SetMarkerColor(1);
    grvn6vn4Ratio_ATLASNpart->SetMarkerStyle(20);
    grvn6vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "Centrality %");
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetTitleOffset(titleOffset);
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetDecimals(2);

    //-- vn{8} / vn{4}
    Double_t ATLASVn8Vn4_xval[16] = {
      2.5,      3.5,      4.5,      7.5,      12.5,
      17.5,     22.5,     27.5,     32.5,     37.5,
      42.5,     47.5,     52.5,     57.5,     62.5,
      67.5};
    Double_t ATLASVn8Vn4_yval[16] = {
      0.6779,      0.9669,      0.9708,      0.9989,      0.9968,
      0.9974,      0.997,       0.9963,      0.9949,      0.9932,
      0.9909,      0.9864,      0.9807,      0.9762,      0.9779,
      0.9823};
    Double_t ATLASVn8Vn4_xerr[16] = {
      0,      0,      0,      0,      0,
      0,      0,      0,      0,      0,
      0,      0,      0,      0,      0,
      0};
    Double_t ATLASVn8Vn4_ystatsys[16] = {
      0.2445,        0.05164,      0.02158,      0.007417,      0.00464,
      0.001148,      0.0002847,    0.001802,     0.004732,      0.002289,
      0.004418,      0.003787,     0.004197,     0.01227,       0.01283,
      0.003727};
    int ATLASVn8Vn4_numpoints = 16;

    grvn8vn4Ratio_ATLASNpart = new TGraphErrors(ATLASVn8Vn4_numpoints, ATLASVn8Vn4_xval, ATLASVn8Vn4_yval, ATLASVn8Vn4_xerr, ATLASVn8Vn4_ystatsys);
    grvn8vn4Ratio_ATLASNpart->SetLineColor(1);
    grvn8vn4Ratio_ATLASNpart->SetMarkerColor(1);
    grvn8vn4Ratio_ATLASNpart->SetMarkerStyle(21);
    grvn8vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "Centrality %");
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetTitleOffset(titleOffset);
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetDecimals(2);

    TLegend * leg64a = new TLegend(0.27, 0.22, 0.69, 0.34);
    leg64a->SetFillStyle(0);
    leg64a->SetBorderSize(0);
    leg64a->AddEntry(grvn6vn4Ratio_Npart,      "5.02 TeV CMS",   "lp");
    leg64a->AddEntry(grvn6vn4Ratio_ATLASNpart, "2.76 TeV ATLAS", "lp");

    TLegend * leg84a = new TLegend(0.27, 0.22, 0.69, 0.34);
    leg84a->SetFillStyle(0);
    leg84a->SetBorderSize(0);
    leg84a->AddEntry(grvn8vn4Ratio_Npart,      "5.02 TeV CMS",   "lp");
    leg84a->AddEntry(grvn8vn4Ratio_ATLASNpart, "2.76 TeV ATLAS", "lp");

    TLine * lone = new TLine(grvn6vn4Ratio_ATLASNpart->GetXaxis()->GetXmin(), 1.0, grvn6vn4Ratio_ATLASNpart->GetXaxis()->GetXmax(), 1.0);
    lone->SetLineStyle(2);
    lone->SetLineWidth(2);

    //-- PAS Fig. 5
    TCanvas * cATLASComp = new TCanvas("cATLASComp", "cATLASComp", 1000, 500);
    cATLASComp->Divide(2,1);

    cATLASComp->cd(1);
    cATLASComp->cd(1)->SetLeftMargin(0.19);
    cATLASComp->cd(1)->SetRightMargin(0.09);
    cATLASComp->cd(1)->SetTopMargin(0.1);
    grvn6vn4Ratio_ATLASNpart->Draw("ap");
    grvn6vn4RatioSys_Npart->Draw("pE2same");
    grvn6vn4Ratio_Npart->Draw("psame");
    lone->Draw("same");
    leg64a->Draw("same");
    latex.DrawLatex(0.19, 0.915, "#bf{CMS} #it{Preliminary}");
    latex.DrawLatex(0.65, 0.915, "PbPb 5.02 TeV");
    latex2.DrawLatex(0.47, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
    latex2.DrawLatex(0.715, 0.76, Form("|#eta| < %.1f", tkEta));

    cATLASComp->cd(2);
    cATLASComp->cd(2)->SetLeftMargin(0.19);
    cATLASComp->cd(2)->SetRightMargin(0.09);
    cATLASComp->cd(2)->SetTopMargin(0.1);
    grvn8vn4Ratio_ATLASNpart->Draw("ap");
    grvn8vn4RatioSys_Npart->Draw("pE2same");
    grvn8vn4Ratio_Npart->Draw("psame");
    lone->Draw("same");
    leg84a->Draw("same");
    latex.DrawLatex(0.19, 0.915, "#bf{CMS} #it{Preliminary}");
    latex.DrawLatex(0.65, 0.915, "PbPb 5.02 TeV");
    latex2.DrawLatex(0.47, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
    latex2.DrawLatex(0.715, 0.76, Form("|#eta| < %.1f", tkEta));


    leg64a->SetTextFont(43);
    leg64a->SetTextSize(23);
    leg84a->SetTextFont(43);
    leg84a->SetTextSize(23);

    cATLASComp->Update();
    cATLASComp->SaveAs("../plots/skew/cSysATLASComp.pdf");

  }

  //-- Save results to file
  fOut = new TFile("PhysicsResults.root", "recreate");
  fOut->cd();
  grVn2Raw->Write("grVn2Raw");
  grVn2RawSys->Write("grVn2RawSys");
  grVn4Raw->Write("grVn4Raw");
  grVn4RawSys->Write("grVn4RawSys");
  grVn6Raw->Write("grVn6Raw");
  grVn6RawSys->Write("grVn6RawSys");
  grVn8Raw->Write("grVn8Raw");
  grVn8RawSys->Write("grVn8RawSys");
  grvn6vn4Ratio->Write("grvn6vn4Ratio");
  grvn6vn4RatioSys->Write("grvn6vn4RatioSys");
  grvn8vn4Ratio->Write("grvn8vn4Ratio");
  grvn8vn4RatioSys->Write("grvn8vn4RatioSys");
  grvn8vn6Ratio->Write("grvn8vn6Ratio");
  grvn8vn6RatioSys->Write("grvn8vn6RatioSys");
  grGamma1Exp->Write("grGamma1Exp");
  grGamma1ExpSys->Write("grGamma1ExpSys");


}

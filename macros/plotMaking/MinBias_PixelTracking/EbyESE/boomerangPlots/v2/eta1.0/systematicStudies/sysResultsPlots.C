#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void sysResultsPlots(){

  TH1D::SetDefaultSumw2();

  int centbin     = 4;
  int norder_     = 2;
  double tkEta    = 1.0;

  bool dosys     = 0;
  bool smoothSys = 1;
  bool compATLAS = 1;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double cumuMin = 0.0;
  double cumuMax = 0.12999;

  double gamm1expMin = -1.0;
  double gamm1expMax = 0.5;
  double ratioMin    = 0.95;
  double ratioMax    = 1.01999;
  double ratioMinVn8Vn6 = 0.994;
  double ratioMaxVn8Vn6 = 1.003;
  double vn46_vn68Min = 0.;
  double vn46_vn68Max = 14.;

  double sysPctMin = 0.00001;
  double sysPctMax = 50.;

  double titleOffset = 1.4;

  TLatex latex;
  TLatex latex2;
  TLatex latex3;

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

  //-- vn{6}, vn{8} splitting prediction
  double ratio_vn46_vn68[NCENT];
  double ratio_vn46_vn68State[NCENT];
  double ratio_vn46_vn68Syse[NCENT];
  TGraphErrors * grvn46_vn68Ratio;
  TGraphErrors * grvn46_vn68RatioSys;

  double trentoPm1_en46_en68[NCENT];
  TGraph * grTrentoPm1_en46_en68;

  double trentoP0_en46_en68[NCENT];
  TGraph * grTrentoP0_en46_en68;

  double trentoP1_en46_en68[NCENT];
  TGraph * grTrentoP1_en46_en68;

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
  TH1D * hVarianceOfMean_Vn46_Vn68;
  TH1D * hVarianceOfMean_Bin[NCENT];

  TFile * fSmoothSys;
  TH1D * SmoothSysTotVn2;
  TH1D * SmoothSysTotVn4;
  TH1D * SmoothSysTotVn6;
  TH1D * SmoothSysTotVn8;
  TH1D * SmoothSysTotVn6Vn4;
  TH1D * SmoothSysTotVn8Vn4;
  TH1D * SmoothSysTotVn8Vn6;
  TH1D * SmoothSysTotVn46_Vn68;
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
  latex3.SetNDC();
  latex3.SetTextFont(43);
  latex3.SetTextSize(26);


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

  if( !dosys ) fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i.root", norder_) );
  else         fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i_dosys.root", norder_) );

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
  hVarianceOfMean_Vn46_Vn68 = (TH1D*) fStat->Get("hVarianceOfMean_Vn46_Vn68");
  for(int icent = 0; icent < NCENT; icent++) hVarianceOfMean_Bin[icent] = (TH1D*) fStat->Get( Form("hVarianceOfMean_Bin_c%i", icent) );

  //-- Smooth Systematics
  fSmoothSys = new TFile("SmoothSysTot.root");
  SmoothSysTotVn2       = (TH1D*) fSmoothSys->Get("SmoothSysTotVn2");
  SmoothSysTotVn4       = (TH1D*) fSmoothSys->Get("SmoothSysTotVn4");
  SmoothSysTotVn6       = (TH1D*) fSmoothSys->Get("SmoothSysTotVn6");
  SmoothSysTotVn8       = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8");
  SmoothSysTotVn6Vn4    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn6Vn4");
  SmoothSysTotVn8Vn4    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8Vn4");
  SmoothSysTotVn8Vn6    = (TH1D*) fSmoothSys->Get("SmoothSysTotVn8Vn6");
  SmoothSysTotVn46_Vn68 = (TH1D*) fSmoothSys->Get("SmoothSysTotVn46_Vn68");
  SmoothSysTotG1E       = (TH1D*) fSmoothSys->Get("SmoothSysTotG1E");

  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    for(int i = 0; i < NITER; i++){

      //-- Get the (re)unfolded histograms
      hUnfold[icent][i] = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );

      //-- Chi squares
      double chi2NDF_Refold = hRefold[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF"); 

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
    
    //-- Truncate, update stat errors, and update colors
    FixUnfold( hUnfold[icent][ii] );
    UpdateUnfoldBinErrors( hUnfold[icent][ii], hVarianceOfMean_Bin[icent]);
    hUnfold[icent][ii]->SetLineColor( centCol[icent] );
    hUnfold[icent][ii]->SetMarkerColor( centCol[icent] );

    //-- Calculate Cumulants
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

    double vn46_vn68;
    if(vn4 == 0 || vn6 == 0 || vn8 == 0) vn46_vn68 = -1000.;
    else                                 vn46_vn68 = (vn4 - vn6)/(vn6 - vn8);
    std::cout<<"vn46_vn68 = "<< vn46_vn68<<std::endl;
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
    double vn46_vn68e = sqrt( hVarianceOfMean_Vn46_Vn68->GetBinContent(icent+1) );

    //-- Sys Error
    double centval = hCentBins.GetBinCenter(icent+1); 
    //int isys = hSysBins.FindBin(centval)-1;

    double sysRelErrTotVn2;
    double sysRelErrTotVn4;
    double sysRelErrTotVn6;
    double sysRelErrTotVn8;
    double sysRelErrTotGamma1Exp;
    double sysRelErrTotVn6Vn4;
    double sysRelErrTotVn8Vn4;
    double sysRelErrTotVn8Vn6;
    double sysRelErrTotVn46_Vn68;

    if( smoothSys ){
      sysRelErrTotVn2       = SmoothSysTotVn2->GetBinContent(icent+1);
      sysRelErrTotVn4       = SmoothSysTotVn4->GetBinContent(icent+1);
      sysRelErrTotVn6       = SmoothSysTotVn6->GetBinContent(icent+1);
      sysRelErrTotVn8       = SmoothSysTotVn8->GetBinContent(icent+1);
      sysRelErrTotGamma1Exp = SmoothSysTotG1E->GetBinContent(icent+1);
      sysRelErrTotVn6Vn4    = SmoothSysTotVn6Vn4->GetBinContent(icent+1);
      sysRelErrTotVn8Vn4    = SmoothSysTotVn8Vn4->GetBinContent(icent+1);
      sysRelErrTotVn8Vn6    = SmoothSysTotVn8Vn6->GetBinContent(icent+1);
      sysRelErrTotVn46_Vn68 = SmoothSysTotVn46_Vn68->GetBinContent(icent+1);
    }
    else{
      /*
      sysRelErrTotVn2       = sysTotalV2_pt0p3_3_eta1p0_Vn2[isys];
      sysRelErrTotVn4       = sysTotalV2_pt0p3_3_eta1p0_Vn4[isys];
      sysRelErrTotVn6       = sysTotalV2_pt0p3_3_eta1p0_Vn6[isys];
      sysRelErrTotVn8       = sysTotalV2_pt0p3_3_eta1p0_Vn8[isys];
      sysRelErrTotGamma1Exp = sysTotalV2_pt0p3_3_eta1p0_Gamma1Exp[isys];
      sysRelErrTotVn6Vn4    = sysTotalV2_pt0p3_3_eta1p0_Vn6Vn4[isys];
      sysRelErrTotVn8Vn4    = sysTotalV2_pt0p3_3_eta1p0_Vn8Vn4[isys];
      sysRelErrTotVn8Vn6    = sysTotalV2_pt0p3_3_eta1p0_Vn8Vn6[isys];
      sysRelErrTotVn46_Vn68 = sysTotalV2_pt0p3_3_eta1p0_Vn46_Vn68[isys];
      */
    }

    double sysErrVn2       = sysRelErrTotVn2 * vn2;
    double sysErrVn4       = sysRelErrTotVn4 * vn4;
    double sysErrVn6       = sysRelErrTotVn6 * vn6;
    double sysErrVn8       = sysRelErrTotVn8 * vn8;
    double sysErrGamma1Exp = sysRelErrTotGamma1Exp * gamma1exp;
    double sysErrVn6Vn4    = sysRelErrTotVn6Vn4 * vn6vn4;
    double sysErrVn8Vn4    = sysRelErrTotVn8Vn4 * vn8vn4;
    double sysErrVn8Vn6    = sysRelErrTotVn8Vn6 * vn8vn6;
    double sysErrVn46_Vn68 = sysRelErrTotVn46_Vn68 * vn46_vn68;

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

    //-- vn{6}, vn{8} splitting prediction 
    ratio_vn46_vn68[icent]      = vn46_vn68;
    ratio_vn46_vn68State[icent] = vn46_vn68e;
    ratio_vn46_vn68Syse[icent]  = sysErrVn46_Vn68; 

  }

  double nullCentErr[NCENT];
  for(int i = 0; i < NCENT; i++) nullCentErr[i] = 0.;

  //-- Initialize cumu graphs

  //-- vn{2}
  grVn2Raw = new TGraphErrors(NCENT, centBinCenter, vn2Raw, nullCentErr, vn2RawState);
  grVn2Raw->SetLineColor(9);
  grVn2Raw->SetMarkerColor(9);
  grVn2Raw->SetMarkerStyle(21);

  grVn2RawSys = new TGraphErrors(NCENT, centBinCenter, vn2Raw, centBinErr, vn2RawSyse);
  grVn2RawSys->SetLineColor(9);
  grVn2RawSys->SetMarkerColor(9);
  grVn2RawSys->SetMarkerStyle(21);
  grVn2RawSys->SetFillColor(17);

  //-- vn{4}
  grVn4Raw = new TGraphErrors(NCENT, centBinCenter, vn4Raw, nullCentErr, vn4RawState);
  grVn4Raw->SetLineColor(kSpring+4);
  grVn4Raw->SetMarkerColor(kSpring+4);
  grVn4Raw->SetMarkerStyle(20);

  grVn4RawSys = new TGraphErrors(NCENT, centBinCenter, vn4Raw, centBinErr, vn4RawSyse);
  grVn4RawSys->SetLineColor(kSpring+4);
  grVn4RawSys->SetMarkerColor(kSpring+4);
  grVn4RawSys->SetMarkerStyle(20);
  grVn4RawSys->SetFillColor(17);

  //-- vn{6}
  grVn6Raw = new TGraphErrors(NCENT, centBinCenter, vn6Raw, nullCentErr, vn6RawState);
  grVn6Raw->SetLineColor(6);
  grVn6Raw->SetMarkerColor(6);
  grVn6Raw->SetMarkerStyle(25);

  grVn6RawSys = new TGraphErrors(NCENT, centBinCenter, vn6Raw, centBinErr, vn6RawSyse);
  grVn6RawSys->SetLineColor(6);
  grVn6RawSys->SetMarkerColor(6);
  grVn6RawSys->SetMarkerStyle(25);
  grVn6RawSys->SetFillColor(17);

  //-- vn{8}
  grVn8Raw = new TGraphErrors(NCENT, centBinCenter, vn8Raw, nullCentErr, vn8RawState);
  grVn8Raw->SetLineColor(kOrange+7);
  grVn8Raw->SetMarkerColor(kOrange+7);
  grVn8Raw->SetMarkerStyle(28);

  grVn8RawSys = new TGraphErrors(NCENT, centBinCenter, vn8Raw, centBinErr, vn8RawSyse);
  grVn8RawSys->SetLineColor(kOrange+7);
  grVn8RawSys->SetMarkerColor(kOrange+7);
  grVn8RawSys->SetMarkerStyle(28);
  grVn8RawSys->SetFillColor(17);

  //-- Gamma1Exp
  grGamma1Exp = new TGraphErrors(NCENT, centBinCenter, gamma1Exp, nullCentErr, gamma1ExpState) ;
  grGamma1Exp->SetLineColor(2);
  grGamma1Exp->SetMarkerColor(2);
  grGamma1Exp->SetMarkerStyle(20);

  grGamma1ExpSys = new TGraphErrors(NCENT, centBinCenter, gamma1Exp, centBinErr, gamma1ExpSyse);
  grGamma1ExpSys->SetLineColor(2);
  grGamma1ExpSys->SetMarkerColor(2);
  grGamma1ExpSys->SetMarkerStyle(20);
  grGamma1ExpSys->SetFillColor(17);

  //-- vn{6} / vn{4}
  grvn6vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn6_vn4, nullCentErr, ratio_vn6_vn4State);
  grvn6vn4Ratio->SetLineColor(4);
  grvn6vn4Ratio->SetMarkerColor(4);
  grvn6vn4Ratio->SetMarkerStyle(21);
  grvn6vn4Ratio->SetMarkerSize(1.2);
  grvn6vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4Ratio->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{6} / #font[12]{v}_{%i}{4}", norder_, norder_) );
  grvn6vn4Ratio->GetYaxis()->SetNdivisions(507);
  grvn6vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn6vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);
  grvn6vn4Ratio->GetYaxis()->SetDecimals(2);

  grvn6vn4RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn6_vn4, centBinErr, ratio_vn6_vn4Syse);
  grvn6vn4RatioSys->SetLineColor(4);
  grvn6vn4RatioSys->SetMarkerColor(4);
  grvn6vn4RatioSys->SetMarkerStyle(21);
  grvn6vn4RatioSys->SetFillColor(17);
  grvn6vn4RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn6vn4RatioSys->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{6} / #font[12]{v}_{%i}{4}", norder_, norder_) );
  grvn6vn4RatioSys->GetYaxis()->SetNdivisions(507);
  grvn6vn4RatioSys->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn6vn4RatioSys->GetYaxis()->SetTitleOffset(titleOffset);
  grvn6vn4RatioSys->GetYaxis()->SetDecimals(2);

  //-- vn{8} / vn{4}  
  grvn8vn4Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn4, nullCentErr, ratio_vn8_vn4State);
  grvn8vn4Ratio->SetLineColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerColor(kGreen+2);
  grvn8vn4Ratio->SetMarkerStyle(21);
  grvn8vn4Ratio->SetMarkerSize(1.2);
  grvn8vn4Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4Ratio->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{8} / #font[12]{v}_{%i}{4}", norder_, norder_) );
  grvn8vn4Ratio->GetYaxis()->SetNdivisions(507);
  grvn8vn4Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn4Ratio->GetYaxis()->SetTitleOffset(titleOffset);
  grvn8vn4Ratio->GetYaxis()->SetDecimals(2);

  grvn8vn4RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn4, centBinErr, ratio_vn8_vn4Syse);
  grvn8vn4RatioSys->SetLineColor(kGreen+2);
  grvn8vn4RatioSys->SetMarkerColor(kGreen+2);
  grvn8vn4RatioSys->SetMarkerStyle(21);
  grvn8vn4RatioSys->SetFillColor(17);
  grvn8vn4RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn4RatioSys->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{8} / #font[12]{v}_{%i}{4}", norder_, norder_) );
  grvn8vn4RatioSys->GetYaxis()->SetNdivisions(507);
  grvn8vn4RatioSys->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn4RatioSys->GetYaxis()->SetTitleOffset(titleOffset);
  grvn8vn4RatioSys->GetYaxis()->SetDecimals(2);

  //-- vn{8} / vn{6}
  grvn8vn6Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn6, nullCentErr, ratio_vn8_vn6State);
  grvn8vn6Ratio->SetLineColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerColor(kViolet-1);
  grvn8vn6Ratio->SetMarkerStyle(21);
  grvn8vn6Ratio->SetMarkerSize(1.2);
  grvn8vn6Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6Ratio->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{8} / #font[12]{v}_{%i}{6}", norder_, norder_) );
  grvn8vn6Ratio->GetYaxis()->SetNdivisions(507);
  grvn8vn6Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grvn8vn6Ratio->GetYaxis()->SetTitleOffset(titleOffset+0.25);
  grvn8vn6Ratio->GetYaxis()->SetDecimals(3);

  grvn8vn6RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn8_vn6, centBinErr, ratio_vn8_vn6Syse);
  grvn8vn6RatioSys->SetLineColor(kViolet-1);
  grvn8vn6RatioSys->SetMarkerColor(kViolet-1);
  grvn8vn6RatioSys->SetMarkerStyle(21);
  grvn8vn6RatioSys->SetFillColor(17);
  grvn8vn6RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn8vn6RatioSys->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{8} / #font[12]{v}_{%i}{6}", norder_, norder_) );
  grvn8vn6RatioSys->GetYaxis()->SetNdivisions(507);
  grvn8vn6RatioSys->GetYaxis()->SetRangeUser(ratioMinVn8Vn6,ratioMaxVn8Vn6);
  grvn8vn6RatioSys->GetYaxis()->SetTitleOffset(titleOffset+0.25);
  grvn8vn6RatioSys->GetYaxis()->SetDecimals(3);

  //-- vn{6}, vn{8} splitting prediction  
  grvn46_vn68Ratio = new TGraphErrors(NCENT, centBinCenter, ratio_vn46_vn68, nullCentErr, ratio_vn46_vn68State);
  grvn46_vn68Ratio->SetLineColor(1);
  grvn46_vn68Ratio->SetMarkerColor(1);
  grvn46_vn68Ratio->SetMarkerStyle(20);
  grvn46_vn68Ratio->GetXaxis()->SetTitle( "Centrality %");
  grvn46_vn68Ratio->GetYaxis()->SetTitle( Form("(#font[12]{v}_{%i}{4} - #font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6} - #font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_) );
  grvn46_vn68Ratio->GetYaxis()->SetTitleOffset(1.45);
  grvn46_vn68Ratio->GetYaxis()->SetDecimals(2);
  grvn46_vn68Ratio->GetYaxis()->SetRangeUser(vn46_vn68Min, vn46_vn68Max);

  grvn46_vn68RatioSys = new TGraphErrors(NCENT, centBinCenter, ratio_vn46_vn68, centBinErr, ratio_vn46_vn68Syse);
  grvn46_vn68RatioSys->SetLineColor(1);
  grvn46_vn68RatioSys->SetMarkerColor(1);
  grvn46_vn68RatioSys->SetMarkerStyle(20);
  grvn46_vn68RatioSys->SetFillColor(17);
  grvn46_vn68RatioSys->GetXaxis()->SetTitle( "Centrality %");
  grvn46_vn68RatioSys->GetYaxis()->SetTitle( Form("(#font[12]{v}_{%i}{4} - #font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6} - #font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_) );
  grvn46_vn68RatioSys->GetYaxis()->SetTitleOffset(1.45);
  grvn46_vn68RatioSys->GetYaxis()->SetDecimals(2);
  grvn46_vn68RatioSys->GetYaxis()->SetRangeUser(vn46_vn68Min, vn46_vn68Max);

  trentoPm1_en46_en68[0]  = -1;
  trentoPm1_en46_en68[1]  = -1;
  trentoPm1_en46_en68[2]  = 12.12903225806452;
  trentoPm1_en46_en68[3]  = 9.827956989247312;
  trentoPm1_en46_en68[4]  = 9.397849462365592;
  trentoPm1_en46_en68[5]  = 8.967741935483872;
  trentoPm1_en46_en68[6]  = 8.623655913978496;
  trentoPm1_en46_en68[7]  = 8.258064516129034;
  trentoPm1_en46_en68[8]  = 8.021505376344088;
  trentoPm1_en46_en68[9]  = 7.67741935483871;
  trentoPm1_en46_en68[10] = 7.376344086021506;
  trentoPm1_en46_en68[11] = 7.118279569892474;
  grTrentoPm1_en46_en68 = new TGraph(NCENT, centBinCenter, trentoPm1_en46_en68);
  grTrentoPm1_en46_en68->SetLineColor(kViolet-1);
  grTrentoPm1_en46_en68->SetMarkerColor(kViolet-1);
  grTrentoPm1_en46_en68->SetMarkerStyle(26);
  grTrentoPm1_en46_en68->GetXaxis()->SetTitle( "Centrality %");
  grTrentoPm1_en46_en68->GetYaxis()->SetTitle( Form("(#font[12]{v}_{%i}{4} - #font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6} - #font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_) );
  grTrentoPm1_en46_en68->GetYaxis()->SetTitleOffset(1.45);
  grTrentoPm1_en46_en68->GetYaxis()->SetDecimals(2);
  grTrentoPm1_en46_en68->GetYaxis()->SetRangeUser(vn46_vn68Min, vn46_vn68Max);

  trentoP0_en46_en68[0]  = -1;
  trentoP0_en46_en68[1]  = -1;
  trentoP0_en46_en68[2]  = 13.075268817204304;
  trentoP0_en46_en68[3]  = 10.451612903225808;
  trentoP0_en46_en68[4]  = 9.505376344086024;
  trentoP0_en46_en68[5]  = 8.9247311827957;
  trentoP0_en46_en68[6]  = 8.451612903225808;
  trentoP0_en46_en68[7]  = 8.150537634408604;
  trentoP0_en46_en68[8]  = 7.827956989247312;
  trentoP0_en46_en68[9]  = 7.526881720430108;
  trentoP0_en46_en68[10] = 7.24731182795699;
  trentoP0_en46_en68[11] = 6.946236559139786;
  grTrentoP0_en46_en68 = new TGraph(NCENT, centBinCenter, trentoP0_en46_en68);
  grTrentoP0_en46_en68->SetLineColor(kCyan+2);
  grTrentoP0_en46_en68->SetMarkerColor(kCyan+2);
  grTrentoP0_en46_en68->SetMarkerStyle(24);
  grTrentoP0_en46_en68->GetXaxis()->SetTitle( "Centrality %");
  grTrentoP0_en46_en68->GetYaxis()->SetTitle( Form("(#font[12]{v}_{%i}{4} - #font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6} - #font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_) );
  grTrentoP0_en46_en68->GetYaxis()->SetTitleOffset(1.45);
  grTrentoP0_en46_en68->GetYaxis()->SetDecimals(2);
  grTrentoP0_en46_en68->GetYaxis()->SetRangeUser(vn46_vn68Min, vn46_vn68Max);

  trentoP1_en46_en68[0]  = -1;
  trentoP1_en46_en68[1]  = -1;
  trentoP1_en46_en68[2]  = 8.752688172043012;
  trentoP1_en46_en68[3]  = 8.623655913978496;
  trentoP1_en46_en68[4]  = 8.344086021505376;
  trentoP1_en46_en68[5]  = 8.38709677419355;
  trentoP1_en46_en68[6]  = 7.56989247311828;
  trentoP1_en46_en68[7]  = 7.440860215053764;
  trentoP1_en46_en68[8]  = 7.096774193548387;
  trentoP1_en46_en68[9]  = 6.795698924731182;
  trentoP1_en46_en68[10] = 6.559139784946238;
  trentoP1_en46_en68[11] = 6.279569892473118;
  grTrentoP1_en46_en68 = new TGraph(NCENT, centBinCenter, trentoP1_en46_en68);
  grTrentoP1_en46_en68->SetLineColor(kOrange+7);
  grTrentoP1_en46_en68->SetMarkerColor(kOrange+7);
  grTrentoP1_en46_en68->SetMarkerStyle(25);
  grTrentoP1_en46_en68->GetXaxis()->SetTitle( "Centrality %");
  grTrentoP1_en46_en68->GetYaxis()->SetTitle( Form("(#font[12]{v}_{%i}{4} - #font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6} - #font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_) );
  grTrentoP1_en46_en68->GetYaxis()->SetTitleOffset(1.45);
  grTrentoP1_en46_en68->GetYaxis()->SetDecimals(2);
  grTrentoP1_en46_en68->GetYaxis()->SetRangeUser(vn46_vn68Min, vn46_vn68Max);

  //-- Theory Curves
  grVn2Theory         = new TGraphErrors("../theoryResults/v2{2}.txt",           "%lg %lg %lg");
  grVn4Theory         = new TGraphErrors("../theoryResults/v2{4}.txt",           "%lg %lg %lg");
  grVn6Theory         = new TGraphErrors("../theoryResults/v2{6}.txt",           "%lg %lg %lg");
  grVn8Theory         = new TGraphErrors("../theoryResults/v2{8}.txt",           "%lg %lg %lg");
  grGamma1ExpTheory   = new TGraphErrors("../theoryResults/skewness.txt",        "%lg %lg %lg");
  grvn6vn4RatioTheory = new TGraphErrors("../theoryResults/v2{6}_over_v2{4}.txt", "%lg %lg %lg");

  formatGraph(grVn2Theory,         "Centrality %", cumuMin,     cumuMax,     Form("#font[12]{v}_{%i}{2}", norder_),                      1,  20, "grVn2Theory");
  formatGraph(grVn4Theory,         "Centrality %", cumuMin,     cumuMax,     Form("#font[12]{v}_{%i}{4}", norder_),                      1,  20, "grVn4Theory");
  formatGraph(grVn6Theory,         "Centrality %", cumuMin,     cumuMax,     Form("#font[12]{v}_{%i}{6}", norder_),                      1,  20, "grVn6Theory");
  formatGraph(grVn8Theory,         "Centrality %", cumuMin,     cumuMax,     Form("#font[12]{v}_{%i}{8}", norder_),                      1,  20, "grVn8Theory");
  formatGraph(grGamma1ExpTheory,   "Centrality %", gamm1expMin, gamm1expMax, "#gamma_{1}^{exp}",                              1, 20, "grGamma1ExpTheory");
  formatGraph(grvn6vn4RatioTheory, "Centrality %", ratioMin,    ratioMax,    Form("#font[12]{v}_{%i}{6} / #font[12]{v}_{%i}{4}", norder_, norder_), kOrange+3, 20, "grvn6vn4RatioTheory");

  grGamma1ExpTheory->SetFillColor(38);
  grGamma1ExpTheory->GetYaxis()->SetDecimals(1);
  //grGamma1ExpTheory->SetFillStyle(3744);

  grvn6vn4RatioTheory->GetYaxis()->SetTitleOffset(titleOffset);
  grvn6vn4RatioTheory->SetFillColor(45);
  grvn6vn4RatioTheory->GetYaxis()->SetNdivisions(507);
  grvn6vn4RatioTheory->GetYaxis()->SetDecimals(2);
  //grvn6vn4RatioTheory->SetFillStyle(3001);


  //-- DRAW!

  // ------------------------ BEGIN PAPER FIGURE 2 ------------------------

  TLegend * legCumu = new TLegend(0.55, 0.27, 0.76, 0.52);
  legCumu->SetBorderSize(0);
  legCumu->SetFillStyle(0);
  legCumu->AddEntry(grVn2RawSys, "#font[12]{m} = 2", "ep");
  legCumu->AddEntry(grVn4RawSys, "#font[12]{m} = 4", "ep");
  legCumu->AddEntry(grVn6RawSys, "#font[12]{m} = 6", "ep");
  legCumu->AddEntry(grVn8RawSys, "#font[12]{m} = 8", "ep");

  double m1[NCENT];
  for(int i = 0; i < NCENT; i++) m1[i] = -1.;
  TGraphErrors * grCumuDummy = new TGraphErrors(NCENT, centBinCenter, m1, nullCentErr, nullCentErr);

  //-- X axis
  grCumuDummy->GetXaxis()->SetTitle( "Centrality %");
  grCumuDummy->GetXaxis()->CenterTitle();
  grCumuDummy->GetXaxis()->SetNdivisions(507);
  grCumuDummy->GetXaxis()->SetLabelFont(43);
  grCumuDummy->GetXaxis()->SetLabelSize(26);
  grCumuDummy->GetXaxis()->SetTitleFont(43);
  grCumuDummy->GetXaxis()->SetTitleSize(35);
  grCumuDummy->GetXaxis()->SetTitleOffset(0.85);

  //-- Y axis
  grCumuDummy->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{#font[12]{m}}", norder_) );
  grCumuDummy->GetYaxis()->CenterTitle();
  grCumuDummy->GetYaxis()->SetRangeUser(cumuMin, cumuMax);
  grCumuDummy->GetYaxis()->SetDecimals(2);
  grCumuDummy->GetYaxis()->SetLabelFont(43);
  grCumuDummy->GetYaxis()->SetLabelSize(26);
  grCumuDummy->GetYaxis()->SetTitleFont(43);
  grCumuDummy->GetYaxis()->SetTitleSize(35);
  grCumuDummy->GetYaxis()->SetTitleOffset(1.3);

  TCanvas * cCumuRaw = new TCanvas("cCumuRaw", "cCumuRaw", 500, 500);
  cCumuRaw->cd();
  cCumuRaw->SetTopMargin(0.1);
  cCumuRaw->SetLeftMargin(0.2);
  cCumuRaw->SetRightMargin(0.1);
  grCumuDummy->Draw("ap");
  grVn2RawSys->Draw("pE2same");
  grVn4RawSys->Draw("psameE2");
  grVn6RawSys->Draw("psameE2");
  grVn8RawSys->Draw("psameE2");
  grVn2Raw->Draw("psame");
  grVn4Raw->Draw("psame");
  grVn6Raw->Draw("psame");  
  grVn8Raw->Draw("psame");
  legCumu->Draw("same");
  legCumu->SetTextFont(43); 
  legCumu->SetTextSize(26);
  latex3.DrawLatex(0.20, 0.913, "#bf{CMS}");
  latex3.DrawLatex(0.385, 0.914, "26 #mub^{-1} (PbPb 5.02 TeV)");
  latex3.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.24, 0.75, Form("|#eta| < %.1f", tkEta));
  latex3.DrawLatex(0.57, 0.54, "#font[12]{v}_{2}{#font[12]{m}}");
  cCumuRaw->SaveAs("../plots/skew/SysCumuRaw.pdf");

  // ------------------------ END PAPER FIGURE 2 ------------------------

  /*
  TCanvas * cCumuVn2 = new TCanvas("cCumuVn2", "cCumuVn2", 500, 500);
  cCumuVn2->cd();
  cCumuVn2->SetTopMargin(0.1);
  cCumuVn2->SetLeftMargin(0.2);
  cCumuVn2->SetRightMargin(0.1);
  grVn2RawSys->GetYaxis()->SetTitle( Form("v_{%i}{2}", norder_) );
  grVn2Raw->GetYaxis()->SetTitle( Form("v_{%i}{2}", norder_) );
  grVn2RawSys->Draw("apE2");
  grVn2Raw->Draw("psame");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.24, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn2->SaveAs("../plots/skew/SysCumuVn2.pdf");

  TCanvas * cCumuVn4 = new TCanvas("cCumuVn4", "cCumuVn4", 500, 500);
  cCumuVn4->cd();
  cCumuVn4->SetTopMargin(0.1);
  cCumuVn4->SetLeftMargin(0.2);
  cCumuVn4->SetRightMargin(0.1);
  grVn4RawSys->GetYaxis()->SetTitle( Form("v_{%i}{4}", norder_) );
  grVn4Raw->GetYaxis()->SetTitle( Form("v_{%i}{4}", norder_) );
  grVn4RawSys->Draw("apE2");
  grVn4Raw->Draw("psame");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.24, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn4->SaveAs("../plots/skew/SysCumuVn4.pdf");

  TCanvas * cCumuVn6 = new TCanvas("cCumuVn6", "cCumuVn6", 500, 500);
  cCumuVn6->cd();
  cCumuVn6->SetTopMargin(0.1);
  cCumuVn6->SetLeftMargin(0.2);
  cCumuVn6->SetRightMargin(0.1);
  grVn6RawSys->GetYaxis()->SetTitle( Form("v_{%i}{6}", norder_) );
  grVn6Raw->GetYaxis()->SetTitle( Form("v_{%i}{6}", norder_) );
  grVn6RawSys->Draw("apE2");
  grVn6Raw->Draw("psame");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.24, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn6->SaveAs("../plots/skew/SysCumuVn6.pdf");

  TCanvas * cCumuVn8 = new TCanvas("cCumuVn8", "cCumuVn8", 500, 500);
  cCumuVn8->cd();
  cCumuVn8->SetTopMargin(0.1);
  cCumuVn8->SetLeftMargin(0.2);
  cCumuVn8->SetRightMargin(0.1);
  grVn8RawSys->GetYaxis()->SetTitle( Form("v_{%i}{8}", norder_) );
  grVn8Raw->GetYaxis()->SetTitle( Form("v_{%i}{8}", norder_) );
  grVn8RawSys->Draw("apE2");
  grVn8Raw->Draw("psame");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.24, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.24, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn8->SaveAs("../plots/skew/SysCumuVn8.pdf");
  */
  /*
  TCanvas * cCumuTheoryComp = new TCanvas("cCumuTheoryComp", "cCumuTheoryComp", 2000, 500);
  cCumuTheoryComp->Divide(4,1);

  cCumuTheoryComp->cd(1);
  grVn2RawSys->Draw("apE2");
  grVn2Raw->Draw("psame");
  grVn2Theory->Draw("E3same");

  cCumuTheoryComp->cd(2);
  grVn4RawSys->Draw("apE2");
  grVn4Raw->Draw("psame");
  grVn4Theory->Draw("E3same");

  cCumuTheoryComp->cd(3);
  grVn6RawSys->Draw("apE2");
  grVn6Raw->Draw("psame");
  grVn6Theory->Draw("E3same");

  cCumuTheoryComp->cd(4);
  grVn8RawSys->Draw("apE2");
  grVn8Raw->Draw("psame");
  grVn8Theory->Draw("E3same");

  cCumuTheoryComp->SaveAs("../plots/skew/cCumuTheoryComp.pdf");
  */

  // ------------------------ BEGIN PAPER FIGURE 4 ------------------------

  TLegend * legg1e = new TLegend(0.22, 0.22, 0.54, 0.375);
  legInit( legg1e );
  legg1e->AddEntry(grGamma1ExpSys,    "#gamma_{1}^{exp}", "ep");
  legg1e->AddEntry(grGamma1ExpTheory, "Hydro",            "f");

  //-- Fig. Gamma1Exp
  double m2[NCENT];
  for(int i = 0; i < NCENT; i++) m2[i] = -100.;
  TGraphErrors * grG1EDummy = new TGraphErrors(NCENT, centBinCenter, m2, nullCentErr, nullCentErr);

  //-- X axis
  grG1EDummy->GetXaxis()->SetTitle( "Centrality %");
  grG1EDummy->GetXaxis()->CenterTitle();
  grG1EDummy->GetXaxis()->SetNdivisions(507);
  grG1EDummy->GetXaxis()->SetLabelFont(43);
  grG1EDummy->GetXaxis()->SetLabelSize(26);
  grG1EDummy->GetXaxis()->SetTitleFont(43);
  grG1EDummy->GetXaxis()->SetTitleSize(35);
  grG1EDummy->GetXaxis()->SetTitleOffset(0.85);

  //-- Y axis
  grG1EDummy->GetYaxis()->SetTitle( "#gamma_{1}^{exp}" );
  grG1EDummy->GetYaxis()->CenterTitle();
  grG1EDummy->GetYaxis()->SetRangeUser(gamm1expMin, gamm1expMax);
  grG1EDummy->GetYaxis()->SetDecimals(1);
  grG1EDummy->GetYaxis()->SetLabelFont(43);
  grG1EDummy->GetYaxis()->SetLabelSize(26);
  grG1EDummy->GetYaxis()->SetTitleFont(43);
  grG1EDummy->GetYaxis()->SetTitleSize(35);
  grG1EDummy->GetYaxis()->SetTitleOffset(1.3);


  TCanvas * cGamma1Exp = new TCanvas("cGamma1Exp", "cGamma1Exp", 500, 500);
  cGamma1Exp->cd();
  cGamma1Exp->SetTopMargin(0.1);
  cGamma1Exp->SetLeftMargin(0.2);
  cGamma1Exp->SetRightMargin(0.1);
  grG1EDummy->Draw("ap");
  grGamma1ExpSys->Draw("pE2same");
  grGamma1ExpTheory->Draw("3same");
  grGamma1ExpTheory->Draw("lXsame");
  grGamma1ExpSys->Draw("pE2same");
  grGamma1Exp->Draw("psame");
  legg1e->Draw("same");
  latex3.DrawLatex(0.20, 0.913, "#bf{CMS}");
  latex3.DrawLatex(0.385, 0.914, "26 #mub^{-1} (PbPb 5.02 TeV)");
  latex3.DrawLatex(0.41, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.68, 0.75, Form("|#eta| < %.1f", tkEta));

  legg1e->SetTextFont(43);
  legg1e->SetTextSize(26);
  cGamma1Exp->Update();

  cGamma1Exp->SaveAs("../plots/skew/SysGamma1Exp.pdf");

  // ------------------------ BEGIN PAPER FIGURE 4 ------------------------


  //-- Fig. Cumu Ratio
  TLegend * leg64 = new TLegend(0.22, 0.21, 0.54, 0.335);
  legInit( leg64 );
  leg64->AddEntry(grvn6vn4Ratio,       "#font[12]{v}_{2}{6} / #font[12]{v}_{2}{4}", "ep");
  leg64->AddEntry(grvn6vn4RatioTheory, "Hydro",                                     "f");

  TLegend * leg84 = new TLegend(0.22, 0.21, 0.54, 0.26);
  legInit( leg84 );
  leg84->AddEntry(grvn8vn4Ratio, "#font[12]{v}_{2}{8} / #font[12]{v}_{2}{4}", "ep");

  TLegend * leg86 = new TLegend(0.25, 0.21, 0.57, 0.26);
  legInit( leg86 );
  leg86->AddEntry(grvn8vn6Ratio, "#font[12]{v}_{2}{8} / #font[12]{v}_{2}{6}", "ep");

  TGraphErrors * grVn64Dummy = new TGraphErrors(NCENT, centBinCenter, m1, nullCentErr, nullCentErr);

  //-- X axis
  grVn64Dummy->GetXaxis()->SetTitle( "Centrality %");
  grVn64Dummy->GetXaxis()->CenterTitle();
  grVn64Dummy->GetXaxis()->SetLimits(0, 61);
  grVn64Dummy->GetXaxis()->SetNdivisions(507);
  grVn64Dummy->GetXaxis()->SetLabelFont(43);
  grVn64Dummy->GetXaxis()->SetLabelSize(26);
  grVn64Dummy->GetXaxis()->SetTitleFont(43);
  grVn64Dummy->GetXaxis()->SetTitleSize(35);
  grVn64Dummy->GetXaxis()->SetTitleOffset(0.85);

  //-- Y axis
  grVn64Dummy->GetYaxis()->SetTitle( "#font[12]{v}_{2}{6} / #font[12]{v}_{2}{4}" );
  grVn64Dummy->GetYaxis()->CenterTitle();
  grVn64Dummy->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  grVn64Dummy->GetYaxis()->SetNdivisions(507);
  grVn64Dummy->GetYaxis()->SetDecimals(2);
  grVn64Dummy->GetYaxis()->SetLabelFont(43);
  grVn64Dummy->GetYaxis()->SetLabelSize(26);
  grVn64Dummy->GetYaxis()->SetTitleFont(43);
  grVn64Dummy->GetYaxis()->SetTitleSize(35);
  grVn64Dummy->GetYaxis()->SetTitleOffset(1.3);

  TGraphErrors * grVn84Dummy = (TGraphErrors*) grVn64Dummy->Clone("grVn84Dummy");
  grVn84Dummy->GetYaxis()->SetTitle( "#font[12]{v}_{2}{8} / #font[12]{v}_{2}{4}" );

  TGraphErrors * grVn86Dummy = (TGraphErrors*) grVn64Dummy->Clone("grVn86Dummy");
  grVn86Dummy->GetYaxis()->SetTitle( "#font[12]{v}_{2}{8} / #font[12]{v}_{2}{6}" );
  grVn86Dummy->GetYaxis()->SetRangeUser(ratioMinVn8Vn6, ratioMaxVn8Vn6);
  grVn86Dummy->GetYaxis()->SetDecimals(3);
  grVn86Dummy->GetYaxis()->SetTitleOffset(1.5);

  TLine * line = new TLine(centbinsDefault[0], 1.0, grVn64Dummy->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLine * line2 = new TLine(centbinsDefault[0], 1.0, grVn84Dummy->GetXaxis()->GetXmax(), 1.0);
  line2->SetLineColor(1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  
  cCumuRatio->cd(1);
  cCumuRatio->cd(1)->SetLeftMargin(0.19);
  cCumuRatio->cd(1)->SetRightMargin(0.07);
  cCumuRatio->cd(1)->SetTopMargin(0.1); 
  grVn64Dummy->Draw("ap");
  grvn6vn4RatioSys->Draw("pE2same");
  grvn6vn4RatioTheory->Draw("3same");
  grvn6vn4RatioTheory->Draw("lXsame");
  grvn6vn4RatioSys->Draw("pE2same");
  grvn6vn4Ratio->Draw("psame");
  line->Draw("same");
  leg64->Draw("same");
  latex3.DrawLatex(0.20, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.58, 0.915, "PbPb 5.02 TeV");
  latex3.DrawLatex(0.413, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.70, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatio->cd(2);
  cCumuRatio->cd(2)->SetLeftMargin(0.19);
  cCumuRatio->cd(2)->SetRightMargin(0.07);
  cCumuRatio->cd(2)->SetTopMargin(0.1);
  grVn84Dummy->Draw("ap");
  grvn8vn4RatioSys->Draw("pE2same");
  grvn8vn4Ratio->Draw("psame");
  line2->Draw("same");
  leg84->Draw("same");
  latex3.DrawLatex(0.20, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.58, 0.915, "PbPb 5.02 TeV");
  latex3.DrawLatex(0.413, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.70, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatio->cd(3);
  cCumuRatio->cd(3)->SetLeftMargin(0.22);
  cCumuRatio->cd(3)->SetRightMargin(0.04);
  cCumuRatio->cd(3)->SetTopMargin(0.1);
  grVn86Dummy->Draw();
  grvn8vn6RatioSys->Draw("pE2same");
  grvn8vn6Ratio->Draw("psame");
  line2->Draw("same");
  leg86->Draw("same");
  latex3.DrawLatex(0.23, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.61, 0.915, "PbPb 5.02 TeV");
  latex3.DrawLatex(0.443, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));

  leg64->SetTextFont(43);
  leg64->SetTextSize(26);

  leg84->SetTextFont(43);
  leg84->SetTextSize(26);

  leg86->SetTextFont(43);
  leg86->SetTextSize(26);

  cCumuRatio->Update();
  cCumuRatio->SaveAs("../plots/skew/SysCumuRatio.pdf");
  /*
  //-- Individuals
  TCanvas * cCumuVn6Vn4 = new TCanvas("cCumuVn6Vn4", "cCumuVn6Vn4", 500, 500);
  cCumuVn6Vn4->cd();
  cCumuVn6Vn4->SetTopMargin(0.1);
  cCumuVn6Vn4->SetLeftMargin(0.2);
  cCumuVn6Vn4->SetRightMargin(0.1);
  grvn6vn4RatioSys->Draw("apE2");
  grvn6vn4Ratio->Draw("psame");
  line->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn6Vn4->SaveAs("../plots/skew/SysCumuVn6Vn4.pdf");

  TCanvas * cCumuVn8Vn4 = new TCanvas("cCumuVn8Vn4", "cCumuVn8Vn4", 500, 500);
  cCumuVn8Vn4->cd();
  cCumuVn8Vn4->SetTopMargin(0.1);
  cCumuVn8Vn4->SetLeftMargin(0.2);
  cCumuVn8Vn4->SetRightMargin(0.1);
  grvn8vn4RatioSys->Draw("apE2");
  grvn8vn4Ratio->Draw("psame");
  line->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn8Vn4->SaveAs("../plots/skew/SysCumuVn8Vn4.pdf");

  TCanvas * cCumuVn8Vn6 = new TCanvas("cCumuVn8Vn6", "cCumuVn8Vn6", 500, 500);
  cCumuVn8Vn6->cd();
  cCumuVn8Vn6->SetTopMargin(0.1);
  cCumuVn8Vn6->SetLeftMargin(0.2);
  cCumuVn8Vn6->SetRightMargin(0.1);
  grvn8vn6RatioSys->Draw("apE2");
  grvn8vn6Ratio->Draw("psame");
  line->Draw("same");
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.48, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));
  cCumuVn8Vn6->SaveAs("../plots/skew/SysCumuVn8Vn6.pdf");
  */
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
    //grvn6vn4RatioSys_Npart->SetFillColorAlpha(46, 0.6);
    grvn6vn4RatioSys_Npart->SetFillColor(46);

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
    //grvn8vn4RatioSys_Npart->SetFillColorAlpha(46, 0.6);
    grvn8vn4RatioSys_Npart->SetFillColor(46);
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
    grvn6vn4Ratio_ATLASNpart->SetMarkerSize(1.2);
    grvn6vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "Centrality %");
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{6} / #font[12]{v}_{%i}{4}", norder_, norder_) );
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
    grvn8vn4Ratio_ATLASNpart->SetMarkerSize(1.2);
    grvn8vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "Centrality %");
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("#font[12]{v}_{%i}{8} / #font[12]{v}_{%i}{4}", norder_, norder_) );
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetTitleOffset(titleOffset);
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetDecimals(2);

    // ------------------------ BEGIN PAPER FIGURE 3 ------------------------

    TLegend * leg64a = new TLegend(0.27, 0.22, 0.69, 0.34);
    leg64a->SetFillStyle(0);
    leg64a->SetBorderSize(0);
    leg64a->AddEntry(grvn6vn4Ratio_Npart,      "CMS",   "ep");
    leg64a->AddEntry(grvn6vn4Ratio_ATLASNpart, "ATLAS", "ep");

    TLegend * leg84a = new TLegend(0.27, 0.22, 0.69, 0.34);
    leg84a->SetFillStyle(0);
    leg84a->SetBorderSize(0);
    leg84a->AddEntry(grvn8vn4Ratio_Npart,      "CMS",   "ep");
    leg84a->AddEntry(grvn8vn4Ratio_ATLASNpart, "ATLAS", "ep");

    TLine * lone = new TLine(grvn6vn4Ratio_ATLASNpart->GetXaxis()->GetXmin(), 1.0, grvn6vn4Ratio_ATLASNpart->GetXaxis()->GetXmax(), 1.0);
    lone->SetLineStyle(2);
    lone->SetLineWidth(2);

    //-- Fig. ATLAS comp
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
    latex.DrawLatex(0.19, 0.915, "#bf{CMS}");
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
    latex.DrawLatex(0.19, 0.915, "#bf{CMS}");
    latex.DrawLatex(0.65, 0.915, "PbPb 5.02 TeV");
    latex2.DrawLatex(0.47, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
    latex2.DrawLatex(0.715, 0.76, Form("|#eta| < %.1f", tkEta));


    leg64a->SetTextFont(43);
    leg64a->SetTextSize(23);
    leg84a->SetTextFont(43);
    leg84a->SetTextSize(23);

    cATLASComp->Update();
    cATLASComp->SaveAs("../plots/skew/cSysATLASComp.pdf");

    // ------------------------ END PAPER FIGURE 3 ------------------------

  }

  TLine * l11 = new TLine(grvn46_vn68RatioSys->GetXaxis()->GetXmin(), 11., grvn46_vn68RatioSys->GetXaxis()->GetXmax(), 11.);
  l11->SetLineColor(1);
  l11->SetLineWidth(2);
  l11->SetLineStyle(2);

  TLegend * leg4668 = new TLegend(0.37, 0.21, 0.89, 0.38);
  legInit(leg4668);
  leg4668->AddEntry(grvn46_vn68Ratio,      Form("(#font[12]{v}_{%i}{4}-#font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6}-#font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_), "ep");
  leg4668->AddEntry(grTrentoP1_en46_en68,  "Trento p = 1",  "ep");
  leg4668->AddEntry(grTrentoP0_en46_en68,  "Trento p = 0",  "ep");
  leg4668->AddEntry(grTrentoPm1_en46_en68, "Trento p = -1", "ep");

  TCanvas * cVn46_Vn68 = new TCanvas("cVn46_Vn68", "cVn46_Vn68", 500, 500);
  cVn46_Vn68->cd();
  cVn46_Vn68->SetTopMargin(0.1);
  cVn46_Vn68->SetLeftMargin(0.2);
  cVn46_Vn68->SetRightMargin(0.1);
  grvn46_vn68RatioSys->Draw("apE2");
  grvn46_vn68Ratio->Draw("psame");
  grTrentoP1_en46_en68->Draw("psame");
  grTrentoP0_en46_en68->Draw("psame");
  grTrentoPm1_en46_en68->Draw("psame");
  l11->Draw("same");
  leg4668->Draw("same");
  latex.DrawLatex(0.20, 0.913, "#bf{CMS}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.49, 0.84, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.725, 0.77, Form("|#eta| < %.1f", tkEta));
  cVn46_Vn68->SaveAs("../plots/skew/SysVn46_Vn68.pdf");

  TLegend * leg46682 = new TLegend(0.37, 0.15, 0.56, 0.32);
  legInit(leg46682);
  leg46682->AddEntry(grvn46_vn68Ratio,      Form("(#font[12]{v}_{%i}{4}-#font[12]{v}_{%i}{6})/(#font[12]{v}_{%i}{6}-#font[12]{v}_{%i}{8})", norder_, norder_, norder_, norder_), "ep");

  //-- Fig. Gamma1Exp with Vn46_Vn68
  TCanvas * cG1eAndVn46_Vn68 = new TCanvas("cG1eAndVn46_Vn68", "cG1eAndVn46_Vn68", 1000, 500);
  cG1eAndVn46_Vn68->Divide(2,1);

  cG1eAndVn46_Vn68->cd(1);
  cG1eAndVn46_Vn68->cd(1)->SetLeftMargin(0.19);
  cG1eAndVn46_Vn68->cd(1)->SetRightMargin(0.09);
  cG1eAndVn46_Vn68->cd(1)->SetTopMargin(0.1);
  grGamma1ExpSys->Draw("apE2");
  grGamma1ExpTheory->Draw("3same");
  grGamma1ExpTheory->Draw("lXsame");
  grGamma1ExpSys->Draw("pE2same");
  grGamma1Exp->Draw("psame");
  legg1e->Draw("same");
  latex.DrawLatex(0.19, 0.915, "#bf{CMS}");
  latex.DrawLatex(0.65, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.47, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.715, 0.76, Form("|#eta| < %.1f", tkEta));

  cG1eAndVn46_Vn68->cd(2);
  cG1eAndVn46_Vn68->cd(2)->SetLeftMargin(0.19);
  cG1eAndVn46_Vn68->cd(2)->SetRightMargin(0.09);
  cG1eAndVn46_Vn68->cd(2)->SetTopMargin(0.1);
  grvn46_vn68RatioSys->Draw("apE2");
  grvn46_vn68Ratio->Draw("psame");
  l11->Draw("same");
  leg46682->Draw("same");
  latex.DrawLatex(0.19, 0.915, "#bf{CMS}");
  latex.DrawLatex(0.65, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.47, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.715, 0.76, Form("|#eta| < %.1f", tkEta));


  leg46682->SetTextFont(43);
  leg46682->SetTextSize(23);

  cG1eAndVn46_Vn68->Update();
  cG1eAndVn46_Vn68->SaveAs("../plots/skew/cSysG1EAndVn46_Vn68.pdf");

  //-- Fig. Cumu Ratio with ATLAS
  grvn6vn4Ratio_ATLASNpart->SetLineColor(1);
  grvn6vn4Ratio_ATLASNpart->SetMarkerColor(1);
  grvn6vn4Ratio_ATLASNpart->SetMarkerStyle(24);

  grvn8vn4Ratio_ATLASNpart->SetLineColor(1);
  grvn8vn4Ratio_ATLASNpart->SetMarkerColor(1);
  grvn8vn4Ratio_ATLASNpart->SetMarkerStyle(25);

  grvn6vn4RatioSys->SetFillColor(17);
  grvn8vn4RatioSys->SetFillColor(17);

  TLegend * leg64c = new TLegend(0.29, 0.20, 0.62, 0.43);
  legInit( leg64c );
  leg64c->AddEntry(grvn6vn4Ratio,            "#font[12]{v}_{2}{6} / #font[12]{v}_{2}{4}", "ep");
  leg64c->AddEntry(grvn6vn4Ratio_ATLASNpart, "ATLAS",      "ep");
  leg64c->AddEntry(grvn6vn4RatioTheory,      "Hydro",      "f");

  TLegend * leg84c = new TLegend(0.29, 0.20, 0.61, 0.33);
  legInit( leg84c );
  leg84c->AddEntry(grvn8vn4Ratio,            "#font[12]{v}_{2}{8} / #font[12]{v}_{2}{4}", "ep");
  leg84c->AddEntry(grvn8vn4Ratio_ATLASNpart, "ATLAS",      "ep");

  TCanvas * cCumuRatioWithATLAS = new TCanvas("cCumuRatioWithATLAS", "cCumuRatioWithATLAS", 1500, 500);
  cCumuRatioWithATLAS->Divide(3,1);

  cCumuRatioWithATLAS->cd(1);
  cCumuRatioWithATLAS->cd(1)->SetLeftMargin(0.19);
  cCumuRatioWithATLAS->cd(1)->SetRightMargin(0.07);
  cCumuRatioWithATLAS->cd(1)->SetTopMargin(0.1);
  grVn64Dummy->Draw("ap");
  grvn6vn4RatioSys->Draw("pE2same");
  grvn6vn4RatioTheory->Draw("3same");
  grvn6vn4RatioTheory->Draw("lXsame");
  grvn6vn4Ratio_ATLASNpart->Draw("psame");
  grvn6vn4RatioSys->Draw("pE2same");
  grvn6vn4Ratio_ATLASNpart->Draw("psame");
  grvn6vn4Ratio->Draw("psame");
  line->Draw("same");
  leg64c->Draw("same");
  latex3.DrawLatex(0.20, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.38, 0.916, "26 #mub^{-1} (PbPb 5.02 TeV)");
  latex3.DrawLatex(0.413, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.70, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatioWithATLAS->cd(2);
  cCumuRatioWithATLAS->cd(2)->SetLeftMargin(0.19);
  cCumuRatioWithATLAS->cd(2)->SetRightMargin(0.07);
  cCumuRatioWithATLAS->cd(2)->SetTopMargin(0.1);
  grVn84Dummy->Draw("ap");
  grvn8vn4RatioSys->Draw("pE2same");
  grvn8vn4Ratio_ATLASNpart->Draw("psame");
  grvn8vn4RatioSys->Draw("pE2same");
  grvn8vn4Ratio_ATLASNpart->Draw("psame");
  grvn8vn4Ratio->Draw("psame");
  line2->Draw("same");
  leg84c->Draw("same");
  latex3.DrawLatex(0.20, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.38, 0.916, "26 #mub^{-1} (PbPb 5.02 TeV)");
  latex3.DrawLatex(0.413, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.70, 0.76, Form("|#eta| < %.1f", tkEta));

  cCumuRatioWithATLAS->cd(3);
  cCumuRatioWithATLAS->cd(3)->SetLeftMargin(0.22);  //0.19
  cCumuRatioWithATLAS->cd(3)->SetRightMargin(0.04); //0.07
  cCumuRatioWithATLAS->cd(3)->SetTopMargin(0.1);
  grVn86Dummy->Draw("ap");
  grvn8vn6RatioSys->Draw("pE2same");
  grvn8vn6Ratio->Draw("psame");
  line2->Draw("same");
  leg86->Draw("same");
  latex3.DrawLatex(0.23, 0.915, "#bf{CMS}");
  latex3.DrawLatex(0.41, 0.916, "26 #mub^{-1} (PbPb 5.02 TeV)");
  latex3.DrawLatex(0.443, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.73, 0.76, Form("|#eta| < %.1f", tkEta));

  leg64c->SetTextFont(43);
  leg64c->SetTextSize(26);

  leg84c->SetTextFont(43);
  leg84c->SetTextSize(26);

  leg86->SetTextFont(43);
  leg86->SetTextSize(26);

  cCumuRatioWithATLAS->Update();
  cCumuRatioWithATLAS->SaveAs("../plots/skew/SysCumuRatioWithATLAS.pdf");

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

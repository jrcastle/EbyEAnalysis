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


void sysClusCompatTune(int n = 2, double e = 1.0){

  const int norder_  = n;
  const double tkEta = e;

  double cumuMin = 0.0;
  double cumuMax = 0.12;

  double vn6vn4Min    = 0.95;
  double vn6vn4Max    = 1.03;
  double vn8vn4Min    = 0.95;
  double vn8vn4Max    = 1.03;
  double vn8vn6Min    = 0.99;
  double vn8vn6Max    = 1.002;
  double vn46_vn68Min = 0.;
  double vn46_vn68Max = 30.;
  double g1eMin       = -1.0;
  double g1eMax       = 0.5;

  double rMin = 0.95;
  double rMax = 1.05;
  double g1rMin = 0.;
  double g1rMax = 0.4;

  TLatex latex;

  //-- Default
  TFile * fAnaDefault;
  TFile * fUnfDefault;
  TH1D * hObsDefault[NCENT];
  TH1D * hMultDefault[NCENT];
  TH1D * hUnfoldDefault[NCENT][NITER];
  TH1D * hRefoldDefault[NCENT][NITER];
  bool stopFoundDefault;
  int iterCutDefault[NCENT];

  TFile * fStatDefault;
  TH1D * hVarianceOfMean_Vn2_Default;
  TH1D * hVarianceOfMean_Vn4_Default;
  TH1D * hVarianceOfMean_Vn6_Default;
  TH1D * hVarianceOfMean_Vn8_Default;
  TH1D * hVarianceOfMean_Vn6Vn4_Default;
  TH1D * hVarianceOfMean_Vn8Vn4_Default;
  TH1D * hVarianceOfMean_Vn8Vn6_Default;
  TH1D * hVarianceOfMean_Gamma1Exp_Default;
  TH1D * hVarianceOfMean_Vn46_Vn68_Default;

  double Vn2Default[NCENT];
  double Vn4Default[NCENT];
  double Vn6Default[NCENT];
  double Vn8Default[NCENT];
  double vn6vn4Default[NCENT];
  double vn8vn4Default[NCENT];
  double vn8vn6Default[NCENT];
  double vn46_vn68Default[NCENT];
  double g1eDefault[NCENT];

  double vn2Default_err[NCENT];
  double vn4Default_err[NCENT];
  double vn6Default_err[NCENT];
  double vn8Default_err[NCENT];
  double vn6vn4Default_err[NCENT];
  double vn8vn4Default_err[NCENT];
  double vn8vn6Default_err[NCENT];
  double vn46_vn68Default_err[NCENT];
  double g1eDefault_err[NCENT];

  TGraphErrors * grVn2Default;
  TGraphErrors * grVn4Default;
  TGraphErrors * grVn6Default;
  TGraphErrors * grVn8Default;
  TGraphErrors * grVn6Vn4Default;
  TGraphErrors * grVn8Vn4Default;
  TGraphErrors * grVn8Vn6Default;
  TGraphErrors * grVn46_Vn68Default;
  TGraphErrors * grG1EDefault;

  //-- 2%
  TFile * fAna2pct;
  TFile * fUnf2pct;
  TH1D * hObs2pct[NCENT];
  TH1D * hMult2pct[NCENT];
  TH1D * hUnfold2pct[NCENT][NITER];
  TH1D * hRefold2pct[NCENT][NITER];
  bool stopFound2pct;
  int iterCut2pct[NCENT];

  TFile * fStat2pct;
  TH1D * hVarianceOfMean_Vn2_2pct;
  TH1D * hVarianceOfMean_Vn4_2pct;
  TH1D * hVarianceOfMean_Vn6_2pct;
  TH1D * hVarianceOfMean_Vn8_2pct;
  TH1D * hVarianceOfMean_Vn6Vn4_2pct;
  TH1D * hVarianceOfMean_Vn8Vn4_2pct;
  TH1D * hVarianceOfMean_Vn8Vn6_2pct;
  TH1D * hVarianceOfMean_Vn46_Vn68_2pct;
  TH1D * hVarianceOfMean_Gamma1Exp_2pct;

  double Vn22pct[NCENT];
  double Vn42pct[NCENT];
  double Vn62pct[NCENT];
  double Vn82pct[NCENT];
  double vn6vn42pct[NCENT];
  double vn8vn42pct[NCENT];
  double vn8vn62pct[NCENT];
  double vn46_vn682pct[NCENT];
  double g1e2pct[NCENT];

  double vn22pct_err[NCENT];
  double vn42pct_err[NCENT];
  double vn62pct_err[NCENT];
  double vn82pct_err[NCENT];
  double vn6vn42pct_err[NCENT];
  double vn8vn42pct_err[NCENT];
  double vn8vn62pct_err[NCENT];
  double vn46_vn682pct_err[NCENT];
  double g1e2pct_err[NCENT];

  double vn22pct_RatioToDefault[NCENT];
  double vn42pct_RatioToDefault[NCENT];
  double vn62pct_RatioToDefault[NCENT];
  double vn82pct_RatioToDefault[NCENT];
  double vn6vn42pct_RatioToDefault[NCENT];
  double vn8vn42pct_RatioToDefault[NCENT];
  double vn8vn62pct_RatioToDefault[NCENT];
  double vn46_vn682pct_RatioToDefault[NCENT];
  double g1e2pct_RatioToDefault[NCENT];

  double vn22pct_RatioToDefault_err[NCENT];
  double vn42pct_RatioToDefault_err[NCENT];
  double vn62pct_RatioToDefault_err[NCENT];
  double vn82pct_RatioToDefault_err[NCENT];
  double vn6vn42pct_RatioToDefault_err[NCENT];
  double vn8vn42pct_RatioToDefault_err[NCENT];
  double vn8vn62pct_RatioToDefault_err[NCENT];
  double vn46_vn682pct_RatioToDefault_err[NCENT];
  double g1e2pct_RatioToDefault_err[NCENT];

  TGraphErrors * grVn22pct;
  TGraphErrors * grVn42pct;
  TGraphErrors * grVn62pct;
  TGraphErrors * grVn82pct;
  TGraphErrors * grVn6Vn42pct;
  TGraphErrors * grVn8Vn42pct;
  TGraphErrors * grVn8Vn62pct;
  TGraphErrors * grVn46_Vn682pct;
  TGraphErrors * grG1E2pct;

  TGraphErrors * grVn22pct_RatioToDefault;
  TGraphErrors * grVn42pct_RatioToDefault;
  TGraphErrors * grVn62pct_RatioToDefault;
  TGraphErrors * grVn82pct_RatioToDefault;
  TGraphErrors * grVn6Vn42pct_RatioToDefault;
  TGraphErrors * grVn8Vn42pct_RatioToDefault;
  TGraphErrors * grVn8Vn62pct_RatioToDefault;
  TGraphErrors * grVn46_Vn682pct_RatioToDefault;
  TGraphErrors * grG1E2pct_RatioToDefault;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get Analyzer files
  fAnaDefault = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fAna2pct    = new TFile( "newCCTune2pct/AnalyzerResults/CastleEbyE.root" );

  //-- Get Unfold files
  fUnfDefault = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnf2pct    = new TFile( Form("newCCTune2pct/UnfoldResults/dataResp/data%i.root", norder_)  );

  //-- Get Statistical Uncertainties
  fStatDefault = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertainties_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_Default       = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn2" );
  hVarianceOfMean_Vn4_Default       = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn4" );
  hVarianceOfMean_Vn6_Default       = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn6" );
  hVarianceOfMean_Vn8_Default       = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn8" );
  hVarianceOfMean_Vn6Vn4_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Vn46_Vn68_Default = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn46_Vn68" );
  hVarianceOfMean_Gamma1Exp_Default = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Gamma1Exp" );

  if( !hVarianceOfMean_Vn2_Default ) {
    std::cout << "WARNING! Statistical resampling procedure not run!\n"
              << "Please run the procedure first and then run this macro\n"
              << "Exiting now..."
              << std::endl;
    exit(0);
  }

  fStat2pct  = new TFile(Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/clusCompatTune/newCCTune2pct/StatisticalUncertainties_v%i.root", norder_, tkEta, norder_));
  hVarianceOfMean_Vn2_2pct       = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn2" );
  hVarianceOfMean_Vn4_2pct       = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn4" );
  hVarianceOfMean_Vn6_2pct       = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn6" );
  hVarianceOfMean_Vn8_2pct       = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn8" );
  hVarianceOfMean_Vn6Vn4_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Vn46_Vn68_2pct = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn46_Vn68" );
  hVarianceOfMean_Gamma1Exp_2pct = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Gamma1Exp" );

  if( !hVarianceOfMean_Vn2_2pct ) {
    std::cout << "WARNING! Statistical resampling procedure not run!\n"
              << "Please run the procedure first and then run this macro\n"
              << "Exiting now..."
              << std::endl;
    exit(0);
  }

  //-- Start Grabbing Histos
  for(int icent = 0; icent < NCENT; icent++){

    //-- Vn Observed 
    hObsDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/hVnFull_c%i", icent) );
    if( !hObsDefault[icent]  ) break;
    hObsDefault[icent]->SetName( Form("hVnFullDefault_c%i", icent) );
    
    hObs2pct[icent] = (TH1D*) fAna2pct->Get( Form("qwebye/hVnFull_c%i", icent) );
    if( !hObs2pct[icent]  ) break;
    hObs2pct[icent]->SetName( Form("hVnFull2pct_c%i", icent) );

    //-- Multiplicity
    hMultDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/Mult_c%i", icent) );
    hMultDefault[icent]->SetLineColor(centCol[icent]);
    hMultDefault[icent]->SetMarkerColor(centCol[icent]);
    hMultDefault[icent]->SetMaximum( 20.*hMultDefault[icent]->GetMaximum() );
    hMultDefault[icent]->GetXaxis()->SetNdivisions(509);
    hMultDefault[icent]->GetXaxis()->SetRange(1, hMultDefault[icent]->FindBin(5500));

    hMult2pct[icent] = (TH1D*) fAna2pct->Get( Form("qwebye/Mult_c%i", icent) );
    hMult2pct[icent]->SetLineColor(centCol[icent]);
    hMult2pct[icent]->SetMarkerColor(centCol[icent]);
    hMult2pct[icent]->SetMaximum( 20.*hMult2pct[icent]->GetMaximum() );
    hMult2pct[icent]->GetXaxis()->SetNdivisions(509);
    hMult2pct[icent]->GetXaxis()->SetRange(1, hMult2pct[icent]->FindBin(5500));

    //-- Reset stopFound boolean
    stopFoundDefault = 0;
    stopFound2pct    = 0;

    //-- Grab un/refolded distns
    for(int i = 0; i < NITER; i++){

      //-- Default
      hUnfoldDefault[icent][i] = (TH1D*) fUnfDefault->Get( Form("hreco%i_c%i", iter[i], icent) );
      if( !hUnfoldDefault[icent][i] ) break;
      hUnfoldDefault[icent][i]->SetName( Form("hrecoDefault%i_c%i", iter[i], icent) );
      hRefoldDefault[icent][i] = (TH1D*) fUnfDefault->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldDefault[icent][i]->SetName( Form("hrefolDefault%i_c%i", iter[i], icent) );

      double chi2Default = hRefoldDefault[icent][i]->Chi2Test(hObsDefault[icent], "UWCHI2/NDF");
      if( chi2Default <= 1.2 && !stopFoundDefault ){
	stopFoundDefault = 1;
	iterCutDefault[icent] = i;
      }
      if( i == NITER-1 && !stopFoundDefault ){
	stopFoundDefault = 1;
        iterCutDefault[icent] = i;
      }

      //-- 2.0%
      hUnfold2pct[icent][i] = (TH1D*) fUnf2pct->Get( Form("hreco%i_c%i", iter[i], icent) );
      if( !hUnfold2pct[icent][i] ) break;
      hUnfold2pct[icent][i]->SetName( Form("hreco2pct%i_c%i", iter[i], icent) );
      hRefold2pct[icent][i] = (TH1D*) fUnf2pct->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold2pct[icent][i]->SetName( Form("hrefol2pct%i_c%i", iter[i], icent) );

      double chi22pct = hRefold2pct[icent][i]->Chi2Test(hObs2pct[icent], "UWCHI2/NDF");
      if( chi22pct <= 1.2 && !stopFound2pct ){
        stopFound2pct = 1;
        iterCut2pct[icent] = i;
      }
      if( i == NITER-1 && !stopFound2pct ){
        stopFound2pct = 1;
        iterCut2pct[icent] = i;
      }

    } //-- End iter loop

  } //-- End cent loop

  if( !hObsDefault[0] || !hObs2pct[0] || !hUnfoldDefault[0][0] || !hUnfold2pct[0][0] ){
    std::cout << "WARNING! Unfolding procedure not run!\n"
              << "Please run the unfolding procedure first and then run this macro\n"
              << "Exiting now..."
              << std::endl;
    exit(0);
  }

  //-- Now that we're done fetching things, let's build some arrays!
  for(int icent = 0; icent < NCENT; icent++){

    //-- Default
    int iDefault = iterCutDefault[icent];
    FixUnfold(hUnfoldDefault[icent][iDefault]);
    EbyECumu cumuDefault(hUnfoldDefault[icent][iDefault]);
    double vn2Default  = cumuDefault.GetCumu_vn2();
    double vn4Default  = cumuDefault.GetCumu_vn4();
    double vn6Default  = cumuDefault.GetCumu_vn6();
    double vn8Default  = cumuDefault.GetCumu_vn8();
    double g1exDefault = cumuDefault.GetGamma1Exp();

    Vn2Default[icent] = vn2Default;
    Vn4Default[icent] = vn4Default;
    Vn6Default[icent] = vn6Default;
    Vn8Default[icent] = vn8Default;

    if(vn4Default == 0 || vn6Default == 0) vn6vn4Default[icent]= -1000.;
    else                                   vn6vn4Default[icent] = vn6Default / vn4Default;
    if(vn4Default == 0 || vn8Default == 0) vn8vn4Default[icent]= -1000.;
    else                                   vn8vn4Default[icent] = vn8Default / vn4Default;
    if(vn6Default == 0 || vn8Default == 0) vn8vn6Default[icent]= -1000.;
    else                                   vn8vn6Default[icent] = vn8Default / vn6Default;
    if(vn4Default == 0 || vn6Default == 0 || vn8Default == 0) vn46_vn68Default[icent]= -1000.;
    else                                                      vn46_vn68Default[icent] = (vn4Default - vn6Default)/(vn6Default - vn8Default);
    g1eDefault[icent] = g1exDefault;

    double vn2Default_staterr       = sqrt( hVarianceOfMean_Vn2_Default->GetBinContent(icent+1) );
    double vn4Default_staterr       = sqrt( hVarianceOfMean_Vn4_Default->GetBinContent(icent+1) );
    double vn6Default_staterr       = sqrt( hVarianceOfMean_Vn6_Default->GetBinContent(icent+1) );
    double vn8Default_staterr       = sqrt( hVarianceOfMean_Vn8_Default->GetBinContent(icent+1) );
    double vn6vn4Default_staterr    = sqrt( hVarianceOfMean_Vn6Vn4_Default->GetBinContent(icent+1) );
    double vn8vn4Default_staterr    = sqrt( hVarianceOfMean_Vn8Vn4_Default->GetBinContent(icent+1) );
    double vn8vn6Default_staterr    = sqrt( hVarianceOfMean_Vn8Vn6_Default->GetBinContent(icent+1) );
    double vn46_vn68Default_staterr = sqrt( hVarianceOfMean_Vn46_Vn68_Default->GetBinContent(icent+1) );
    double g1eDefault_staterr       = sqrt( hVarianceOfMean_Gamma1Exp_Default->GetBinContent(icent+1) );

    vn2Default_err[icent]       = vn2Default_staterr;
    vn4Default_err[icent]       = vn4Default_staterr;
    vn6Default_err[icent]       = vn6Default_staterr;
    vn8Default_err[icent]       = vn8Default_staterr;
    vn6vn4Default_err[icent]    = vn6vn4Default_staterr;
    vn8vn4Default_err[icent]    = vn8vn4Default_staterr;
    vn8vn6Default_err[icent]    = vn8vn6Default_staterr;
    vn46_vn68Default_err[icent] = vn46_vn68Default_staterr;
    g1eDefault_err[icent]       = g1eDefault_staterr;

    //-- 2%
    //-- Raw ////////////////////////////////////////////START HERE
    int i2pct = iterCut2pct[icent];
    FixUnfold(hUnfold2pct[icent][i2pct]);
    EbyECumu cumu2pct(hUnfold2pct[icent][i2pct]);
    double vn22pct  = cumu2pct.GetCumu_vn2();
    double vn42pct  = cumu2pct.GetCumu_vn4();
    double vn62pct  = cumu2pct.GetCumu_vn6();
    double vn82pct  = cumu2pct.GetCumu_vn8();
    double g1ex2pct = cumu2pct.GetGamma1Exp();

    Vn22pct[icent] = vn22pct;
    Vn42pct[icent] = vn42pct;
    Vn62pct[icent] = vn62pct;
    Vn82pct[icent] = vn82pct;

    if(vn42pct == 0 || vn62pct == 0) vn6vn42pct[icent] = -1000.;
    else                             vn6vn42pct[icent] = vn62pct / vn42pct;
    if(vn42pct == 0 || vn82pct == 0) vn8vn42pct[icent] = -1000.;
    else                             vn8vn42pct[icent] = vn82pct / vn42pct;
    if(vn62pct == 0 || vn82pct == 0) vn8vn62pct[icent] = -1000.;
    else                             vn8vn62pct[icent] = vn82pct / vn62pct;
    if(vn42pct == 0 || vn62pct == 0 || vn82pct == 0) vn46_vn682pct[icent]= -1000.;
    else                                             vn46_vn682pct[icent] = (vn42pct - vn62pct)/(vn62pct - vn82pct);
    g1e2pct[icent] = g1ex2pct;

    double vn22pct_staterr       = sqrt( hVarianceOfMean_Vn2_2pct->GetBinContent(icent+1) );
    double vn42pct_staterr       = sqrt( hVarianceOfMean_Vn4_2pct->GetBinContent(icent+1) );
    double vn62pct_staterr       = sqrt( hVarianceOfMean_Vn6_2pct->GetBinContent(icent+1) );
    double vn82pct_staterr       = sqrt( hVarianceOfMean_Vn8_2pct->GetBinContent(icent+1) );
    double vn6vn42pct_staterr    = sqrt( hVarianceOfMean_Vn6Vn4_2pct->GetBinContent(icent+1) );
    double vn8vn42pct_staterr    = sqrt( hVarianceOfMean_Vn8Vn4_2pct->GetBinContent(icent+1) );
    double vn8vn62pct_staterr    = sqrt( hVarianceOfMean_Vn8Vn6_2pct->GetBinContent(icent+1) );
    double vn46_vn682pct_staterr = sqrt( hVarianceOfMean_Vn46_Vn68_2pct->GetBinContent(icent+1) );
    double g1e2pct_staterr       = sqrt( hVarianceOfMean_Gamma1Exp_2pct->GetBinContent(icent+1) );
 
    vn22pct_err[icent]       = vn22pct_staterr;
    vn42pct_err[icent]       = vn42pct_staterr;
    vn62pct_err[icent]       = vn62pct_staterr;
    vn82pct_err[icent]       = vn82pct_staterr;
    vn6vn42pct_err[icent]    = vn6vn42pct_staterr;
    vn8vn42pct_err[icent]    = vn8vn42pct_staterr;
    vn8vn62pct_err[icent]    = vn8vn62pct_staterr;
    vn46_vn682pct_err[icent] = vn46_vn682pct_staterr;
    g1e2pct_err[icent]       = g1e2pct_staterr;

    //-- Ratio to default
    if(vn22pct == 0 || vn2Default == 0) vn22pct_RatioToDefault[icent] = 0.;
    else                                vn22pct_RatioToDefault[icent] = vn22pct / vn2Default;
    if(vn42pct == 0 || vn4Default == 0) vn42pct_RatioToDefault[icent] = 0.;
    else                                vn42pct_RatioToDefault[icent] = vn42pct/ vn4Default;
    if(vn62pct == 0 || vn6Default == 0) vn62pct_RatioToDefault[icent] = 0.;
    else                                vn62pct_RatioToDefault[icent] = vn62pct/ vn6Default;
    if(vn82pct == 0 || vn8Default == 0) vn82pct_RatioToDefault[icent] = 0.;
    else                                vn82pct_RatioToDefault[icent] = vn82pct/ vn8Default;

    if(vn42pct == 0 || vn62pct == 0 || vn4Default == 0 || vn6Default == 0) vn6vn42pct_RatioToDefault[icent] = 0.;
    else                                                                   vn6vn42pct_RatioToDefault[icent] = vn6vn42pct[icent] / vn6vn4Default[icent];
    if(vn42pct == 0 || vn82pct == 0 || vn4Default == 0 || vn8Default == 0) vn8vn42pct_RatioToDefault[icent] = 0.;
    else                                                                   vn8vn42pct_RatioToDefault[icent] = vn8vn42pct[icent] / vn8vn4Default[icent];
    if(vn62pct == 0 || vn82pct == 0 || vn6Default == 0 || vn8Default == 0) vn8vn62pct_RatioToDefault[icent] = 0.;
    else                                                                   vn8vn62pct_RatioToDefault[icent] = vn8vn62pct[icent] / vn8vn6Default[icent];
    if(vn46_vn682pct[icent] <= 0 || vn46_vn68Default[icent] <= 0) vn46_vn682pct_RatioToDefault[icent] = 0.;
    else                                                          vn46_vn682pct_RatioToDefault[icent] = vn46_vn682pct[icent] / vn46_vn68Default[icent];
    if(g1ex2pct== -10000. || g1exDefault == -10000.) g1e2pct_RatioToDefault[icent] = 0.;
    else                                                          g1e2pct_RatioToDefault[icent]= fabs(1.0 - g1e2pct[icent]/ g1eDefault[icent]);

    vn22pct_RatioToDefault_err[icent]       = sqrt( pow(vn22pct_err[icent]/Vn2Default[icent],2) + pow(Vn22pct[icent]*vn2Default_err[icent]/Vn2Default[icent]/Vn2Default[icent],2) );
    vn42pct_RatioToDefault_err[icent]       = sqrt( pow(vn42pct_err[icent]/Vn4Default[icent],2) + pow(Vn42pct[icent]*vn4Default_err[icent]/Vn4Default[icent]/Vn4Default[icent],2) );
    vn62pct_RatioToDefault_err[icent]       = sqrt( pow(vn62pct_err[icent]/Vn6Default[icent],2) + pow(Vn62pct[icent]*vn6Default_err[icent]/Vn6Default[icent]/Vn6Default[icent],2) );
    vn82pct_RatioToDefault_err[icent]       = sqrt( pow(vn82pct_err[icent]/Vn8Default[icent],2) + pow(Vn82pct[icent]*vn8Default_err[icent]/Vn8Default[icent]/Vn8Default[icent],2) );
    vn6vn42pct_RatioToDefault_err[icent]    = sqrt( pow(vn6vn42pct_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn42pct[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn42pct_RatioToDefault_err[icent]    = sqrt( pow(vn8vn42pct_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn42pct[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn62pct_RatioToDefault_err[icent]    = sqrt( pow(vn8vn62pct_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn62pct[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    vn46_vn682pct_RatioToDefault_err[icent] = sqrt( pow(vn46_vn682pct_err[icent]/vn46_vn68Default[icent],2) + pow(vn46_vn682pct[icent]*vn46_vn68Default_err[icent]/vn46_vn68Default[icent]/vn46_vn68Default[icent],2) );
    g1e2pct_RatioToDefault_err[icent]       = sqrt( pow(g1e2pct_err[icent]/g1eDefault[icent],2) + pow(g1e2pct[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

  } //-- End cent loop

  //-- TGraph time!
  double c_err[NCENT];
  for(int i = 0; i < NCENT; i++) c_err[i] = 0.;

  //-- Default
  grVn2Default       = new TGraphErrors(NCENT, centBinCenter, Vn2Default,       c_err, vn2Default_err);
  grVn4Default       = new TGraphErrors(NCENT, centBinCenter, Vn4Default,       c_err, vn4Default_err);
  grVn6Default       = new TGraphErrors(NCENT, centBinCenter, Vn6Default,       c_err, vn6Default_err);
  grVn8Default       = new TGraphErrors(NCENT, centBinCenter, Vn8Default,       c_err, vn8Default_err);
  grVn6Vn4Default    = new TGraphErrors(NCENT, centBinCenter, vn6vn4Default,    c_err, vn6vn4Default_err);
  grVn8Vn4Default    = new TGraphErrors(NCENT, centBinCenter, vn8vn4Default,    c_err, vn8vn4Default_err);
  grVn8Vn6Default    = new TGraphErrors(NCENT, centBinCenter, vn8vn6Default,    c_err, vn8vn6Default_err);
  grVn46_Vn68Default = new TGraphErrors(NCENT, centBinCenter, vn46_vn68Default, c_err, vn46_vn68Default_err);
  grG1EDefault       = new TGraphErrors(NCENT, centBinCenter, g1eDefault,       c_err, g1eDefault_err);

  formatGraph(grVn2Default,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{2}", norder_),                                                                  9,         21, "grVn2Default");
  formatGraph(grVn4Default,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{4}", norder_),                                                                  kSpring+4, 21, "grVn4Default");
  formatGraph(grVn6Default,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{6}", norder_),                                                                  6,         21, "grVn6Default");
  formatGraph(grVn8Default,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{8}", norder_),                                                                  kOrange+7, 21, "grVn8Default");
  formatGraph(grVn6Vn4Default,    "Centrality %", vn6vn4Min,    vn6vn4Max,    Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_),                                               4,         21, "grVn6Vn4Default");
  formatGraph(grVn8Vn4Default,    "Centrality %", vn8vn4Min,    vn8vn4Max,    Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_),                                               kGreen+2,  34, "grVn8Vn4Default");
  formatGraph(grVn8Vn6Default,    "Centrality %", vn8vn6Min,    vn8vn6Max,    Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_),                                               kViolet-1, 33, "grVn8Vn6Default");
  formatGraph(grVn46_Vn68Default, "Centrality %", vn46_vn68Min, vn46_vn68Max, Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   22, "grVn46_Vn68Default");
  formatGraph(grG1EDefault,       "Centrality %", g1eMin,       g1eMax,       "#gamma_{1}^{exp}",                                                                          2,         20, "grG1EDefault");

  //-- 2.0%
  grVn22pct       = new TGraphErrors(NCENT, centBinCenter, Vn22pct,       c_err, vn22pct_err);
  grVn42pct       = new TGraphErrors(NCENT, centBinCenter, Vn42pct,       c_err, vn42pct_err);
  grVn62pct       = new TGraphErrors(NCENT, centBinCenter, Vn62pct,       c_err, vn62pct_err);
  grVn82pct       = new TGraphErrors(NCENT, centBinCenter, Vn82pct,       c_err, vn82pct_err);
  grVn6Vn42pct    = new TGraphErrors(NCENT, centBinCenter, vn6vn42pct,    c_err, vn6vn42pct_err);
  grVn8Vn42pct    = new TGraphErrors(NCENT, centBinCenter, vn8vn42pct,    c_err, vn8vn42pct_err);
  grVn8Vn62pct    = new TGraphErrors(NCENT, centBinCenter, vn8vn62pct,    c_err, vn8vn62pct_err);
  grVn46_Vn682pct = new TGraphErrors(NCENT, centBinCenter, vn46_vn682pct, c_err, vn46_vn682pct_err);
  grG1E2pct       = new TGraphErrors(NCENT, centBinCenter, g1e2pct,       c_err, g1e2pct_err);

  formatGraph(grVn22pct,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{2}", norder_),                                                                  9,         25, "grVn22pct");
  formatGraph(grVn42pct,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{4}", norder_),                                                                  kSpring+4, 25, "grVn42pct");
  formatGraph(grVn62pct,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{6}", norder_),                                                                  6,         25, "grVn62pct");
  formatGraph(grVn82pct,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{8}", norder_),                                                                  kOrange+7, 25, "grVn82pct");
  formatGraph(grVn6Vn42pct,    "Centrality %", vn6vn4Min,    vn6vn4Max,    Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_),                                               4,         25, "grVn6Vn42pct");
  formatGraph(grVn8Vn42pct,    "Centrality %", vn8vn4Min,    vn8vn4Max,    Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_),                                               kGreen+2,  28, "grVn8Vn42pct");
  formatGraph(grVn8Vn62pct,    "Centrality %", vn8vn6Min,    vn8vn6Max,    Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_),                                               kViolet-1, 27, "grVn8Vn62pct");
  formatGraph(grVn46_Vn682pct, "Centrality %", vn46_vn68Min, vn46_vn68Max, Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   26, "grVn46_Vn682pct");
  formatGraph(grG1E2pct,       "Centrality %", g1eMin,       g1eMax,       "#gamma_{1}^{exp}",                                                                          2,         24, "grG1E2pct");

  //-- Ratio to Default
  grVn22pct_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn22pct_RatioToDefault,       c_err, vn22pct_RatioToDefault_err);
  grVn42pct_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn42pct_RatioToDefault,       c_err, vn42pct_RatioToDefault_err);
  grVn62pct_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn62pct_RatioToDefault,       c_err, vn62pct_RatioToDefault_err);
  grVn82pct_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn82pct_RatioToDefault,       c_err, vn82pct_RatioToDefault_err);
  grVn6Vn42pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn6vn42pct_RatioToDefault,    c_err, vn6vn42pct_RatioToDefault_err);
  grVn8Vn42pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn8vn42pct_RatioToDefault,    c_err, vn8vn42pct_RatioToDefault_err);
  grVn8Vn62pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn8vn62pct_RatioToDefault,    c_err, vn8vn62pct_RatioToDefault_err);
  grVn46_Vn682pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn46_vn682pct_RatioToDefault, c_err, vn46_vn682pct_RatioToDefault_err);
  grG1E2pct_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, g1e2pct_RatioToDefault,       c_err, g1e2pct_RatioToDefault_err);

  formatGraph(grVn22pct_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{2} Ratio", norder_),                                                            9,         24, "grVn22pct_RatioToDefault");
  formatGraph(grVn42pct_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{4} Ratio", norder_),                                                            kSpring+4, 24, "grVn42pct_RatioToDefault");
  formatGraph(grVn62pct_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{6} Ratio", norder_),                                                            6,         24, "grVn62pct_RatioToDefault");
  formatGraph(grVn82pct_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{8} Ratio", norder_),                                                            kOrange+7, 24, "grVn82pct_RatioToDefault");
  formatGraph(grVn6Vn42pct_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{6}/v_{%i}{4} Ratio", norder_, norder_),                                         4,         24, "grVn6Vn42pct_RatioToDefault");
  formatGraph(grVn8Vn42pct_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{8}/v_{%i}{4} Ratio", norder_, norder_),                                         kGreen+2,  24, "grVn8Vn42pct_RatioToDefault");
  formatGraph(grVn8Vn62pct_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{8}/v_{%i}{6} Ratio", norder_, norder_),                                         kViolet-1, 24, "grVn8Vn62pct_RatioToDefault");
  formatGraph(grVn46_Vn682pct_RatioToDefault, "Centrality %", 0.85,   1.15,   Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   24, "grVn46_Vn682pct_RatioToDefault");
  formatGraph(grG1E2pct_RatioToDefault,       "Centrality %", g1rMin, g1rMax, "|1 - #gamma_{1}^{exp} Ratio|",                                                              2,         24, "grG1E2pct_RatioToDefault");

  //-- DRAW
  TLine * lOne = new TLine(0, 1.0, grVn2Default->GetXaxis()->GetXmax(), 1.0);
  lOne->SetLineStyle(2);
  lOne->SetLineWidth(2);

  //-- cumu
  TLegend * legvn2  = new TLegend(0.68, 0.22, 0.97, 0.39);
  legvn2->SetBorderSize(0);
  legvn2->SetFillStyle(0);
  legvn2->AddEntry(grVn2Default, "Default", "lp");
  legvn2->AddEntry(grVn22pct,    "2%",      "lp");

  TLegend * legvn4  = new TLegend(0.68, 0.22, 0.97, 0.39);
  legvn4->SetBorderSize(0);
  legvn4->SetFillStyle(0);
  legvn4->AddEntry(grVn4Default, "Default", "lp");
  legvn4->AddEntry(grVn42pct,    "2%",      "lp");

  TLegend * legvn6  = new TLegend(0.68, 0.22, 0.97, 0.39);
  legvn6->SetBorderSize(0);
  legvn6->SetFillStyle(0);
  legvn6->AddEntry(grVn6Default, "Default", "lp");
  legvn6->AddEntry(grVn62pct,    "2%",      "lp");

  TLegend * legvn8  = new TLegend(0.68, 0.22, 0.97, 0.39);
  legvn8->SetBorderSize(0);
  legvn8->SetFillStyle(0);
  legvn8->AddEntry(grVn8Default, "Default", "lp");
  legvn8->AddEntry(grVn82pct,    "2%",      "lp");

  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 1000, 1000);
  cCumu->Divide(2,2);

  cCumu->cd(1);
  grVn2Default->Draw("ap");
  grVn22pct->Draw("psame");
  legvn2->Draw("same");

  cCumu->cd(2);
  grVn4Default->Draw("ap");
  grVn42pct->Draw("psame");
  legvn4->Draw("same");

  cCumu->cd(3);
  grVn6Default->Draw("ap");
  grVn62pct->Draw("psame");
  legvn6->Draw("same");

  cCumu->cd(4);
  grVn8Default->Draw("ap");
  grVn82pct->Draw("psame");
  legvn8->Draw("same");

  cCumu->SaveAs("cCumu.pdf");

  //-- Cumu Ratio to Default
  TCanvas * cCumuRatToDef = new TCanvas("cCumuRatToDef", "cCumuRatToDef", 1000, 1000);
  cCumuRatToDef->Divide(2,2);

  cCumuRatToDef->cd(1);
  grVn22pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(2);
  grVn42pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(3);
  grVn62pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(4);
  grVn82pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->SaveAs("../../plots/systematicStudies/cSysNewCC_CumuCent.pdf");


  //-- cumuRatio
  TLegend * legvn6vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Default, "Default", "lp");
  legvn6vn4->AddEntry(grVn6Vn42pct,    "2%",      "lp");

  TLegend * legvn8vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Default, "Default", "lp");
  legvn8vn4->AddEntry(grVn8Vn42pct,    "2%",      "lp");

  TLegend * legvn8vn6 = new TLegend(0.4520, 0.1973, 0.94, 0.3703);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Default, "Default", "lp");
  legvn8vn6->AddEntry(grVn8Vn62pct,    "2%",      "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);

  cCumuRatio->cd(1);
  grVn6Vn4Default->Draw("ap");
  grVn6Vn42pct->Draw("psame");
  legvn6vn4->Draw("same");

  cCumuRatio->cd(2);
  grVn8Vn4Default->Draw("ap");
  grVn8Vn42pct->Draw("psame");
  legvn8vn4->Draw("same");

  cCumuRatio->cd(3);
  grVn8Vn6Default->Draw("ap");
  grVn8Vn62pct->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- cumuRatio To Default
  TCanvas * cCumuRatio_RatToDef = new TCanvas("cCumuRatio_RatToDef", "cCumuRatio_RatToDef", 1500, 500);
  cCumuRatio_RatToDef->Divide(3,1);

  cCumuRatio_RatToDef->cd(1);
  cCumuRatio_RatToDef->cd(1)->SetLeftMargin(0.2);
  grVn6Vn42pct_RatioToDefault->GetYaxis()->SetTitleOffset(1.6);
  grVn6Vn42pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->cd(2);
  grVn8Vn42pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->cd(3);
  grVn8Vn62pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->SaveAs("../../plots/systematicStudies/cSysNewCC_CumuRatioCent.pdf");

  //-- g1e
  TLegend * legG1E = new TLegend(0.4538, 0.7421, 0.9529, 0.9172);
  legG1E->SetBorderSize(0);
  legG1E->SetFillStyle(0);
  legG1E->AddEntry(grG1EDefault, "Default", "lp");
  legG1E->AddEntry(grG1E2pct,    "2%",      "lp");

  TCanvas* cG1e = new TCanvas("cG1e", "cG1e", 500, 500);
  cG1e->cd();
  grG1EDefault->Draw("ap");
  grG1E2pct->Draw("psame");
  legG1E->Draw("same");
  cG1e->SaveAs("cG1e.pdf");

  //-- G1E ratio to default
  TCanvas* cG1eR = new TCanvas("cG1eR", "cG1eR", 500, 500);
  cG1eR->cd();
  grG1E2pct_RatioToDefault->Draw("ap");
  //lOne->Draw("same");
  cG1eR->SaveAs("../../plots/systematicStudies/cSysNewCC_G1ECent.pdf");

  //-- vn46_vn68 
  TCanvas* cvn46_vn68 = new TCanvas("cvn46_vn68", "cvn46_vn68", 500, 500);
  cvn46_vn68->cd();
  grVn46_Vn682pct_RatioToDefault->Draw("ap");
  lOne->Draw("same");
  cvn46_vn68->SaveAs("../../plots/systematicStudies/cSysNewCC_Vn46_Vn68.pdf");

  //----------------------------------------------------------------------------------------------------
  //-- Save plots for smoothing
  TFile * fSave = new TFile("SysNewCC.root", "recreate");
  fSave->cd();
  grVn22pct_RatioToDefault->Write("grVn22pct_RatioToDefault");
  grVn42pct_RatioToDefault->Write("grVn42pct_RatioToDefault");
  grVn62pct_RatioToDefault->Write("grVn62pct_RatioToDefault");
  grVn82pct_RatioToDefault->Write("grVn82pct_RatioToDefault");
  grVn6Vn42pct_RatioToDefault->Write("grVn6Vn42pct_RatioToDefault");
  grVn8Vn42pct_RatioToDefault->Write("grVn8Vn42pct_RatioToDefault");
  grVn8Vn62pct_RatioToDefault->Write("grVn8Vn62pct_RatioToDefault");
  grVn46_Vn682pct_RatioToDefault->Write("grVn46_Vn682pct_RatioToDefault");
  grG1E2pct_RatioToDefault->Write("grG1E2pct_RatioToDefault");

  //-- Save the final unfolded distns from the 2% study
  for(int icent = 0; icent < NCENT; icent++){
    int i = iterCut2pct[icent];
    std::cout<<i<<std::endl;
    hUnfold2pct[icent][i]->SetLineColor(1);
    hUnfold2pct[icent][i]->SetMarkerColor(1);
    hUnfold2pct[icent][i]->Write( Form("hFinalUnfold_SysNewCC_c%i", icent) );
  }


  //-- Mult Distns
  /*
  TCanvas * cMult = new TCanvas("cMult", "cMult", 1000, 1000);
  cMult->Divide(3,2);
  cMult->cd(1)->SetLogy();
  cMult->cd(2)->SetLogy();
  cMult->cd(3)->SetLogy();
  cMult->cd(4)->SetLogy();
  cMult->cd(5)->SetLogy();

  cMult->cd(1);
  hMultDefault[0]->Draw();
  latex.DrawLatex(0.5, 0.88, "Default Tune");
  cMult->cd(2);
  hMult0p5pct[0]->Draw();
  latex.DrawLatex(0.5, 0.88, "0.5% Tune");
  cMult->cd(3);
  hMult1pct[0]->Draw();
  latex.DrawLatex(0.5, 0.88, "1% Tune");
  cMult->cd(4);
  hMult2pct[0]->Draw();
  latex.DrawLatex(0.5, 0.88, "2% Tune");
  cMult->cd(5);
  hMult3pct[0]->Draw();
  latex.DrawLatex(0.5, 0.88, "3% Tune");

  for( int icent = 1; icent < NCENT; icent++){
    cMult->cd(1);
    hMultDefault[icent]->Draw("same");
    cMult->cd(2);
    hMult0p5pct[icent]->Draw("same");
    cMult->cd(3);
    hMult1pct[icent]->Draw("same");
    cMult->cd(4);
    hMult2pct[icent]->Draw("same");
    cMult->cd(5);
    hMult3pct[icent]->Draw("same");
  }
  */

}

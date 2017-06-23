#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void sysPtRes(){

  const int norder_ = 2;
  const double tkEta = 1.0;

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

  double rMin = 0.98;
  double rMax = 1.02;
  double g1rMin = 0.7;
  double g1rMax = 1.3;

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
  TFile * fAnaSmearPt5;
  TFile * fUnfSmearPt5;
  TH1D * hObsSmearPt5[NCENT];
  TH1D * hMultSmearPt5[NCENT];
  TH1D * hUnfoldSmearPt5[NCENT][NITER];
  TH1D * hRefoldSmearPt5[NCENT][NITER];
  bool stopFoundSmearPt5;
  int iterCutSmearPt5[NCENT];

  TFile * fStatSmearPt5;
  TH1D * hVarianceOfMean_Vn2_SmearPt5;
  TH1D * hVarianceOfMean_Vn4_SmearPt5;
  TH1D * hVarianceOfMean_Vn6_SmearPt5;
  TH1D * hVarianceOfMean_Vn8_SmearPt5;
  TH1D * hVarianceOfMean_Vn6Vn4_SmearPt5;
  TH1D * hVarianceOfMean_Vn8Vn4_SmearPt5;
  TH1D * hVarianceOfMean_Vn8Vn6_SmearPt5;
  TH1D * hVarianceOfMean_Vn46_Vn68_SmearPt5;
  TH1D * hVarianceOfMean_Gamma1Exp_SmearPt5;

  double Vn2SmearPt5[NCENT];
  double Vn4SmearPt5[NCENT];
  double Vn6SmearPt5[NCENT];
  double Vn8SmearPt5[NCENT];
  double vn6vn4SmearPt5[NCENT];
  double vn8vn4SmearPt5[NCENT];
  double vn8vn6SmearPt5[NCENT];
  double vn46_vn68SmearPt5[NCENT];
  double g1eSmearPt5[NCENT];

  double vn2SmearPt5_err[NCENT];
  double vn4SmearPt5_err[NCENT];
  double vn6SmearPt5_err[NCENT];
  double vn8SmearPt5_err[NCENT];
  double vn6vn4SmearPt5_err[NCENT];
  double vn8vn4SmearPt5_err[NCENT];
  double vn8vn6SmearPt5_err[NCENT];
  double vn46_vn68SmearPt5_err[NCENT];
  double g1eSmearPt5_err[NCENT];

  double vn2SmearPt5_RatioToDefault[NCENT];
  double vn4SmearPt5_RatioToDefault[NCENT];
  double vn6SmearPt5_RatioToDefault[NCENT];
  double vn8SmearPt5_RatioToDefault[NCENT];
  double vn6vn4SmearPt5_RatioToDefault[NCENT];
  double vn8vn4SmearPt5_RatioToDefault[NCENT];
  double vn8vn6SmearPt5_RatioToDefault[NCENT];
  double vn46_vn68SmearPt5_RatioToDefault[NCENT];
  double g1eSmearPt5_RatioToDefault[NCENT];

  double vn2SmearPt5_RatioToDefault_err[NCENT];
  double vn4SmearPt5_RatioToDefault_err[NCENT];
  double vn6SmearPt5_RatioToDefault_err[NCENT];
  double vn8SmearPt5_RatioToDefault_err[NCENT];
  double vn6vn4SmearPt5_RatioToDefault_err[NCENT];
  double vn8vn4SmearPt5_RatioToDefault_err[NCENT];
  double vn8vn6SmearPt5_RatioToDefault_err[NCENT];
  double vn46_vn68SmearPt5_RatioToDefault_err[NCENT];
  double g1eSmearPt5_RatioToDefault_err[NCENT];

  TGraphErrors * grVn2SmearPt5;
  TGraphErrors * grVn4SmearPt5;
  TGraphErrors * grVn6SmearPt5;
  TGraphErrors * grVn8SmearPt5;
  TGraphErrors * grVn6Vn4SmearPt5;
  TGraphErrors * grVn8Vn4SmearPt5;
  TGraphErrors * grVn8Vn6SmearPt5;
  TGraphErrors * grVn46_Vn68SmearPt5;
  TGraphErrors * grG1ESmearPt5;

  TGraphErrors * grVn2SmearPt5_RatioToDefault;
  TGraphErrors * grVn4SmearPt5_RatioToDefault;
  TGraphErrors * grVn6SmearPt5_RatioToDefault;
  TGraphErrors * grVn8SmearPt5_RatioToDefault;
  TGraphErrors * grVn6Vn4SmearPt5_RatioToDefault;
  TGraphErrors * grVn8Vn4SmearPt5_RatioToDefault;
  TGraphErrors * grVn8Vn6SmearPt5_RatioToDefault;
  TGraphErrors * grVn46_Vn68SmearPt5_RatioToDefault;
  TGraphErrors * grG1ESmearPt5_RatioToDefault;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get Analyzer files
  fAnaDefault = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fAnaSmearPt5    = new TFile( "AnalyzerResults/CastleEbyE.root" );

  //-- Get Unfold files
  fUnfDefault = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnfSmearPt5    = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_)  );

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

  fStatSmearPt5  = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/ptResolution/StatisticalUncertainties_v%i.root", norder_, tkEta, norder_) );
  hVarianceOfMean_Vn2_SmearPt5       = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn2" );
  hVarianceOfMean_Vn4_SmearPt5       = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn4" );
  hVarianceOfMean_Vn6_SmearPt5       = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn6" );
  hVarianceOfMean_Vn8_SmearPt5       = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn8" );
  hVarianceOfMean_Vn6Vn4_SmearPt5    = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_SmearPt5    = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_SmearPt5    = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Vn46_Vn68_SmearPt5 = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Vn46_Vn68" );
  hVarianceOfMean_Gamma1Exp_SmearPt5 = (TH1D*) fStatSmearPt5->Get( "hVarianceOfMean_Gamma1Exp" );

  //-- Start Grabbing Histos
  for(int icent = 0; icent < NCENT; icent++){

    //-- Vn Observed 
    hObsDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsDefault[icent]->SetName( Form("hVnFullDefault_c%i", icent) );
    hObsSmearPt5[icent] = (TH1D*) fAnaSmearPt5->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsSmearPt5[icent]->SetName( Form("hVnFullSmearPt5_c%i", icent) );

    //-- Multiplicity
    hMultDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/Mult_c%i", icent) );
    hMultDefault[icent]->SetLineColor(centCol[icent]);
    hMultDefault[icent]->SetMarkerColor(centCol[icent]);
    hMultDefault[icent]->SetMaximum( 20.*hMultDefault[icent]->GetMaximum() );
    hMultDefault[icent]->GetXaxis()->SetNdivisions(509);
    hMultDefault[icent]->GetXaxis()->SetRange(1, hMultDefault[icent]->FindBin(5500));

    hMultSmearPt5[icent] = (TH1D*) fAnaSmearPt5->Get( Form("qwebye/Mult_c%i", icent) );
    hMultSmearPt5[icent]->SetLineColor(centCol[icent]);
    hMultSmearPt5[icent]->SetMarkerColor(centCol[icent]);
    hMultSmearPt5[icent]->SetMaximum( 20.*hMultSmearPt5[icent]->GetMaximum() );
    hMultSmearPt5[icent]->GetXaxis()->SetNdivisions(509);
    hMultSmearPt5[icent]->GetXaxis()->SetRange(1, hMultSmearPt5[icent]->FindBin(5500));

    //-- Reset stopFound boolean
    stopFoundDefault = 0;
    stopFoundSmearPt5    = 0;

    //-- Grab un/refolded distns
    for(int i = 0; i < NITER; i++){

      //-- Default
      hUnfoldDefault[icent][i] = (TH1D*) fUnfDefault->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldDefault[icent][i]->SetName( Form("hrecoDefault%i_c%i", iter[i], icent) );
      hRefoldDefault[icent][i] = (TH1D*) fUnfDefault->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldDefault[icent][i]->SetName( Form("hrefolDefault%i_c%i", iter[i], icent) );

      double chi2Default = hRefoldDefault[icent][i]->Chi2Test(hObsDefault[icent], "CHI2/NDF");
      if( chi2Default <= 1.2 && !stopFoundDefault ){
	stopFoundDefault = 1;
	iterCutDefault[icent] = i;
      }
      if( i == NITER-1 && !stopFoundDefault ){
	stopFoundDefault = 1;
        iterCutDefault[icent] = i;
      }

      //-- 2.0%
      hUnfoldSmearPt5[icent][i] = (TH1D*) fUnfSmearPt5->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldSmearPt5[icent][i]->SetName( Form("hrecoSmearPt5%i_c%i", iter[i], icent) );
      hRefoldSmearPt5[icent][i] = (TH1D*) fUnfSmearPt5->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldSmearPt5[icent][i]->SetName( Form("hrefolSmearPt5%i_c%i", iter[i], icent) );

      double chi2SmearPt5 = hRefoldSmearPt5[icent][i]->Chi2Test(hObsSmearPt5[icent], "CHI2/NDF");
      if( chi2SmearPt5 <= 1.2 && !stopFoundSmearPt5 ){
        stopFoundSmearPt5 = 1;
        iterCutSmearPt5[icent] = i;
      }
      if( i == NITER-1 && !stopFoundSmearPt5 ){
        stopFoundSmearPt5 = 1;
        iterCutSmearPt5[icent] = i;
      }

    } //-- End iter loop

  } //-- End cent loop


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
    int iSmearPt5 = iterCutSmearPt5[icent];
    FixUnfold(hUnfoldSmearPt5[icent][iSmearPt5]);
    EbyECumu cumuSmearPt5(hUnfoldSmearPt5[icent][iSmearPt5]);
    double vn2SmearPt5  = cumuSmearPt5.GetCumu_vn2();
    double vn4SmearPt5  = cumuSmearPt5.GetCumu_vn4();
    double vn6SmearPt5  = cumuSmearPt5.GetCumu_vn6();
    double vn8SmearPt5  = cumuSmearPt5.GetCumu_vn8();
    double g1exSmearPt5 = cumuSmearPt5.GetGamma1Exp();

    Vn2SmearPt5[icent] = vn2SmearPt5;
    Vn4SmearPt5[icent] = vn4SmearPt5;
    Vn6SmearPt5[icent] = vn6SmearPt5;
    Vn8SmearPt5[icent] = vn8SmearPt5;

    if(vn4SmearPt5 == 0 || vn6SmearPt5 == 0) vn6vn4SmearPt5[icent] = -1000.;
    else                             vn6vn4SmearPt5[icent] = vn6SmearPt5 / vn4SmearPt5;
    if(vn4SmearPt5 == 0 || vn8SmearPt5 == 0) vn8vn4SmearPt5[icent] = -1000.;
    else                             vn8vn4SmearPt5[icent] = vn8SmearPt5 / vn4SmearPt5;
    if(vn6SmearPt5 == 0 || vn8SmearPt5 == 0) vn8vn6SmearPt5[icent] = -1000.;
    else                             vn8vn6SmearPt5[icent] = vn8SmearPt5 / vn6SmearPt5;
    if(vn4SmearPt5 == 0 || vn6SmearPt5 == 0 || vn8SmearPt5 == 0) vn46_vn68SmearPt5[icent]= -1000.;
    else                                             vn46_vn68SmearPt5[icent] = (vn4SmearPt5 - vn6SmearPt5)/(vn6SmearPt5 - vn8SmearPt5);
    g1eSmearPt5[icent] = g1exSmearPt5;

    double vn2SmearPt5_staterr       = sqrt( hVarianceOfMean_Vn2_SmearPt5->GetBinContent(icent+1) );
    double vn4SmearPt5_staterr       = sqrt( hVarianceOfMean_Vn4_SmearPt5->GetBinContent(icent+1) );
    double vn6SmearPt5_staterr       = sqrt( hVarianceOfMean_Vn6_SmearPt5->GetBinContent(icent+1) );
    double vn8SmearPt5_staterr       = sqrt( hVarianceOfMean_Vn8_SmearPt5->GetBinContent(icent+1) );
    double vn6vn4SmearPt5_staterr    = sqrt( hVarianceOfMean_Vn6Vn4_SmearPt5->GetBinContent(icent+1) );
    double vn8vn4SmearPt5_staterr    = sqrt( hVarianceOfMean_Vn8Vn4_SmearPt5->GetBinContent(icent+1) );
    double vn8vn6SmearPt5_staterr    = sqrt( hVarianceOfMean_Vn8Vn6_SmearPt5->GetBinContent(icent+1) );
    double vn46_vn68SmearPt5_staterr = sqrt( hVarianceOfMean_Vn46_Vn68_SmearPt5->GetBinContent(icent+1) );
    double g1eSmearPt5_staterr       = sqrt( hVarianceOfMean_Gamma1Exp_SmearPt5->GetBinContent(icent+1) );
 
    vn2SmearPt5_err[icent]       = vn2SmearPt5_staterr;
    vn4SmearPt5_err[icent]       = vn4SmearPt5_staterr;
    vn6SmearPt5_err[icent]       = vn6SmearPt5_staterr;
    vn8SmearPt5_err[icent]       = vn8SmearPt5_staterr;
    vn6vn4SmearPt5_err[icent]    = vn6vn4SmearPt5_staterr;
    vn8vn4SmearPt5_err[icent]    = vn8vn4SmearPt5_staterr;
    vn8vn6SmearPt5_err[icent]    = vn8vn6SmearPt5_staterr;
    vn46_vn68SmearPt5_err[icent] = vn46_vn68SmearPt5_staterr;
    g1eSmearPt5_err[icent]       = g1eSmearPt5_staterr;

    //-- Ratio to default
    if(vn2SmearPt5 == 0 || vn2Default == 0) vn2SmearPt5_RatioToDefault[icent] = 0.;
    else                                vn2SmearPt5_RatioToDefault[icent] = vn2SmearPt5 / vn2Default;
    if(vn4SmearPt5 == 0 || vn4Default == 0) vn4SmearPt5_RatioToDefault[icent] = 0.;
    else                                vn4SmearPt5_RatioToDefault[icent] = vn4SmearPt5/ vn4Default;
    if(vn6SmearPt5 == 0 || vn6Default == 0) vn6SmearPt5_RatioToDefault[icent] = 0.;
    else                                vn6SmearPt5_RatioToDefault[icent] = vn6SmearPt5/ vn6Default;
    if(vn8SmearPt5 == 0 || vn8Default == 0) vn8SmearPt5_RatioToDefault[icent] = 0.;
    else                                vn8SmearPt5_RatioToDefault[icent] = vn8SmearPt5/ vn8Default;

    if(vn4SmearPt5 == 0 || vn6SmearPt5 == 0 || vn4Default == 0 || vn6Default == 0) vn6vn4SmearPt5_RatioToDefault[icent] = 0.;
    else                                                                           vn6vn4SmearPt5_RatioToDefault[icent] = vn6vn4SmearPt5[icent] / vn6vn4Default[icent];
    if(vn4SmearPt5 == 0 || vn8SmearPt5 == 0 || vn4Default == 0 || vn8Default == 0) vn8vn4SmearPt5_RatioToDefault[icent] = 0.;
    else                                                                           vn8vn4SmearPt5_RatioToDefault[icent] = vn8vn4SmearPt5[icent] / vn8vn4Default[icent];
    if(vn6SmearPt5 == 0 || vn8SmearPt5 == 0 || vn6Default == 0 || vn8Default == 0) vn8vn6SmearPt5_RatioToDefault[icent] = 0.;
    else                                                                           vn8vn6SmearPt5_RatioToDefault[icent] = vn8vn6SmearPt5[icent] / vn8vn6Default[icent];
    if(vn46_vn68SmearPt5[icent] <= 0 || vn46_vn68Default[icent] <= 0) vn46_vn68SmearPt5_RatioToDefault[icent] = 0.;
    else                                                              vn46_vn68SmearPt5_RatioToDefault[icent] = vn46_vn68SmearPt5[icent] / vn46_vn68Default[icent];
    if(g1exSmearPt5== -10000. || g1exDefault == -10000.) g1eSmearPt5_RatioToDefault[icent] = 0.;
    else                                                 g1eSmearPt5_RatioToDefault[icent]= g1eSmearPt5[icent]/ g1eDefault[icent];

    vn2SmearPt5_RatioToDefault_err[icent]       = sqrt( pow(vn2SmearPt5_err[icent]/Vn2Default[icent],2) + pow(Vn2SmearPt5[icent]*vn2Default_err[icent]/Vn2Default[icent]/Vn2Default[icent],2) );
    vn4SmearPt5_RatioToDefault_err[icent]       = sqrt( pow(vn4SmearPt5_err[icent]/Vn4Default[icent],2) + pow(Vn4SmearPt5[icent]*vn4Default_err[icent]/Vn4Default[icent]/Vn4Default[icent],2) );
    vn6SmearPt5_RatioToDefault_err[icent]       = sqrt( pow(vn6SmearPt5_err[icent]/Vn6Default[icent],2) + pow(Vn6SmearPt5[icent]*vn6Default_err[icent]/Vn6Default[icent]/Vn6Default[icent],2) );
    vn8SmearPt5_RatioToDefault_err[icent]       = sqrt( pow(vn8SmearPt5_err[icent]/Vn8Default[icent],2) + pow(Vn8SmearPt5[icent]*vn8Default_err[icent]/Vn8Default[icent]/Vn8Default[icent],2) );
    vn6vn4SmearPt5_RatioToDefault_err[icent]    = sqrt( pow(vn6vn4SmearPt5_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn4SmearPt5[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn4SmearPt5_RatioToDefault_err[icent]    = sqrt( pow(vn8vn4SmearPt5_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn4SmearPt5[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn6SmearPt5_RatioToDefault_err[icent]    = sqrt( pow(vn8vn6SmearPt5_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn6SmearPt5[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    vn46_vn68SmearPt5_RatioToDefault_err[icent] = sqrt( pow(vn46_vn68SmearPt5_err[icent]/vn46_vn68Default[icent],2) + pow(vn46_vn68SmearPt5[icent]*vn46_vn68Default_err[icent]/vn46_vn68Default[icent]/vn46_vn68Default[icent],2) );
    g1eSmearPt5_RatioToDefault_err[icent]       = sqrt( pow(g1eSmearPt5_err[icent]/g1eDefault[icent],2) + pow(g1eSmearPt5[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

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
  formatGraph(grVn8Vn4Default,    "Centrality %", vn8vn4Min,    vn8vn4Max,    Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_),                                               kGreen+2,  21, "grVn8Vn4Default");
  formatGraph(grVn8Vn6Default,    "Centrality %", vn8vn6Min,    vn8vn6Max,    Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_),                                               kViolet-1, 21, "grVn8Vn6Default");
  formatGraph(grVn46_Vn68Default, "Centrality %", vn46_vn68Min, vn46_vn68Max, Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   21, "grVn46_Vn68Default");
  formatGraph(grG1EDefault,       "Centrality %", g1eMin,       g1eMax,       "#gamma_{1}^{exp}",                                                                          2,         21, "grG1EDefault");

  //-- 2.0%
  grVn2SmearPt5       = new TGraphErrors(NCENT, centBinCenter, Vn2SmearPt5,       c_err, vn2SmearPt5_err);
  grVn4SmearPt5       = new TGraphErrors(NCENT, centBinCenter, Vn4SmearPt5,       c_err, vn4SmearPt5_err);
  grVn6SmearPt5       = new TGraphErrors(NCENT, centBinCenter, Vn6SmearPt5,       c_err, vn6SmearPt5_err);
  grVn8SmearPt5       = new TGraphErrors(NCENT, centBinCenter, Vn8SmearPt5,       c_err, vn8SmearPt5_err);
  grVn6Vn4SmearPt5    = new TGraphErrors(NCENT, centBinCenter, vn6vn4SmearPt5,    c_err, vn6vn4SmearPt5_err);
  grVn8Vn4SmearPt5    = new TGraphErrors(NCENT, centBinCenter, vn8vn4SmearPt5,    c_err, vn8vn4SmearPt5_err);
  grVn8Vn6SmearPt5    = new TGraphErrors(NCENT, centBinCenter, vn8vn6SmearPt5,    c_err, vn8vn6SmearPt5_err);
  grVn46_Vn68SmearPt5 = new TGraphErrors(NCENT, centBinCenter, vn46_vn68SmearPt5, c_err, vn46_vn68SmearPt5_err);
  grG1ESmearPt5       = new TGraphErrors(NCENT, centBinCenter, g1eSmearPt5,       c_err, g1eSmearPt5_err);

  formatGraph(grVn2SmearPt5,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{2}", norder_),                                                                  9,         24, "grVn2SmearPt5");
  formatGraph(grVn4SmearPt5,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{4}", norder_),                                                                  kSpring+4, 24, "grVn4SmearPt5");
  formatGraph(grVn6SmearPt5,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{6}", norder_),                                                                  6,         24, "grVn6SmearPt5");
  formatGraph(grVn8SmearPt5,       "Centrality %", cumuMin,      cumuMax,      Form("v_{%i}{8}", norder_),                                                                  kOrange+7, 24, "grVn8SmearPt5");
  formatGraph(grVn6Vn4SmearPt5,    "Centrality %", vn6vn4Min,    vn6vn4Max,    Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_),                                               4,         24, "grVn6Vn4SmearPt5");
  formatGraph(grVn8Vn4SmearPt5,    "Centrality %", vn8vn4Min,    vn8vn4Max,    Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_),                                               kGreen+2,  24, "grVn8Vn4SmearPt5");
  formatGraph(grVn8Vn6SmearPt5,    "Centrality %", vn8vn6Min,    vn8vn6Max,    Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_),                                               kViolet-1, 34, "grVn8Vn6SmearPt5");
  formatGraph(grVn46_Vn68SmearPt5, "Centrality %", vn46_vn68Min, vn46_vn68Max, Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   24, "grVn46_Vn68SmearPt5");
  formatGraph(grG1ESmearPt5,       "Centrality %", g1eMin,       g1eMax,       "#gamma_{1}^{exp}",                                                                          2,         24, "grG1ESmearPt5");

  //-- Ratio to Default
  grVn2SmearPt5_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn2SmearPt5_RatioToDefault,       c_err, vn2SmearPt5_RatioToDefault_err);
  grVn4SmearPt5_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn4SmearPt5_RatioToDefault,       c_err, vn4SmearPt5_RatioToDefault_err);
  grVn6SmearPt5_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn6SmearPt5_RatioToDefault,       c_err, vn6SmearPt5_RatioToDefault_err);
  grVn8SmearPt5_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, vn8SmearPt5_RatioToDefault,       c_err, vn8SmearPt5_RatioToDefault_err);
  grVn6Vn4SmearPt5_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn6vn4SmearPt5_RatioToDefault,    c_err, vn6vn4SmearPt5_RatioToDefault_err);
  grVn8Vn4SmearPt5_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn8vn4SmearPt5_RatioToDefault,    c_err, vn8vn4SmearPt5_RatioToDefault_err);
  grVn8Vn6SmearPt5_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, vn8vn6SmearPt5_RatioToDefault,    c_err, vn8vn6SmearPt5_RatioToDefault_err);
  grVn46_Vn68SmearPt5_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn46_vn68SmearPt5_RatioToDefault, c_err, vn46_vn68SmearPt5_RatioToDefault_err);
  grG1ESmearPt5_RatioToDefault       = new TGraphErrors(NCENT, centBinCenter, g1eSmearPt5_RatioToDefault,       c_err, g1eSmearPt5_RatioToDefault_err);

  formatGraph(grVn2SmearPt5_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{2} Ratio", norder_),                                                            9,         24, "grVn2SmearPt5_RatioToDefault");
  formatGraph(grVn4SmearPt5_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{4} Ratio", norder_),                                                            kSpring+4, 24, "grVn4SmearPt5_RatioToDefault");
  formatGraph(grVn6SmearPt5_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{6} Ratio", norder_),                                                            6,         24, "grVn6SmearPt5_RatioToDefault");
  formatGraph(grVn8SmearPt5_RatioToDefault,       "Centrality %", rMin,   rMax,   Form("v_{%i}{8} Ratio", norder_),                                                            kOrange+7, 24, "grVn8SmearPt5_RatioToDefault");
  formatGraph(grVn6Vn4SmearPt5_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{6}/v_{%i}{4} Ratio", norder_, norder_),                                         4,         24, "grVn6Vn4SmearPt5_RatioToDefault");
  formatGraph(grVn8Vn4SmearPt5_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{8}/v_{%i}{4} Ratio", norder_, norder_),                                         kGreen+2,  24, "grVn8Vn4SmearPt5_RatioToDefault");
  formatGraph(grVn8Vn6SmearPt5_RatioToDefault,    "Centrality %", rMin,   rMax,   Form("v_{%i}{8}/v_{%i}{6} Ratio", norder_, norder_),                                         kViolet-1, 24, "grVn8Vn6SmearPt5_RatioToDefault");
  formatGraph(grVn46_Vn68SmearPt5_RatioToDefault, "Centrality %", 0.85,   1.15,   Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8})", norder_, norder_, norder_, norder_), kGray+2,   24, "grVn46_Vn68SmearPt5_RatioToDefault");
  formatGraph(grG1ESmearPt5_RatioToDefault,       "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",                                                                    2,         24, "grG1ESmearPt5_RatioToDefault");

  //-- DRAW
  TLine * lOne = new TLine(0, 1.0, grVn2Default->GetXaxis()->GetXmax(), 1.0);
  lOne->SetLineStyle(2);
  lOne->SetLineWidth(2);

  //-- cumu
  TLegend * legvn2  = new TLegend(0.47, 0.23, 0.96, 0.41);
  legvn2->SetBorderSize(0);
  legvn2->SetFillStyle(0);
  legvn2->AddEntry(grVn2Default,  "Default",        "lp");
  legvn2->AddEntry(grVn2SmearPt5, "5% p_{T} smear", "lp");

  TLegend * legvn4  = new TLegend(0.47, 0.23, 0.96, 0.41);
  legvn4->SetBorderSize(0);
  legvn4->SetFillStyle(0);
  legvn4->AddEntry(grVn4Default,  "Default",        "lp");
  legvn4->AddEntry(grVn4SmearPt5, "5% p_{T} smear", "lp");

  TLegend * legvn6  = new TLegend(0.47, 0.23, 0.96, 0.41);
  legvn6->SetBorderSize(0);
  legvn6->SetFillStyle(0);
  legvn6->AddEntry(grVn6Default,  "Default",        "lp");
  legvn6->AddEntry(grVn6SmearPt5, "5% p_{T} smear", "lp");

  TLegend * legvn8  = new TLegend(0.47, 0.23, 0.96, 0.41);
  legvn8->SetBorderSize(0);
  legvn8->SetFillStyle(0);
  legvn8->AddEntry(grVn8Default,  "Default",        "lp");
  legvn8->AddEntry(grVn8SmearPt5, "5% p_{T} smear", "lp");

  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 1000, 1000);
  cCumu->Divide(2,2);

  cCumu->cd(1);
  grVn2Default->Draw("ap");
  grVn2SmearPt5->Draw("psame");
  legvn2->Draw("same");

  cCumu->cd(2);
  grVn4Default->Draw("ap");
  grVn4SmearPt5->Draw("psame");
  legvn4->Draw("same");

  cCumu->cd(3);
  grVn6Default->Draw("ap");
  grVn6SmearPt5->Draw("psame");
  legvn6->Draw("same");

  cCumu->cd(4);
  grVn8Default->Draw("ap");
  grVn8SmearPt5->Draw("psame");
  legvn8->Draw("same");

  cCumu->SaveAs("cCumu.pdf");

  //-- Cumu Ratio to Default
  TCanvas * cCumuRatToDef = new TCanvas("cCumuRatToDef", "cCumuRatToDef", 1000, 1000);
  cCumuRatToDef->Divide(2,2);

  double mar = 0.2;
  double offs = 1.6;

  grVn2SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);
  grVn4SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);
  grVn6SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);
  grVn8SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);

  cCumuRatToDef->cd(1)->SetLeftMargin(mar);
  grVn2SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(2)->SetLeftMargin(mar);
  grVn4SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(3)->SetLeftMargin(mar);
  grVn6SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->cd(4)->SetLeftMargin(mar);
  grVn8SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatToDef->SaveAs("cSysPtRes_CumuCent.pdf");


  //-- cumuRatio
  TLegend * legvn6vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Default,  "Default",        "lp");
  legvn6vn4->AddEntry(grVn6Vn4SmearPt5, "5% p_{T} smear", "lp");

  TLegend * legvn8vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Default,  "Default",        "lp");
  legvn8vn4->AddEntry(grVn8Vn4SmearPt5, "5% p_{T} smear", "lp");

  TLegend * legvn8vn6 = new TLegend(0.4520, 0.1973, 0.94, 0.3703);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Default,  "Default",        "lp");
  legvn8vn6->AddEntry(grVn8Vn6SmearPt5, "5% p_{T} smear", "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);

  cCumuRatio->cd(1);
  grVn6Vn4Default->Draw("ap");
  grVn6Vn4SmearPt5->Draw("psame");
  legvn6vn4->Draw("same");

  cCumuRatio->cd(2);
  grVn8Vn4Default->Draw("ap");
  grVn8Vn4SmearPt5->Draw("psame");
  legvn8vn4->Draw("same");

  cCumuRatio->cd(3);
  grVn8Vn6Default->Draw("ap");
  grVn8Vn6SmearPt5->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- cumuRatio To Default
  TCanvas * cCumuRatio_RatToDef = new TCanvas("cCumuRatio_RatToDef", "cCumuRatio_RatToDef", 1500, 500);
  cCumuRatio_RatToDef->Divide(3,1);

  grVn6Vn4SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);
  grVn8Vn4SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);
  grVn8Vn6SmearPt5_RatioToDefault->GetYaxis()->SetTitleOffset(offs);

  cCumuRatio_RatToDef->cd(1)->SetLeftMargin(mar);
  grVn6Vn4SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->cd(2)->SetLeftMargin(mar);
  grVn8Vn4SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->cd(3)->SetLeftMargin(mar);
  grVn8Vn6SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");

  cCumuRatio_RatToDef->SaveAs("cSysPtRes_CumuRatioCent.pdf");

  //-- g1e
  TLegend * legG1E = new TLegend(0.4538, 0.7421, 0.9529, 0.9172);
  legG1E->SetBorderSize(0);
  legG1E->SetFillStyle(0);
  legG1E->AddEntry(grG1EDefault,  "Default",        "lp");
  legG1E->AddEntry(grG1ESmearPt5, "5% p_{T} smear", "lp");

  TCanvas* cG1e = new TCanvas("cG1e", "cG1e", 500, 500);
  cG1e->cd();
  grG1EDefault->Draw("ap");
  grG1ESmearPt5->Draw("psame");
  legG1E->Draw("same");
  cG1e->SaveAs("cG1e.pdf");

  //-- G1E ratio to default
  TCanvas* cG1eR = new TCanvas("cG1eR", "cG1eR", 500, 500);
  cG1eR->cd();
  grG1ESmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");
  cG1eR->SaveAs("cSysPtRes_G1ECent.pdf");

  //-- vn46_vn68 
  TCanvas* cvn46_vn68 = new TCanvas("cvn46_vn68", "cvn46_vn68", 500, 500);
  cvn46_vn68->cd();
  grVn46_Vn68SmearPt5_RatioToDefault->Draw("ap");
  lOne->Draw("same");
  cvn46_vn68->SaveAs("cSysPtRes_Vn46_Vn68.pdf");

  //----------------------------------------------------------------------------------------------------
  //-- Save plots for smoothing
  TFile * fSave = new TFile("SysPtRes.root", "recreate");
  fSave->cd();
  grVn2SmearPt5_RatioToDefault->Write("grVn2SmearPt5_RatioToDefault");
  grVn4SmearPt5_RatioToDefault->Write("grVn4SmearPt5_RatioToDefault");
  grVn6SmearPt5_RatioToDefault->Write("grVn6SmearPt5_RatioToDefault");
  grVn8SmearPt5_RatioToDefault->Write("grVn8SmearPt5_RatioToDefault");
  grVn6Vn4SmearPt5_RatioToDefault->Write("grVn6Vn4SmearPt5_RatioToDefault");
  grVn8Vn4SmearPt5_RatioToDefault->Write("grVn8Vn4SmearPt5_RatioToDefault");
  grVn8Vn6SmearPt5_RatioToDefault->Write("grVn8Vn6SmearPt5_RatioToDefault");
  grVn46_Vn68SmearPt5_RatioToDefault->Write("grVn46_Vn68SmearPt5_RatioToDefault");
  grG1ESmearPt5_RatioToDefault->Write("grG1ESmearPt5_RatioToDefault");

  //-- Save the final unfolded distns from the 5% p_{T} smear study
  for(int icent = 0; icent < NCENT; icent++){
    int i = iterCutSmearPt5[icent];
    std::cout<<i<<std::endl;
    hUnfoldSmearPt5[icent][i]->SetLineColor(1);
    hUnfoldSmearPt5[icent][i]->SetMarkerColor(1);
    hUnfoldSmearPt5[icent][i]->Write( Form("hFinalUnfold_SysPtRes_c%i", icent) );
  }

}

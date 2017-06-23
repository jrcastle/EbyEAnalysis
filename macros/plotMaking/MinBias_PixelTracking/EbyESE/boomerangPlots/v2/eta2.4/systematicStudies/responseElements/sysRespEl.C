#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace ebyese;

void sysRespEl(){

  int norder_ = 2;

  double ratioMin = 0.9;
  double ratioMax = 1.1;

  double ratioG1Min = 0.;
  double ratioG1Max = 1.;

  double ratioCumuRatioMin = 0.95;
  double ratioCumuRatioMax = 1.05;

  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Standard Unfolding
  TFile * fUnf;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- DoSys Unfolding
  TFile * fUnfDoSys;
  TH1D * hUnfoldDoSys[NCENT][NITER];
  TH1D * hRefoldDoSys[NCENT][NITER];

  //-- Iter Cutoffs
  bool iterStopFound[NCENT];
  int iterCutoff[NCENT];

  bool iterStopFoundDoSys[NCENT];
  int iterCutoffDoSys[NCENT];

  TLatex latex;

  //-- Systematics
  double vn2DoSys_RatioToNominal[NCENT];
  double vn4DoSys_RatioToNominal[NCENT];
  double vn6DoSys_RatioToNominal[NCENT];
  double vn8DoSys_RatioToNominal[NCENT];
  double gamma1expDoSys_RatioToNominal[NCENT];
  double vn6vn4DoSys_RatioToNominal[NCENT];
  double vn8vn4DoSys_RatioToNominal[NCENT];
  double vn8vn6DoSys_RatioToNominal[NCENT];

  double vn2DoSys_RatioToNominal_staterr[NCENT];
  double vn4DoSys_RatioToNominal_staterr[NCENT];
  double vn6DoSys_RatioToNominal_staterr[NCENT];
  double vn8DoSys_RatioToNominal_staterr[NCENT];
  double gamma1expDoSys_RatioToNominal_staterr[NCENT];
  double vn6vn4DoSys_RatioToNominal_staterr[NCENT];
  double vn8vn4DoSys_RatioToNominal_staterr[NCENT];
  double vn8vn6DoSys_RatioToNominal_staterr[NCENT];

  double vn2DoSys_PctDiffToNominal[NCENT];
  double vn4DoSys_PctDiffToNominal[NCENT];
  double vn6DoSys_PctDiffToNominal[NCENT];
  double vn8DoSys_PctDiffToNominal[NCENT];
  double gamma1expDoSys_PctDiffToNominal[NCENT];
  double vn6vn4DoSys_PctDiffToNominal[NCENT];
  double vn8vn4DoSys_PctDiffToNominal[NCENT];
  double vn8vn6DoSys_PctDiffToNominal[NCENT];

  //-- Systematic Performance Plots
  TGraphErrors * grVn2DoSys_RatioToNominal;
  TGraphErrors * grVn4DoSys_RatioToNominal;
  TGraphErrors * grVn6DoSys_RatioToNominal;
  TGraphErrors * grVn8DoSys_RatioToNominal;
  TGraphErrors * grGamma1ExpDoSys_RatioToNominal;
  TGraphErrors * grVn6Vn4DoSys_RatioToNominal;
  TGraphErrors * grVn8Vn4DoSys_RatioToNominal;
  TGraphErrors * grVn8Vn6DoSys_RatioToNominal;

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

  TFile * fStat_DoSys;
  TH1D * hVarianceOfMean_Vn2_DoSys;
  TH1D * hVarianceOfMean_Vn4_DoSys;
  TH1D * hVarianceOfMean_Vn6_DoSys;
  TH1D * hVarianceOfMean_Vn8_DoSys;
  TH1D * hVarianceOfMean_Gamma1Exp_DoSys;
  TH1D * hVarianceOfMean_Vn6Vn4_DoSys;
  TH1D * hVarianceOfMean_Vn8Vn4_DoSys;
  TH1D * hVarianceOfMean_Vn8Vn6_DoSys;

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
  fUnfDoSys = new TFile( Form("../../UnfoldResults/dataResp/data%iGauss.root", norder_) );

  fRelErr   = new TFile("relErrorResponse.root", "recreate");
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

  fStat_DoSys = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/responseElements/StatUncertRespEl_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn2_DoSys       = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn2_RespEl" );
  hVarianceOfMean_Vn4_DoSys       = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn4_RespEl" );
  hVarianceOfMean_Vn6_DoSys       = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn6_RespEl" );
  hVarianceOfMean_Vn8_DoSys       = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn8_RespEl" );
  hVarianceOfMean_Gamma1Exp_DoSys = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Gamma1Exp_RespEl" );
  hVarianceOfMean_Vn6Vn4_DoSys    = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn6Vn4_RespEl" );
  hVarianceOfMean_Vn8Vn4_DoSys    = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn8Vn4_RespEl" );
  hVarianceOfMean_Vn8Vn6_DoSys    = (TH1D*) fStat_DoSys->Get( "hVarianceOfMean_Vn8Vn6_RespEl" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    //-- Initialize iteration cutoff values/booleans
    iterStopFound[icent] = 0;
    iterCutoff[icent]    = 0;

    iterStopFoundDoSys[icent] = 0;
    iterCutoffDoSys[icent]    = 0;

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

      //-- Get the DoSys unfolded histograms
      hUnfoldDoSys[icent][i] = (TH1D*) fUnfDoSys->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldDoSys[icent][i]->SetLineColor(col[i]);
      hUnfoldDoSys[icent][i]->SetMarkerColor(col[i]);

      //-- Get the DoSys refolded histograms
      hRefoldDoSys[icent][i] = (TH1D*) fUnfDoSys->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldDoSys[icent][i]->SetLineColor(col[i]);
      hRefoldDoSys[icent][i]->SetMarkerColor(col[i]);

      //-- Chi squares
      double chi2NDF_Refold       = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      double chi2NDF_Refold_DoSys = hRefoldDoSys[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      std::cout << "chi2NDF_Refold_DoSys  = " << chi2NDF_Refold_DoSys << std::endl;

      //-- Normal unfolding chi2 check
      if( chi2NDF_Refold < 1.2 && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }

      //-- DoSys unfolding chi2 check
      if( chi2NDF_Refold_DoSys < 1.2 && !iterStopFoundDoSys[icent] ){
        iterStopFoundDoSys[icent] = 1;
        iterCutoffDoSys[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFoundDoSys[icent] ){
	iterStopFoundDoSys[icent] = 1;
        iterCutoffDoSys[icent]    = i;
      }

    } //-- End iter loop

  } //-- End cent loop

  //-- Figure out the response matrix uncertainty component
  for(int icent = 0; icent < NCENT; icent++){

    std::cout << Form("=================== Cent %i ===================", icent) << std::endl;

    int i      = iterCutoff[icent];
    int iDoSys = iterCutoffDoSys[icent];

    FixUnfold( hUnfold[icent][i] );
    FixUnfold( hUnfoldDoSys[icent][iDoSys] );

    EbyECumu cumu( hUnfold[icent][i] );
    EbyECumu cumuDoSys( hUnfoldDoSys[icent][iDoSys] );

    double vn2        = cumu.GetCumu_vn2();
    double vn4        = cumu.GetCumu_vn4();
    double vn6        = cumu.GetCumu_vn6();
    double vn8        = cumu.GetCumu_vn8();
    double vn2DoSys   = cumuDoSys.GetCumu_vn2();
    double vn4DoSys   = cumuDoSys.GetCumu_vn4();
    double vn6DoSys   = cumuDoSys.GetCumu_vn6();
    double vn8DoSys   = cumuDoSys.GetCumu_vn8();

    double vn2_staterr        = sqrt( hVarianceOfMean_Vn2_Nominal->GetBinContent(icent+1) );
    double vn4_staterr        = sqrt( hVarianceOfMean_Vn4_Nominal->GetBinContent(icent+1) );
    double vn6_staterr        = sqrt( hVarianceOfMean_Vn6_Nominal->GetBinContent(icent+1) );
    double vn8_staterr        = sqrt( hVarianceOfMean_Vn8_Nominal->GetBinContent(icent+1) );
    double vn2DoSys_staterr   = sqrt( hVarianceOfMean_Vn2_DoSys->GetBinContent(icent+1) );
    double vn4DoSys_staterr   = sqrt( hVarianceOfMean_Vn4_DoSys->GetBinContent(icent+1) );
    double vn6DoSys_staterr   = sqrt( hVarianceOfMean_Vn6_DoSys->GetBinContent(icent+1) );
    double vn8DoSys_staterr   = sqrt( hVarianceOfMean_Vn8_DoSys->GetBinContent(icent+1) );

    double gamma1exp      = cumu.GetGamma1Exp();
    double gamma1expDoSys = cumuDoSys.GetGamma1Exp();

    double gamma1exp_staterr      = sqrt( hVarianceOfMean_Gamma1Exp_Nominal->GetBinContent(icent+1) );
    double gamma1expDoSys_staterr = sqrt( hVarianceOfMean_Gamma1Exp_DoSys->GetBinContent(icent+1) );

    double vn6vn4;
    double vn6vn4DoSys;
    if( vn6 == 0 || vn4 == 0 ) vn6vn4 = 0;
    else                       vn6vn4 = vn6 / vn4;
    if( vn6DoSys == 0 || vn4DoSys == 0) vn6vn4DoSys = 0;
    else                                vn6vn4DoSys = vn6DoSys / vn4DoSys;

    double vn8vn4;
    double vn8vn4DoSys;
    if( vn8 == 0 || vn4 == 0 ) vn8vn4 = 0;
    else                       vn8vn4 = vn8 / vn4;
    if( vn8DoSys == 0 || vn4DoSys == 0) vn8vn4DoSys = 0;
    else                                vn8vn4DoSys = vn8DoSys / vn4DoSys;

    double vn8vn6;
    double vn8vn6DoSys;
    if( vn8 == 0 || vn6 == 0 ) vn8vn6 = 0;
    else                       vn8vn6 = vn8 / vn6;
    if( vn8DoSys == 0 || vn6DoSys == 0) vn8vn6DoSys = 0;
    else                                vn8vn6DoSys = vn8DoSys / vn6DoSys;

    double vn6vn4_staterr      = sqrt( hVarianceOfMean_Vn6Vn4_Nominal->GetBinContent(icent+1) );
    double vn6vn4DoSys_staterr = sqrt( hVarianceOfMean_Vn6Vn4_DoSys->GetBinContent(icent+1) );
    double vn8vn4_staterr      = sqrt( hVarianceOfMean_Vn8Vn4_Nominal->GetBinContent(icent+1) );
    double vn8vn4DoSys_staterr = sqrt( hVarianceOfMean_Vn8Vn4_DoSys->GetBinContent(icent+1) );
    double vn8vn6_staterr      = sqrt( hVarianceOfMean_Vn8Vn6_Nominal->GetBinContent(icent+1) );
    double vn8vn6DoSys_staterr = sqrt( hVarianceOfMean_Vn8Vn6_DoSys->GetBinContent(icent+1) );

    //-- Calculate ratios
    if( vn2 == 0 || vn2DoSys == 0 ) vn2DoSys_RatioToNominal[icent] = 0;
    else                            vn2DoSys_RatioToNominal[icent] = vn2DoSys / vn2;

    if( vn4 == 0 || vn4DoSys == 0 ) vn4DoSys_RatioToNominal[icent] = 0;
    else                            vn4DoSys_RatioToNominal[icent] = vn4DoSys / vn4;

    if( vn6 == 0 || vn6DoSys == 0 ) vn6DoSys_RatioToNominal[icent] = 0;
    else                            vn6DoSys_RatioToNominal[icent] = vn6DoSys / vn6;

    if( vn8 == 0 || vn8DoSys == 0 ) vn8DoSys_RatioToNominal[icent] = 0;
    else                            vn8DoSys_RatioToNominal[icent] = vn8DoSys / vn8;

    if( gamma1exp == 0 || gamma1expDoSys == 0 ) gamma1expDoSys_RatioToNominal[icent] = 0;
    else                                        gamma1expDoSys_RatioToNominal[icent] = gamma1expDoSys / gamma1exp;

    if( vn6vn4 == 0 || vn6vn4DoSys == 0 ) vn6vn4DoSys_RatioToNominal[icent] = 0;
    else                                  vn6vn4DoSys_RatioToNominal[icent] = vn6vn4DoSys / vn6vn4;

    if( vn8vn4 == 0 || vn8vn4DoSys == 0 ) vn8vn4DoSys_RatioToNominal[icent] = 0;
    else                                  vn8vn4DoSys_RatioToNominal[icent] = vn8vn4DoSys / vn8vn4;

    if( vn8vn6 == 0 || vn8vn6DoSys == 0 ) vn8vn6DoSys_RatioToNominal[icent] = 0;
    else                                  vn8vn6DoSys_RatioToNominal[icent] = vn8vn6DoSys / vn8vn6;

    //-- Ratio errors
    vn2DoSys_RatioToNominal_staterr[icent]       = sqrt( pow( vn2DoSys_staterr/vn2,2) + pow(vn2DoSys*vn2_staterr/vn2/vn2,2) );
    vn4DoSys_RatioToNominal_staterr[icent]       = sqrt( pow( vn4DoSys_staterr/vn4,2) + pow(vn4DoSys*vn4_staterr/vn4/vn4,2) );
    vn6DoSys_RatioToNominal_staterr[icent]       = sqrt( pow( vn6DoSys_staterr/vn6,2) + pow(vn6DoSys*vn6_staterr/vn6/vn6,2) );
    vn8DoSys_RatioToNominal_staterr[icent]       = sqrt( pow( vn8DoSys_staterr/vn8,2) + pow(vn8DoSys*vn8_staterr/vn8/vn8,2) );
    gamma1expDoSys_RatioToNominal_staterr[icent] = sqrt( pow( gamma1expDoSys_staterr/gamma1exp,2) + pow(gamma1expDoSys*gamma1exp_staterr/gamma1exp/gamma1exp,2) );
    vn6vn4DoSys_RatioToNominal_staterr[icent]    = sqrt( pow( vn6vn4DoSys_staterr/vn6vn4,2) + pow(vn6vn4DoSys*vn6vn4_staterr/vn6vn4/vn6vn4,2) );
    vn8vn4DoSys_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn4DoSys_staterr/vn8vn4,2) + pow(vn8vn4DoSys*vn8vn4_staterr/vn8vn4/vn8vn4,2) );
    vn8vn6DoSys_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn6DoSys_staterr/vn8vn6,2) + pow(vn8vn6DoSys*vn8vn6_staterr/vn8vn6/vn8vn6,2) );

    //-- Calculate pct difference relative to nominal
    if( vn2 == 0 || vn2DoSys == 0 || compatibleWithOne(vn2DoSys_RatioToNominal[icent], vn2DoSys_RatioToNominal_staterr[icent]) ) vn2DoSys_PctDiffToNominal[icent] = 0;
    else                            vn2DoSys_PctDiffToNominal[icent] = fabs( vn2DoSys - vn2 ) / fabs( vn2 );

    if( vn4 == 0 || vn4DoSys == 0 || compatibleWithOne(vn4DoSys_RatioToNominal[icent], vn4DoSys_RatioToNominal_staterr[icent])) vn4DoSys_PctDiffToNominal[icent] = 0;
    else                            vn4DoSys_PctDiffToNominal[icent] = fabs( vn4DoSys - vn4 ) / fabs( vn4 );

    if( vn6 == 0 || vn6DoSys == 0 || compatibleWithOne(vn6DoSys_RatioToNominal[icent], vn6DoSys_RatioToNominal_staterr[icent])) vn6DoSys_PctDiffToNominal[icent] = 0;
    else                            vn6DoSys_PctDiffToNominal[icent] = fabs( vn6DoSys - vn6 ) / fabs( vn6 );

    if( vn8 == 0 || vn8DoSys == 0 || compatibleWithOne(vn8DoSys_RatioToNominal[icent], vn8DoSys_RatioToNominal_staterr[icent])) vn8DoSys_PctDiffToNominal[icent] = 0;
    else                            vn8DoSys_PctDiffToNominal[icent] = fabs( vn8DoSys - vn8 ) / fabs( vn8 );

    if( gamma1exp == 0 || gamma1expDoSys == 0 || compatibleWithOne(gamma1expDoSys_RatioToNominal[icent], gamma1expDoSys_RatioToNominal_staterr[icent]) ) gamma1expDoSys_PctDiffToNominal[icent] = 0;
    else                                        gamma1expDoSys_PctDiffToNominal[icent] = fabs( 1. - fabs( gamma1expDoSys - gamma1exp ) / fabs( gamma1exp ) );

    if( vn6vn4 == 0 || vn6vn4DoSys == 0 || compatibleWithOne(vn6vn4DoSys_RatioToNominal[icent], vn6vn4DoSys_RatioToNominal_staterr[icent]) ) vn6vn4DoSys_PctDiffToNominal[icent] = 0;
    else                                  vn6vn4DoSys_PctDiffToNominal[icent] = fabs( vn6vn4DoSys - vn6vn4 ) / fabs( vn6vn4 );

    if( vn8vn4 == 0 || vn8vn4DoSys == 0 || compatibleWithOne(vn8vn4DoSys_RatioToNominal[icent], vn8vn4DoSys_RatioToNominal_staterr[icent]) ) vn8vn4DoSys_PctDiffToNominal[icent] = 0;
    else                                  vn8vn4DoSys_PctDiffToNominal[icent] = fabs( vn8vn4DoSys - vn8vn4 ) / fabs( vn8vn4 );

    if( vn8vn6 == 0 || vn8vn6DoSys == 0 || compatibleWithOne(vn8vn6DoSys_RatioToNominal[icent], vn8vn6DoSys_RatioToNominal_staterr[icent]) ) vn8vn6DoSys_PctDiffToNominal[icent] = 0;
    else                                  vn8vn6DoSys_PctDiffToNominal[icent] = fabs( vn8vn6DoSys - vn8vn6 ) / fabs( vn8vn6 );

  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT];
  for(int i = 0; i < NCENT; i++) cErr[i] = 0;
  grVn2DoSys_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn2DoSys_RatioToNominal, cErr, vn2DoSys_RatioToNominal_staterr);
  formatGraph(grVn2DoSys_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} [RespErr] / [Nominal]", norder_), 1, 24, "grVn2DoSys_RatioToNominal");
  grVn4DoSys_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn4DoSys_RatioToNominal, cErr, vn4DoSys_RatioToNominal_staterr);
  formatGraph(grVn4DoSys_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} [RespErr] / [Nominal]", norder_), kSpring+4, 25, "grVn4DoSys_RatioToNominal");
  grVn6DoSys_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn6DoSys_RatioToNominal, cErr, vn6DoSys_RatioToNominal_staterr);
  formatGraph(grVn6DoSys_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} [RespErr] / [Nominal]", norder_), 6, 28, "grVn6DoSys_RatioToNominal");
  grVn8DoSys_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn8DoSys_RatioToNominal, cErr, vn8DoSys_RatioToNominal_staterr);
  formatGraph(grVn8DoSys_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} [RespErr] / [Nominal]", norder_), kOrange+7, 27, "grVn8DoSys_RatioToNominal");
  grGamma1ExpDoSys_RatioToNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expDoSys_RatioToNominal, cErr, gamma1expDoSys_RatioToNominal_staterr);
  formatGraph(grGamma1ExpDoSys_RatioToNominal, "Centrality %", ratioG1Min, ratioG1Max, "|1-#gamma_{1}^{exp} [RespErr] / [Nominal]|", 2, 20, "grGamma1ExpDoSys_RatioToNominal");
  grVn6Vn4DoSys_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4DoSys_RatioToNominal, cErr, vn6vn4DoSys_RatioToNominal_staterr);
  formatGraph(grVn6Vn4DoSys_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{6}/v_{%i}{4} [RespErr] / [Nominal]", norder_, norder_), 4, 21, "grVn6Vn4DoSys_RatioToNominal");
  grVn8Vn4DoSys_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4DoSys_RatioToNominal, cErr, vn8vn4DoSys_RatioToNominal_staterr);
  formatGraph(grVn8Vn4DoSys_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{4} [RespErr] / [Nominal]", norder_, norder_), kGreen+2, 34, "grVn8Vn4DoSys_RatioToNominal");
  grVn8Vn6DoSys_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6DoSys_RatioToNominal, cErr, vn8vn6DoSys_RatioToNominal_staterr);
  formatGraph(grVn8Vn6DoSys_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{6} [RespErr] / [Nominal]", norder_, norder_), kViolet-1, 33, "grVn8Vn6DoSys_RatioToNominal");


  TLine * line = new TLine(grVn2DoSys_RatioToNominal->GetXaxis()->GetXmin(), 1.0, grVn2DoSys_RatioToNominal->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);


  TCanvas * cCumuSys = new TCanvas("cCumuSys", "cCumuSys", 1000, 1000);
  cCumuSys->Divide(2,2);
  cCumuSys->cd(1);
  grVn2DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(2);
  grVn4DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(3);
  grVn6DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(4);
  grVn8DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->SaveAs("../../plots/systematicStudies/cSysResp_CumuCent.pdf");

  TCanvas * cGamma1ExpSys = new TCanvas("cGamma1ExpSys", "cGamma1ExpSys", 500, 500);
  cGamma1ExpSys->cd();
  grGamma1ExpDoSys_RatioToNominal->Draw("ap");
  //line->Draw("same");
  cGamma1ExpSys->SaveAs("../../plots/systematicStudies/cSysResp_Gamma1ExpCent.pdf");

  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  grVn6Vn4DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6DoSys_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/cSysResp_CumuRatioCent.pdf");

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){

    std::cout<<"---------------"<<std::endl;

    relErrVn2->SetBinContent(icent+1, vn2DoSys_PctDiffToNominal[icent]);
    relErrVn4->SetBinContent(icent+1, vn4DoSys_PctDiffToNominal[icent]);
    relErrVn6->SetBinContent(icent+1, vn6DoSys_PctDiffToNominal[icent]);
    relErrVn8->SetBinContent(icent+1, vn6DoSys_PctDiffToNominal[icent]);
    relErrGamma1Exp->SetBinContent(icent+1, gamma1expDoSys_PctDiffToNominal[icent]);
    relErrVn6Vn4->SetBinContent(icent+1, vn6vn4DoSys_PctDiffToNominal[icent]);
    relErrVn8Vn4->SetBinContent(icent+1, vn8vn4DoSys_PctDiffToNominal[icent]);
    relErrVn8Vn6->SetBinContent(icent+1, vn8vn6DoSys_PctDiffToNominal[icent]);

  }

  fRelErr->Write();

}

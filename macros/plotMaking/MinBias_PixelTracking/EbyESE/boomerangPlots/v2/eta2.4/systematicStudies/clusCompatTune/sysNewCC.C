#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace ebyese;

void sysNewCC(){

  int norder_ = 2;

  double ratioMin = 0.98;
  double ratioMax = 1.02;

  double ratioG1Min = 0.0;
  double ratioG1Max = 1.0;

  double ratioCumuRatioMin = 0.98;
  double ratioCumuRatioMax = 1.02;

  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Standard Unfolding
  TFile * fUnf;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- NewCC Unfolding
  TFile * fUnfNewCC;
  TH1D * hUnfoldNewCC[NCENT][NITER];
  TH1D * hRefoldNewCC[NCENT][NITER];

  //-- Iter Cutoffs
  bool iterStopFound[NCENT];
  int iterCutoff[NCENT];

  bool iterStopFoundNewCC[NCENT];
  int iterCutoffNewCC[NCENT];

  TLatex latex;

  //-- Systematics
  double vn2NewCC_RatioToNominal[NCENT];
  double vn4NewCC_RatioToNominal[NCENT];
  double vn6NewCC_RatioToNominal[NCENT];
  double vn8NewCC_RatioToNominal[NCENT];
  double gamma1expNewCC_RatioToNominal[NCENT];
  double vn6vn4NewCC_RatioToNominal[NCENT];
  double vn8vn4NewCC_RatioToNominal[NCENT];
  double vn8vn6NewCC_RatioToNominal[NCENT];

  double vn2NewCC_RatioToNominal_staterr[NCENT];
  double vn4NewCC_RatioToNominal_staterr[NCENT];
  double vn6NewCC_RatioToNominal_staterr[NCENT];
  double vn8NewCC_RatioToNominal_staterr[NCENT];
  double gamma1expNewCC_RatioToNominal_staterr[NCENT];
  double vn6vn4NewCC_RatioToNominal_staterr[NCENT];
  double vn8vn4NewCC_RatioToNominal_staterr[NCENT];
  double vn8vn6NewCC_RatioToNominal_staterr[NCENT];

  double vn2NewCC_PctDiffToNominal[NCENT];
  double vn4NewCC_PctDiffToNominal[NCENT];
  double vn6NewCC_PctDiffToNominal[NCENT];
  double vn8NewCC_PctDiffToNominal[NCENT];
  double gamma1expNewCC_PctDiffToNominal[NCENT];
  double vn6vn4NewCC_PctDiffToNominal[NCENT];
  double vn8vn4NewCC_PctDiffToNominal[NCENT];
  double vn8vn6NewCC_PctDiffToNominal[NCENT];

  //-- Systematic Performance Plots
  TGraphErrors * grVn2NewCC_RatioToNominal;
  TGraphErrors * grVn4NewCC_RatioToNominal;
  TGraphErrors * grVn6NewCC_RatioToNominal;
  TGraphErrors * grVn8NewCC_RatioToNominal;
  TGraphErrors * grGamma1ExpNewCC_RatioToNominal;
  TGraphErrors * grVn6Vn4NewCC_RatioToNominal;
  TGraphErrors * grVn8Vn4NewCC_RatioToNominal;
  TGraphErrors * grVn8Vn6NewCC_RatioToNominal;

  //-- Compare raw values
  double vn2NEWCC[NCENT];
  double vn4NEWCC[NCENT];
  double vn6NEWCC[NCENT];
  double vn8NEWCC[NCENT];
  double g1eNEWCC[NCENT];
  double vn6vn4NEWCC[NCENT];
  double vn8vn4NEWCC[NCENT];
  double vn8vn6NEWCC[NCENT];

  double vn2NEWCCe[NCENT];
  double vn4NEWCCe[NCENT];
  double vn6NEWCCe[NCENT];
  double vn8NEWCCe[NCENT];
  double g1eNEWCCe[NCENT];
  double vn6vn4NEWCCe[NCENT];
  double vn8vn4NEWCCe[NCENT];
  double vn8vn6NEWCCe[NCENT];

  TGraphErrors * grVn2NewCC;
  TGraphErrors * grVn4NewCC;
  TGraphErrors * grVn6NewCC;
  TGraphErrors * grVn8NewCC;
  TGraphErrors * grG1eNewCC;
  TGraphErrors * grVn6Vn4NewCC;
  TGraphErrors * grVn8Vn4NewCC;
  TGraphErrors * grVn8Vn6NewCC;

  double vn2Nominal[NCENT];
  double vn4Nominal[NCENT];
  double vn6Nominal[NCENT];
  double vn8Nominal[NCENT];
  double g1eNominal[NCENT];
  double vn6vn4Nominal[NCENT];
  double vn8vn4Nominal[NCENT];
  double vn8vn6Nominal[NCENT];

  double vn2Nominale[NCENT];
  double vn4Nominale[NCENT];
  double vn6Nominale[NCENT];
  double vn8Nominale[NCENT];
  double g1eNominale[NCENT];
  double vn6vn4Nominale[NCENT];
  double vn8vn4Nominale[NCENT];
  double vn8vn6Nominale[NCENT];

  TGraphErrors * grVn2Nominal;
  TGraphErrors * grVn4Nominal;
  TGraphErrors * grVn6Nominal;
  TGraphErrors * grVn8Nominal;
  TGraphErrors * grG1eNominal;
  TGraphErrors * grVn6Vn4Nominal;
  TGraphErrors * grVn8Vn4Nominal;
  TGraphErrors * grVn8Vn6Nominal;


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

  TFile * fStat_NewCC;
  TH1D * hVarianceOfMean_Vn2_NewCC;
  TH1D * hVarianceOfMean_Vn4_NewCC;
  TH1D * hVarianceOfMean_Vn6_NewCC;
  TH1D * hVarianceOfMean_Vn8_NewCC;
  TH1D * hVarianceOfMean_Gamma1Exp_NewCC;
  TH1D * hVarianceOfMean_Vn6Vn4_NewCC;
  TH1D * hVarianceOfMean_Vn8Vn4_NewCC;
  TH1D * hVarianceOfMean_Vn8Vn6_NewCC;

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
  fUnfNewCC = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorNewCC.root", "recreate");
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

  fStat_NewCC = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/clusCompatTune/StatUncertNewCC_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn2_NewCC       = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn2_NewCC" );
  hVarianceOfMean_Vn4_NewCC       = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn4_NewCC" );
  hVarianceOfMean_Vn6_NewCC       = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn6_NewCC" );
  hVarianceOfMean_Vn8_NewCC       = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn8_NewCC" );
  hVarianceOfMean_Gamma1Exp_NewCC = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Gamma1Exp_NewCC" );
  hVarianceOfMean_Vn6Vn4_NewCC    = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn6Vn4_NewCC" );
  hVarianceOfMean_Vn8Vn4_NewCC    = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn8Vn4_NewCC" );
  hVarianceOfMean_Vn8Vn6_NewCC    = (TH1D*) fStat_NewCC->Get( "hVarianceOfMean_Vn8Vn6_NewCC" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    //-- Initialize iteration cutoff values/booleans
    iterStopFound[icent] = 0;
    iterCutoff[icent]    = 0;

    iterStopFoundNewCC[icent] = 0;
    iterCutoffNewCC[icent]    = 0;

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

      //-- Get the NewCC unfolded histograms
      hUnfoldNewCC[icent][i] = (TH1D*) fUnfNewCC->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldNewCC[icent][i]->SetLineColor(col[i]);
      hUnfoldNewCC[icent][i]->SetMarkerColor(col[i]);

      //-- Get the NewCC refolded histograms
      hRefoldNewCC[icent][i] = (TH1D*) fUnfNewCC->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefoldNewCC[icent][i]->SetLineColor(col[i]);
      hRefoldNewCC[icent][i]->SetMarkerColor(col[i]);


      //-- Chi squares
      double chi2NDF_Refold       = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      double chi2NDF_Refold_NewCC = hRefoldNewCC[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");


      //-- Normal unfolding chi2 check
      if( chi2NDF_Refold < 1.2 && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFound[icent] ){
	iterStopFound[icent] = 1;
	iterCutoff[icent]    = i;
      }

      //-- NewCC unfolding chi2 check
      if( chi2NDF_Refold_NewCC < 1.2 && !iterStopFoundNewCC[icent] ){
        iterStopFoundNewCC[icent] = 1;
        iterCutoffNewCC[icent]    = i;
      }
      if( i == (NITER - 1) && !iterStopFoundNewCC[icent] ){
	iterStopFoundNewCC[icent] = 1;
        iterCutoffNewCC[icent]    = i;
      }

    } //-- End iter loop

  } //-- End cent loop

  //-- Figure out the response matrix uncertainty component
  for(int icent = 0; icent < NCENT; icent++){

    std::cout << Form("=================== Cent %i ===================", icent) << std::endl;

    int i      = iterCutoff[icent];
    int iNewCC = iterCutoffNewCC[icent];

    FixUnfold( hUnfold[icent][i] );
    FixUnfold( hUnfoldNewCC[icent][iNewCC] );

    EbyECumu cumu( hUnfold[icent][i] );
    EbyECumu cumuNewCC( hUnfoldNewCC[icent][iNewCC] );

    double vn2        = cumu.GetCumu_vn2();
    double vn4        = cumu.GetCumu_vn4();
    double vn6        = cumu.GetCumu_vn6();
    double vn8        = cumu.GetCumu_vn8();
    double vn2NewCC   = cumuNewCC.GetCumu_vn2();
    double vn4NewCC   = cumuNewCC.GetCumu_vn4();
    double vn6NewCC   = cumuNewCC.GetCumu_vn6();
    double vn8NewCC   = cumuNewCC.GetCumu_vn8();

    std::cout << "vn2 = " << vn2 << "\tvn2NewCC = " << vn2NewCC << "\tratio = " << vn2NewCC / vn2 << std::endl;

    double vn2_staterr         = sqrt( hVarianceOfMean_Vn2_Nominal->GetBinContent(icent+1) );
    double vn4_staterr         = sqrt( hVarianceOfMean_Vn4_Nominal->GetBinContent(icent+1) );
    double vn6_staterr         = sqrt( hVarianceOfMean_Vn6_Nominal->GetBinContent(icent+1) );
    double vn8_staterr         = sqrt( hVarianceOfMean_Vn8_Nominal->GetBinContent(icent+1) );
    double vn2NewCC_staterr   = sqrt( hVarianceOfMean_Vn2_NewCC->GetBinContent(icent+1) );
    double vn4NewCC_staterr   = sqrt( hVarianceOfMean_Vn4_NewCC->GetBinContent(icent+1) );
    double vn6NewCC_staterr   = sqrt( hVarianceOfMean_Vn6_NewCC->GetBinContent(icent+1) );
    double vn8NewCC_staterr   = sqrt( hVarianceOfMean_Vn8_NewCC->GetBinContent(icent+1) );

    double gamma1exp      = cumu.GetGamma1Exp();
    double gamma1expNewCC = cumuNewCC.GetGamma1Exp();

    double gamma1exp_staterr       = sqrt( hVarianceOfMean_Gamma1Exp_Nominal->GetBinContent(icent+1) );
    double gamma1expNewCC_staterr = sqrt( hVarianceOfMean_Gamma1Exp_NewCC->GetBinContent(icent+1) );

    double vn6vn4;
    double vn6vn4NewCC;
    if( vn6 == 0 || vn4 == 0 ) vn6vn4 = 0;
    else                       vn6vn4 = vn6 / vn4;
    if( vn6NewCC == 0 || vn4NewCC == 0) vn6vn4NewCC = 0;
    else                                  vn6vn4NewCC = vn6NewCC / vn4NewCC;

    double vn8vn4;
    double vn8vn4NewCC;
    if( vn8 == 0 || vn4 == 0 ) vn8vn4 = 0;
    else                       vn8vn4 = vn8 / vn4;
    if( vn8NewCC == 0 || vn4NewCC == 0) vn8vn4NewCC = 0;
    else                                  vn8vn4NewCC = vn8NewCC / vn4NewCC;

    double vn8vn6;
    double vn8vn6NewCC;
    if( vn8 == 0 || vn6 == 0 ) vn8vn6 = 0;
    else                       vn8vn6 = vn8 / vn6;
    if( vn8NewCC == 0 || vn6NewCC == 0) vn8vn6NewCC = 0;
    else                                  vn8vn6NewCC = vn8NewCC / vn6NewCC;

    double vn6vn4_staterr       = sqrt( hVarianceOfMean_Vn6Vn4_Nominal->GetBinContent(icent+1) );
    double vn6vn4NewCC_staterr = sqrt( hVarianceOfMean_Vn6Vn4_NewCC->GetBinContent(icent+1) );
    double vn8vn4_staterr       = sqrt( hVarianceOfMean_Vn8Vn4_Nominal->GetBinContent(icent+1) );
    double vn8vn4NewCC_staterr = sqrt( hVarianceOfMean_Vn8Vn4_NewCC->GetBinContent(icent+1) );
    double vn8vn6_staterr       = sqrt( hVarianceOfMean_Vn8Vn6_Nominal->GetBinContent(icent+1) );
    double vn8vn6NewCC_staterr = sqrt( hVarianceOfMean_Vn8Vn6_NewCC->GetBinContent(icent+1) );

    //-- Calculate ratios
    if( vn2 == 0 || vn2NewCC == 0 ) vn2NewCC_RatioToNominal[icent] = 0;
    else                             vn2NewCC_RatioToNominal[icent] = vn2NewCC / vn2;

    if( vn4 == 0 || vn4NewCC == 0 ) vn4NewCC_RatioToNominal[icent] = 0;
    else                             vn4NewCC_RatioToNominal[icent] = vn4NewCC / vn4;

    if( vn6 == 0 || vn6NewCC == 0 ) vn6NewCC_RatioToNominal[icent] = 0;
    else                             vn6NewCC_RatioToNominal[icent] = vn6NewCC / vn6;

    if( vn8 == 0 || vn8NewCC == 0 ) vn8NewCC_RatioToNominal[icent] = 0;
    else                             vn8NewCC_RatioToNominal[icent] = vn8NewCC / vn8;

    if( gamma1exp == 0 || gamma1expNewCC == 0 ) gamma1expNewCC_RatioToNominal[icent] = 0;
    else                                         gamma1expNewCC_RatioToNominal[icent] = fabs(1. - gamma1expNewCC / gamma1exp );

    if( vn6vn4 == 0 || vn6vn4NewCC == 0 ) vn6vn4NewCC_RatioToNominal[icent] = 0;
    else                                   vn6vn4NewCC_RatioToNominal[icent] = vn6vn4NewCC / vn6vn4;

    if( vn8vn4 == 0 || vn8vn4NewCC == 0 ) vn8vn4NewCC_RatioToNominal[icent] = 0;
    else                                   vn8vn4NewCC_RatioToNominal[icent] = vn8vn4NewCC / vn8vn4;

    if( vn8vn6 == 0 || vn8vn6NewCC == 0 ) vn8vn6NewCC_RatioToNominal[icent] = 0;
    else                                   vn8vn6NewCC_RatioToNominal[icent] = vn8vn6NewCC / vn8vn6;

    //-- Raw Values
    if(vn2 == 0) vn2 = -1.;
    if(vn4 == 0) vn4 = -1.;
    if(vn6 == 0) vn6 = -1.;
    if(vn8 == 0) vn8 = -1.;
    if(vn2NewCC== 0) vn2NewCC = -1;
    if(vn4NewCC == 0) vn4NewCC = -1;
    if(vn6NewCC == 0) vn6NewCC = -1;
    if(vn8NewCC == 0) vn8NewCC = -1;

    vn2Nominal[icent] = vn2;
    vn4Nominal[icent] = vn4;
    vn6Nominal[icent] = vn6;
    vn8Nominal[icent] = vn8;
    g1eNominal[icent] = gamma1exp;
    vn6vn4Nominal[icent] = vn6vn4;
    vn8vn4Nominal[icent] = vn8vn4;
    vn8vn6Nominal[icent] = vn8vn6;

    vn2Nominale[icent] = vn2_staterr;
    vn4Nominale[icent] = vn4_staterr;
    vn6Nominale[icent] = vn6_staterr;
    vn8Nominale[icent] = vn8_staterr;
    g1eNominale[icent] = gamma1exp_staterr;
    vn6vn4Nominale[icent] = vn6vn4_staterr;
    vn8vn4Nominale[icent] = vn8vn4_staterr;
    vn8vn6Nominale[icent] = vn8vn6_staterr;


    vn2NEWCC[icent] = vn2NewCC;
    vn4NEWCC[icent] = vn4NewCC;
    vn6NEWCC[icent] = vn6NewCC;
    vn8NEWCC[icent] = vn8NewCC;
    g1eNEWCC[icent] = gamma1expNewCC;
    vn6vn4NEWCC[icent] = vn6vn4NewCC;
    vn8vn4NEWCC[icent] = vn8vn4NewCC;
    vn8vn6NEWCC[icent] = vn8vn6NewCC;
    std::cout<< vn8vn6NewCC << "\t" << vn8vn6 << std::endl;
    vn2NEWCCe[icent] = vn2NewCC_staterr;
    vn4NEWCCe[icent] = vn4NewCC_staterr;
    vn6NEWCCe[icent] = vn6NewCC_staterr;
    vn8NEWCCe[icent] = vn8NewCC_staterr;
    g1eNEWCCe[icent] = gamma1expNewCC_staterr;
    vn6vn4NEWCCe[icent] = vn6vn4NewCC_staterr;
    vn8vn4NEWCCe[icent] = vn8vn4NewCC_staterr;
    vn8vn6NEWCCe[icent] = vn8vn6NewCC_staterr;

    //-- Ratio errors
    vn2NewCC_RatioToNominal_staterr[icent]       = sqrt( pow( vn2NewCC_staterr/vn2,2) + pow(vn2NewCC*vn2_staterr/vn2/vn2,2) );
    vn4NewCC_RatioToNominal_staterr[icent]       = sqrt( pow( vn4NewCC_staterr/vn4,2) + pow(vn4NewCC*vn4_staterr/vn4/vn4,2) );
    vn6NewCC_RatioToNominal_staterr[icent]       = sqrt( pow( vn6NewCC_staterr/vn6,2) + pow(vn6NewCC*vn6_staterr/vn6/vn6,2) );
    vn8NewCC_RatioToNominal_staterr[icent]       = sqrt( pow( vn8NewCC_staterr/vn8,2) + pow(vn8NewCC*vn8_staterr/vn8/vn8,2) );
    gamma1expNewCC_RatioToNominal_staterr[icent] = sqrt( pow( gamma1expNewCC_staterr/gamma1exp,2) + pow(gamma1expNewCC*gamma1exp_staterr/gamma1exp/gamma1exp,2) );
    vn6vn4NewCC_RatioToNominal_staterr[icent]    = sqrt( pow( vn6vn4NewCC_staterr/vn6vn4,2) + pow(vn6vn4NewCC*vn6vn4_staterr/vn6vn4/vn6vn4,2) );
    vn8vn4NewCC_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn4NewCC_staterr/vn8vn4,2) + pow(vn8vn4NewCC*vn8vn4_staterr/vn8vn4/vn8vn4,2) );
    vn8vn6NewCC_RatioToNominal_staterr[icent]    = sqrt( pow( vn8vn6NewCC_staterr/vn8vn6,2) + pow(vn8vn6NewCC*vn8vn6_staterr/vn8vn6/vn8vn6,2) );

    //-- Calculate pct difference relative to nominal
    if( vn2 == 0 || vn2NewCC == 0 || compatibleWithOne(vn2NewCC_RatioToNominal[icent], vn2NewCC_RatioToNominal_staterr[icent]) ) vn2NewCC_PctDiffToNominal[icent] = 0;
    else                            vn2NewCC_PctDiffToNominal[icent] = fabs( vn2NewCC - vn2 ) / fabs( vn2 );

    if( vn4 == 0 || vn4NewCC == 0 || compatibleWithOne(vn4NewCC_RatioToNominal[icent], vn4NewCC_RatioToNominal_staterr[icent])) vn4NewCC_PctDiffToNominal[icent] = 0;
    else                            vn4NewCC_PctDiffToNominal[icent] = fabs( vn4NewCC - vn4 ) / fabs( vn4 );

    if( vn6 == 0 || vn6NewCC == 0 || compatibleWithOne(vn6NewCC_RatioToNominal[icent], vn6NewCC_RatioToNominal_staterr[icent])) vn6NewCC_PctDiffToNominal[icent] = 0;
    else                            vn6NewCC_PctDiffToNominal[icent] = fabs( vn6NewCC - vn6 ) / fabs( vn6 );

    if( vn8 == 0 || vn8NewCC == 0 || compatibleWithOne(vn8NewCC_RatioToNominal[icent], vn8NewCC_RatioToNominal_staterr[icent])) vn8NewCC_PctDiffToNominal[icent] = 0;
    else                            vn8NewCC_PctDiffToNominal[icent] = fabs( vn8NewCC - vn8 ) / fabs( vn8 );

    if( gamma1exp == 0 || gamma1expNewCC == 0 || compatibleWithOne(gamma1expNewCC_RatioToNominal[icent], gamma1expNewCC_RatioToNominal_staterr[icent]) ) gamma1expNewCC_PctDiffToNominal[icent] = 0;
    else                                        gamma1expNewCC_PctDiffToNominal[icent] = fabs( gamma1expNewCC - gamma1exp ) / fabs( gamma1exp );

    if( vn6vn4 == 0 || vn6vn4NewCC == 0 || compatibleWithOne(vn6vn4NewCC_RatioToNominal[icent], vn6vn4NewCC_RatioToNominal_staterr[icent]) ) vn6vn4NewCC_PctDiffToNominal[icent] = 0;
    else                                  vn6vn4NewCC_PctDiffToNominal[icent] = fabs( vn6vn4NewCC - vn6vn4 ) / fabs( vn6vn4 );

    if( vn8vn4 == 0 || vn8vn4NewCC == 0 || compatibleWithOne(vn8vn4NewCC_RatioToNominal[icent], vn8vn4NewCC_RatioToNominal_staterr[icent]) ) vn8vn4NewCC_PctDiffToNominal[icent] = 0;
    else                                  vn8vn4NewCC_PctDiffToNominal[icent] = fabs( vn8vn4NewCC - vn8vn4 ) / fabs( vn8vn4 );

    if( vn8vn6 == 0 || vn8vn6NewCC == 0 || compatibleWithOne(vn8vn6NewCC_RatioToNominal[icent], vn8vn6NewCC_RatioToNominal_staterr[icent]) ) vn8vn6NewCC_PctDiffToNominal[icent] = 0;
    else                                  vn8vn6NewCC_PctDiffToNominal[icent] = fabs( vn8vn6NewCC - vn8vn6 ) / fabs( vn8vn6 );

  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //-- Make sweet, sweet TGraphErrors
  double cErr[NCENT];
  for(int i = 0; i < NCENT; i++) cErr[i] = 0;
  grVn2NewCC_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn2NewCC_RatioToNominal, cErr, vn2NewCC_RatioToNominal_staterr);
  formatGraph(grVn2NewCC_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} [NewCC] / [Nominal]", norder_), 1, 24, "grVn2NewCC_RatioToNominal");
  grVn4NewCC_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn4NewCC_RatioToNominal, cErr, vn4NewCC_RatioToNominal_staterr);
  formatGraph(grVn4NewCC_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} [NewCC] / [Nominal]", norder_), kSpring+4, 25, "grVn4NewCC_RatioToNominal");
  grVn6NewCC_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn6NewCC_RatioToNominal, cErr, vn6NewCC_RatioToNominal_staterr);
  formatGraph(grVn6NewCC_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} [NewCC] / [Nominal]", norder_), 6, 28, "grVn6NewCC_RatioToNominal");
  grVn8NewCC_RatioToNominal       = new TGraphErrors(NCENT, centBinCenter, vn8NewCC_RatioToNominal, cErr, vn8NewCC_RatioToNominal_staterr);
  formatGraph(grVn8NewCC_RatioToNominal, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} [NewCC] / [Nominal]", norder_), kOrange+7, 27, "grVn8NewCC_RatioToNominal");
  grGamma1ExpNewCC_RatioToNominal = new TGraphErrors(NCENT, centBinCenter, gamma1expNewCC_RatioToNominal, cErr, gamma1expNewCC_RatioToNominal_staterr);
  formatGraph(grGamma1ExpNewCC_RatioToNominal, "Centrality %", ratioG1Min, ratioG1Max, "|1-#gamma_{1}^{exp} [NewCC] / [Nominal]|", 2, 20, "grGamma1ExpNewCC_RatioToNominal");
  grVn6Vn4NewCC_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn6vn4NewCC_RatioToNominal, cErr, vn6vn4NewCC_RatioToNominal_staterr);
  formatGraph(grVn6Vn4NewCC_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{6}/v_{%i}{4} [NewCC] / [Nominal]", norder_, norder_), 4, 21, "grVn6Vn4NewCC_RatioToNominal");
  grVn8Vn4NewCC_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn4NewCC_RatioToNominal, cErr, vn8vn4NewCC_RatioToNominal_staterr);
  formatGraph(grVn8Vn4NewCC_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{4} [NewCC] / [Nominal]", norder_, norder_), kGreen+2, 34, "grVn8Vn4NewCC_RatioToNominal");
  grVn8Vn6NewCC_RatioToNominal    = new TGraphErrors(NCENT, centBinCenter, vn8vn6NewCC_RatioToNominal, cErr, vn8vn6NewCC_RatioToNominal_staterr);
  formatGraph(grVn8Vn6NewCC_RatioToNominal, "Centrality %", ratioCumuRatioMin, ratioCumuRatioMax, Form("v_{%i}{8}/v_{%i}{6} [NewCC] / [Nominal]", norder_, norder_), kViolet-1, 33, "grVn8Vn6NewCC_RatioToNominal");

  //-- Raw Values
  grVn2Nominal    = new TGraphErrors(NCENT, centBinCenter, vn2Nominal,    cErr, vn2Nominale);
  grVn4Nominal    = new TGraphErrors(NCENT, centBinCenter, vn4Nominal,    cErr, vn4Nominale);
  grVn6Nominal    = new TGraphErrors(NCENT, centBinCenter, vn6Nominal,    cErr, vn6Nominale);
  grVn8Nominal    = new TGraphErrors(NCENT, centBinCenter, vn8Nominal,    cErr, vn8Nominale);
  grG1eNominal    = new TGraphErrors(NCENT, centBinCenter, g1eNominal,    cErr, g1eNominale);
  grVn6Vn4Nominal = new TGraphErrors(NCENT, centBinCenter, vn6vn4Nominal, cErr, vn6vn4Nominale);
  grVn8Vn4Nominal = new TGraphErrors(NCENT, centBinCenter, vn8vn4Nominal, cErr, vn8vn4Nominale);
  grVn8Vn6Nominal = new TGraphErrors(NCENT, centBinCenter, vn8vn6Nominal, cErr, vn8vn6Nominale);

  formatGraph(grVn2Nominal,    "Centrality %", 0.,    0.15,  Form("v_{%i}{2}", norder_),                    1,         20, "grVn2Nominal");
  formatGraph(grVn4Nominal,    "Centrality %", 0.,    0.15,  Form("v_{%i}{4}", norder_),                    kSpring+4, 21, "grVn4Nominal");
  formatGraph(grVn6Nominal,    "Centrality %", 0.,    0.15,  Form("v_{%i}{6}", norder_),                    6,         22, "grVn6Nominal");
  formatGraph(grVn8Nominal,    "Centrality %", 0.,    0.15,  Form("v_{%i}{8}", norder_),                    kOrange+7, 23, "grVn8Nominal");
  formatGraph(grG1eNominal,    "Centrality %", -1.0,  0.5,   "#gamma_{1}^{exp}",                            2,         20, "grG1eNominal");
  formatGraph(grVn6Vn4Nominal, "Centrality %", 0.95,  1.03,  Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         21, "grVn6Vn4Nominal");
  formatGraph(grVn8Vn4Nominal, "Centrality %", 0.95,  1.03,  Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  34, "grVn8Vn4Nominal");
  formatGraph(grVn8Vn6Nominal, "Centrality %", 0.993, 1.003, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 20, "grVn8Vn6Nominal");

  grVn2NewCC    = new TGraphErrors(NCENT, centBinCenter, vn2NEWCC,    cErr, vn2NEWCCe);
  grVn4NewCC    = new TGraphErrors(NCENT, centBinCenter, vn4NEWCC,    cErr, vn4NEWCCe);
  grVn6NewCC    = new TGraphErrors(NCENT, centBinCenter, vn6NEWCC,    cErr, vn6NEWCCe);
  grVn8NewCC    = new TGraphErrors(NCENT, centBinCenter, vn8NEWCC,    cErr, vn8NEWCCe);
  grG1eNewCC    = new TGraphErrors(NCENT, centBinCenter, g1eNEWCC,    cErr, g1eNEWCCe);
  grVn6Vn4NewCC = new TGraphErrors(NCENT, centBinCenter, vn6vn4NEWCC, cErr, vn6vn4NEWCCe);
  grVn8Vn4NewCC = new TGraphErrors(NCENT, centBinCenter, vn8vn4NEWCC, cErr, vn8vn4NEWCCe);
  grVn8Vn6NewCC = new TGraphErrors(NCENT, centBinCenter, vn8vn6NEWCC, cErr, vn8vn6NEWCCe);

  formatGraph(grVn2NewCC,    "Centrality %", 0.,    0.15,  Form("v_{%i}{2}", norder_),                    1,         24, "grVn2NewCC");
  formatGraph(grVn4NewCC,    "Centrality %", 0.,    0.15,  Form("v_{%i}{4}", norder_),                    kSpring+4, 25, "grVn4NewCC");
  formatGraph(grVn6NewCC,    "Centrality %", 0.,    0.15,  Form("v_{%i}{6}", norder_),                    6,         26, "grVn6NewCC");
  formatGraph(grVn8NewCC,    "Centrality %", 0.,    0.15,  Form("v_{%i}{8}", norder_),                    kOrange+7, 32, "grVn8NewCC");
  formatGraph(grG1eNewCC,    "Centrality %", -1.0,  0.5,   "#gamma_{1}^{exp}",                            2,         24, "grG1eNewCC");
  formatGraph(grVn6Vn4NewCC, "Centrality %", 0.95,  1.03,  Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         25, "grVn6Vn4NewCC");
  formatGraph(grVn8Vn4NewCC, "Centrality %", 0.95,  1.03,  Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  28, "grVn8Vn4NewCC");
  formatGraph(grVn8Vn6NewCC, "Centrality %", 0.993, 1.003, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 24, "grVn8Vn6NewCC");


  TLine * line = new TLine(grVn2NewCC_RatioToNominal->GetXaxis()->GetXmin(), 1.0, grVn2NewCC_RatioToNominal->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);


  TCanvas * cCumuSys = new TCanvas("cCumuSys", "cCumuSys", 1000, 1000);
  cCumuSys->Divide(2,2);
  cCumuSys->cd(1);
  grVn2NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(2);
  grVn4NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(3);
  grVn6NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->cd(4);
  grVn8NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuSys->SaveAs("../../plots/systematicStudies/cSysNewCC_CumuCent.pdf");

  TCanvas * cGamma1ExpSys = new TCanvas("cGamma1ExpSys", "cGamma1ExpSys", 500, 500);
  cGamma1ExpSys->cd();
  grGamma1ExpNewCC_RatioToNominal->Draw("ap");
  //line->Draw("same");
  cGamma1ExpSys->SaveAs("../../plots/systematicStudies/cSysNewCC_Gamma1ExpCent.pdf");

  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  grVn6Vn4NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6NewCC_RatioToNominal->Draw("ap");
  line->Draw("same");
  cCumuRatioSys->SaveAs("../../plots/systematicStudies/cSysNewCC_CumuRatioCent.pdf");

  //-- Raw Values -------------------------------
  TLegend * legvn2 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn2->SetBorderSize(0);
  legvn2->SetFillStyle(0);
  legvn2->AddEntry(grVn2Nominal, "Old Tune", "lp");
  legvn2->AddEntry(grVn2NewCC,   "New Tune", "lp");

  TLegend * legvn4 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn4->SetBorderSize(0);
  legvn4->SetFillStyle(0);
  legvn4->AddEntry(grVn4Nominal, "Old Tune", "lp");
  legvn4->AddEntry(grVn4NewCC,   "New Tune", "lp");

  TLegend * legvn6 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn6->SetBorderSize(0);
  legvn6->SetFillStyle(0);
  legvn6->AddEntry(grVn6Nominal, "Old Tune", "lp");
  legvn6->AddEntry(grVn6NewCC,   "New Tune", "lp");

  TLegend * legvn8 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn8->SetBorderSize(0);
  legvn8->SetFillStyle(0);
  legvn8->AddEntry(grVn8Nominal, "Old Tune", "lp");
  legvn8->AddEntry(grVn8NewCC,   "New Tune", "lp");


  TCanvas * cCumuRaw = new TCanvas("cCumuRaw", "cCumuRaw", 1000, 1000);
  cCumuRaw->Divide(2,2);
  cCumuRaw->cd(1);
  grVn2NewCC->Draw("ap");
  grVn2Nominal->Draw("psame");
  legvn2->Draw("same");
  cCumuRaw->cd(2);
  grVn4NewCC->Draw("ap");
  grVn4Nominal->Draw("psame");
  legvn4->Draw("same");
  cCumuRaw->cd(3);
  grVn6NewCC->Draw("ap");
  grVn6Nominal->Draw("psame");
  legvn6->Draw("same");
  cCumuRaw->cd(4);
  grVn8NewCC->Draw("ap");
  grVn8Nominal->Draw("psame");
  legvn8->Draw("same");
  cCumuRaw->SaveAs("../../plots/systematicStudies/cRawNewCC_CumuCent.pdf");

  TLegend * legG1e = new TLegend(0.6, 0.8, 0.9, 0.9);
  legG1e->SetBorderSize(0);
  legG1e->SetFillStyle(0);
  legG1e->AddEntry(grG1eNominal, "Old Tune", "lp");
  legG1e->AddEntry(grG1eNewCC,   "New Tune", "lp");

  TCanvas * cGamma1ExpRaw = new TCanvas("cGamma1ExpRaw", "cGamma1ExpRaw", 500, 500);
  cGamma1ExpRaw->cd();
  grG1eNewCC->Draw("ap");
  grG1eNominal->Draw("psame");
  legG1e->Draw("same");
  cGamma1ExpRaw->SaveAs("../../plots/systematicStudies/cRawNewCC_Gamma1ExpCent.pdf");

  TLegend * legvn6vn4 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Nominal, "Old Tune", "lp");
  legvn6vn4->AddEntry(grVn6Vn4NewCC,   "New Tune", "lp");

  TLegend * legvn8vn4 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Nominal, "Old Tune", "lp");
  legvn8vn4->AddEntry(grVn8Vn4NewCC,   "New Tune", "lp");

  TLegend * legvn8vn6 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Nominal, "Old Tune", "lp");
  legvn8vn6->AddEntry(grVn8Vn6NewCC,   "New Tune", "lp");

  TCanvas * cCumuRatioRaw = new TCanvas("cCumuRatioRaw", "cCumuRatioRaw", 1500, 500);
  cCumuRatioRaw->Divide(3,1);
  cCumuRatioRaw->cd(1);
  grVn6Vn4NewCC->Draw("ap");
  grVn6Vn4Nominal->Draw("psame");
  legvn6vn4->Draw("same");
  cCumuRatioRaw->cd(2);
  grVn8Vn4NewCC->Draw("ap");
  grVn8Vn4Nominal->Draw("psame");
  legvn8vn4->Draw("same");
  cCumuRatioRaw->cd(3);
  grVn8Vn6NewCC->Draw("ap");
  grVn8Vn6Nominal->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatioRaw->SaveAs("../../plots/systematicStudies/cRawNewCC_CumuRatioCent.pdf");

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){

    std::cout<<"---------------"<<std::endl;

    relErrVn2->SetBinContent(icent+1, vn2NewCC_PctDiffToNominal[icent]);
    relErrVn4->SetBinContent(icent+1, vn4NewCC_PctDiffToNominal[icent]);
    relErrVn6->SetBinContent(icent+1, vn6NewCC_PctDiffToNominal[icent]);
    relErrVn8->SetBinContent(icent+1, vn6NewCC_PctDiffToNominal[icent]);
    relErrGamma1Exp->SetBinContent(icent+1, gamma1expNewCC_PctDiffToNominal[icent]);
    relErrVn6Vn4->SetBinContent(icent+1, vn6vn4NewCC_PctDiffToNominal[icent]);
    relErrVn8Vn4->SetBinContent(icent+1, vn8vn4NewCC_PctDiffToNominal[icent]);
    relErrVn8Vn6->SetBinContent(icent+1, vn8vn6NewCC_PctDiffToNominal[icent]);

  }

  fRelErr->Write();

}

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;


void sysClusCompatTune(){

  const int norder_ = 2;

  double vn6vn4Min = 0.95;
  double vn6vn4Max = 1.03;
  double vn8vn4Min = 0.95;
  double vn8vn4Max = 1.03;
  double vn8vn6Min = 0.99;
  double vn8vn6Max = 1.002;
  double g1eMin    = -1.0;
  double g1eMax    = 0.5;

  double rMin = 0.95;
  double rMax = 1.05;
  double g1rMin = 0.;
  double g1rMax = 5;

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
  TH1D * hVarianceOfMean_Vn6Vn4_Default;
  TH1D * hVarianceOfMean_Vn8Vn4_Default;
  TH1D * hVarianceOfMean_Vn8Vn6_Default;
  TH1D * hVarianceOfMean_Gamma1Exp_Default;

  double vn6vn4Default[NCENT];
  double vn8vn4Default[NCENT];
  double vn8vn6Default[NCENT];
  double g1eDefault[NCENT];

  double vn6vn4Default_err[NCENT];
  double vn8vn4Default_err[NCENT];
  double vn8vn6Default_err[NCENT];
  double g1eDefault_err[NCENT];

  TGraphErrors * grVn6Vn4Default;
  TGraphErrors * grVn8Vn4Default;
  TGraphErrors * grVn8Vn6Default;
  TGraphErrors * grG1EDefault;

  //-- 0.5%
  TFile * fAna0p5pct;
  TFile * fUnf0p5pct;
  TH1D * hObs0p5pct[NCENT];
  TH1D * hMult0p5pct[NCENT];
  TH1D * hUnfold0p5pct[NCENT][NITER];
  TH1D * hRefold0p5pct[NCENT][NITER];
  bool stopFound0p5pct;
  int iterCut0p5pct[NCENT];

  TFile * fStat0p5pct;
  TH1D * hVarianceOfMean_Vn6Vn4_0p5pct;
  TH1D * hVarianceOfMean_Vn8Vn4_0p5pct;
  TH1D * hVarianceOfMean_Vn8Vn6_0p5pct;
  TH1D * hVarianceOfMean_Gamma1Exp_0p5pct;

  double vn6vn40p5pct[NCENT];
  double vn8vn40p5pct[NCENT];
  double vn8vn60p5pct[NCENT];
  double g1e0p5pct[NCENT];

  double vn6vn40p5pct_err[NCENT];
  double vn8vn40p5pct_err[NCENT];
  double vn8vn60p5pct_err[NCENT];
  double g1e0p5pct_err[NCENT];

  double vn6vn40p5pct_RatioToDefault[NCENT];
  double vn8vn40p5pct_RatioToDefault[NCENT];
  double vn8vn60p5pct_RatioToDefault[NCENT];
  double g1e0p5pct_RatioToDefault[NCENT];

  double vn6vn40p5pct_RatioToDefault_err[NCENT];
  double vn8vn40p5pct_RatioToDefault_err[NCENT];
  double vn8vn60p5pct_RatioToDefault_err[NCENT];
  double g1e0p5pct_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn40p5pct;
  TGraphErrors * grVn8Vn40p5pct;
  TGraphErrors * grVn8Vn60p5pct;
  TGraphErrors * grG1E0p5pct;

  TGraphErrors * grVn6Vn40p5pct_RatioToDefault;
  TGraphErrors * grVn8Vn40p5pct_RatioToDefault;
  TGraphErrors * grVn8Vn60p5pct_RatioToDefault;
  TGraphErrors * grG1E0p5pct_RatioToDefault;

  //-- 1%
  TFile * fAna1pct;
  TFile * fUnf1pct;
  TH1D * hObs1pct[NCENT];
  TH1D * hMult1pct[NCENT];
  TH1D * hUnfold1pct[NCENT][NITER];
  TH1D * hRefold1pct[NCENT][NITER];
  bool stopFound1pct;
  int iterCut1pct[NCENT];

  TFile * fStat1pct;
  TH1D * hVarianceOfMean_Vn6Vn4_1pct;
  TH1D * hVarianceOfMean_Vn8Vn4_1pct;
  TH1D * hVarianceOfMean_Vn8Vn6_1pct;
  TH1D * hVarianceOfMean_Gamma1Exp_1pct;

  double vn6vn41pct[NCENT];
  double vn8vn41pct[NCENT];
  double vn8vn61pct[NCENT];
  double g1e1pct[NCENT];

  double vn6vn41pct_err[NCENT];
  double vn8vn41pct_err[NCENT];
  double vn8vn61pct_err[NCENT];
  double g1e1pct_err[NCENT];

  double vn6vn41pct_RatioToDefault[NCENT];
  double vn8vn41pct_RatioToDefault[NCENT];
  double vn8vn61pct_RatioToDefault[NCENT];
  double g1e1pct_RatioToDefault[NCENT];

  double vn6vn41pct_RatioToDefault_err[NCENT];
  double vn8vn41pct_RatioToDefault_err[NCENT];
  double vn8vn61pct_RatioToDefault_err[NCENT];
  double g1e1pct_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn41pct;
  TGraphErrors * grVn8Vn41pct;
  TGraphErrors * grVn8Vn61pct;
  TGraphErrors * grG1E1pct;

  TGraphErrors * grVn6Vn41pct_RatioToDefault;
  TGraphErrors * grVn8Vn41pct_RatioToDefault;
  TGraphErrors * grVn8Vn61pct_RatioToDefault;
  TGraphErrors * grG1E1pct_RatioToDefault;

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
  TH1D * hVarianceOfMean_Vn6Vn4_2pct;
  TH1D * hVarianceOfMean_Vn8Vn4_2pct;
  TH1D * hVarianceOfMean_Vn8Vn6_2pct;
  TH1D * hVarianceOfMean_Gamma1Exp_2pct;

  double vn6vn42pct[NCENT];
  double vn8vn42pct[NCENT];
  double vn8vn62pct[NCENT];
  double g1e2pct[NCENT];

  double vn6vn42pct_err[NCENT];
  double vn8vn42pct_err[NCENT];
  double vn8vn62pct_err[NCENT];
  double g1e2pct_err[NCENT];

  double vn6vn42pct_RatioToDefault[NCENT];
  double vn8vn42pct_RatioToDefault[NCENT];
  double vn8vn62pct_RatioToDefault[NCENT];
  double g1e2pct_RatioToDefault[NCENT];

  double vn6vn42pct_RatioToDefault_err[NCENT];
  double vn8vn42pct_RatioToDefault_err[NCENT];
  double vn8vn62pct_RatioToDefault_err[NCENT];
  double g1e2pct_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn42pct;
  TGraphErrors * grVn8Vn42pct;
  TGraphErrors * grVn8Vn62pct;
  TGraphErrors * grG1E2pct;

  TGraphErrors * grVn6Vn42pct_RatioToDefault;
  TGraphErrors * grVn8Vn42pct_RatioToDefault;
  TGraphErrors * grVn8Vn62pct_RatioToDefault;
  TGraphErrors * grG1E2pct_RatioToDefault;

  //-- 3%
  TFile * fAna3pct;
  TFile * fUnf3pct;
  TH1D * hObs3pct[NCENT];
  TH1D * hMult3pct[NCENT];
  TH1D * hUnfold3pct[NCENT][NITER];
  TH1D * hRefold3pct[NCENT][NITER];
  bool stopFound3pct;
  int iterCut3pct[NCENT];

  TFile * fStat3pct;
  TH1D * hVarianceOfMean_Vn6Vn4_3pct;
  TH1D * hVarianceOfMean_Vn8Vn4_3pct;
  TH1D * hVarianceOfMean_Vn8Vn6_3pct;
  TH1D * hVarianceOfMean_Gamma1Exp_3pct;

  double vn6vn43pct[NCENT];
  double vn8vn43pct[NCENT];
  double vn8vn63pct[NCENT];
  double g1e3pct[NCENT];

  double vn6vn43pct_err[NCENT];
  double vn8vn43pct_err[NCENT];
  double vn8vn63pct_err[NCENT];
  double g1e3pct_err[NCENT];

  double vn6vn43pct_RatioToDefault[NCENT];
  double vn8vn43pct_RatioToDefault[NCENT];
  double vn8vn63pct_RatioToDefault[NCENT];
  double g1e3pct_RatioToDefault[NCENT];

  double vn6vn43pct_RatioToDefault_err[NCENT];
  double vn8vn43pct_RatioToDefault_err[NCENT];
  double vn8vn63pct_RatioToDefault_err[NCENT];
  double g1e3pct_RatioToDefault_err[NCENT];

  TGraphErrors * grVn6Vn43pct;
  TGraphErrors * grVn8Vn43pct;
  TGraphErrors * grVn8Vn63pct;
  TGraphErrors * grG1E3pct;

  TGraphErrors * grVn6Vn43pct_RatioToDefault;
  TGraphErrors * grVn8Vn43pct_RatioToDefault;
  TGraphErrors * grVn8Vn63pct_RatioToDefault;
  TGraphErrors * grG1E3pct_RatioToDefault;


  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get Analyzer files
  fAnaDefault = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fAna0p5pct  = new TFile( "newCCTune0p5pct/AnalyzerResults/CastleEbyE.root" );
  fAna1pct    = new TFile( "newCCTune1pct/AnalyzerResults/CastleEbyE.root" );
  fAna2pct    = new TFile( "newCCTune2pct/AnalyzerResults/CastleEbyE.root" );
  fAna3pct    = new TFile( "newCCTune3pct/AnalyzerResults/CastleEbyE.root" );

  //-- Get Unfold files
  fUnfDefault = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnf0p5pct  = new TFile( Form("newCCTune0p5pct/UnfoldResults/dataResp/data%i.root", norder_)  );
  fUnf1pct    = new TFile( Form("newCCTune1pct/UnfoldResults/dataResp/data%i.root", norder_)  );
  fUnf2pct    = new TFile( Form("newCCTune2pct/UnfoldResults/dataResp/data%i.root", norder_)  );
  fUnf3pct    = new TFile( Form("newCCTune3pct/UnfoldResults/dataResp/data%i.root", norder_)  );

  //-- Get Statistical Uncertainties
  fStatDefault = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_Default    = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_Default = (TH1D*) fStatDefault->Get( "hVarianceOfMean_Gamma1Exp" );

  fStat0p5pct  = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/clusCompatTune/newCCTune0p5pct/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_0p5pct    = (TH1D*) fStat0p5pct->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_0p5pct    = (TH1D*) fStat0p5pct->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_0p5pct    = (TH1D*) fStat0p5pct->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_0p5pct = (TH1D*) fStat0p5pct->Get( "hVarianceOfMean_Gamma1Exp" );

  fStat1pct  = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/clusCompatTune/newCCTune1pct/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_1pct    = (TH1D*) fStat1pct->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_1pct    = (TH1D*) fStat1pct->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_1pct    = (TH1D*) fStat1pct->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_1pct = (TH1D*) fStat1pct->Get( "hVarianceOfMean_Gamma1Exp" );

  fStat2pct  = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/clusCompatTune/newCCTune2pct/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_2pct    = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_2pct = (TH1D*) fStat2pct->Get( "hVarianceOfMean_Gamma1Exp" );

  fStat3pct  = new TFile( Form("../../../../statErrorHandle/v%i/eta2.4/systematicStudies/clusCompatTune/newCCTune3pct/StatisticalUncertainties_v%i.root", norder_, norder_) );
  hVarianceOfMean_Vn6Vn4_3pct    = (TH1D*) fStat3pct->Get( "hVarianceOfMean_Vn6Vn4" );
  hVarianceOfMean_Vn8Vn4_3pct    = (TH1D*) fStat3pct->Get( "hVarianceOfMean_Vn8Vn4" );
  hVarianceOfMean_Vn8Vn6_3pct    = (TH1D*) fStat3pct->Get( "hVarianceOfMean_Vn8Vn6" );
  hVarianceOfMean_Gamma1Exp_3pct = (TH1D*) fStat3pct->Get( "hVarianceOfMean_Gamma1Exp" );


  //-- Start Grabbing Histos
  for(int icent = 0; icent < NCENT; icent++){

    //-- Vn Observed 
    hObsDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObsDefault[icent]->SetName( Form("hVnFullDefault_c%i", icent) );
    hObs0p5pct[icent]  = (TH1D*) fAna0p5pct->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs0p5pct[icent]->SetName( Form("hVnFull0p5pct_c%i", icent) );
    hObs1pct[icent] = (TH1D*) fAna1pct->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs1pct[icent]->SetName( Form("hVnFull1pct_c%i", icent) );
    hObs2pct[icent] = (TH1D*) fAna2pct->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs2pct[icent]->SetName( Form("hVnFull2pct_c%i", icent) );
    hObs3pct[icent] = (TH1D*) fAna3pct->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs3pct[icent]->SetName( Form("hVnFull3pct_c%i", icent) );

    //-- Multiplicity
    hMultDefault[icent] = (TH1D*) fAnaDefault->Get( Form("qwebye/Mult_c%i", icent) );
    hMultDefault[icent]->SetLineColor(centCol[icent]);
    hMultDefault[icent]->SetMarkerColor(centCol[icent]);
    hMultDefault[icent]->SetMaximum( 20.*hMultDefault[icent]->GetMaximum() );
    hMultDefault[icent]->GetXaxis()->SetNdivisions(509);
    hMultDefault[icent]->GetXaxis()->SetRange(1, hMultDefault[icent]->FindBin(5500));

    hMult0p5pct[icent] = (TH1D*) fAna0p5pct->Get( Form("qwebye/Mult_c%i", icent) );
    hMult0p5pct[icent]->SetLineColor(centCol[icent]);
    hMult0p5pct[icent]->SetMarkerColor(centCol[icent]);
    hMult0p5pct[icent]->SetMaximum( 20.*hMult0p5pct[icent]->GetMaximum() );
    hMult0p5pct[icent]->GetXaxis()->SetNdivisions(509);
    hMult0p5pct[icent]->GetXaxis()->SetRange(1, hMult0p5pct[icent]->FindBin(5500));

    hMult1pct[icent] = (TH1D*) fAna1pct->Get( Form("qwebye/Mult_c%i", icent) );
    hMult1pct[icent]->SetLineColor(centCol[icent]);
    hMult1pct[icent]->SetMarkerColor(centCol[icent]);
    hMult1pct[icent]->SetMaximum( 20.*hMult1pct[icent]->GetMaximum() );
    hMult1pct[icent]->GetXaxis()->SetNdivisions(509);
    hMult1pct[icent]->GetXaxis()->SetRange(1, hMult1pct[icent]->FindBin(5500));

    hMult2pct[icent] = (TH1D*) fAna2pct->Get( Form("qwebye/Mult_c%i", icent) );
    hMult2pct[icent]->SetLineColor(centCol[icent]);
    hMult2pct[icent]->SetMarkerColor(centCol[icent]);
    hMult2pct[icent]->SetMaximum( 20.*hMult2pct[icent]->GetMaximum() );
    hMult2pct[icent]->GetXaxis()->SetNdivisions(509);
    hMult2pct[icent]->GetXaxis()->SetRange(1, hMult2pct[icent]->FindBin(5500));

    hMult3pct[icent] = (TH1D*) fAna3pct->Get( Form("qwebye/Mult_c%i", icent) );
    hMult3pct[icent]->SetLineColor(centCol[icent]);
    hMult3pct[icent]->SetMarkerColor(centCol[icent]);
    hMult3pct[icent]->SetMaximum( 20.*hMult3pct[icent]->GetMaximum() );
    hMult3pct[icent]->GetXaxis()->SetNdivisions(509);
    hMult3pct[icent]->GetXaxis()->SetRange(1, hMult3pct[icent]->FindBin(5500));

    //-- Reset stopFound boolean
    stopFoundDefault = 0;
    stopFound0p5pct  = 0;
    stopFound1pct    = 0;
    stopFound2pct    = 0;
    stopFound3pct    = 0;

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

      //-- 0.5%
      hUnfold0p5pct[icent][i] = (TH1D*) fUnf0p5pct->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold0p5pct[icent][i]->SetName( Form("hreco0p5pct%i_c%i", iter[i], icent) );
      hRefold0p5pct[icent][i] = (TH1D*) fUnf0p5pct->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold0p5pct[icent][i]->SetName( Form("hrefol0p5pct%i_c%i", iter[i], icent) );

      double chi20p5pct = hRefold0p5pct[icent][i]->Chi2Test(hObs0p5pct[icent], "CHI2/NDF");
      if( chi20p5pct <= 1.2 && !stopFound0p5pct ){
        stopFound0p5pct = 1;
        iterCut0p5pct[icent] = i;
      } 
      if( i == NITER-1 && !stopFound0p5pct ){
        stopFound0p5pct = 1;
        iterCut0p5pct[icent] = i;
      }

      //-- 1.0%
      hUnfold1pct[icent][i] = (TH1D*) fUnf1pct->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold1pct[icent][i]->SetName( Form("hreco1pct%i_c%i", iter[i], icent) );
      hRefold1pct[icent][i] = (TH1D*) fUnf1pct->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold1pct[icent][i]->SetName( Form("hrefol1pct%i_c%i", iter[i], icent) );

      double chi21pct = hRefold1pct[icent][i]->Chi2Test(hObs1pct[icent], "CHI2/NDF");
      if( chi21pct <= 1.2 && !stopFound1pct ){
        stopFound1pct = 1;
        iterCut1pct[icent] = i;
      } 
      if( i == NITER-1 && !stopFound1pct ){
        stopFound1pct = 1;
        iterCut1pct[icent] = i;
      } 

      //-- 2.0%
      hUnfold2pct[icent][i] = (TH1D*) fUnf2pct->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold2pct[icent][i]->SetName( Form("hreco2pct%i_c%i", iter[i], icent) );
      hRefold2pct[icent][i] = (TH1D*) fUnf2pct->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold2pct[icent][i]->SetName( Form("hrefol2pct%i_c%i", iter[i], icent) );

      double chi22pct = hRefold2pct[icent][i]->Chi2Test(hObs2pct[icent], "CHI2/NDF");
      if( chi22pct <= 1.2 && !stopFound2pct ){
        stopFound2pct = 1;
        iterCut2pct[icent] = i;
      }
      if( i == NITER-1 && !stopFound2pct ){
        stopFound2pct = 1;
        iterCut2pct[icent] = i;
      }

      //-- 3.0%
      hUnfold3pct[icent][i] = (TH1D*) fUnf3pct->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold3pct[icent][i]->SetName( Form("hreco3pct%i_c%i", iter[i], icent) );
      hRefold3pct[icent][i] = (TH1D*) fUnf3pct->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold3pct[icent][i]->SetName( Form("hrefol3pct%i_c%i", iter[i], icent) );

      double chi23pct = hRefold3pct[icent][i]->Chi2Test(hObs3pct[icent], "CHI2/NDF");
      if( chi23pct <= 1.2 && !stopFound3pct ){
        stopFound3pct = 1;
        iterCut3pct[icent] = i;
      }
      if( i == NITER-1 && !stopFound3pct ){
        stopFound3pct = 1;
        iterCut3pct[icent] = i;
      }

    } //-- End iter loop

  } //-- End cent loop


  //-- Now that we're done fetching things, let's build some arrays!
  for(int icent = 0; icent < NCENT; icent++){

    //-- Default
    int iDefault = iterCutDefault[icent];
    FixUnfold(hUnfoldDefault[icent][iDefault]);
    EbyECumu cumuDefault(hUnfoldDefault[icent][iDefault]);
    double vn4Default  = cumuDefault.GetCumu_vn4();
    double vn6Default  = cumuDefault.GetCumu_vn6();
    double vn8Default  = cumuDefault.GetCumu_vn8();
    double g1exDefault = cumuDefault.GetGamma1Exp();

    if(vn4Default == 0 || vn6Default == 0) vn6vn4Default[icent]= -1000.;
    else                                             vn6vn4Default[icent] = vn6Default / vn4Default;
    if(vn4Default == 0 || vn8Default == 0) vn8vn4Default[icent]= -1000.;
    else                                             vn8vn4Default[icent] = vn8Default / vn4Default;
    if(vn6Default == 0 || vn8Default == 0) vn8vn6Default[icent]= -1000.;
    else                                             vn8vn6Default[icent] = vn8Default / vn6Default;
    g1eDefault[icent] = g1exDefault;

    double vn6vn4Default_staterr = sqrt( hVarianceOfMean_Vn6Vn4_Default->GetBinContent(icent+1) );
    double vn8vn4Default_staterr = sqrt( hVarianceOfMean_Vn8Vn4_Default->GetBinContent(icent+1) );
    double vn8vn6Default_staterr = sqrt( hVarianceOfMean_Vn8Vn6_Default->GetBinContent(icent+1) );
    double g1eDefault_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_Default->GetBinContent(icent+1) );

    vn6vn4Default_err[icent] = vn6vn4Default_staterr;
    vn8vn4Default_err[icent] = vn8vn4Default_staterr;
    vn8vn6Default_err[icent] = vn8vn6Default_staterr;
    g1eDefault_err[icent]    = g1eDefault_staterr;

    //-- 0.5%
    //-- Raw
    int i0p5pct = iterCut0p5pct[icent];
    FixUnfold(hUnfold0p5pct[icent][i0p5pct]);
    EbyECumu cumu0p5pct(hUnfold0p5pct[icent][i0p5pct]);
    double vn40p5pct  = cumu0p5pct.GetCumu_vn4();
    double vn60p5pct  = cumu0p5pct.GetCumu_vn6();
    double vn80p5pct  = cumu0p5pct.GetCumu_vn8();
    double g1ex0p5pct = cumu0p5pct.GetGamma1Exp();

    if(vn40p5pct == 0 || vn60p5pct == 0) vn6vn40p5pct[icent] = -1000.;
    else                                           vn6vn40p5pct[icent] = vn60p5pct / vn40p5pct;
    if(vn40p5pct == 0 || vn80p5pct == 0) vn8vn40p5pct[icent] = -1000.;
    else                                           vn8vn40p5pct[icent] = vn80p5pct / vn40p5pct;
    if(vn60p5pct == 0 || vn80p5pct == 0) vn8vn60p5pct[icent] = -1000.;
    else                                           vn8vn60p5pct[icent] = vn80p5pct / vn60p5pct;
    g1e0p5pct[icent] = g1ex0p5pct;

    double vn6vn40p5pct_staterr = sqrt( hVarianceOfMean_Vn6Vn4_0p5pct->GetBinContent(icent+1) );
    double vn8vn40p5pct_staterr = sqrt( hVarianceOfMean_Vn8Vn4_0p5pct->GetBinContent(icent+1) );
    double vn8vn60p5pct_staterr = sqrt( hVarianceOfMean_Vn8Vn6_0p5pct->GetBinContent(icent+1) );
    double g1e0p5pct_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_0p5pct->GetBinContent(icent+1) );

    vn6vn40p5pct_err[icent] = vn6vn40p5pct_staterr;
    vn8vn40p5pct_err[icent] = vn8vn40p5pct_staterr;
    vn8vn60p5pct_err[icent] = vn8vn60p5pct_staterr;
    g1e0p5pct_err[icent]    = g1e0p5pct_staterr;

    //-- Ratio to default
    if(vn40p5pct == 0 || vn60p5pct == 0 || vn4Default == 0 || vn6Default == 0) vn6vn40p5pct_RatioToDefault[icent] = -1000.;
    else                                                                                           vn6vn40p5pct_RatioToDefault[icent] = vn6vn40p5pct[icent] / vn6vn4Default[icent];
    if(vn40p5pct == 0 || vn80p5pct == 0 || vn4Default == 0 || vn8Default == 0) vn8vn40p5pct_RatioToDefault[icent] = -1000.;
    else                                                                                           vn8vn40p5pct_RatioToDefault[icent] = vn8vn40p5pct[icent] / vn8vn4Default[icent];
    if(vn60p5pct == 0 || vn80p5pct == 0 || vn6Default == 0 || vn8Default == 0) vn8vn60p5pct_RatioToDefault[icent] = -1000.;
    else                                                                                           vn8vn60p5pct_RatioToDefault[icent] = vn8vn60p5pct[icent] / vn8vn6Default[icent];
    if( g1ex0p5pct == -10000. || g1exDefault == -10000.) g1e0p5pct_RatioToDefault[icent] = -1000.;
    else                                                           g1e0p5pct_RatioToDefault[icent] = g1e0p5pct[icent] / g1eDefault[icent];

    vn6vn40p5pct_RatioToDefault_err[icent] = sqrt( pow(vn6vn40p5pct_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn40p5pct[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn40p5pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn40p5pct_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn40p5pct[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn60p5pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn60p5pct_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn60p5pct[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    g1e0p5pct_RatioToDefault_err[icent]    = sqrt( pow(g1e0p5pct_err[icent]/g1eDefault[icent],2) + pow(g1e0p5pct[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

    //-- 1%
    //-- Raw
    int i1pct = iterCut1pct[icent];
    FixUnfold(hUnfold1pct[icent][i1pct]);
    EbyECumu cumu1pct(hUnfold1pct[icent][i1pct]);
    double vn41pct  = cumu1pct.GetCumu_vn4();
    double vn61pct  = cumu1pct.GetCumu_vn6();
    double vn81pct  = cumu1pct.GetCumu_vn8();
    double g1ex1pct = cumu1pct.GetGamma1Exp();

    if(vn41pct == 0 || vn61pct == 0) vn6vn41pct[icent] = -1000.;
    else                                             vn6vn41pct[icent] = vn61pct / vn41pct;
    if(vn41pct == 0 || vn81pct == 0) vn8vn41pct[icent] = -1000.;
    else                                             vn8vn41pct[icent] = vn81pct / vn41pct;
    if(vn61pct == 0 || vn81pct == 0) vn8vn61pct[icent] = -1000.;
    else                                             vn8vn61pct[icent] = vn81pct / vn61pct;
    g1e1pct[icent] = g1ex1pct;

    double vn6vn41pct_staterr = sqrt( hVarianceOfMean_Vn6Vn4_1pct->GetBinContent(icent+1) );
    double vn8vn41pct_staterr = sqrt( hVarianceOfMean_Vn8Vn4_1pct->GetBinContent(icent+1) );
    double vn8vn61pct_staterr = sqrt( hVarianceOfMean_Vn8Vn6_1pct->GetBinContent(icent+1) );
    double g1e1pct_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_1pct->GetBinContent(icent+1) );

    vn6vn41pct_err[icent] = vn6vn41pct_staterr;
    vn8vn41pct_err[icent] = vn8vn41pct_staterr;
    vn8vn61pct_err[icent] = vn8vn61pct_staterr;
    g1e1pct_err[icent]    = g1e1pct_staterr;

    //-- Ratio to default
    if(vn41pct == 0 || vn61pct == 0 || vn4Default == 0 || vn6Default == 0) vn6vn41pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn6vn41pct_RatioToDefault[icent] = vn6vn41pct[icent] / vn6vn4Default[icent];
    if(vn41pct == 0 || vn81pct == 0 || vn4Default == 0 || vn8Default == 0) vn8vn41pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn41pct_RatioToDefault[icent] = vn8vn41pct[icent] / vn8vn4Default[icent];
    if(vn61pct == 0 || vn81pct == 0 || vn6Default == 0 || vn8Default == 0) vn8vn61pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn61pct_RatioToDefault[icent] = vn8vn61pct[icent] / vn8vn6Default[icent];
    if(g1ex1pct== -10000. || g1exDefault == -10000.) g1e1pct_RatioToDefault[icent] = -1000.;
    else                                                          g1e1pct_RatioToDefault[icent]= g1e1pct[icent]/ g1eDefault[icent];

    vn6vn41pct_RatioToDefault_err[icent] = sqrt( pow(vn6vn41pct_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn41pct[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn41pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn41pct_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn41pct[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn61pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn61pct_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn61pct[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    g1e1pct_RatioToDefault_err[icent]    = sqrt( pow(g1e1pct_err[icent]/g1eDefault[icent],2) + pow(g1e1pct[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

    //-- 2%
    //-- Raw
    int i2pct = iterCut2pct[icent];
    FixUnfold(hUnfold2pct[icent][i2pct]);
    EbyECumu cumu2pct(hUnfold2pct[icent][i2pct]);
    double vn42pct  = cumu2pct.GetCumu_vn4();
    double vn62pct  = cumu2pct.GetCumu_vn6();
    double vn82pct  = cumu2pct.GetCumu_vn8();
    double g1ex2pct = cumu2pct.GetGamma1Exp();

    if(vn42pct == 0 || vn62pct == 0) vn6vn42pct[icent] = -1000.;
    else                                             vn6vn42pct[icent] = vn62pct / vn42pct;
    if(vn42pct == 0 || vn82pct == 0) vn8vn42pct[icent] = -1000.;
    else                                             vn8vn42pct[icent] = vn82pct / vn42pct;
    if(vn62pct == 0 || vn82pct == 0) vn8vn62pct[icent] = -1000.;
    else                                             vn8vn62pct[icent] = vn82pct / vn62pct;
    g1e2pct[icent] = g1ex2pct;

    double vn6vn42pct_staterr = sqrt( hVarianceOfMean_Vn6Vn4_2pct->GetBinContent(icent+1) );
    double vn8vn42pct_staterr = sqrt( hVarianceOfMean_Vn8Vn4_2pct->GetBinContent(icent+1) );
    double vn8vn62pct_staterr = sqrt( hVarianceOfMean_Vn8Vn6_2pct->GetBinContent(icent+1) );
    double g1e2pct_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_2pct->GetBinContent(icent+1) );

    vn6vn42pct_err[icent] = vn6vn42pct_staterr;
    vn8vn42pct_err[icent] = vn8vn42pct_staterr;
    vn8vn62pct_err[icent] = vn8vn62pct_staterr;
    g1e2pct_err[icent]    = g1e2pct_staterr;

    //-- Ratio to default
    if(vn42pct == 0 || vn62pct == 0 || vn4Default == 0 || vn6Default == 0) vn6vn42pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn6vn42pct_RatioToDefault[icent] = vn6vn42pct[icent] / vn6vn4Default[icent];
    if(vn42pct == 0 || vn82pct == 0 || vn4Default == 0 || vn8Default == 0) vn8vn42pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn42pct_RatioToDefault[icent] = vn8vn42pct[icent] / vn8vn4Default[icent];
    if(vn62pct == 0 || vn82pct == 0 || vn6Default == 0 || vn8Default == 0) vn8vn62pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn62pct_RatioToDefault[icent] = vn8vn62pct[icent] / vn8vn6Default[icent];
    if(g1ex2pct== -10000. || g1exDefault == -10000.) g1e2pct_RatioToDefault[icent] = -1000.;
    else                                                          g1e2pct_RatioToDefault[icent]= g1e2pct[icent]/ g1eDefault[icent];

    vn6vn42pct_RatioToDefault_err[icent] = sqrt( pow(vn6vn42pct_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn42pct[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn42pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn42pct_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn42pct[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn62pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn62pct_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn62pct[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    g1e2pct_RatioToDefault_err[icent]    = sqrt( pow(g1e2pct_err[icent]/g1eDefault[icent],2) + pow(g1e2pct[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

    //-- 1%
    //-- Raw
    int i3pct = iterCut3pct[icent];
    FixUnfold(hUnfold3pct[icent][i3pct]);
    EbyECumu cumu3pct(hUnfold3pct[icent][i3pct]);
    double vn43pct  = cumu3pct.GetCumu_vn4();
    double vn63pct  = cumu3pct.GetCumu_vn6();
    double vn83pct  = cumu3pct.GetCumu_vn8();
    double g1ex3pct = cumu3pct.GetGamma1Exp();

    if(vn43pct == 0 || vn63pct == 0) vn6vn43pct[icent] = -1000.;
    else                                             vn6vn43pct[icent] = vn63pct / vn43pct;
    if(vn43pct == 0 || vn83pct == 0) vn8vn43pct[icent] = -1000.;
    else                                             vn8vn43pct[icent] = vn83pct / vn43pct;
    if(vn63pct == 0 || vn83pct == 0) vn8vn63pct[icent] = -1000.;
    else                                             vn8vn63pct[icent] = vn83pct / vn63pct;
    g1e3pct[icent] = g1ex3pct;

    double vn6vn43pct_staterr = sqrt( hVarianceOfMean_Vn6Vn4_3pct->GetBinContent(icent+1) );
    double vn8vn43pct_staterr = sqrt( hVarianceOfMean_Vn8Vn4_3pct->GetBinContent(icent+1) );
    double vn8vn63pct_staterr = sqrt( hVarianceOfMean_Vn8Vn6_3pct->GetBinContent(icent+1) );
    double g1e3pct_staterr    = sqrt( hVarianceOfMean_Gamma1Exp_3pct->GetBinContent(icent+1) );

    vn6vn43pct_err[icent] = vn6vn43pct_staterr;
    vn8vn43pct_err[icent] = vn8vn43pct_staterr;
    vn8vn63pct_err[icent] = vn8vn63pct_staterr;
    g1e3pct_err[icent]    = g1e3pct_staterr;

    //-- Ratio to default
    if(vn43pct == 0 || vn63pct == 0 || vn4Default == 0 || vn6Default == 0) vn6vn43pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn6vn43pct_RatioToDefault[icent] = vn6vn43pct[icent] / vn6vn4Default[icent];
    if(vn43pct == 0 || vn83pct == 0 || vn4Default == 0 || vn8Default == 0) vn8vn43pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn43pct_RatioToDefault[icent] = vn8vn43pct[icent] / vn8vn4Default[icent];
    if(vn63pct == 0 || vn83pct == 0 || vn6Default == 0 || vn8Default == 0) vn8vn63pct_RatioToDefault[icent] = -1000.;
    else                                                                                             vn8vn63pct_RatioToDefault[icent] = vn8vn63pct[icent] / vn8vn6Default[icent];
    if(g1ex3pct== -10000. || g1exDefault == -10000.) g1e3pct_RatioToDefault[icent] = -1000.;
    else                                                          g1e3pct_RatioToDefault[icent]= g1e3pct[icent]/ g1eDefault[icent];

    vn6vn43pct_RatioToDefault_err[icent] = sqrt( pow(vn6vn43pct_err[icent]/vn6vn4Default[icent],2) + pow(vn6vn43pct[icent]*vn6vn4Default_err[icent]/vn6vn4Default[icent]/vn6vn4Default[icent],2) );
    vn8vn43pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn43pct_err[icent]/vn8vn4Default[icent],2) + pow(vn8vn43pct[icent]*vn8vn4Default_err[icent]/vn8vn4Default[icent]/vn8vn4Default[icent],2) );
    vn8vn63pct_RatioToDefault_err[icent] = sqrt( pow(vn8vn63pct_err[icent]/vn8vn6Default[icent],2) + pow(vn8vn63pct[icent]*vn8vn6Default_err[icent]/vn8vn6Default[icent]/vn8vn6Default[icent],2) );
    g1e3pct_RatioToDefault_err[icent]    = sqrt( pow(g1e3pct_err[icent]/g1eDefault[icent],2) + pow(g1e3pct[icent]*g1eDefault_err[icent]/g1eDefault[icent]/g1eDefault[icent],2) );

  } //-- End cent loop

  //-- TGraph time!
  double c_err[NCENT];
  for(int i = 0; i < NCENT; i++) c_err[i] = 0.;

  //-- Default
  grVn6Vn4Default = new TGraphErrors(NCENT, centBinCenter, vn6vn4Default, c_err, vn6vn4Default_err);
  grVn8Vn4Default = new TGraphErrors(NCENT, centBinCenter, vn8vn4Default, c_err, vn8vn4Default_err);
  grVn8Vn6Default = new TGraphErrors(NCENT, centBinCenter, vn8vn6Default, c_err, vn8vn6Default_err);
  grG1EDefault    = new TGraphErrors(NCENT, centBinCenter, g1eDefault,    c_err, g1eDefault_err);

  formatGraph(grVn6Vn4Default, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         21, "grVn6Vn4Default");
  formatGraph(grVn8Vn4Default, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  34, "grVn8Vn4Default");
  formatGraph(grVn8Vn6Default, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 33, "grVn8Vn6Default");
  formatGraph(grG1EDefault,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         20, "grG1EDefault");

  //-- 0.5%
  grVn6Vn40p5pct = new TGraphErrors(NCENT, centBinCenter, vn6vn40p5pct, c_err, vn6vn40p5pct_err);
  grVn8Vn40p5pct = new TGraphErrors(NCENT, centBinCenter, vn8vn40p5pct, c_err, vn8vn40p5pct_err);
  grVn8Vn60p5pct = new TGraphErrors(NCENT, centBinCenter, vn8vn60p5pct, c_err, vn8vn60p5pct_err);
  grG1E0p5pct    = new TGraphErrors(NCENT, centBinCenter, g1e0p5pct,    c_err, g1e0p5pct_err);

  formatGraph(grVn6Vn40p5pct, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         26, "grVn6Vn40p5pct");
  formatGraph(grVn8Vn40p5pct, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  26, "grVn8Vn40p5pct");
  formatGraph(grVn8Vn60p5pct, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 26, "grVn8Vn60p5pct");
  formatGraph(grG1E0p5pct,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         26, "grG1E0p5pct");

  //-- Ratio to Default
  grVn6Vn40p5pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn40p5pct_RatioToDefault, c_err, vn6vn40p5pct_RatioToDefault_err);
  grVn8Vn40p5pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn40p5pct_RatioToDefault, c_err, vn8vn40p5pct_RatioToDefault_err);
  grVn8Vn60p5pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn60p5pct_RatioToDefault, c_err, vn8vn60p5pct_RatioToDefault_err);
  grG1E0p5pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1e0p5pct_RatioToDefault,    c_err, g1e0p5pct_RatioToDefault_err);

  formatGraph(grVn6Vn40p5pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         26, "grVn6Vn40p5pct_RatioToDefault");
  formatGraph(grVn8Vn40p5pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  26, "grVn8Vn40p5pct_RatioToDefault");
  formatGraph(grVn8Vn60p5pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 26, "grVn8Vn60p5pct_RatioToDefault");
  formatGraph(grG1E0p5pct_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         26, "grG1E0p5pct_RatioToDefault");

  //-- 1.0%
  grVn6Vn41pct = new TGraphErrors(NCENT, centBinCenter, vn6vn41pct, c_err, vn6vn41pct_err);
  grVn8Vn41pct = new TGraphErrors(NCENT, centBinCenter, vn8vn41pct, c_err, vn8vn41pct_err);
  grVn8Vn61pct = new TGraphErrors(NCENT, centBinCenter, vn8vn61pct, c_err, vn8vn61pct_err);
  grG1E1pct    = new TGraphErrors(NCENT, centBinCenter, g1e1pct,    c_err, g1e1pct_err);

  formatGraph(grVn6Vn41pct, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         25, "grVn6Vn41pct");
  formatGraph(grVn8Vn41pct, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  25, "grVn8Vn41pct");
  formatGraph(grVn8Vn61pct, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 25, "grVn8Vn61pct");
  formatGraph(grG1E1pct,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         25, "grG1E1pct");

  //-- Ratio to Default
  grVn6Vn41pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn41pct_RatioToDefault, c_err, vn6vn41pct_RatioToDefault_err);
  grVn8Vn41pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn41pct_RatioToDefault, c_err, vn8vn41pct_RatioToDefault_err);
  grVn8Vn61pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn61pct_RatioToDefault, c_err, vn8vn61pct_RatioToDefault_err);
  grG1E1pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1e1pct_RatioToDefault,    c_err, g1e1pct_RatioToDefault_err);

  formatGraph(grVn6Vn41pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         25, "grVn6Vn41pct_RatioToDefault");
  formatGraph(grVn8Vn41pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  25, "grVn8Vn41pct_RatioToDefault");
  formatGraph(grVn8Vn61pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 25, "grVn8Vn61pct_RatioToDefault");
  formatGraph(grG1E1pct_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         25, "grG1E1pct_RatioToDefault");

  //-- 2.0%
  grVn6Vn42pct = new TGraphErrors(NCENT, centBinCenter, vn6vn42pct, c_err, vn6vn42pct_err);
  grVn8Vn42pct = new TGraphErrors(NCENT, centBinCenter, vn8vn42pct, c_err, vn8vn42pct_err);
  grVn8Vn62pct = new TGraphErrors(NCENT, centBinCenter, vn8vn62pct, c_err, vn8vn62pct_err);
  grG1E2pct    = new TGraphErrors(NCENT, centBinCenter, g1e2pct,    c_err, g1e2pct_err);

  formatGraph(grVn6Vn42pct, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         24, "grVn6Vn42pct");
  formatGraph(grVn8Vn42pct, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  24, "grVn8Vn42pct");
  formatGraph(grVn8Vn62pct, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 24, "grVn8Vn62pct");
  formatGraph(grG1E2pct,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         24, "grG1E2pct");

  //-- Ratio to Default
  grVn6Vn42pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn42pct_RatioToDefault, c_err, vn6vn42pct_RatioToDefault_err);
  grVn8Vn42pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn42pct_RatioToDefault, c_err, vn8vn42pct_RatioToDefault_err);
  grVn8Vn62pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn62pct_RatioToDefault, c_err, vn8vn62pct_RatioToDefault_err);
  grG1E2pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1e2pct_RatioToDefault,    c_err, g1e2pct_RatioToDefault_err);

  formatGraph(grVn6Vn42pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         24, "grVn6Vn42pct_RatioToDefault");
  formatGraph(grVn8Vn42pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  24, "grVn8Vn42pct_RatioToDefault");
  formatGraph(grVn8Vn62pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 24, "grVn8Vn62pct_RatioToDefault");
  formatGraph(grG1E2pct_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         24, "grG1E2pct_RatioToDefault");

  //-- 3.0%
  grVn6Vn43pct = new TGraphErrors(NCENT, centBinCenter, vn6vn43pct, c_err, vn6vn43pct_err);
  grVn8Vn43pct = new TGraphErrors(NCENT, centBinCenter, vn8vn43pct, c_err, vn8vn43pct_err);
  grVn8Vn63pct = new TGraphErrors(NCENT, centBinCenter, vn8vn63pct, c_err, vn8vn63pct_err);
  grG1E3pct    = new TGraphErrors(NCENT, centBinCenter, g1e3pct,    c_err, g1e3pct_err);

  formatGraph(grVn6Vn43pct, "Centrality %", vn6vn4Min, vn6vn4Max, "v_{2}{6}/v_{2}{4}", 4,         32, "grVn6Vn43pct");
  formatGraph(grVn8Vn43pct, "Centrality %", vn8vn4Min, vn8vn4Max, "v_{2}{8}/v_{2}{4}", kGreen+2,  32, "grVn8Vn43pct");
  formatGraph(grVn8Vn63pct, "Centrality %", vn8vn6Min, vn8vn6Max, "v_{2}{8}/v_{2}{6}", kViolet-1, 32, "grVn8Vn63pct");
  formatGraph(grG1E3pct,    "Centrality %", g1eMin,    g1eMax,    "#gamma_{1}^{exp}",  2,         32, "grG1E3pct");

  //-- Ratio to Default
  grVn6Vn43pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn6vn43pct_RatioToDefault, c_err, vn6vn43pct_RatioToDefault_err);
  grVn8Vn43pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn43pct_RatioToDefault, c_err, vn8vn43pct_RatioToDefault_err);
  grVn8Vn63pct_RatioToDefault = new TGraphErrors(NCENT, centBinCenter, vn8vn63pct_RatioToDefault, c_err, vn8vn63pct_RatioToDefault_err);
  grG1E3pct_RatioToDefault    = new TGraphErrors(NCENT, centBinCenter, g1e3pct_RatioToDefault,    c_err, g1e3pct_RatioToDefault_err);

  formatGraph(grVn6Vn43pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{6}/v_{2}{4} Ratio", 4,         32, "grVn6Vn43pct_RatioToDefault");
  formatGraph(grVn8Vn43pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{4} Ratio", kGreen+2,  32, "grVn8Vn43pct_RatioToDefault");
  formatGraph(grVn8Vn63pct_RatioToDefault, "Centrality %", rMin,   rMax,   "v_{2}{8}/v_{2}{6} Ratio", kViolet-1, 32, "grVn8Vn63pct_RatioToDefault");
  formatGraph(grG1E3pct_RatioToDefault,    "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio",  2,         32, "grG1E3pct_RatioToDefault");

  //-- DRAW

  //-- cumuRatio
  TLegend * legvn6vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Default, "Default", "lp");
  legvn6vn4->AddEntry(grVn6Vn40p5pct,  "0.5%",    "lp");
  legvn6vn4->AddEntry(grVn6Vn41pct,    "1%",      "lp");
  legvn6vn4->AddEntry(grVn6Vn42pct,    "2%",      "lp");
  legvn6vn4->AddEntry(grVn6Vn43pct,    "3%",      "lp");

  TLegend * legvn8vn4 = new TLegend(0.4538, 0.7421, 0.94, 0.9172);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Default, "Default", "lp");
  legvn8vn4->AddEntry(grVn8Vn40p5pct,  "0.5%",    "lp");
  legvn8vn4->AddEntry(grVn8Vn41pct,    "1%",      "lp");
  legvn8vn4->AddEntry(grVn8Vn42pct,    "2%",      "lp");
  legvn8vn4->AddEntry(grVn8Vn43pct,    "3%",      "lp");

  TLegend * legvn8vn6 = new TLegend(0.4520, 0.1973, 0.94, 0.3703);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Default, "Default", "lp");
  legvn8vn6->AddEntry(grVn8Vn60p5pct,  "0.5%",    "lp");
  legvn8vn6->AddEntry(grVn8Vn61pct,    "1%",      "lp");
  legvn8vn6->AddEntry(grVn8Vn62pct,    "2%",      "lp");
  legvn8vn6->AddEntry(grVn8Vn63pct,    "3%",      "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);

  cCumuRatio->cd(1);
  grVn6Vn4Default->Draw("ap");
  grVn6Vn40p5pct->Draw("psame");
  grVn6Vn41pct->Draw("psame");
  grVn6Vn42pct->Draw("psame");
  grVn6Vn43pct->Draw("psame");
  legvn6vn4->Draw("same");

  cCumuRatio->cd(2);
  grVn8Vn4Default->Draw("ap");
  grVn8Vn40p5pct->Draw("psame");
  grVn8Vn41pct->Draw("psame");
  grVn8Vn42pct->Draw("psame");
  grVn8Vn43pct->Draw("psame");
  legvn8vn4->Draw("same");

  cCumuRatio->cd(3);
  grVn8Vn6Default->Draw("ap");
  grVn8Vn60p5pct->Draw("psame");
  grVn8Vn61pct->Draw("psame");
  grVn8Vn62pct->Draw("psame");
  grVn8Vn63pct->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- cumuRatio To Default
  TCanvas * cCumuRatio_RatToDef = new TCanvas("cCumuRatio_RatToDef", "cCumuRatio_RatToDef", 1500, 500);
  cCumuRatio_RatToDef->Divide(3,1);

  cCumuRatio_RatToDef->cd(1);
  grVn6Vn40p5pct_RatioToDefault->Draw("ap");
  grVn6Vn41pct_RatioToDefault->Draw("psame");
  grVn6Vn42pct_RatioToDefault->Draw("psame");
  grVn6Vn43pct_RatioToDefault->Draw("psame");

  cCumuRatio_RatToDef->cd(2);
  grVn8Vn40p5pct_RatioToDefault->Draw("ap");
  grVn8Vn41pct_RatioToDefault->Draw("psame");
  grVn8Vn42pct_RatioToDefault->Draw("psame");
  grVn8Vn43pct_RatioToDefault->Draw("psame");

  cCumuRatio_RatToDef->cd(3);
  grVn8Vn60p5pct_RatioToDefault->Draw("ap");
  grVn8Vn61pct_RatioToDefault->Draw("psame");
  grVn8Vn62pct_RatioToDefault->Draw("psame");
  grVn8Vn63pct_RatioToDefault->Draw("psame");
  cCumuRatio_RatToDef->SaveAs("cCumuRatio_RatToDef.pdf");

  //-- g1e
  TLegend * legG1E = new TLegend(0.4538, 0.7421, 0.9529, 0.9172);
  legG1E->SetBorderSize(0);
  legG1E->SetFillStyle(0);
  legG1E->AddEntry(grG1EDefault, "Default", "lp");
  legG1E->AddEntry(grG1E0p5pct,  "0.5%",    "lp");
  legG1E->AddEntry(grG1E1pct,    "1%",      "lp");
  legG1E->AddEntry(grG1E2pct,    "2%",      "lp");
  legG1E->AddEntry(grG1E3pct,    "3%",      "lp");

  TCanvas* cG1e = new TCanvas("cG1e", "cG1e", 1000, 500);
  cG1e->Divide(2,1);

  cG1e->cd(1);
  grG1EDefault->Draw("ap");
  grG1E0p5pct->Draw("psame");
  grG1E1pct->Draw("psame");
  grG1E2pct->Draw("psame");
  grG1E3pct->Draw("psame");
  legG1E->Draw("same");

  cG1e->cd(2);
  grG1E0p5pct_RatioToDefault->Draw("ap");
  grG1E1pct_RatioToDefault->Draw("psame");
  grG1E2pct_RatioToDefault->Draw("psame");
  grG1E3pct_RatioToDefault->Draw("psame");
  cG1e->SaveAs("cG1e.pdf");

  //-- Mult Distns
  TLegend * legCent = new TLegend(0.0642, 0.2027, 0.9935, 0.9420);
  legCent->SetFillStyle(0);
  legCent->SetBorderSize(0);
  legCent->SetNColumns(2);
  legCent->AddEntry(hMultDefault[0], Form("Cent %i-%i%s",cent_min[0], cent_max[0], "%"),"lp");

  TCanvas * cMult = new TCanvas("cMult", "cMult", 1500, 1000);
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

    legCent->AddEntry(hMultDefault[icent], Form("Cent %i-%i%s", cent_min[icent], cent_max[icent], "%"), "lp");

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

  cMult->cd(6);
  legCent->Draw("same");
  cMult->SaveAs("cMult.pdf");

}

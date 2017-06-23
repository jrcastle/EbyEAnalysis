#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace ebyese;

void sysChi2Cutoff(){

  int norder_     = 2;
  double tkEta    = 1.0;

  bool dosys     = 0;
  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double ratioGamma1ExpMin = 0.6;
  double ratioGamma1ExpMax = 1.4;
  double ratioVn6Vn4Min    = 0.9;
  double ratioVn6Vn4Max    = 1.1;
  double ratioVn8Vn4Min    = 0.9;
  double ratioVn8Vn4Max    = 1.1;
  double ratioVn8Vn6Min    = 0.9;
  double ratioVn8Vn6Max    = 1.1;

  double chi2Min = 0.1;
  double chi2Max = 10000;

  double ratioMin = 0.98;
  double ratioMax = 1.02;

  TLatex latex;

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Unfolding output
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  //-- Chi2 Refold vs iter
  TGraphErrors * grChi2NDFvsIter[NCENT];
  double chi2NDF[NCENT][NITER];

  //-- Stat Uncert Output
  TFile * fStat;
  TH1D * hVarianceOfMean_Vn2[NCHI2];
  TH1D * hVarianceOfMean_Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn6[NCHI2];
  TH1D * hVarianceOfMean_Vn8[NCHI2];
  TH1D * hVarianceOfMean_Gamma1Exp[NCHI2];
  TH1D * hVarianceOfMean_Vn6Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn8Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn8Vn6[NCHI2];
  TH1D * hVarianceOfMean_Vn46_Vn68[NCHI2];

  //-- Cumu vs chi2 Cutoff
  double unfold_Vn2VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn4VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn6VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn8VSChi2Cut[NCENT][NCHI2];
  double unfold_Gamma1ExpVSChi2Cut[NCENT][NCHI2];
  double unfold_Vn6Vn4VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn8Vn4VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn8Vn6VSChi2Cut[NCENT][NCHI2];
  double unfold_Vn46_Vn68VSChi2Cut[NCENT][NCHI2];

  double unfold_Vn2VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn4VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn6VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn8VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Gamma1ExpVSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn6Vn4VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn8Vn4VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn8Vn6VSChi2Cut_staterr[NCENT][NCHI2];
  double unfold_Vn46_Vn68VSChi2Cut_staterr[NCENT][NCHI2];

  //-- Ratio to the Chi2=1 cutoff scenario
  double unfoldVn2_RatioToChi21[NCENT][NCHI2];
  double unfoldVn4_RatioToChi21[NCENT][NCHI2];
  double unfoldVn6_RatioToChi21[NCENT][NCHI2];
  double unfoldVn8_RatioToChi21[NCENT][NCHI2];
  double unfoldGamma1Exp_RatioToChi21[NCENT][NCHI2];
  double unfoldVn6Vn4_RatioToChi21[NCENT][NCHI2];
  double unfoldVn8Vn4_RatioToChi21[NCENT][NCHI2];
  double unfoldVn8Vn6_RatioToChi21[NCENT][NCHI2];
  double unfoldVn46_Vn68_RatioToChi21[NCENT][NCHI2];

  double unfoldVn2_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn4_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn6_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn8_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldGamma1Exp_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn6Vn4_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn8Vn4_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn8Vn6_RatioToChi21_staterr[NCENT][NCHI2];
  double unfoldVn46_Vn68_RatioToChi21_staterr[NCENT][NCHI2];

  bool iterCutoff[NCHI2];

  //-- Ratio chi2 = 2 / chi2 = 1 vs cent
  double ratioChi2_21_vn2[NCENT];
  double ratioChi2_21_vn4[NCENT];
  double ratioChi2_21_vn6[NCENT];
  double ratioChi2_21_vn8[NCENT];
  double ratioChi2_21_gamma1exp[NCENT];
  double ratioChi2_21_vn6vn4[NCENT];
  double ratioChi2_21_vn8vn4[NCENT];
  double ratioChi2_21_vn8vn6[NCENT];
  double ratioChi2_21_vn46_vn68[NCENT];

  double ratioChi2_21_vn2_staterr[NCENT];
  double ratioChi2_21_vn4_staterr[NCENT];
  double ratioChi2_21_vn6_staterr[NCENT];
  double ratioChi2_21_vn8_staterr[NCENT];
  double ratioChi2_21_gamma1exp_staterr[NCENT];
  double ratioChi2_21_vn6vn4_staterr[NCENT];
  double ratioChi2_21_vn8vn4_staterr[NCENT];
  double ratioChi2_21_vn8vn6_staterr[NCENT];
  double ratioChi2_21_vn46_vn68_staterr[NCENT];

  TGraphErrors * grRatioChi2_21_vn2;
  TGraphErrors * grRatioChi2_21_vn4;
  TGraphErrors * grRatioChi2_21_vn6;
  TGraphErrors * grRatioChi2_21_vn8;
  TGraphErrors * grRatioChi2_21_gamma1exp;
  TGraphErrors * grRatioChi2_21_vn6vn4;
  TGraphErrors * grRatioChi2_21_vn8vn4;
  TGraphErrors * grRatioChi2_21_vn8vn6;
  TGraphErrors * grRatioChi2_21_vn46_vn68;

  int finalIteration[NCENT][NCHI2];

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  //gErrorIgnoreLevel = kWarning;

  //-- Get the Analyzer output file
  fAna = new TFile( "../../AnalyzerResults/CastleEbyE.root" );

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
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/gaussResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/gaussResp/data%i_dosys.root", norder_) );
  }
  if( studTResp ){
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/studTResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/studTResp/data%i_dosys.root", norder_) );
  }
  if( dataResp ){
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/dataResp/data%i_dosys.root", norder_) );
  }

  //-- Grab the statistical errors
  fStat = new TFile( Form("../../../../statErrorHandle/v%i/eta%.1f/systematicStudies/chi2Cutoff/StatUncertChi2Cutoff_v%i.root", norder_, tkEta, norder_) );
  for(int ichi = 0; ichi < NCHI2; ichi++){
    hVarianceOfMean_Vn2[ichi]       = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn2_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn4[ichi]       = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn4_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn6[ichi]       = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn6_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn8[ichi]       = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn8_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Gamma1Exp[ichi] = (TH1D*) fStat->Get( Form("hVarianceOfMean_Gamma1Exp_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn6Vn4[ichi]    = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn6Vn4_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn8Vn4[ichi]    = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn8Vn4_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn8Vn6[ichi]    = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn8Vn6_chi2Cut%.1f", chi2Cutoff[ichi]) );
    hVarianceOfMean_Vn46_Vn68[ichi] = (TH1D*) fStat->Get( Form("hVarianceOfMean_Vn46_Vn68_chi2Cut%.1f", chi2Cutoff[ichi]) );

    for(int icent = 0; icent < NCENT; icent++){
      unfold_Vn2VSChi2Cut_staterr[icent][ichi]       = sqrt( hVarianceOfMean_Vn2[ichi]->GetBinContent(icent+1) );
      unfold_Vn4VSChi2Cut_staterr[icent][ichi]       = sqrt( hVarianceOfMean_Vn4[ichi]->GetBinContent(icent+1) );
      unfold_Vn6VSChi2Cut_staterr[icent][ichi]       = sqrt( hVarianceOfMean_Vn6[ichi]->GetBinContent(icent+1) );
      unfold_Vn8VSChi2Cut_staterr[icent][ichi]       = sqrt( hVarianceOfMean_Vn8[ichi]->GetBinContent(icent+1) );
      unfold_Gamma1ExpVSChi2Cut_staterr[icent][ichi] = sqrt( hVarianceOfMean_Gamma1Exp[ichi]->GetBinContent(icent+1) );
      unfold_Vn6Vn4VSChi2Cut_staterr[icent][ichi]    = sqrt( hVarianceOfMean_Vn6Vn4[ichi]->GetBinContent(icent+1) );
      unfold_Vn8Vn4VSChi2Cut_staterr[icent][ichi]    = sqrt( hVarianceOfMean_Vn8Vn4[ichi]->GetBinContent(icent+1) );
      unfold_Vn8Vn6VSChi2Cut_staterr[icent][ichi]    = sqrt( hVarianceOfMean_Vn8Vn6[ichi]->GetBinContent(icent+1) );
      unfold_Vn46_Vn68VSChi2Cut_staterr[icent][ichi] = sqrt( hVarianceOfMean_Vn46_Vn68[ichi]->GetBinContent(icent+1) );
    }
  }

  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    //-- Reset the flags that say whether or no the chi2 cutoff has been reached
    for(int ichi = 0; ichi < NCHI2; ichi++) iterCutoff[ichi] = 0;

    //-- Loop over iterations
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

      //-- Calculate the cumulants of each unfolded histogram
      FixUnfold( hUnfold[icent][i] );
      EbyECumu cumu(hUnfold[icent][i]);
      double vn2       = cumu.GetCumu_vn2();
      double vn4       = cumu.GetCumu_vn4();
      double vn6       = cumu.GetCumu_vn6();
      double vn8       = cumu.GetCumu_vn8();
      double gamma1exp = cumu.GetGamma1Exp();
      double vn6vn4;
      if(vn4 == 0) vn6vn4 = 0;
      else         vn6vn4 = vn6 / vn4;
      double vn8vn4;
      if(vn4 == 0) vn8vn4 = 0;
      else         vn8vn4 = vn8 / vn4;
      double vn8vn6;
      if(vn6 == 0) vn8vn6 = 0;
      else         vn8vn6 = vn8 / vn6;
      double vn46_vn68;
      if(vn4 == 0 || vn6 == 0 || vn8 == 0) vn46_vn68 = 0;
      else                                 vn46_vn68 = (vn4 - vn6)/(vn6 - vn8);

      //-- Calculate the refold chi2
      double chi2NDF_Refold = hRefold[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF");
      chi2NDF[icent][i] = chi2NDF_Refold;

      //-- loop over the scenarios and choose the iteration for each
      for(int ichi = 0; ichi < NCHI2; ichi++){

	if( chi2NDF_Refold <= chi2Cutoff[ichi] && !iterCutoff[ichi] ){
	  unfold_Vn2VSChi2Cut[icent][ichi]         = vn2;
	  unfold_Vn4VSChi2Cut[icent][ichi]         = vn4;
	  unfold_Vn6VSChi2Cut[icent][ichi]         = vn6;
	  unfold_Vn8VSChi2Cut[icent][ichi]         = vn8;
	  unfold_Gamma1ExpVSChi2Cut[icent][ichi]   = gamma1exp;
	  unfold_Vn6Vn4VSChi2Cut[icent][ichi]      = vn6vn4;
	  unfold_Vn8Vn4VSChi2Cut[icent][ichi]      = vn8vn4;
	  unfold_Vn8Vn6VSChi2Cut[icent][ichi]      = vn8vn6;
	  unfold_Vn46_Vn68VSChi2Cut[icent][ichi]   = vn46_vn68;

	  iterCutoff[ichi] = 1;
	  finalIteration[icent][ichi] = i;
	}
	if( i == NITER-1 && !iterCutoff[ichi] ){
	  unfold_Vn2VSChi2Cut[icent][ichi]         = vn2;
          unfold_Vn4VSChi2Cut[icent][ichi]         = vn4;
          unfold_Vn6VSChi2Cut[icent][ichi]         = vn6;
	  unfold_Vn8VSChi2Cut[icent][ichi]         = vn8;
          unfold_Gamma1ExpVSChi2Cut[icent][ichi]   = gamma1exp;
          unfold_Vn6Vn4VSChi2Cut[icent][ichi]      = vn6vn4;
	  unfold_Vn8Vn4VSChi2Cut[icent][ichi]      = vn8vn4;
          unfold_Vn8Vn6VSChi2Cut[icent][ichi]      = vn8vn6;
	  unfold_Vn46_Vn68VSChi2Cut[icent][ichi]   = vn46_vn68;

	  iterCutoff[ichi] = 1;
	  finalIteration[icent][ichi] = i;
	}

      } //-- End chi2 scenario loop

    } //-- End unfold iteration loop

    //-- Ratio to the Chi2=1 cutoff scenario
    for(int ichi = 0; ichi < NCHI2; ichi++){

      double vn2_i  = unfold_Vn2VSChi2Cut[icent][ichi];
      double vn2_0  = unfold_Vn2VSChi2Cut[icent][0];
      double vn4_i  = unfold_Vn4VSChi2Cut[icent][ichi];
      double vn4_0  = unfold_Vn4VSChi2Cut[icent][0];
      double vn6_i  = unfold_Vn6VSChi2Cut[icent][ichi];
      double vn6_0  = unfold_Vn6VSChi2Cut[icent][0];
      double vn8_i  = unfold_Vn8VSChi2Cut[icent][ichi];
      double vn8_0  = unfold_Vn8VSChi2Cut[icent][0];

      double gamma1exp_i = unfold_Gamma1ExpVSChi2Cut[icent][ichi];
      double gamma1exp_0 = unfold_Gamma1ExpVSChi2Cut[icent][0];
      double vn6vn4_i    = unfold_Vn6Vn4VSChi2Cut[icent][ichi];
      double vn6vn4_0    = unfold_Vn6Vn4VSChi2Cut[icent][0];
      double vn8vn4_i    = unfold_Vn8Vn4VSChi2Cut[icent][ichi];
      double vn8vn4_0    = unfold_Vn8Vn4VSChi2Cut[icent][0];
      double vn8vn6_i    = unfold_Vn8Vn6VSChi2Cut[icent][ichi];
      double vn8vn6_0    = unfold_Vn8Vn6VSChi2Cut[icent][0];
      double vn46_vn68_i = unfold_Vn46_Vn68VSChi2Cut[icent][ichi];
      double vn46_vn68_0 = unfold_Vn46_Vn68VSChi2Cut[icent][0];

      double vn2Rat       = vn2_i / vn2_0;
      double vn4Rat       = vn4_i / vn4_0;
      double vn6Rat       = vn6_i / vn6_0;
      double vn8Rat       = vn8_i / vn8_0;
      double g1expRat     = gamma1exp_i / gamma1exp_0;
      double vn6vn4Rat    = vn6vn4_i / vn6vn4_0;
      double vn8vn4Rat    = vn8vn4_i / vn8vn4_0;
      double vn8vn6Rat    = vn8vn6_i / vn8vn6_0;
      double vn46_vn68Rat = vn46_vn68_i / vn46_vn68_0;

      double vn2_i_se  = unfold_Vn2VSChi2Cut_staterr[icent][ichi];
      double vn2_0_se  = unfold_Vn2VSChi2Cut_staterr[icent][0];
      double vn4_i_se  = unfold_Vn4VSChi2Cut_staterr[icent][ichi];
      double vn4_0_se  = unfold_Vn4VSChi2Cut_staterr[icent][0];
      double vn6_i_se  = unfold_Vn6VSChi2Cut_staterr[icent][ichi];
      double vn6_0_se  = unfold_Vn6VSChi2Cut_staterr[icent][0];
      double vn8_i_se  = unfold_Vn8VSChi2Cut_staterr[icent][ichi];
      double vn8_0_se  = unfold_Vn8VSChi2Cut_staterr[icent][0];

      double gamma1exp_i_se = unfold_Gamma1ExpVSChi2Cut_staterr[icent][ichi];
      double gamma1exp_0_se = unfold_Gamma1ExpVSChi2Cut_staterr[icent][0];
      double vn6vn4_i_se    = unfold_Vn6Vn4VSChi2Cut_staterr[icent][ichi];
      double vn6vn4_0_se    = unfold_Vn6Vn4VSChi2Cut_staterr[icent][0];
      double vn8vn4_i_se    = unfold_Vn8Vn4VSChi2Cut_staterr[icent][ichi];
      double vn8vn4_0_se    = unfold_Vn8Vn4VSChi2Cut_staterr[icent][0];
      double vn8vn6_i_se    = unfold_Vn8Vn6VSChi2Cut_staterr[icent][ichi];
      double vn8vn6_0_se    = unfold_Vn8Vn6VSChi2Cut_staterr[icent][0];
      double vn46_vn68_i_se = unfold_Vn46_Vn68VSChi2Cut_staterr[icent][ichi];
      double vn46_vn68_0_se = unfold_Vn46_Vn68VSChi2Cut_staterr[icent][0];

      double vn2Rat_staterr       = sqrt( pow(vn2_i_se/vn2_0, 2) + pow( vn2_i*vn2_0_se/vn2_0/vn2_0, 2) );
      double vn4Rat_staterr       = sqrt( pow(vn4_i_se/vn4_0, 2) + pow( vn4_i*vn4_0_se/vn4_0/vn4_0, 2) );
      double vn6Rat_staterr       = sqrt( pow(vn6_i_se/vn6_0, 2) + pow( vn6_i*vn6_0_se/vn6_0/vn6_0, 2) );
      double vn8Rat_staterr       = sqrt( pow(vn8_i_se/vn8_0, 2) + pow( vn8_i*vn8_0_se/vn8_0/vn8_0, 2) );
      double g1expRat_staterr     = sqrt( pow(gamma1exp_i_se/gamma1exp_0, 2) + pow( gamma1exp_i*gamma1exp_0_se/gamma1exp_0/gamma1exp_0, 2) );
      double vn6vn4Rat_staterr    = sqrt( pow(vn6vn4_i_se/vn6vn4_0, 2) + pow( vn6vn4_i*vn6vn4_0_se/vn6vn4_0/vn6vn4_0, 2) );
      double vn8vn4Rat_staterr    = sqrt( pow(vn8vn4_i_se/vn8vn4_0, 2) + pow( vn8vn4_i*vn8vn4_0_se/vn8vn4_0/vn8vn4_0, 2) );
      double vn8vn6Rat_staterr    = sqrt( pow(vn8vn6_i_se/vn8vn6_0, 2) + pow( vn8vn6_i*vn8vn6_0_se/vn8vn6_0/vn8vn6_0, 2) );
      double vn46_vn68Rat_staterr = sqrt( pow(vn46_vn68_i_se/vn46_vn68_0, 2) + pow( vn46_vn68_i*vn46_vn68_0_se/vn46_vn68_0/vn46_vn68_0, 2) );

      unfoldVn2_RatioToChi21[icent][ichi]        = vn2Rat;
      unfoldVn4_RatioToChi21[icent][ichi]        = vn4Rat;
      unfoldVn6_RatioToChi21[icent][ichi]        = vn6Rat;
      unfoldVn8_RatioToChi21[icent][ichi]        = vn8Rat;
      unfoldGamma1Exp_RatioToChi21[icent][ichi]  = g1expRat;
      unfoldVn6Vn4_RatioToChi21[icent][ichi]     = vn6vn4Rat;
      unfoldVn8Vn4_RatioToChi21[icent][ichi]     = vn8vn4Rat;
      unfoldVn8Vn6_RatioToChi21[icent][ichi]     = vn8vn6Rat;
      unfoldVn46_Vn68_RatioToChi21[icent][ichi]  = vn46_vn68Rat;

      unfoldVn2_RatioToChi21_staterr[icent][ichi]       = vn2Rat_staterr;
      unfoldVn4_RatioToChi21_staterr[icent][ichi]       = vn4Rat_staterr;
      unfoldVn6_RatioToChi21_staterr[icent][ichi]       = vn6Rat_staterr;
      unfoldVn8_RatioToChi21_staterr[icent][ichi]       = vn8Rat_staterr;
      unfoldGamma1Exp_RatioToChi21_staterr[icent][ichi] = g1expRat_staterr;
      unfoldVn6Vn4_RatioToChi21_staterr[icent][ichi]    = vn6vn4Rat_staterr;
      unfoldVn8Vn4_RatioToChi21_staterr[icent][ichi]    = vn8vn4Rat_staterr;
      unfoldVn8Vn6_RatioToChi21_staterr[icent][ichi]    = vn8vn6Rat_staterr;
      unfoldVn46_Vn68_RatioToChi21_staterr[icent][ichi] = vn46_vn68Rat_staterr;

    } //-- End chi2 loop for ratios


    //-- Set up all the wonderful graphs
    grChi2NDFvsIter[icent] = new TGraphErrors(NITER, diter, chi2NDF[icent], iterErr, iterErr); 
    formatGraph(grChi2NDFvsIter[icent], "Iteration", chi2Min, chi2Max, "#chi^{2}/NDF", 1, 20, Form("grChi2NDFvsIter_c%i", icent) );

  } //-- End cent loop
  //gErrorIgnoreLevel = kError;

  //----------------------------------------------------------------------------------------------------
  //-- Chi2/NDF big plot

  TLine * l1 = new TLine(diter[0], 1.2, diter[NITER-1], 1.2);
  l1->SetLineColor(1);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);

  TCanvas * cChi2NDF_Big   = new TCanvas( "cChi2NDF_Big_%s", "cChi2NDF_Big_%s", 1500, 2000);
  cChi2NDF_Big->Divide(3,4);

  for(int icent = 0; icent < NCENT; icent++){
    cChi2NDF_Big->cd(icent+1);
    cChi2NDF_Big->cd(icent+1)->SetLogx();
    grChi2NDFvsIter[icent]->Draw("alp");
    formatGraph(grChi2NDFvsIter[icent], "Iteration", chi2Min, chi2Max, "#chi^{2}/NDF", 1, 20, Form("grChi2NDFvsIter_c%i", icent) );
    cChi2NDF_Big->cd(icent+1)->SetLogy();
    latex.DrawLatex(0.18, 0.88, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    l1->Draw("same");
  }

  cChi2NDF_Big->Update();
  cChi2NDF_Big->SaveAs( "../../plots/systematicStudies/cChi2NDF_Big.pdf" );


  //----------------------------------------------------------------------------------------------------
  //-- Ratio chi2 = 2 / chi2 = 1 vs cent
  double cErr[NCENT];
  for(int icent = 0; icent < NCENT; icent++){
    ratioChi2_21_vn2[icent]        = unfoldVn2_RatioToChi21[icent][3];
    ratioChi2_21_vn4[icent]        = unfoldVn4_RatioToChi21[icent][3];
    ratioChi2_21_vn6[icent]        = unfoldVn6_RatioToChi21[icent][3];
    ratioChi2_21_vn8[icent]        = unfoldVn8_RatioToChi21[icent][3];
    ratioChi2_21_gamma1exp[icent]  = unfoldGamma1Exp_RatioToChi21[icent][3];
    ratioChi2_21_vn6vn4[icent]     = unfoldVn6Vn4_RatioToChi21[icent][3];
    ratioChi2_21_vn8vn4[icent]     = unfoldVn8Vn4_RatioToChi21[icent][3];
    ratioChi2_21_vn8vn6[icent]     = unfoldVn8Vn6_RatioToChi21[icent][3];
    ratioChi2_21_vn46_vn68[icent]  = unfoldVn46_Vn68_RatioToChi21[icent][3];

    std::cout<<ratioChi2_21_gamma1exp[icent]<<std::endl;

    ratioChi2_21_vn2_staterr[icent]        = unfoldVn2_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn4_staterr[icent]        = unfoldVn4_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn6_staterr[icent]        = unfoldVn6_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn8_staterr[icent]        = unfoldVn8_RatioToChi21_staterr[icent][3];
    ratioChi2_21_gamma1exp_staterr[icent]  = unfoldGamma1Exp_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn6vn4_staterr[icent]     = unfoldVn6Vn4_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn8vn4_staterr[icent]     = unfoldVn8Vn4_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn8vn6_staterr[icent]     = unfoldVn8Vn6_RatioToChi21_staterr[icent][3];
    ratioChi2_21_vn46_vn68_staterr[icent]  = unfoldVn46_Vn68_RatioToChi21_staterr[icent][3];

    cErr[icent] = 0;
  }

  grRatioChi2_21_vn2 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn2, cErr, ratioChi2_21_vn2_staterr);
  formatGraph(grRatioChi2_21_vn2, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{2} [#chi^{2}=2]/[#chi^{2}=1]", norder_), 1, 24, "grRatioChi2_21_vn2");
  grRatioChi2_21_vn4 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn4, cErr, ratioChi2_21_vn4_staterr);
  formatGraph(grRatioChi2_21_vn4, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{4} [#chi^{2}=2]/[#chi^{2}=1]", norder_), kSpring+4, 25, "grRatioChi2_21_vn4");
  grRatioChi2_21_vn6 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn6, cErr, ratioChi2_21_vn6_staterr);
  formatGraph(grRatioChi2_21_vn6, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6} [#chi^{2}=2]/[#chi^{2}=1]", norder_), 6, 28, "grRatioChi2_21_vn6");
  grRatioChi2_21_vn8 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn8, cErr, ratioChi2_21_vn8_staterr);
  formatGraph(grRatioChi2_21_vn8, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8} [#chi^{2}=2]/[#chi^{2}=1]", norder_), kOrange+7, 27, "grRatioChi2_21_vn8");

  grRatioChi2_21_gamma1exp = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_gamma1exp, cErr, ratioChi2_21_gamma1exp_staterr);
  formatGraph(grRatioChi2_21_gamma1exp, "Centrality %", ratioGamma1ExpMin, ratioGamma1ExpMax, "#gamma_{1}^{exp} [#chi^{2}=2]/[#chi^{2}=1]", 2, 20, "grRatioChi2_21_gamma1exp");
  grRatioChi2_21_vn6vn4 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn6vn4, cErr, ratioChi2_21_vn6vn4_staterr);
  formatGraph(grRatioChi2_21_vn6vn4, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{6}/v_{%i}{4} [#chi^{2}=2]/[#chi^{2}=1]", norder_, norder_), 4, 21, "grRatioChi2_21_vn6vn4");
  grRatioChi2_21_vn8vn4 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn8vn4, cErr, ratioChi2_21_vn8vn4_staterr);
  formatGraph(grRatioChi2_21_vn8vn4, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8}/v_{%i}{4} [#chi^{2}=2]/[#chi^{2}=1]", norder_, norder_), kGreen+2, 34, "grRatioChi2_21_vn8vn4");
  grRatioChi2_21_vn8vn6 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn8vn6, cErr, ratioChi2_21_vn8vn6_staterr);
  formatGraph(grRatioChi2_21_vn8vn6, "Centrality %", ratioMin, ratioMax, Form("v_{%i}{8}/v_{%i}{6} [#chi^{2}=2]/[#chi^{2}=1]", norder_, norder_), kViolet-1, 33, "grRatioChi2_21_vn8vn6");
  grRatioChi2_21_vn46_vn68 = new TGraphErrors(NCENT, centBinCenter, ratioChi2_21_vn46_vn68, cErr, ratioChi2_21_vn46_vn68_staterr);
  formatGraph(grRatioChi2_21_vn46_vn68, "Centrality %", 0.80, 1.20, Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) [#chi^{2}=2]/[#chi^{2}=1]", norder_, norder_, norder_, norder_), kGray+2, 22, "grRatioChi2_21_vn46_vn68");

  TLegend * leg2 = new TLegend(0.1946, 0.1995, 0.3452, 0.4157);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(grRatioChi2_21_vn2, Form("v_{%i}{%i}", norder_, 2), "lp");
  leg2->AddEntry(grRatioChi2_21_vn4, Form("v_{%i}{%i}", norder_, 4), "lp");
  leg2->AddEntry(grRatioChi2_21_vn6, Form("v_{%i}{%i}", norder_, 6), "lp");

  TLine * line = new TLine(grRatioChi2_21_vn2->GetXaxis()->GetXmin(), 1.0, grRatioChi2_21_vn2->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * cCumuCentDep = new TCanvas("cCumuCentDep","cCumuCentDep", 1000, 1000);
  cCumuCentDep->Divide(2,2);
  
  double mar = 0.2;
  double offs = 1.6;

  grRatioChi2_21_vn2->GetYaxis()->SetTitleOffset(offs);
  grRatioChi2_21_vn4->GetYaxis()->SetTitleOffset(offs);
  grRatioChi2_21_vn6->GetYaxis()->SetTitleOffset(offs);
  grRatioChi2_21_vn8->GetYaxis()->SetTitleOffset(offs);

  cCumuCentDep->cd(1)->SetLeftMargin(mar);
  grRatioChi2_21_vn2->Draw("ap");
  line->Draw("same");
  cCumuCentDep->cd(2)->SetLeftMargin(mar);
  grRatioChi2_21_vn4->Draw("ap");
  line->Draw("same");
  cCumuCentDep->cd(3)->SetLeftMargin(mar);
  grRatioChi2_21_vn6->Draw("ap");
  line->Draw("same");
  cCumuCentDep->cd(4)->SetLeftMargin(mar);
  grRatioChi2_21_vn8->Draw("ap");
  line->Draw("same");
  cCumuCentDep->SaveAs("../../plots/systematicStudies/cSysChi2Cut_CumuCent.pdf" );

  TCanvas * cGamma1ExpCentDep = new TCanvas("cGamma1ExpCentDep","cGamma1ExpCentDep", 500, 500);
  cGamma1ExpCentDep->cd();
  grRatioChi2_21_gamma1exp->Draw("ap");
  line->Draw("same");
  cGamma1ExpCentDep->SaveAs("../../plots/systematicStudies/cSysChi2Cut_Gamma1ExpCent.pdf" );

  TCanvas * cvnCumuRatioCentDep = new TCanvas("cvnCumuRatioCentDep","cvnCumuRatioCentDep", 1500, 500);
  cvnCumuRatioCentDep->Divide(3,1);

  grRatioChi2_21_vn6vn4->GetYaxis()->SetTitleOffset(offs);
  grRatioChi2_21_vn8vn4->GetYaxis()->SetTitleOffset(offs);
  grRatioChi2_21_vn8vn6->GetYaxis()->SetTitleOffset(offs);

  cvnCumuRatioCentDep->cd(1)->SetLeftMargin(mar);
  grRatioChi2_21_vn6vn4->Draw("ap");
  line->Draw("same");
  cvnCumuRatioCentDep->cd(2)->SetLeftMargin(mar);
  grRatioChi2_21_vn8vn4->Draw("ap");
  line->Draw("same");
  cvnCumuRatioCentDep->cd(3)->SetLeftMargin(mar);
  grRatioChi2_21_vn8vn6->Draw("ap");
  line->Draw("same");
  cvnCumuRatioCentDep->SaveAs("../../plots/systematicStudies/cSysChi2Cut_CumuRatioCent.pdf" );

  TCanvas * cVn46_Vn68 = new TCanvas("cVn46_Vn68", "cVn46_Vn68", 500, 500);
  cVn46_Vn68->cd();
  grRatioChi2_21_vn46_vn68->Draw("ap");
  line->Draw("same");
  cVn46_Vn68->SaveAs("../../plots/systematicStudies/cSysChi2Cut_Vn46_Vn68.pdf");

  //---------------------------------------------------------------------------------------------------- 
  //-- Save plots for smoothing
  TFile * fSave = new TFile("SysReg.root", "recreate");
  fSave->cd();
  grRatioChi2_21_vn2->Write("grRatioChi2_21_vn2");
  grRatioChi2_21_vn4->Write("grRatioChi2_21_vn4");
  grRatioChi2_21_vn6->Write("grRatioChi2_21_vn6");
  grRatioChi2_21_vn8->Write("grRatioChi2_21_vn8");
  grRatioChi2_21_vn6vn4->Write("grRatioChi2_21_vn6vn4");
  grRatioChi2_21_vn8vn4->Write("grRatioChi2_21_vn8vn4");
  grRatioChi2_21_vn8vn6->Write("grRatioChi2_21_vn8vn6");
  grRatioChi2_21_vn46_vn68->Write("grRatioChi2_21_vn46_vn68");
  grRatioChi2_21_gamma1exp->Write("grRatioChi2_21_gamma1exp");

  //-- Save the unfolded distns for when the cutoff is chi2=2.
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIteration[icent][3];
    std::cout<<i<<std::endl;
    hUnfold[icent][i]->SetLineColor(1);
    hUnfold[icent][i]->SetMarkerColor(1);
    hUnfold[icent][i]->Write( Form("hFinalUnfold_SysReg_c%i", icent) );
  }


  
}

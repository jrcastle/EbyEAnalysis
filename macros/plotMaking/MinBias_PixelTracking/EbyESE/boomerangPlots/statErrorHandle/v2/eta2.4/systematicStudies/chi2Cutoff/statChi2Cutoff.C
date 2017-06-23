#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void statChi2Cutoff(){

  int norder_ = 2;

  //-- Analyzer Output
  TFile * fAna[NSPLIT];
  TH1D * hObs[NCENT][NSPLIT];

  //-- Unfolding output
  TFile * fUnf[NSPLIT];
  TH1D * hUnfold[NCENT][NSPLIT][NITER];
  TH1D * hRefold[NCENT][NSPLIT][NITER];

  bool iterCutoff[NCHI2][NCENT][NSPLIT];
  bool cutoffFound[NCHI2];

  //-- Stat Uncert Output
  TFile * fOut;
  TH1D * hVarianceOfMean_Vn2[NCHI2];
  TH1D * hVarianceOfMean_Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn6[NCHI2];
  TH1D * hVarianceOfMean_Vn8[NCHI2];
  TH1D * hVarianceOfMean_Gamma1Exp[NCHI2];
  TH1D * hVarianceOfMean_Vn6Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn8Vn4[NCHI2];
  TH1D * hVarianceOfMean_Vn8Vn6[NCHI2];

  //
  // MAIN
  //
  setTDRStyle();
  //-- Turn off warning messages for the chi2 test:
  gErrorIgnoreLevel = kWarning;

  //-- Set Up Output files
  fOut = new TFile( Form("StatUncertChi2Cutoff_v%i.root", norder_), "recreate");

  for(int ichi = 0; ichi < NCHI2; ichi++){

    fOut->cd();
    hVarianceOfMean_Vn2[ichi]    = new TH1D( Form("hVarianceOfMean_Vn2_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn2_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn2[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn2[ichi]->GetXaxis()->SetTitle("#sigma^{2}");
    
    hVarianceOfMean_Vn4[ichi]    = new TH1D( Form("hVarianceOfMean_Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn4[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn4[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Vn6[ichi]    = new TH1D( Form("hVarianceOfMean_Vn6_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn6_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn6[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn6[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Vn8[ichi]    = new TH1D( Form("hVarianceOfMean_Vn8_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn8_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn8[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn8[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Gamma1Exp[ichi] = new TH1D( Form("hVarianceOfMean_Gamma1Exp_chi2Cut%.1f", chi2Cutoff[ichi]), Form("hVarianceOfMean_Gamma1Exp_chi2Cut%.1f", chi2Cutoff[ichi]), NCENT, centbinsDefault);
    hVarianceOfMean_Gamma1Exp[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Gamma1Exp[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Vn6Vn4[ichi]    = new TH1D( Form("hVarianceOfMean_Vn6Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn6Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn6Vn4[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn6Vn4[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Vn8Vn4[ichi]    = new TH1D( Form("hVarianceOfMean_Vn8Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn8Vn4_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn8Vn4[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn8Vn4[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

    hVarianceOfMean_Vn8Vn6[ichi]    = new TH1D( Form("hVarianceOfMean_Vn8Vn6_chi2Cut%.1f", chi2Cutoff[ichi]),    Form("hVarianceOfMean_Vn8Vn6_chi2Cut%.1f", chi2Cutoff[ichi]),    NCENT, centbinsDefault);
    hVarianceOfMean_Vn8Vn6[ichi]->GetXaxis()->SetTitle("Centrality %");
    hVarianceOfMean_Vn8Vn6[ichi]->GetXaxis()->SetTitle("#sigma^{2}");

  }

  for(int iS = 0; iS < NSPLIT; iS++){

    //-- Grab the analyzer file
    fAna[iS] = new TFile( Form("../../AnalyzerResults/CastleEbyE_Split%i.root", iS) );

    //-- Grab the unfolding file
    fUnf[iS] = new TFile( Form("../../UnfoldResults/dataResp/data%i_Split%i.root", norder_, iS) );

    for(int icent = 0; icent < NCENT; icent++){

      //-- Grab the observed histo
      hObs[icent][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_c%i", icent) );

      //-- Reset the flags that say whether or no the chi2 cutoff has been reached
      for(int ichi = 0; ichi < NCHI2; ichi++) cutoffFound[ichi] = 0;

      for(int i = 0; i < NITER; i++){

	//-- Grab the un(re)folded histos
	hUnfold[icent][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hreco%i_c%i",   iter[i], icent) );
	hRefold[icent][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hrefold%i_c%i", iter[i], icent) );


	double chi2Refold = hRefold[icent][iS][i]->Chi2Test(hObs[icent][iS], "CHI2/NDF");

	//-- loop over the scenarios and choose the iteration for each
	for(int ichi = 0; ichi < NCHI2; ichi++){

	  if( chi2Refold <= chi2Cutoff[ichi] && !cutoffFound[ichi] ){
	    iterCutoff[ichi][icent][iS] = i;
	    cutoffFound[ichi] = 1;
	  }
	  if( i == NITER-1 && !cutoffFound[ichi] ){
	    iterCutoff[ichi][icent][iS] = i;
	    cutoffFound[ichi] = 1;
	  }

	} //-- End chi2 scenario loop	

      } //-- End iter loop
    } //-- End cent loop
  } //-- End split loop


  //-- Now that we have the distributions, let's figure out the scatter on the reported quantities


  //--  First, get the sample mean of the split data
  double sampleMean_Vn2[NCHI2][NCENT];
  double sampleMean_Vn4[NCHI2][NCENT];
  double sampleMean_Vn6[NCHI2][NCENT];
  double sampleMean_Vn8[NCHI2][NCENT];
  double sampleMean_Gamma1Exp[NCHI2][NCENT];
  double sampleMean_Vn6Vn4[NCHI2][NCENT];
  double sampleMean_Vn8Vn4[NCHI2][NCENT];
  double sampleMean_Vn8Vn6[NCHI2][NCENT];

  double NSUCCESS_VN2[NCHI2][NCENT];
  double NSUCCESS_VN4[NCHI2][NCENT];
  double NSUCCESS_VN6[NCHI2][NCENT];
  double NSUCCESS_VN8[NCHI2][NCENT];
  double NSUCCESS_G1E[NCHI2][NCENT];
  double NSUCCESS_VN6VN4[NCHI2][NCENT];
  double NSUCCESS_VN8VN4[NCHI2][NCENT];
  double NSUCCESS_VN8VN6[NCHI2][NCENT];

  for(int icent = 0; icent < NCENT; icent++){
    for(int ichi = 0; ichi < NCHI2; ichi++){

      sampleMean_Vn2[ichi][icent]       = 0.;
      sampleMean_Vn4[ichi][icent]       = 0.;
      sampleMean_Vn6[ichi][icent]       = 0.;
      sampleMean_Vn8[ichi][icent]       = 0.;
      sampleMean_Gamma1Exp[ichi][icent] = 0.;
      sampleMean_Vn6Vn4[ichi][icent]    = 0.;
      sampleMean_Vn8Vn4[ichi][icent]    = 0.;
      sampleMean_Vn8Vn6[ichi][icent]    = 0.;

      NSUCCESS_VN2[ichi][icent]    = (double) NSPLIT;
      NSUCCESS_VN4[ichi][icent]    = (double) NSPLIT;
      NSUCCESS_VN6[ichi][icent]    = (double) NSPLIT;
      NSUCCESS_VN8[ichi][icent]    = (double) NSPLIT;
      NSUCCESS_G1E[ichi][icent]    = (double) NSPLIT;
      NSUCCESS_VN6VN4[ichi][icent] = (double) NSPLIT;
      NSUCCESS_VN8VN4[ichi][icent] = (double) NSPLIT;
      NSUCCESS_VN8VN6[ichi][icent] = (double) NSPLIT;

      for(int iS = 0; iS < NSPLIT; iS++){

	int iter = iterCutoff[ichi][icent][iS];
	FixUnfold( hUnfold[icent][iS][iter] );
	EbyECumu cumu(hUnfold[icent][iS][iter]);
	double vn2       = cumu.GetCumu_vn2();
	double vn4       = cumu.GetCumu_vn4();
	double vn6       = cumu.GetCumu_vn6();
	double vn8       = cumu.GetCumu_vn8();
	double gamma1exp = cumu.GetGamma1Exp();


	if( vn2 == 0 ) NSUCCESS_VN2[ichi][icent]   -= 1;
	else           sampleMean_Vn2[ichi][icent] += vn2;

	if( vn4 == 0 ) NSUCCESS_VN4[ichi][icent]   -= 1;
	else           sampleMean_Vn4[ichi][icent] += vn4;

	if( vn6 == 0 ) NSUCCESS_VN6[ichi][icent]   -= 1;
	else           sampleMean_Vn6[ichi][icent] += vn6;

	if( vn8 == 0 ) NSUCCESS_VN8[ichi][icent]   -= 1;
	else           sampleMean_Vn8[ichi][icent] += vn8;

	if( gamma1exp == 0 ) NSUCCESS_G1E[ichi][icent] -= 1;
	else                 sampleMean_Gamma1Exp[ichi][icent] += gamma1exp;

	if( vn6 == 0 || vn4 == 0 ) NSUCCESS_VN6VN4[ichi][icent]   -= 1;
	else                       sampleMean_Vn6Vn4[ichi][icent] += (vn6 / vn4);

	if( vn8 == 0 || vn4 == 0 ) NSUCCESS_VN8VN4[ichi][icent]   -= 1;
	else                       sampleMean_Vn8Vn4[ichi][icent] += (vn8 / vn4);

	if( vn8 == 0 || vn6 == 0 ) NSUCCESS_VN8VN6[ichi][icent]   -= 1;
	else                       sampleMean_Vn8Vn6[ichi][icent] += (vn8 / vn6);



      }

      if( NSUCCESS_VN2[ichi][icent] > 0 )    sampleMean_Vn2[ichi][icent]       /= NSUCCESS_VN2[ichi][icent];
      if( NSUCCESS_VN4[ichi][icent] > 0 )    sampleMean_Vn4[ichi][icent]       /= NSUCCESS_VN4[ichi][icent];
      if( NSUCCESS_VN6[ichi][icent] > 0 )    sampleMean_Vn6[ichi][icent]       /= NSUCCESS_VN6[ichi][icent];
      if( NSUCCESS_VN8[ichi][icent] > 0 )    sampleMean_Vn8[ichi][icent]       /= NSUCCESS_VN8[ichi][icent];
      if( NSUCCESS_G1E[ichi][icent] > 0 )    sampleMean_Gamma1Exp[ichi][icent] /= NSUCCESS_G1E[ichi][icent];
      if( NSUCCESS_VN6VN4[ichi][icent] > 0 ) sampleMean_Vn6Vn4[ichi][icent]    /= NSUCCESS_VN6VN4[ichi][icent];
      if( NSUCCESS_VN8VN4[ichi][icent] > 0 ) sampleMean_Vn8Vn4[ichi][icent]    /= NSUCCESS_VN8VN4[ichi][icent];
      if( NSUCCESS_VN8VN6[ichi][icent] > 0 ) sampleMean_Vn8Vn6[ichi][icent]    /= NSUCCESS_VN8VN6[ichi][icent];

    }
  }

  //-- Now get the sample variance
  double sampleVariance_Vn2[NCHI2][NCENT];
  double sampleVariance_Vn4[NCHI2][NCENT];
  double sampleVariance_Vn6[NCHI2][NCENT];
  double sampleVariance_Vn8[NCHI2][NCENT];
  double sampleVariance_Gamma1Exp[NCHI2][NCENT];
  double sampleVariance_Vn6Vn4[NCHI2][NCENT];
  double sampleVariance_Vn8Vn4[NCHI2][NCENT];
  double sampleVariance_Vn8Vn6[NCHI2][NCENT];

  double varianceOfMean_Vn2[NCHI2][NCENT];
  double varianceOfMean_Vn4[NCHI2][NCENT];
  double varianceOfMean_Vn6[NCHI2][NCENT];
  double varianceOfMean_Vn8[NCHI2][NCENT];
  double varianceOfMean_Gamma1Exp[NCHI2][NCENT];
  double varianceOfMean_Vn6Vn4[NCHI2][NCENT];
  double varianceOfMean_Vn8Vn4[NCHI2][NCENT];
  double varianceOfMean_Vn8Vn6[NCHI2][NCENT];

  for(int icent = 0; icent < NCENT; icent++){
    for(int ichi = 0; ichi < NCHI2; ichi++){

      sampleVariance_Vn2[ichi][icent]       = 0.;
      sampleVariance_Vn4[ichi][icent]       = 0.;
      sampleVariance_Vn6[ichi][icent]       = 0.;
      sampleVariance_Vn8[ichi][icent]       = 0.;
      sampleVariance_Gamma1Exp[ichi][icent] = 0.;
      sampleVariance_Vn6Vn4[ichi][icent]    = 0.;
      sampleVariance_Vn8Vn4[ichi][icent]    = 0.;
      sampleVariance_Vn8Vn6[ichi][icent]    = 0.;

      for(int iS = 0; iS < NSPLIT; iS++){

	int iter = iterCutoff[ichi][icent][iS];
	EbyECumu cumu(hUnfold[icent][iS][iter]);
	double vn2       = cumu.GetCumu_vn2();
	double vn4       = cumu.GetCumu_vn4();
	double vn6       = cumu.GetCumu_vn6();
	double vn8       = cumu.GetCumu_vn8();
	double gamma1exp = cumu.GetGamma1Exp();

	if( vn2 != 0 )             sampleVariance_Vn2[ichi][icent]       += pow( vn2 - sampleMean_Vn2[ichi][icent], 2);
	if( vn4 != 0 )             sampleVariance_Vn4[ichi][icent]       += pow( vn4 - sampleMean_Vn4[ichi][icent], 2);
	if( vn6 != 0 )             sampleVariance_Vn6[ichi][icent]       += pow( vn6 - sampleMean_Vn6[ichi][icent], 2);
	if( vn8 != 0 )             sampleVariance_Vn8[ichi][icent]       += pow( vn8 - sampleMean_Vn8[ichi][icent], 2);
	if( gamma1exp != 0 )       sampleVariance_Gamma1Exp[ichi][icent] += pow(gamma1exp - sampleMean_Gamma1Exp[ichi][icent], 2);
	if( vn6 != 0 && vn4 != 0 ) sampleVariance_Vn6Vn4[ichi][icent]    += pow( (vn6/vn4) - sampleMean_Vn6Vn4[ichi][icent], 2);
	if( vn8 != 0 && vn4 != 0 ) sampleVariance_Vn8Vn4[ichi][icent]    += pow( (vn8/vn4) - sampleMean_Vn8Vn4[ichi][icent], 2);
	if( vn8 != 0 && vn6 != 0 ) sampleVariance_Vn8Vn6[ichi][icent]    += pow( (vn8/vn6) - sampleMean_Vn8Vn6[ichi][icent], 2);
    }

      if( NSUCCESS_VN2[ichi][icent] > 1 )    sampleVariance_Vn2[ichi][icent]       /= (NSUCCESS_VN2[ichi][icent]-1);
      if( NSUCCESS_VN4[ichi][icent] > 1 )    sampleVariance_Vn4[ichi][icent]       /= (NSUCCESS_VN4[ichi][icent]-1);
      if( NSUCCESS_VN6[ichi][icent] > 1 )    sampleVariance_Vn6[ichi][icent]       /= (NSUCCESS_VN6[ichi][icent]-1);
      if( NSUCCESS_VN8[ichi][icent] > 1 )    sampleVariance_Vn8[ichi][icent]       /= (NSUCCESS_VN8[ichi][icent]-1);
      if( NSUCCESS_G1E[ichi][icent] > 1 )    sampleVariance_Gamma1Exp[ichi][icent] /= (NSUCCESS_G1E[ichi][icent]-1);
      if( NSUCCESS_VN6VN4[ichi][icent] > 1 ) sampleVariance_Vn6Vn4[ichi][icent]    /= (NSUCCESS_VN6VN4[ichi][icent]-1);
      if( NSUCCESS_VN8VN4[ichi][icent] > 1 ) sampleVariance_Vn8Vn4[ichi][icent]    /= (NSUCCESS_VN8VN4[ichi][icent]-1);
      if( NSUCCESS_VN8VN6[ichi][icent] > 1 ) sampleVariance_Vn8Vn6[ichi][icent]    /= (NSUCCESS_VN8VN6[ichi][icent]-1);

      if( NSUCCESS_VN2[ichi][icent] > 0 )    varianceOfMean_Vn2[ichi][icent]       = sampleVariance_Vn2[ichi][icent] / NSUCCESS_VN2[ichi][icent];
      if( NSUCCESS_VN4[ichi][icent] > 0 )    varianceOfMean_Vn4[ichi][icent]       = sampleVariance_Vn4[ichi][icent] / NSUCCESS_VN4[ichi][icent];
      if( NSUCCESS_VN6[ichi][icent] > 0 )    varianceOfMean_Vn6[ichi][icent]       = sampleVariance_Vn6[ichi][icent] / NSUCCESS_VN6[ichi][icent];
      if( NSUCCESS_VN8[ichi][icent] > 0 )    varianceOfMean_Vn8[ichi][icent]       = sampleVariance_Vn8[ichi][icent] / NSUCCESS_VN8[ichi][icent];
      if( NSUCCESS_G1E[ichi][icent] > 0 )    varianceOfMean_Gamma1Exp[ichi][icent] = sampleVariance_Gamma1Exp[ichi][icent] / NSUCCESS_G1E[ichi][icent];
      if( NSUCCESS_VN6VN4[ichi][icent] > 0 ) varianceOfMean_Vn6Vn4[ichi][icent]    = sampleVariance_Vn6Vn4[ichi][icent] / NSUCCESS_VN6VN4[ichi][icent];
      if( NSUCCESS_VN8VN4[ichi][icent] > 0 ) varianceOfMean_Vn8Vn4[ichi][icent]    = sampleVariance_Vn8Vn4[ichi][icent] / NSUCCESS_VN8VN4[ichi][icent];
      if( NSUCCESS_VN8VN6[ichi][icent] > 0 ) varianceOfMean_Vn8Vn6[ichi][icent]    = sampleVariance_Vn8Vn6[ichi][icent] / NSUCCESS_VN8VN6[ichi][icent];

      hVarianceOfMean_Vn2[ichi]->SetBinContent(icent+1, varianceOfMean_Vn2[ichi][icent]);
      hVarianceOfMean_Vn4[ichi]->SetBinContent(icent+1, varianceOfMean_Vn4[ichi][icent]);
      hVarianceOfMean_Vn6[ichi]->SetBinContent(icent+1, varianceOfMean_Vn6[ichi][icent]);
      hVarianceOfMean_Vn8[ichi]->SetBinContent(icent+1, varianceOfMean_Vn8[ichi][icent]);
      hVarianceOfMean_Gamma1Exp[ichi]->SetBinContent(icent+1, varianceOfMean_Gamma1Exp[ichi][icent]);
      hVarianceOfMean_Vn6Vn4[ichi]->SetBinContent(icent+1, varianceOfMean_Vn6Vn4[ichi][icent]);
      hVarianceOfMean_Vn8Vn4[ichi]->SetBinContent(icent+1, varianceOfMean_Vn8Vn4[ichi][icent]);
      hVarianceOfMean_Vn8Vn6[ichi]->SetBinContent(icent+1, varianceOfMean_Vn8Vn6[ichi][icent]);

    }
  }

  fOut->Write();

}

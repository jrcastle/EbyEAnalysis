#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TError.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void statUncertNewCC(){

  int norder_ = 2;

  //-- Analyzer Output
  TFile * fAna[NSPLIT];
  TH1D * hObs[NCENT][NSPLIT];

  //-- Unfolding output
  TFile * fUnf[NSPLIT];
  TH1D * hUnfold[NCENT][NSPLIT][NITER];
  TH1D * hRefold[NCENT][NSPLIT][NITER];

  bool iterCutoff[NCENT][NSPLIT];

  //-- Stat Uncert Output
  TFile * fOut;
  TH1D * hVarianceOfMean_Vn2;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
  TH1D * hVarianceOfMean_Gamma1Exp;
  TH1D * hVarianceOfMean_Vn6Vn4;
  TH1D * hVarianceOfMean_Vn8Vn4;
  TH1D * hVarianceOfMean_Vn8Vn6;

  //
  // MAIN
  //
  setTDRStyle();
  //-- Turn off warning messages for the chi2 test:
  gErrorIgnoreLevel = kWarning;

  //-- Set Up Output files
  fOut = new TFile( Form("StatUncertNewCC_v%i.root", norder_), "recreate");

  fOut->cd();
  hVarianceOfMean_Vn2    = new TH1D("hVarianceOfMean_Vn2_NewCC",    "hVarianceOfMean_Vn2_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn2->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn2->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn4    = new TH1D("hVarianceOfMean_Vn4_NewCC",    "hVarianceOfMean_Vn4_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn4->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn4->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn6    = new TH1D("hVarianceOfMean_Vn6_NewCC",    "hVarianceOfMean_Vn6_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn6->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn6->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn8    = new TH1D("hVarianceOfMean_Vn8_NewCC",    "hVarianceOfMean_Vn8_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn8->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn8->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Gamma1Exp = new TH1D("hVarianceOfMean_Gamma1Exp_NewCC", "hVarianceOfMean_Gamma1Exp_NewCC", NCENT, centbinsDefault);
  hVarianceOfMean_Gamma1Exp->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Gamma1Exp->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn6Vn4    = new TH1D("hVarianceOfMean_Vn6Vn4_NewCC",    "hVarianceOfMean_Vn6Vn4_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn6Vn4->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn6Vn4->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn8Vn4    = new TH1D("hVarianceOfMean_Vn8Vn4_NewCC",    "hVarianceOfMean_Vn8Vn4_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn8Vn4->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn8Vn4->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn8Vn6    = new TH1D("hVarianceOfMean_Vn8Vn6_NewCC",    "hVarianceOfMean_Vn8Vn6_NewCC",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn8Vn6->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn8Vn6->GetXaxis()->SetTitle("#sigma^{2}");

  for(int iS = 0; iS < NSPLIT; iS++){

    //-- Grab the analyzer file
    fAna[iS] = new TFile( Form("AnalyzerResults/CastleEbyE_Split%i.root", iS) );

    //-- Grab the unfolding file
    fUnf[iS] = new TFile( Form("UnfoldResults/dataResp/data%i_Split%i.root", norder_, iS) );

    for(int icent = 0; icent < NCENT; icent++){

      //-- Grab the observed histo
      hObs[icent][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_c%i", icent) );
      iterCutoff[icent][iS] = 0;

      for(int i = 0; i < NITER; i++){

	//-- Grab the un(re)folded histos
	hUnfold[icent][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hreco%i_c%i",   iter[i], icent) );
	hRefold[icent][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hrefold%i_c%i", iter[i], icent) );


	double chi2Refold = hRefold[icent][iS][i]->Chi2Test(hObs[icent][iS], "CHI2/NDF");
	if( chi2Refold < 1.2 ){
	  iterCutoff[icent][iS] = i;
	  break;
	}
	if( i == NITER-1 ){
	  iterCutoff[icent][iS] = i;
          break;
	}

      } //-- End iter loop
    } //-- End cent loop
  } //-- End split loop


  //-- Now that we have the distributions, let's figure out the scatter on the reported quantities


  //--  First, get the sample mean of the split data
  double sampleMean_Vn2[NCENT];
  double sampleMean_Vn4[NCENT];
  double sampleMean_Vn6[NCENT];
  double sampleMean_Vn8[NCENT];
  double sampleMean_Gamma1Exp[NCENT];
  double sampleMean_Vn6Vn4[NCENT];
  double sampleMean_Vn8Vn4[NCENT];
  double sampleMean_Vn8Vn6[NCENT];

  double NSUCCESS_VN2[NCENT];
  double NSUCCESS_VN4[NCENT];
  double NSUCCESS_VN6[NCENT];
  double NSUCCESS_VN8[NCENT];
  double NSUCCESS_G1E[NCENT];
  double NSUCCESS_VN6VN4[NCENT];
  double NSUCCESS_VN8VN4[NCENT];
  double NSUCCESS_VN8VN6[NCENT];

  for(int icent = 0; icent < NCENT; icent++){

    sampleMean_Vn2[icent]       = 0.;
    sampleMean_Vn4[icent]       = 0.;
    sampleMean_Vn6[icent]       = 0.;
    sampleMean_Vn8[icent]       = 0.;
    sampleMean_Gamma1Exp[icent] = 0.;
    sampleMean_Vn6Vn4[icent]    = 0.;
    sampleMean_Vn8Vn4[icent]    = 0.;
    sampleMean_Vn8Vn6[icent]    = 0.;

    NSUCCESS_VN2[icent]    = (double) NSPLIT;
    NSUCCESS_VN4[icent]    = (double) NSPLIT;
    NSUCCESS_VN6[icent]    = (double) NSPLIT;
    NSUCCESS_VN8[icent]    = (double) NSPLIT;
    NSUCCESS_G1E[icent]    = (double) NSPLIT;
    NSUCCESS_VN6VN4[icent] = (double) NSPLIT;
    NSUCCESS_VN8VN4[icent] = (double) NSPLIT;
    NSUCCESS_VN8VN6[icent] = (double) NSPLIT;

    for(int iS = 0; iS < NSPLIT; iS++){

      int iter = iterCutoff[icent][iS];
      FixUnfold( hUnfold[icent][iS][iter] );
      EbyECumu cumu(hUnfold[icent][iS][iter]);
      double vn2       = cumu.GetCumu_vn2();
      double vn4       = cumu.GetCumu_vn4();
      double vn6       = cumu.GetCumu_vn6();
      double vn8       = cumu.GetCumu_vn8();
      double gamma1exp = cumu.GetGamma1Exp();


      if( vn2 == 0 ) NSUCCESS_VN2[icent]   -= 1;
      else           sampleMean_Vn2[icent] += vn2;

      if( vn4 == 0 ) NSUCCESS_VN4[icent]   -= 1;
      else           sampleMean_Vn4[icent] += vn4;

      if( vn6 == 0 ) NSUCCESS_VN6[icent]   -= 1;
      else           sampleMean_Vn6[icent] += vn6;

      if( vn8 == 0 ) NSUCCESS_VN8[icent]   -= 1;
      else           sampleMean_Vn8[icent] += vn8;

      if( gamma1exp == 0 ) NSUCCESS_G1E[icent] -= 1;
      else                 sampleMean_Gamma1Exp[icent] += gamma1exp;

      if( vn6 == 0 || vn4 == 0 ) NSUCCESS_VN6VN4[icent]   -= 1;
      else                       sampleMean_Vn6Vn4[icent] += (vn6 / vn4);

      if( vn8 == 0 || vn4 == 0 ) NSUCCESS_VN8VN4[icent]   -= 1;
      else                       sampleMean_Vn8Vn4[icent] += (vn8 / vn4);

      if( vn8 == 0 || vn6 == 0 ) NSUCCESS_VN8VN6[icent]   -= 1;
      else                       sampleMean_Vn8Vn6[icent] += (vn8 / vn6);



    }

    if( NSUCCESS_VN2[icent] > 0 )    sampleMean_Vn2[icent]       /= NSUCCESS_VN2[icent];
    if( NSUCCESS_VN4[icent] > 0 )    sampleMean_Vn4[icent]       /= NSUCCESS_VN4[icent];
    if( NSUCCESS_VN6[icent] > 0 )    sampleMean_Vn6[icent]       /= NSUCCESS_VN6[icent];
    if( NSUCCESS_VN8[icent] > 0 )    sampleMean_Vn8[icent]       /= NSUCCESS_VN8[icent];
    if( NSUCCESS_G1E[icent] > 0 )    sampleMean_Gamma1Exp[icent] /= NSUCCESS_G1E[icent];
    if( NSUCCESS_VN6VN4[icent] > 0 ) sampleMean_Vn6Vn4[icent]    /= NSUCCESS_VN6VN4[icent];
    if( NSUCCESS_VN8VN4[icent] > 0 ) sampleMean_Vn8Vn4[icent]    /= NSUCCESS_VN8VN4[icent];
    if( NSUCCESS_VN8VN6[icent] > 0 ) sampleMean_Vn8Vn6[icent]    /= NSUCCESS_VN8VN6[icent];

  }

  //-- Now get the sample variance
  double sampleVariance_Vn2[NCENT];
  double sampleVariance_Vn4[NCENT];
  double sampleVariance_Vn6[NCENT];
  double sampleVariance_Vn8[NCENT];
  double sampleVariance_Gamma1Exp[NCENT];
  double sampleVariance_Vn6Vn4[NCENT];
  double sampleVariance_Vn8Vn4[NCENT];
  double sampleVariance_Vn8Vn6[NCENT];

  double varianceOfMean_Vn2[NCENT];
  double varianceOfMean_Vn4[NCENT];
  double varianceOfMean_Vn6[NCENT];
  double varianceOfMean_Vn8[NCENT];
  double varianceOfMean_Gamma1Exp[NCENT];
  double varianceOfMean_Vn6Vn4[NCENT];
  double varianceOfMean_Vn8Vn4[NCENT];
  double varianceOfMean_Vn8Vn6[NCENT];

  for(int icent = 0; icent < NCENT; icent++){

    sampleVariance_Vn2[icent]       = 0.;
    sampleVariance_Vn4[icent]       = 0.;
    sampleVariance_Vn6[icent]       = 0.;
    sampleVariance_Vn8[icent]       = 0.;
    sampleVariance_Gamma1Exp[icent] = 0.;
    sampleVariance_Vn6Vn4[icent]    = 0.;
    sampleVariance_Vn8Vn4[icent]    = 0.;
    sampleVariance_Vn8Vn6[icent]    = 0.;

    for(int iS = 0; iS < NSPLIT; iS++){

      int iter = iterCutoff[icent][iS];
      EbyECumu cumu(hUnfold[icent][iS][iter]);
      double vn2       = cumu.GetCumu_vn2();
      double vn4       = cumu.GetCumu_vn4();
      double vn6       = cumu.GetCumu_vn6();
      double vn8       = cumu.GetCumu_vn8();
      double gamma1exp = cumu.GetGamma1Exp();

      if( vn2 != 0 )             sampleVariance_Vn2[icent]       += pow( vn2 - sampleMean_Vn2[icent], 2);
      if( vn4 != 0 )             sampleVariance_Vn4[icent]       += pow( vn4 - sampleMean_Vn4[icent], 2);
      if( vn6 != 0 )             sampleVariance_Vn6[icent]       += pow( vn6 - sampleMean_Vn6[icent], 2);
      if( vn8 != 0 )             sampleVariance_Vn8[icent]       += pow( vn8 - sampleMean_Vn8[icent], 2);
      if( gamma1exp != 0 )       sampleVariance_Gamma1Exp[icent] += pow(gamma1exp - sampleMean_Gamma1Exp[icent], 2);
      if( vn6 != 0 && vn4 != 0 ) sampleVariance_Vn6Vn4[icent]    += pow( (vn6/vn4) - sampleMean_Vn6Vn4[icent], 2);
      if( vn8 != 0 && vn4 != 0 ) sampleVariance_Vn8Vn4[icent]    += pow( (vn8/vn4) - sampleMean_Vn8Vn4[icent], 2);
      if( vn8 != 0 && vn6 != 0 ) sampleVariance_Vn8Vn6[icent]    += pow( (vn8/vn6) - sampleMean_Vn8Vn6[icent], 2);
    }

    if( NSUCCESS_VN2[icent] > 1 )    sampleVariance_Vn2[icent]       /= (NSUCCESS_VN2[icent]-1);
    if( NSUCCESS_VN4[icent] > 1 )    sampleVariance_Vn4[icent]       /= (NSUCCESS_VN4[icent]-1);
    if( NSUCCESS_VN6[icent] > 1 )    sampleVariance_Vn6[icent]       /= (NSUCCESS_VN6[icent]-1);
    if( NSUCCESS_VN8[icent] > 1 )    sampleVariance_Vn8[icent]       /= (NSUCCESS_VN8[icent]-1);
    if( NSUCCESS_G1E[icent] > 1 )    sampleVariance_Gamma1Exp[icent] /= (NSUCCESS_G1E[icent]-1);
    if( NSUCCESS_VN6VN4[icent] > 1 ) sampleVariance_Vn6Vn4[icent]    /= (NSUCCESS_VN6VN4[icent]-1);
    if( NSUCCESS_VN8VN4[icent] > 1 ) sampleVariance_Vn8Vn4[icent]    /= (NSUCCESS_VN8VN4[icent]-1);
    if( NSUCCESS_VN8VN6[icent] > 1 ) sampleVariance_Vn8Vn6[icent]    /= (NSUCCESS_VN8VN6[icent]-1);

    if( NSUCCESS_VN2[icent] > 0 )    varianceOfMean_Vn2[icent]       = sampleVariance_Vn2[icent] / NSUCCESS_VN2[icent];
    if( NSUCCESS_VN4[icent] > 0 )    varianceOfMean_Vn4[icent]       = sampleVariance_Vn4[icent] / NSUCCESS_VN4[icent];
    if( NSUCCESS_VN6[icent] > 0 )    varianceOfMean_Vn6[icent]       = sampleVariance_Vn6[icent] / NSUCCESS_VN6[icent];
    if( NSUCCESS_VN8[icent] > 0 )    varianceOfMean_Vn8[icent]       = sampleVariance_Vn8[icent] / NSUCCESS_VN8[icent];
    if( NSUCCESS_G1E[icent] > 0 )    varianceOfMean_Gamma1Exp[icent] = sampleVariance_Gamma1Exp[icent] / NSUCCESS_G1E[icent];
    if( NSUCCESS_VN6VN4[icent] > 0 ) varianceOfMean_Vn6Vn4[icent]    = sampleVariance_Vn6Vn4[icent] / NSUCCESS_VN6VN4[icent];
    if( NSUCCESS_VN8VN4[icent] > 0 ) varianceOfMean_Vn8Vn4[icent]    = sampleVariance_Vn8Vn4[icent] / NSUCCESS_VN8VN4[icent];
    if( NSUCCESS_VN8VN6[icent] > 0 ) varianceOfMean_Vn8Vn6[icent]    = sampleVariance_Vn8Vn6[icent] / NSUCCESS_VN8VN6[icent];

    hVarianceOfMean_Vn2->SetBinContent(icent+1, varianceOfMean_Vn2[icent]);
    hVarianceOfMean_Vn4->SetBinContent(icent+1, varianceOfMean_Vn4[icent]);
    hVarianceOfMean_Vn6->SetBinContent(icent+1, varianceOfMean_Vn6[icent]);
    hVarianceOfMean_Vn8->SetBinContent(icent+1, varianceOfMean_Vn8[icent]);
    hVarianceOfMean_Gamma1Exp->SetBinContent(icent+1, varianceOfMean_Gamma1Exp[icent]);
    hVarianceOfMean_Vn6Vn4->SetBinContent(icent+1, varianceOfMean_Vn6Vn4[icent]);
    hVarianceOfMean_Vn8Vn4->SetBinContent(icent+1, varianceOfMean_Vn8Vn4[icent]);
    hVarianceOfMean_Vn8Vn6->SetBinContent(icent+1, varianceOfMean_Vn8Vn6[icent]);

  }

  //-- Spit out values so I can quickly assess things
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn2> = " << sampleMean_Vn2[icent] << " +/- " << sqrt( varianceOfMean_Vn2[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn4> = " << sampleMean_Vn4[icent] << " +/- " << sqrt( varianceOfMean_Vn4[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn6> = " << sampleMean_Vn6[icent] << " +/- " << sqrt( varianceOfMean_Vn6[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn8> = " << sampleMean_Vn8[icent] << " +/- " << sqrt( varianceOfMean_Vn8[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Gamma1Exp> = " << sampleMean_Gamma1Exp[icent] << " +/- " << sqrt( varianceOfMean_Gamma1Exp[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn6 / Vn4> = " << sampleMean_Vn6Vn4[icent] << " +/- " << sqrt( varianceOfMean_Vn6Vn4[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn8 / Vn4> = " << sampleMean_Vn8Vn4[icent] << " +/- " << sqrt( varianceOfMean_Vn8Vn4[icent] ) << std::endl;
  std::cout << "\n\n" << std::endl;
  for(int icent = 0; icent < NCENT; icent++) std::cout << Form("Cent Bin = %i", icent) << "\t<Vn8 / Vn6> = " << sampleMean_Vn8Vn6[icent] << " +/- " << sqrt( varianceOfMean_Vn8Vn6[icent] ) << std::endl;

  fOut->Write();

}

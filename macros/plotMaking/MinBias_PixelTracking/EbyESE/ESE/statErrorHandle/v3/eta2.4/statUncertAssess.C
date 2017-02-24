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

void statUncertAssess(){

  int norder_   = 3;

  //-- Analyzer Output
  TFile * fAna[NSPLIT];
  TH1D * hObs[NCENT][NEPSymm][NQN][NSPLIT];

  //-- Unfolding output
  TFile * fUnf[NSPLIT];
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NSPLIT][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NSPLIT][NITER];

  bool iterCutoff[NCENT][NEPSymm][NQN][NSPLIT];

  //-- Stat Uncert Output
  TFile * fOut;
  TH1D * hVarianceOfMean_RMSVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_MeanVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_StDevVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_RelFluctVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn2[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn6[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Gamma1Exp[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn6Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8Vn6[NEPSymm][NQN];

  //
  // MAIN
  //
  setTDRStyle();
  //-- Turn off warning messages for the chi2 test:
  gErrorIgnoreLevel = kWarning;

  //-- Set Up Output files
  fOut = new TFile( Form("StatisticalUncertainties_v%i.root", norder_), "recreate");

  for(int iEP = 0; iEP < NEPSymm; iEP++){
    if( iEP != EPSymmBin ) continue;
    for(int iqn = 0; iqn < NQN; iqn++){

      fOut->cd();
      hVarianceOfMean_RMSVn[iEP][iqn]   = new TH1D(Form("hVarianceOfMean_RMSVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),   Form("hVarianceOfMean_RMSVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),   NCENT, centbinsDefault);
      hVarianceOfMean_RMSVn[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_RMSVn[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_MeanVn[iEP][iqn]  = new TH1D(Form("hVarianceOfMean_MeanVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),  Form("hVarianceOfMean_MeanVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),  NCENT, centbinsDefault);
      hVarianceOfMean_MeanVn[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_MeanVn[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_StDevVn[iEP][iqn] = new TH1D(Form("hVarianceOfMean_StDevVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), Form("hVarianceOfMean_StDevVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), NCENT, centbinsDefault);
      hVarianceOfMean_StDevVn[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_StDevVn[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_RelFluctVn[iEP][iqn] = new TH1D(Form("hVarianceOfMean_RelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), Form("hVarianceOfMean_RelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), NCENT, centbinsDefault);
      hVarianceOfMean_RelFluctVn[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_RelFluctVn[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn2[iEP][iqn]     = new TH1D(Form("hVarianceOfMean_Vn2_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn2_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn2[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn2[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn4[iEP][iqn]     = new TH1D(Form("hVarianceOfMean_Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn4[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn4[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn6[iEP][iqn]     = new TH1D(Form("hVarianceOfMean_Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn6[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn6[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn8[iEP][iqn]     = new TH1D(Form("hVarianceOfMean_Vn8_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn8_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn8[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn8[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Gamma1Exp[iEP][iqn] = new TH1D(Form("hVarianceOfMean_Gamma1Exp_%s_qbin%i", EPSymmNames[iEP].data(), iqn), Form("hVarianceOfMean_Gamma1Exp_%s_qbin%i", EPSymmNames[iEP].data(), iqn), NCENT, centbinsDefault);
      hVarianceOfMean_Gamma1Exp[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Gamma1Exp[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn6Vn4[iEP][iqn]    = new TH1D(Form("hVarianceOfMean_Vn6Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn6Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn6Vn4[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn6Vn4[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");
      
      hVarianceOfMean_Vn8Vn4[iEP][iqn]    = new TH1D(Form("hVarianceOfMean_Vn8Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn8Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn8Vn4[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn8Vn4[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

      hVarianceOfMean_Vn8Vn6[iEP][iqn]    = new TH1D(Form("hVarianceOfMean_Vn8Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("hVarianceOfMean_Vn8Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      hVarianceOfMean_Vn8Vn6[iEP][iqn]->GetXaxis()->SetTitle("Centrality %");
      hVarianceOfMean_Vn8Vn6[iEP][iqn]->GetXaxis()->SetTitle("#sigma^{2}");

    }

  }

  for(int iS = 0; iS < NSPLIT; iS++){

    //-- Grab the analyzer file
    fAna[iS] = new TFile( Form("AnalyzerResults/CastleEbyE_Split%i.root", iS) );

    //-- Grab the unfolding file
    fUnf[iS] = new TFile( Form("UnfoldResults/dataResp/data%i_Split%i.root", norder_, iS) );
    
    for(int icent = 0; icent < NCENT; icent++){
      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;
	for(int iqn = 0; iqn < NQN; iqn++){
	  
	  //-- Grab the observed histo
	  hObs[icent][iEP][iqn][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	  iterCutoff[icent][iEP][iqn][iS] = 0;

	  for(int i = 0; i < NITER; i++){

	    //-- Grab the un(re)folded histos
	    hUnfold[icent][iEP][iqn][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	    hRefold[icent][iEP][iqn][iS][i] = (TH1D*) fUnf[iS]->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	    double chi2Refold = hRefold[icent][iEP][iqn][iS][i]->Chi2Test(hObs[icent][iEP][iqn][iS], "CHI2/NDF");
	    if( chi2Refold < 1.2 ){
	      iterCutoff[icent][iEP][iqn][iS] = i;
	      break;
	    }
	    if( i == NITER-1 ){
	      iterCutoff[icent][iEP][iqn][iS] = i;
	      break;
	    }

	  } //-- End iter loop
	} //-- End qn loop
      } //-- End cent loop
    } //-- End EP loop
  } //-- End split loop


  //-- Now that we have the distributions, let's figure out the scatter on the reported quantities


  //--  First, get the sample mean of the split data
  double sampleMean_RMSVn[NCENT][NEPSymm][NQN];
  double sampleMean_MeanVn[NCENT][NEPSymm][NQN];
  double sampleMean_StDevVn[NCENT][NEPSymm][NQN];
  double sampleMean_RelFluctVn[NCENT][NEPSymm][NQN];
  double sampleMean_Vn2[NCENT][NEPSymm][NQN];
  double sampleMean_Vn4[NCENT][NEPSymm][NQN];
  double sampleMean_Vn6[NCENT][NEPSymm][NQN];
  double sampleMean_Vn8[NCENT][NEPSymm][NQN];
  double sampleMean_Gamma1Exp[NCENT][NEPSymm][NQN];
  double sampleMean_Vn6Vn4[NCENT][NEPSymm][NQN];
  double sampleMean_Vn8Vn4[NCENT][NEPSymm][NQN];
  double sampleMean_Vn8Vn6[NCENT][NEPSymm][NQN];

  double NSUCCESS_RMSVN[NCENT][NEPSymm][NQN];
  double NSUCCESS_MEANVN[NCENT][NEPSymm][NQN];
  double NSUCCESS_STDEVVN[NCENT][NEPSymm][NQN];
  double NSUCCESS_RELFLUCTVN[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN2[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN4[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN6[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN8[NCENT][NEPSymm][NQN];
  double NSUCCESS_G1E[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN6VN4[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN8VN4[NCENT][NEPSymm][NQN];
  double NSUCCESS_VN8VN6[NCENT][NEPSymm][NQN];

  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	sampleMean_RMSVn[icent][iEP][iqn]      = 0.;
	sampleMean_MeanVn[icent][iEP][iqn]     = 0.;
	sampleMean_StDevVn[icent][iEP][iqn]    = 0.;
	sampleMean_RelFluctVn[icent][iEP][iqn] = 0.;
	sampleMean_Vn2[icent][iEP][iqn]        = 0.;
	sampleMean_Vn4[icent][iEP][iqn]        = 0.;
	sampleMean_Vn6[icent][iEP][iqn]        = 0.;
	sampleMean_Vn8[icent][iEP][iqn]        = 0.;
	sampleMean_Gamma1Exp[icent][iEP][iqn]  = 0.;
	sampleMean_Vn6Vn4[icent][iEP][iqn]     = 0.;
	sampleMean_Vn8Vn4[icent][iEP][iqn]     = 0.;
	sampleMean_Vn8Vn6[icent][iEP][iqn]     = 0.;

	NSUCCESS_RMSVN[icent][iEP][iqn]      = (double) NSPLIT;
	NSUCCESS_MEANVN[icent][iEP][iqn]     = (double) NSPLIT;
	NSUCCESS_STDEVVN[icent][iEP][iqn]    = (double) NSPLIT;
	NSUCCESS_RELFLUCTVN[icent][iEP][iqn] = (double) NSPLIT;
	NSUCCESS_VN2[icent][iEP][iqn]        = (double) NSPLIT;
	NSUCCESS_VN4[icent][iEP][iqn]        = (double) NSPLIT;
	NSUCCESS_VN6[icent][iEP][iqn]        = (double) NSPLIT;
	NSUCCESS_VN8[icent][iEP][iqn]        = (double) NSPLIT;
	NSUCCESS_G1E[icent][iEP][iqn]        = (double) NSPLIT;
	NSUCCESS_VN6VN4[icent][iEP][iqn]     = (double) NSPLIT;
	NSUCCESS_VN8VN4[icent][iEP][iqn]     = (double) NSPLIT;
	NSUCCESS_VN8VN6[icent][iEP][iqn]     = (double) NSPLIT;

	for(int iS = 0; iS < NSPLIT; iS++){
	  
	  int iter = iterCutoff[icent][iEP][iqn][iS];
	  FixUnfold( hUnfold[icent][iEP][iqn][iS][iter] );
	  EbyECumu cumu(hUnfold[icent][iEP][iqn][iS][iter]);
	  double vn2        = cumu.GetCumu_vn2();
	  double vn4        = cumu.GetCumu_vn4();
	  double vn6        = cumu.GetCumu_vn6();
	  double vn8        = cumu.GetCumu_vn8();
	  double gamma1exp  = cumu.GetGamma1Exp();
	  double meanvn     = hUnfold[icent][iEP][iqn][iS][iter]->GetMean();
	  double stdevvn    = hUnfold[icent][iEP][iqn][iS][iter]->GetRMS();
	  double rmsvn      = sqrt( meanvn*meanvn + stdevvn*stdevvn );
	  double relfluctvn = stdevvn / meanvn;

	  sampleMean_RMSVn[icent][iEP][iqn]      += rmsvn;
	  sampleMean_MeanVn[icent][iEP][iqn]     += meanvn;
	  sampleMean_StDevVn[icent][iEP][iqn]    += stdevvn;
	  sampleMean_RelFluctVn[icent][iEP][iqn] += relfluctvn;

	  if( vn2 == 0 ) NSUCCESS_VN2[icent][iEP][iqn]   -= 1;
	  else           sampleMean_Vn2[icent][iEP][iqn] += vn2;

	  if( vn4 == 0 ) NSUCCESS_VN4[icent][iEP][iqn]   -= 1;
	  else           sampleMean_Vn4[icent][iEP][iqn] += vn4;

	  if( vn6 == 0 ) NSUCCESS_VN6[icent][iEP][iqn]   -= 1;
	  else           sampleMean_Vn6[icent][iEP][iqn] += vn6;

	  if( vn8 == 0 ) NSUCCESS_VN8[icent][iEP][iqn]   -= 1;
	  else           sampleMean_Vn8[icent][iEP][iqn] += vn8;

	  if( gamma1exp == 0 ) NSUCCESS_G1E[icent][iEP][iqn] -= 1;
	  else                 sampleMean_Gamma1Exp[icent][iEP][iqn] += gamma1exp;

	  if( vn6 == 0 || vn4 == 0 ) NSUCCESS_VN6VN4[icent][iEP][iqn]   -= 1;
	  else                       sampleMean_Vn6Vn4[icent][iEP][iqn] += (vn6 / vn4);

	  if( vn8 == 0 || vn4 == 0 ) NSUCCESS_VN8VN4[icent][iEP][iqn]   -= 1;
	  else                       sampleMean_Vn8Vn4[icent][iEP][iqn] += (vn8 / vn4);

	  if( vn8 == 0 || vn6 == 0 ) NSUCCESS_VN8VN6[icent][iEP][iqn]   -= 1;
	  else                       sampleMean_Vn8Vn6[icent][iEP][iqn] += (vn8 / vn6);

	} //-- End split loop
	if( NSUCCESS_RMSVN[icent][iEP][iqn] > 0 )      sampleMean_RMSVn[icent][iEP][iqn] /= NSUCCESS_RMSVN[icent][iEP][iqn];
	if( NSUCCESS_MEANVN[icent][iEP][iqn] > 0 )     sampleMean_MeanVn[icent][iEP][iqn] /= NSUCCESS_MEANVN[icent][iEP][iqn];
	if( NSUCCESS_STDEVVN[icent][iEP][iqn] > 0 )    sampleMean_StDevVn[icent][iEP][iqn] /= NSUCCESS_STDEVVN[icent][iEP][iqn];
	if( NSUCCESS_RELFLUCTVN[icent][iEP][iqn] > 0 ) sampleMean_RelFluctVn[icent][iEP][iqn] /= NSUCCESS_RELFLUCTVN[icent][iEP][iqn];
	if( NSUCCESS_VN2[icent][iEP][iqn] > 0 )        sampleMean_Vn2[icent][iEP][iqn] /= NSUCCESS_VN2[icent][iEP][iqn];
	if( NSUCCESS_VN4[icent][iEP][iqn] > 0 )        sampleMean_Vn4[icent][iEP][iqn] /= NSUCCESS_VN4[icent][iEP][iqn];
	if( NSUCCESS_VN6[icent][iEP][iqn] > 0 )        sampleMean_Vn6[icent][iEP][iqn] /= NSUCCESS_VN6[icent][iEP][iqn];
	if( NSUCCESS_VN8[icent][iEP][iqn] > 0 )        sampleMean_Vn8[icent][iEP][iqn] /= NSUCCESS_VN8[icent][iEP][iqn];
	if( NSUCCESS_G1E[icent][iEP][iqn] > 0 )        sampleMean_Gamma1Exp[icent][iEP][iqn] /= NSUCCESS_G1E[icent][iEP][iqn];
	if( NSUCCESS_VN6VN4[icent][iEP][iqn] > 0 )     sampleMean_Vn6Vn4[icent][iEP][iqn] /= NSUCCESS_VN6VN4[icent][iEP][iqn];
	if( NSUCCESS_VN8VN4[icent][iEP][iqn] > 0 )     sampleMean_Vn8Vn4[icent][iEP][iqn] /= NSUCCESS_VN8VN4[icent][iEP][iqn];
	if( NSUCCESS_VN8VN6[icent][iEP][iqn] > 0 )     sampleMean_Vn8Vn6[icent][iEP][iqn] /= NSUCCESS_VN8VN6[icent][iEP][iqn];

      } //-- End  qn loop
    } //-- End EP loop
  } //-- End cent loop

  //-- Now get the sample variance
  double sampleVariance_RMSVn[NCENT][NEPSymm][NQN];
  double sampleVariance_MeanVn[NCENT][NEPSymm][NQN];
  double sampleVariance_StDevVn[NCENT][NEPSymm][NQN];
  double sampleVariance_RelFluctVn[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn2[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn4[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn6[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn8[NCENT][NEPSymm][NQN];
  double sampleVariance_Gamma1Exp[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn6Vn4[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn8Vn4[NCENT][NEPSymm][NQN];
  double sampleVariance_Vn8Vn6[NCENT][NEPSymm][NQN];

  double varianceOfMean_RMSVn[NCENT][NEPSymm][NQN];
  double varianceOfMean_MeanVn[NCENT][NEPSymm][NQN];
  double varianceOfMean_StDevVn[NCENT][NEPSymm][NQN];
  double varianceOfMean_RelFluctVn[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn2[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn4[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn6[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn8[NCENT][NEPSymm][NQN];
  double varianceOfMean_Gamma1Exp[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn6Vn4[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn8Vn4[NCENT][NEPSymm][NQN];
  double varianceOfMean_Vn8Vn6[NCENT][NEPSymm][NQN];

  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	sampleVariance_RMSVn[icent][iEP][iqn]      = 0.;
	sampleVariance_MeanVn[icent][iEP][iqn]     = 0.;
	sampleVariance_StDevVn[icent][iEP][iqn]    = 0.;
	sampleVariance_RelFluctVn[icent][iEP][iqn] = 0.;
	sampleVariance_Vn2[icent][iEP][iqn]        = 0.;
	sampleVariance_Vn4[icent][iEP][iqn]        = 0.;
	sampleVariance_Vn6[icent][iEP][iqn]        = 0.;
	sampleVariance_Vn8[icent][iEP][iqn]        = 0.;
	sampleVariance_Gamma1Exp[icent][iEP][iqn]  = 0.;
	sampleVariance_Vn6Vn4[icent][iEP][iqn]     = 0.;
	sampleVariance_Vn8Vn4[icent][iEP][iqn]     = 0.;
	sampleVariance_Vn8Vn6[icent][iEP][iqn]     = 0.;

	for(int iS = 0; iS < NSPLIT; iS++){

	  int iter = iterCutoff[icent][iEP][iqn][iS];
	  EbyECumu cumu(hUnfold[icent][iEP][iqn][iS][iter]);
	  double vn2        = cumu.GetCumu_vn2();
	  double vn4        = cumu.GetCumu_vn4();
	  double vn6        = cumu.GetCumu_vn6();
	  double vn8        = cumu.GetCumu_vn8();
	  double gamma1exp  = cumu.GetGamma1Exp();
	  double meanvn     = hUnfold[icent][iEP][iqn][iS][iter]->GetMean();
          double stdevvn    = hUnfold[icent][iEP][iqn][iS][iter]->GetRMS();
          double rmsvn      = sqrt( meanvn*meanvn + stdevvn*stdevvn );
	  double relfluctvn = stdevvn / meanvn;

	  sampleVariance_RMSVn[icent][iEP][iqn]      += pow( rmsvn - sampleMean_RMSVn[icent][iEP][iqn], 2);
	  sampleVariance_MeanVn[icent][iEP][iqn]     += pow( meanvn - sampleMean_MeanVn[icent][iEP][iqn], 2);
	  sampleVariance_StDevVn[icent][iEP][iqn]    += pow( stdevvn - sampleMean_StDevVn[icent][iEP][iqn], 2);
	  sampleVariance_RelFluctVn[icent][iEP][iqn] += pow( relfluctvn - sampleMean_RelFluctVn[icent][iEP][iqn], 2);

	  if( vn2 != 0 )             sampleVariance_Vn2[icent][iEP][iqn]       += pow( vn2 - sampleMean_Vn2[icent][iEP][iqn], 2);
	  if( vn4 != 0 )             sampleVariance_Vn4[icent][iEP][iqn]       += pow( vn4 - sampleMean_Vn4[icent][iEP][iqn], 2);
	  if( vn6 != 0 )             sampleVariance_Vn6[icent][iEP][iqn]       += pow( vn6 - sampleMean_Vn6[icent][iEP][iqn], 2);
	  if( vn8 != 0 )             sampleVariance_Vn8[icent][iEP][iqn]       += pow( vn8 - sampleMean_Vn8[icent][iEP][iqn], 2);
	  if( gamma1exp != 0 )       sampleVariance_Gamma1Exp[icent][iEP][iqn] += pow(gamma1exp - sampleMean_Gamma1Exp[icent][iEP][iqn], 2);
	  if( vn6 != 0 && vn4 != 0 ) sampleVariance_Vn6Vn4[icent][iEP][iqn]    += pow( (vn6/vn4) - sampleMean_Vn6Vn4[icent][iEP][iqn], 2);
	  if( vn8 != 0 && vn4 != 0 ) sampleVariance_Vn8Vn4[icent][iEP][iqn]    += pow( (vn8/vn4) - sampleMean_Vn8Vn4[icent][iEP][iqn], 2);
	  if( vn8 != 0 && vn6 != 0 ) sampleVariance_Vn8Vn6[icent][iEP][iqn]    += pow( (vn8/vn6) - sampleMean_Vn8Vn6[icent][iEP][iqn], 2);

	}

	if( NSUCCESS_RMSVN[icent][iEP][iqn] > 1 )      sampleVariance_RMSVn[icent][iEP][iqn] /= (NSUCCESS_RMSVN[icent][iEP][iqn]-1);
	if( NSUCCESS_MEANVN[icent][iEP][iqn] > 1 )     sampleVariance_MeanVn[icent][iEP][iqn] /= (NSUCCESS_MEANVN[icent][iEP][iqn]-1);
	if( NSUCCESS_STDEVVN[icent][iEP][iqn] > 1 )    sampleVariance_StDevVn[icent][iEP][iqn] /= (NSUCCESS_STDEVVN[icent][iEP][iqn]-1);
	if( NSUCCESS_RELFLUCTVN[icent][iEP][iqn] > 1 ) sampleVariance_RelFluctVn[icent][iEP][iqn] /= (NSUCCESS_RELFLUCTVN[icent][iEP][iqn]-1);
	if( NSUCCESS_VN2[icent][iEP][iqn] > 1 )        sampleVariance_Vn2[icent][iEP][iqn] /= (NSUCCESS_VN2[icent][iEP][iqn]-1);
	if( NSUCCESS_VN4[icent][iEP][iqn] > 1 )        sampleVariance_Vn4[icent][iEP][iqn] /= (NSUCCESS_VN4[icent][iEP][iqn]-1);
	if( NSUCCESS_VN6[icent][iEP][iqn] > 1 )        sampleVariance_Vn6[icent][iEP][iqn] /= (NSUCCESS_VN6[icent][iEP][iqn]-1);
	if( NSUCCESS_VN8[icent][iEP][iqn] > 1 )        sampleVariance_Vn8[icent][iEP][iqn] /= (NSUCCESS_VN8[icent][iEP][iqn]-1);
	if( NSUCCESS_G1E[icent][iEP][iqn] > 1 )        sampleVariance_Gamma1Exp[icent][iEP][iqn] /= (NSUCCESS_G1E[icent][iEP][iqn]-1);
	if( NSUCCESS_VN6VN4[icent][iEP][iqn] > 1 )     sampleVariance_Vn6Vn4[icent][iEP][iqn] /= (NSUCCESS_VN6VN4[icent][iEP][iqn]-1);
	if( NSUCCESS_VN8VN4[icent][iEP][iqn] > 1 )     sampleVariance_Vn8Vn4[icent][iEP][iqn] /= (NSUCCESS_VN8VN4[icent][iEP][iqn]-1);
	if( NSUCCESS_VN8VN6[icent][iEP][iqn] > 1 )     sampleVariance_Vn8Vn6[icent][iEP][iqn] /= (NSUCCESS_VN8VN6[icent][iEP][iqn]-1);

	if( NSUCCESS_RMSVN[icent][iEP][iqn] > 0 )      varianceOfMean_RMSVn[icent][iEP][iqn]      = sampleVariance_RMSVn[icent][iEP][iqn] / NSUCCESS_RMSVN[icent][iEP][iqn];
	if( NSUCCESS_MEANVN[icent][iEP][iqn] > 0 )     varianceOfMean_MeanVn[icent][iEP][iqn]     = sampleVariance_MeanVn[icent][iEP][iqn] / NSUCCESS_MEANVN[icent][iEP][iqn];
	if( NSUCCESS_STDEVVN[icent][iEP][iqn] > 0 )    varianceOfMean_StDevVn[icent][iEP][iqn]    = sampleVariance_StDevVn[icent][iEP][iqn] / NSUCCESS_STDEVVN[icent][iEP][iqn];
	if( NSUCCESS_RELFLUCTVN[icent][iEP][iqn] > 0 ) varianceOfMean_RelFluctVn[icent][iEP][iqn] = sampleVariance_RelFluctVn[icent][iEP][iqn] / NSUCCESS_RELFLUCTVN[icent][iEP][iqn];
	if( NSUCCESS_VN2[icent][iEP][iqn] > 0 )        varianceOfMean_Vn2[icent][iEP][iqn]        = sampleVariance_Vn2[icent][iEP][iqn] / NSUCCESS_VN2[icent][iEP][iqn];
	if( NSUCCESS_VN4[icent][iEP][iqn] > 0 )        varianceOfMean_Vn4[icent][iEP][iqn]        = sampleVariance_Vn4[icent][iEP][iqn] / NSUCCESS_VN4[icent][iEP][iqn];
	if( NSUCCESS_VN6[icent][iEP][iqn] > 0 )        varianceOfMean_Vn6[icent][iEP][iqn]        = sampleVariance_Vn6[icent][iEP][iqn] / NSUCCESS_VN6[icent][iEP][iqn];
	if( NSUCCESS_VN8[icent][iEP][iqn] > 0 )        varianceOfMean_Vn8[icent][iEP][iqn]        = sampleVariance_Vn8[icent][iEP][iqn] / NSUCCESS_VN8[icent][iEP][iqn];
	if( NSUCCESS_G1E[icent][iEP][iqn] > 0 )        varianceOfMean_Gamma1Exp[icent][iEP][iqn]  = sampleVariance_Gamma1Exp[icent][iEP][iqn] / NSUCCESS_G1E[icent][iEP][iqn];
	if( NSUCCESS_VN6VN4[icent][iEP][iqn] > 0 )     varianceOfMean_Vn6Vn4[icent][iEP][iqn]     = sampleVariance_Vn6Vn4[icent][iEP][iqn] / NSUCCESS_VN6VN4[icent][iEP][iqn];
	if( NSUCCESS_VN8VN4[icent][iEP][iqn] > 0 )     varianceOfMean_Vn8Vn4[icent][iEP][iqn]     = sampleVariance_Vn8Vn4[icent][iEP][iqn] / NSUCCESS_VN8VN4[icent][iEP][iqn];
	if( NSUCCESS_VN8VN6[icent][iEP][iqn] > 0 )     varianceOfMean_Vn8Vn6[icent][iEP][iqn]     = sampleVariance_Vn8Vn6[icent][iEP][iqn] / NSUCCESS_VN8VN6[icent][iEP][iqn];

	hVarianceOfMean_RMSVn[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_RMSVn[icent][iEP][iqn]);
	hVarianceOfMean_MeanVn[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_MeanVn[icent][iEP][iqn]);
	hVarianceOfMean_StDevVn[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_StDevVn[icent][iEP][iqn]);
	hVarianceOfMean_RelFluctVn[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_RelFluctVn[icent][iEP][iqn]);
	hVarianceOfMean_Vn2[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn2[icent][iEP][iqn]);
	hVarianceOfMean_Vn4[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn4[icent][iEP][iqn]);
	hVarianceOfMean_Vn6[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn6[icent][iEP][iqn]);
	hVarianceOfMean_Vn8[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn8[icent][iEP][iqn]);
	hVarianceOfMean_Gamma1Exp[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Gamma1Exp[icent][iEP][iqn]);
	hVarianceOfMean_Vn6Vn4[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn6Vn4[icent][iEP][iqn]);
	hVarianceOfMean_Vn8Vn4[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn8Vn4[icent][iEP][iqn]);
	hVarianceOfMean_Vn8Vn6[iEP][iqn]->SetBinContent(icent+1, varianceOfMean_Vn8Vn6[icent][iEP][iqn]);

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

  fOut->Write();

}

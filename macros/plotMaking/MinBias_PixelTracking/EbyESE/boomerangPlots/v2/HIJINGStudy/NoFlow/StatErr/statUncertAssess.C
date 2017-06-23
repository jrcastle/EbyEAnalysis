#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void statUncertAssess(){

  int norder_ = 2;
  const int c = 0;

  //-- Analyzer Output
  TFile * fAna[NSPLIT];
  TH1D * hObs[NCENT][NSPLIT];
  TH1D * hV2True[NCENT][NSPLIT];

  //-- Unfolding output
  TFile * fUnf[NSPLIT];
  TH1D * hUnfold[NCENT][NSPLIT][NITER];
  TH1D * hRefold[NCENT][NSPLIT][NITER];

  int iterCutoff[NCENT][NSPLIT];

  //-- Stat Uncert Output
  TFile * fOut;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
  TH1D * hVarianceOfMean_Vn4True;
  TH1D * hVarianceOfMean_Vn6True;
  TH1D * hVarianceOfMean_Vn8True;

  //
  // MAIN
  //
  setTDRStyle();
  //-- Turn off warning messages for the chi2 test:
  gErrorIgnoreLevel = kWarning;

  //-- Set Up Output files
  fOut = new TFile( Form("StatisticalUncertainties_v%i.root", norder_), "recreate");

  fOut->cd();
  hVarianceOfMean_Vn4     = new TH1D("hVarianceOfMean_Vn4",    "hVarianceOfMean_Vn4",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn4->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn4->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn6     = new TH1D("hVarianceOfMean_Vn6",    "hVarianceOfMean_Vn6",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn6->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn6->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn8     = new TH1D("hVarianceOfMean_Vn8",    "hVarianceOfMean_Vn8",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn8->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn8->GetXaxis()->SetTitle("#sigma^{2}");
  /*
  hVarianceOfMean_Vn4True = new TH1D("hVarianceOfMean_Vn4True",    "hVarianceOfMean_Vn4True",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn4True->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn4True->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn6True = new TH1D("hVarianceOfMean_Vn6True",    "hVarianceOfMean_Vn6True",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn6True->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn6True->GetXaxis()->SetTitle("#sigma^{2}");

  hVarianceOfMean_Vn8True = new TH1D("hVarianceOfMean_Vn8True",    "hVarianceOfMean_Vn8True",    NCENT, centbinsDefault);
  hVarianceOfMean_Vn8True->GetXaxis()->SetTitle("Centrality %");
  hVarianceOfMean_Vn8True->GetXaxis()->SetTitle("#sigma^{2}");
  */
  for(int iS = 0; iS < NSPLIT; iS++){

    //-- Grab the analyzer file
    fAna[iS] = new TFile( Form("AnalyzerResults/CastleEbyE_Split%i.root", iS) );

    //-- Grab the unfolding file
    fUnf[iS] = new TFile( Form("UnfoldResults/dataResp/data%i_Split%i.root", norder_, iS) );

    for(int icent = 0; icent < NCENT; icent++){

      if(icent != c) continue;

      //-- Grab the observed and truth histo
      hObs[icent][iS]    = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_c%i", icent) );
      //hV2True[icent][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hV2True_c%i", icent) );
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
  double sampleMean_Vn4[NCENT];
  double sampleMean_Vn6[NCENT];
  double sampleMean_Vn8[NCENT];
  //double sampleMean_Vn4True[NCENT];
  //double sampleMean_Vn6True[NCENT];
  //double sampleMean_Vn8True[NCENT];

  double NSUCCESS_VN4[NCENT];
  double NSUCCESS_VN6[NCENT];
  double NSUCCESS_VN8[NCENT];
  //double NSUCCESS_VN4TRUE[NCENT];
  //double NSUCCESS_VN6TRUE[NCENT];
  //double NSUCCESS_VN8TRUE[NCENT];


  for(int icent = 0; icent < NCENT; icent++){

    if(icent != c) continue;

    sampleMean_Vn4[icent]     = 0.;
    sampleMean_Vn6[icent]     = 0.;
    sampleMean_Vn8[icent]     = 0.;
    //sampleMean_Vn4True[icent] = 0.;
    //sampleMean_Vn6True[icent] = 0.;
    //sampleMean_Vn8True[icent] = 0.;

    NSUCCESS_VN4[icent]     = (double) NSPLIT;
    NSUCCESS_VN6[icent]     = (double) NSPLIT;
    NSUCCESS_VN8[icent]     = (double) NSPLIT;
    //NSUCCESS_VN4TRUE[icent] = (double) NSPLIT;
    //NSUCCESS_VN6TRUE[icent] = (double) NSPLIT;
    //NSUCCESS_VN8TRUE[icent] = (double) NSPLIT;

    for(int iS = 0; iS < NSPLIT; iS++){

      int iter = iterCutoff[icent][iS];
      FixUnfold( hUnfold[icent][iS][iter] );
      EbyECumu cumu(hUnfold[icent][iS][iter]);
      double vn4 = cumu.GetCumu_vn4();
      double vn6 = cumu.GetCumu_vn6();
      double vn8 = cumu.GetCumu_vn8();

      //EbyECumu cumut(hV2True[icent][iS]);
      //double vn4t = cumut.GetCumu_vn4();
      //double vn6t = cumut.GetCumu_vn6();
      //double vn8t = cumut.GetCumu_vn8();

      if( vn4 == 0 ) NSUCCESS_VN4[icent]   -= 1;
      else           sampleMean_Vn4[icent] += vn4;
      if( vn6 == 0 ) NSUCCESS_VN6[icent]   -= 1;
      else           sampleMean_Vn6[icent] += vn6;
      if( vn8 == 0 ) NSUCCESS_VN8[icent]   -= 1;
      else           sampleMean_Vn8[icent] += vn8;

      //if( vn4t == 0 ) NSUCCESS_VN4TRUE[icent]   -= 1;
      //else            sampleMean_Vn4True[icent] += vn4;
      //if( vn6t == 0 ) NSUCCESS_VN6TRUE[icent]   -= 1;
      //else            sampleMean_Vn6True[icent] += vn6;
      //if( vn8t == 0 ) NSUCCESS_VN8TRUE[icent]   -= 1;
      //else            sampleMean_Vn8True[icent] += vn8;

    }

    if( NSUCCESS_VN4[icent] > 0 )     sampleMean_Vn4[icent] /= NSUCCESS_VN4[icent];
    if( NSUCCESS_VN6[icent] > 0 )     sampleMean_Vn6[icent] /= NSUCCESS_VN6[icent];
    if( NSUCCESS_VN8[icent] > 0 )     sampleMean_Vn8[icent] /= NSUCCESS_VN8[icent];
    //if( NSUCCESS_VN4TRUE[icent] > 0 ) sampleMean_Vn4True[icent] /= NSUCCESS_VN4TRUE[icent];
    //if( NSUCCESS_VN6TRUE[icent] > 0 ) sampleMean_Vn6True[icent] /= NSUCCESS_VN6TRUE[icent];
    //if( NSUCCESS_VN8TRUE[icent] > 0 ) sampleMean_Vn8True[icent] /= NSUCCESS_VN8TRUE[icent];

  }

  //-- Now get the sample variance
  double sampleVariance_Vn4[NCENT];
  double sampleVariance_Vn6[NCENT];
  double sampleVariance_Vn8[NCENT];
  //double sampleVariance_Vn4True[NCENT];
  //double sampleVariance_Vn6True[NCENT];
  //double sampleVariance_Vn8True[NCENT];

  double varianceOfMean_Vn4[NCENT];
  double varianceOfMean_Vn6[NCENT];
  double varianceOfMean_Vn8[NCENT];
  //double varianceOfMean_Vn4True[NCENT];
  //double varianceOfMean_Vn6True[NCENT];
  //double varianceOfMean_Vn8True[NCENT];

  for(int icent = 0; icent < NCENT; icent++){

    if(icent != c) continue;

    sampleVariance_Vn4[icent]     = 0.;
    sampleVariance_Vn6[icent]     = 0.;
    sampleVariance_Vn8[icent]     = 0.;
    //sampleVariance_Vn4True[icent] = 0.;
    //sampleVariance_Vn6True[icent] = 0.;
    //sampleVariance_Vn8True[icent] = 0.;

    for(int iS = 0; iS < NSPLIT; iS++){

      int iter = iterCutoff[icent][iS];
      EbyECumu cumu(hUnfold[icent][iS][iter]);
      double vn4       = cumu.GetCumu_vn4();
      double vn6       = cumu.GetCumu_vn6();
      double vn8       = cumu.GetCumu_vn8();

      //EbyECumu cumut(hV2True[icent][iS]);
      //double vn4t = cumut.GetCumu_vn4();
      //double vn6t = cumut.GetCumu_vn6();
      //double vn8t = cumut.GetCumu_vn8();

      if( vn4 != 0 )  sampleVariance_Vn4[icent]     += pow( vn4 - sampleMean_Vn4[icent], 2);
      if( vn6 != 0 )  sampleVariance_Vn6[icent]     += pow( vn6 - sampleMean_Vn6[icent], 2);
      if( vn8 != 0 )  sampleVariance_Vn8[icent]     += pow( vn8 - sampleMean_Vn8[icent], 2);
      //if( vn4t != 0 ) sampleVariance_Vn4True[icent] += pow( vn4t - sampleMean_Vn4True[icent], 2);
      //if( vn6t != 0 ) sampleVariance_Vn6True[icent] += pow( vn6t - sampleMean_Vn6True[icent], 2);
      //if( vn8t != 0 ) sampleVariance_Vn8True[icent] += pow( vn8t - sampleMean_Vn8True[icent], 2);

    }

    if( NSUCCESS_VN4[icent] > 1 )     sampleVariance_Vn4[icent]     /= (NSUCCESS_VN4[icent]-1);
    if( NSUCCESS_VN6[icent] > 1 )     sampleVariance_Vn6[icent]     /= (NSUCCESS_VN6[icent]-1);
    if( NSUCCESS_VN8[icent] > 1 )     sampleVariance_Vn8[icent]     /= (NSUCCESS_VN8[icent]-1);
    //if( NSUCCESS_VN4TRUE[icent] > 1 ) sampleVariance_Vn4True[icent] /= (NSUCCESS_VN4TRUE[icent]-1);
    //if( NSUCCESS_VN6TRUE[icent] > 1 ) sampleVariance_Vn6True[icent] /= (NSUCCESS_VN6TRUE[icent]-1);
    //if( NSUCCESS_VN8TRUE[icent] > 1 ) sampleVariance_Vn8True[icent] /= (NSUCCESS_VN8TRUE[icent]-1);

    if( NSUCCESS_VN4[icent] > 0 )       varianceOfMean_Vn4[icent]     = sampleVariance_Vn4[icent] / NSUCCESS_VN4[icent];
    if( NSUCCESS_VN6[icent] > 0 )       varianceOfMean_Vn6[icent]     = sampleVariance_Vn6[icent] / NSUCCESS_VN6[icent];
    if( NSUCCESS_VN8[icent] > 0 )       varianceOfMean_Vn8[icent]     = sampleVariance_Vn8[icent] / NSUCCESS_VN8[icent];
    //if( NSUCCESS_VN4TRUE[icent] > 0 )   varianceOfMean_Vn4True[icent] = sampleVariance_Vn4True[icent] / NSUCCESS_VN4TRUE[icent];
    //if( NSUCCESS_VN6TRUE[icent] > 0 )   varianceOfMean_Vn6True[icent] = sampleVariance_Vn6True[icent] / NSUCCESS_VN6TRUE[icent];
    //if( NSUCCESS_VN8TRUE[icent] > 0 )   varianceOfMean_Vn8True[icent] = sampleVariance_Vn8True[icent] / NSUCCESS_VN8TRUE[icent];

    hVarianceOfMean_Vn4->SetBinContent(icent+1, varianceOfMean_Vn4[icent]);
    hVarianceOfMean_Vn6->SetBinContent(icent+1, varianceOfMean_Vn6[icent]);
    hVarianceOfMean_Vn8->SetBinContent(icent+1, varianceOfMean_Vn8[icent]);
    //hVarianceOfMean_Vn4True->SetBinContent(icent+1, varianceOfMean_Vn4True[icent]);
    //hVarianceOfMean_Vn6True->SetBinContent(icent+1, varianceOfMean_Vn6True[icent]);
    //hVarianceOfMean_Vn8True->SetBinContent(icent+1, varianceOfMean_Vn8True[icent]);

  }

  fOut->Write();
  std::cout << "DONE!" << std::endl;

}

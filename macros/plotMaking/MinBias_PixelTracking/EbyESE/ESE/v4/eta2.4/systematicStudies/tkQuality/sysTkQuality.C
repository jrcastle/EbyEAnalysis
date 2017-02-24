#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;
using namespace hi;

void sysTkQuality(){

  int norder_ = 4;

  //-- Loose
  TFile * fAnaLoose;
  TFile * fUnfLoose;
  TH1D * hObsLoose[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldLoose[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldLoose[NCENT][NEPSymm][NQN][NITER];

  double rmsVnLoose[NCENT][NEPSymm][NQN];
  double meanVnLoose[NCENT][NEPSymm][NQN];
  double stDevVnLoose[NCENT][NEPSymm][NQN];
  double relFluctVnLoose[NCENT][NEPSymm][NQN];
  double vn2Loose[NCENT][NEPSymm][NQN];
  double vn4Loose[NCENT][NEPSymm][NQN];
  double vn6Loose[NCENT][NEPSymm][NQN];
  double vn8Loose[NCENT][NEPSymm][NQN];
  double gamma1expLoose[NCENT][NEPSymm][NQN];
  double vn6vn4Loose[NCENT][NEPSymm][NQN];
  double vn8vn4Loose[NCENT][NEPSymm][NQN];
  double vn8vn6Loose[NCENT][NEPSymm][NQN];

  double rmsLoose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double meanLoose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double stDevLoose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double relFluctLoose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn2Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn4Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn6Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double gamma1expLoose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn6vn4Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8vn4Loose_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8vn6Loose_RatiotoNominal[NCENT][NEPSymm][NQN];

  bool iterCutoffLoose[NCENT][NEPSymm][NQN];

  //-- Tight
  TFile * fAnaTight;
  TFile * fUnfTight;
  TH1D * hObsTight[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldTight[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldTight[NCENT][NEPSymm][NQN][NITER];

  double rmsVnTight[NCENT][NEPSymm][NQN];
  double meanVnTight[NCENT][NEPSymm][NQN];
  double stDevVnTight[NCENT][NEPSymm][NQN];
  double relFluctVnTight[NCENT][NEPSymm][NQN];
  double vn2Tight[NCENT][NEPSymm][NQN];
  double vn4Tight[NCENT][NEPSymm][NQN];
  double vn6Tight[NCENT][NEPSymm][NQN];
  double vn8Tight[NCENT][NEPSymm][NQN];
  double gamma1expTight[NCENT][NEPSymm][NQN];
  double vn6vn4Tight[NCENT][NEPSymm][NQN];
  double vn8vn4Tight[NCENT][NEPSymm][NQN];
  double vn8vn6Tight[NCENT][NEPSymm][NQN];

  double rmsTight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double meanTight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double stDevTight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double relFluctTight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn2Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn4Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn6Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double gamma1expTight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn6vn4Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8vn4Tight_RatiotoNominal[NCENT][NEPSymm][NQN];
  double vn8vn6Tight_RatiotoNominal[NCENT][NEPSymm][NQN];

  bool iterCutoffTight[NCENT][NEPSymm][NQN];

  //-- Nominal
  TFile * fAnaNominal;
  TFile * fUnfNominal;
  TH1D * hObsNominal[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldNominal[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldNominal[NCENT][NEPSymm][NQN][NITER];

  double rmsVnNominal[NCENT][NEPSymm][NQN];
  double meanVnNominal[NCENT][NEPSymm][NQN];
  double stDevVnNominal[NCENT][NEPSymm][NQN];
  double relFluctVnNominal[NCENT][NEPSymm][NQN];
  double vn2Nominal[NCENT][NEPSymm][NQN];
  double vn4Nominal[NCENT][NEPSymm][NQN];
  double vn6Nominal[NCENT][NEPSymm][NQN];
  double vn8Nominal[NCENT][NEPSymm][NQN];
  double gamma1expNominal[NCENT][NEPSymm][NQN];
  double vn6vn4Nominal[NCENT][NEPSymm][NQN];
  double vn8vn4Nominal[NCENT][NEPSymm][NQN];
  double vn8vn6Nominal[NCENT][NEPSymm][NQN];

  bool iterCutoffNominal[NCENT][NEPSymm][NQN];

  //-- RelErrors
  TFile * fRelErr;
  TH1D * relErrRMSVn[NEPSymm][NQN];
  TH1D * relErrMeanVn[NEPSymm][NQN];
  TH1D * relErrStDevVn[NEPSymm][NQN];
  TH1D * relErrRelFluctVn[NEPSymm][NQN];
  TH1D * relErrVn2[NEPSymm][NQN];
  TH1D * relErrVn4[NEPSymm][NQN];
  TH1D * relErrVn6[NEPSymm][NQN];
  TH1D * relErrVn8[NEPSymm][NQN];
  TH1D * relErrGamma1Exp[NEPSymm][NQN];
  TH1D * relErrVn6Vn4[NEPSymm][NQN];
  TH1D * relErrVn8Vn4[NEPSymm][NQN];
  TH1D * relErrVn8Vn6[NEPSymm][NQN];

  TLatex latex;

  //
  //-- MAIN
  // 

  setTDRStyle();
  latex.SetNDC();

  //-- Grab files
  fAnaLoose    = new TFile( "loose/AnalyzerResults/CastleEbyE.root" );
  fUnfLoose    = new TFile( Form("loose/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaTight = new TFile( "tight/AnalyzerResults/CastleEbyE.root" );
  fUnfTight = new TFile( Form("tight/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaNominal   = new TFile( "../../AnalyzerResults/CastleEbyE.root"  );
  fUnfNominal   = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorTkQuality.root", "recreate");
  for(int iEP = 0; iEP < NEPSymm; iEP++){
    if( iEP != EPSymmBin ) continue;
    for(int iqn = 0; iqn < NQN; iqn++){
      fRelErr->cd();
      relErrRMSVn[iEP][iqn]      = new TH1D(Form("relErrRMSVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),      Form("relErrRMSVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),      NCENT, centbinsDefault);
      relErrMeanVn[iEP][iqn]     = new TH1D(Form("relErrMeanVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     Form("relErrMeanVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     NCENT, centbinsDefault);
      relErrStDevVn[iEP][iqn]    = new TH1D(Form("relErrStDevVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    Form("relErrStDevVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn),    NCENT, centbinsDefault);
      relErrRelFluctVn[iEP][iqn] = new TH1D(Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn), NCENT, centbinsDefault);
      relErrVn2[iEP][iqn]        = new TH1D(Form("relErrVn2_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        Form("relErrVn2_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        NCENT, centbinsDefault);
      relErrVn4[iEP][iqn]        = new TH1D(Form("relErrVn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        Form("relErrVn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        NCENT, centbinsDefault);
      relErrVn6[iEP][iqn]        = new TH1D(Form("relErrVn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        Form("relErrVn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        NCENT, centbinsDefault);
      relErrVn8[iEP][iqn]        = new TH1D(Form("relErrVn8_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        Form("relErrVn8_%s_qbin%i", EPSymmNames[iEP].data(), iqn),        NCENT, centbinsDefault);
      relErrGamma1Exp[iEP][iqn]  = new TH1D(Form("relErrGamma1Exp_%s_qbin%i", EPSymmNames[iEP].data(), iqn),  Form("relErrGamma1Exp_%s_qbin%i", EPSymmNames[iEP].data(), iqn),  NCENT, centbinsDefault);
      relErrVn6Vn4[iEP][iqn]     = new TH1D(Form("relErrVn6Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     Form("relErrVn6Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     NCENT, centbinsDefault);
      relErrVn8Vn4[iEP][iqn]     = new TH1D(Form("relErrVn8Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     Form("relErrVn8Vn4_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     NCENT, centbinsDefault);
      relErrVn8Vn6[iEP][iqn]     = new TH1D(Form("relErrVn8Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     Form("relErrVn8Vn6_%s_qbin%i", EPSymmNames[iEP].data(), iqn),     NCENT, centbinsDefault);
    }
  }

  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	//-- -------------------- Loose --------------------
	hObsLoose[icent][iEP][iqn] = (TH1D*) fAnaLoose->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsLoose[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Loose", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffLoose[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldLoose[icent][iEP][iqn][i] = (TH1D*) fUnfLoose->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldLoose[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Loose", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldLoose[icent][iEP][iqn][i] = (TH1D*) fUnfLoose->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldLoose[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Loose", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldLoose[icent][iEP][iqn][i]->Chi2Test(hObsLoose[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffLoose[icent][iEP][iqn] = i;
	
	    FixUnfold( hUnfoldLoose[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldLoose[icent][iEP][iqn][i]);
	    vn2Loose[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Loose[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Loose[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Loose[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expLoose[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Loose[icent][iEP][iqn] == 0 ){
	      vn6vn4Loose[icent][iEP][iqn] = 0;
	      vn8vn4Loose[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Loose[icent][iEP][iqn] = vn6Loose[icent][iEP][iqn] / vn4Loose[icent][iEP][iqn];
	      vn8vn4Loose[icent][iEP][iqn] = vn8Loose[icent][iEP][iqn] / vn4Loose[icent][iEP][iqn];
	    }
	    if( vn6Loose[icent][iEP][iqn] == 0 ) vn8vn6Loose[icent][iEP][iqn] = 0;
	    else                                vn8vn6Loose[icent][iEP][iqn] = vn8Loose[icent][iEP][iqn] / vn6Loose[icent][iEP][iqn];

	    double mean = hUnfoldLoose[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldLoose[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnLoose[icent][iEP][iqn]   = rms;
            meanVnLoose[icent][iEP][iqn]  = mean;
            stDevVnLoose[icent][iEP][iqn] = stdev;
	    relFluctVnLoose[icent][iEP][iqn] = relfluct;

	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffLoose[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldLoose[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldLoose[icent][iEP][iqn][i]);
	    vn2Loose[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Loose[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Loose[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Loose[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expLoose[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Loose[icent][iEP][iqn] == 0 ){
	      vn6vn4Loose[icent][iEP][iqn] = 0;
	      vn8vn4Loose[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Loose[icent][iEP][iqn] = vn6Loose[icent][iEP][iqn] / vn4Loose[icent][iEP][iqn];
	      vn8vn4Loose[icent][iEP][iqn] = vn8Loose[icent][iEP][iqn] / vn4Loose[icent][iEP][iqn];
	    }
	    if( vn6Loose[icent][iEP][iqn] == 0 ) vn8vn6Loose[icent][iEP][iqn] = 0;
	    else                                vn8vn6Loose[icent][iEP][iqn] = vn8Loose[icent][iEP][iqn] / vn6Loose[icent][iEP][iqn];

	    double mean = hUnfoldLoose[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldLoose[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct= stdev / mean;
	    rmsVnLoose[icent][iEP][iqn]   = rms;
	    meanVnLoose[icent][iEP][iqn]  = mean;
	    stDevVnLoose[icent][iEP][iqn] = stdev;
	    relFluctVnLoose[icent][iEP][iqn] = relfluct;
	    break;
	  }
 
	} //-- End iter loop

	//-- -------------------- Tight --------------------
	hObsTight[icent][iEP][iqn] = (TH1D*) fAnaTight->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsTight[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Tight", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffTight[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldTight[icent][iEP][iqn][i] = (TH1D*) fUnfTight->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldTight[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Tight", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldTight[icent][iEP][iqn][i] = (TH1D*) fUnfTight->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldTight[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Tight", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldTight[icent][iEP][iqn][i]->Chi2Test(hObsTight[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffTight[icent][iEP][iqn] = i;
	    
	    FixUnfold( hUnfoldTight[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldTight[icent][iEP][iqn][i]);
	    vn2Tight[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Tight[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Tight[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Tight[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expTight[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Tight[icent][iEP][iqn] == 0 ){
	      vn6vn4Tight[icent][iEP][iqn] = 0;
	      vn8vn4Tight[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Tight[icent][iEP][iqn] = vn6Tight[icent][iEP][iqn] / vn4Tight[icent][iEP][iqn];
	      vn8vn4Tight[icent][iEP][iqn] = vn8Tight[icent][iEP][iqn] / vn4Tight[icent][iEP][iqn];
	    }
	    if( vn6Tight[icent][iEP][iqn] == 0 ) vn8vn6Tight[icent][iEP][iqn] = 0;
	    else                                vn8vn6Tight[icent][iEP][iqn] = vn8Tight[icent][iEP][iqn] / vn6Tight[icent][iEP][iqn];

	    double mean = hUnfoldTight[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldTight[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnTight[icent][iEP][iqn]   = rms;
            meanVnTight[icent][iEP][iqn]  = mean;
            stDevVnTight[icent][iEP][iqn] = stdev;
	    relFluctVnTight[icent][iEP][iqn] = relfluct;
	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffTight[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldTight[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldTight[icent][iEP][iqn][i]);
	    vn2Tight[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Tight[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Tight[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Tight[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expTight[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Tight[icent][iEP][iqn] == 0 ){
	      vn6vn4Tight[icent][iEP][iqn] = 0;
	      vn8vn4Tight[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Tight[icent][iEP][iqn] = vn6Tight[icent][iEP][iqn] / vn4Tight[icent][iEP][iqn];
	      vn8vn4Tight[icent][iEP][iqn] = vn8Tight[icent][iEP][iqn] / vn4Tight[icent][iEP][iqn];
	    }
	    if( vn6Tight[icent][iEP][iqn] == 0 ) vn8vn6Tight[icent][iEP][iqn] = 0;
	    else                                vn8vn6Tight[icent][iEP][iqn] = vn8Tight[icent][iEP][iqn] / vn6Tight[icent][iEP][iqn];

	    double mean = hUnfoldTight[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldTight[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
	    rmsVnTight[icent][iEP][iqn]   = rms;
	    meanVnTight[icent][iEP][iqn]  = mean;
	    stDevVnTight[icent][iEP][iqn] = stdev;
	    relFluctVnTight[icent][iEP][iqn] = relfluct;

	    break;
	  }
	} //-- End iter loop

	//-- -------------------- Nominal --------------------
	hObsNominal[icent][iEP][iqn] = (TH1D*) fAnaNominal->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsNominal[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Nominal", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffNominal[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldNominal[icent][iEP][iqn][i] = (TH1D*) fUnfNominal->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldNominal[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Nominal", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldNominal[icent][iEP][iqn][i] = (TH1D*) fUnfNominal->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldNominal[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Nominal", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldNominal[icent][iEP][iqn][i]->Chi2Test(hObsNominal[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffNominal[icent][iEP][iqn] = i;
	    
	    FixUnfold( hUnfoldNominal[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldNominal[icent][iEP][iqn][i]);
	    vn2Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expNominal[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Nominal[icent][iEP][iqn] == 0 ){
	      vn6vn4Nominal[icent][iEP][iqn] = 0;
	      vn8vn4Nominal[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Nominal[icent][iEP][iqn] = vn6Nominal[icent][iEP][iqn] / vn4Nominal[icent][iEP][iqn];
	      vn8vn4Nominal[icent][iEP][iqn] = vn8Nominal[icent][iEP][iqn] / vn4Nominal[icent][iEP][iqn];
	    }
	    if( vn6Nominal[icent][iEP][iqn] == 0 ) vn8vn6Nominal[icent][iEP][iqn] = 0;
	    else                                vn8vn6Nominal[icent][iEP][iqn] = vn8Nominal[icent][iEP][iqn] / vn6Nominal[icent][iEP][iqn];

	    double mean = hUnfoldNominal[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldNominal[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnNominal[icent][iEP][iqn]   = rms;
            meanVnNominal[icent][iEP][iqn]  = mean;
            stDevVnNominal[icent][iEP][iqn] = stdev;
	    relFluctVnNominal[icent][iEP][iqn] = relfluct;

	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffNominal[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldNominal[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldNominal[icent][iEP][iqn][i]);
	    vn2Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Nominal[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expNominal[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Nominal[icent][iEP][iqn] == 0 ){
	      vn6vn4Nominal[icent][iEP][iqn] = 0;
	      vn8vn4Nominal[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Nominal[icent][iEP][iqn] = vn6Nominal[icent][iEP][iqn] / vn4Nominal[icent][iEP][iqn];
	      vn8vn4Nominal[icent][iEP][iqn] = vn8Nominal[icent][iEP][iqn] / vn4Nominal[icent][iEP][iqn];
	    }
	    if( vn6Nominal[icent][iEP][iqn] == 0 ) vn8vn6Nominal[icent][iEP][iqn] = 0;
	    else                                vn8vn6Nominal[icent][iEP][iqn] = vn8Nominal[icent][iEP][iqn] / vn6Nominal[icent][iEP][iqn];

	    double mean = hUnfoldNominal[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldNominal[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct= stdev/ mean;
	    rmsVnNominal[icent][iEP][iqn]   = rms;
	    meanVnNominal[icent][iEP][iqn]  = mean;
	    stDevVnNominal[icent][iEP][iqn] = stdev;
	    relFluctVnNominal[icent][iEP][iqn] = relfluct;

	    break;
	  }
 
	} //-- End iter loop

      }  //-- End qn loop
    } //-- End EP loop
  } //-- End Cent loop

  //-- Calculate ratios to the Nominal case
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	double rmsvnL      = rmsVnLoose[icent][iEP][iqn];
	double meanvnL     = meanVnLoose[icent][iEP][iqn];
	double stdevvnL    = stDevVnLoose[icent][iEP][iqn];
	double relfluctvnL = relFluctVnLoose[icent][iEP][iqn];
	double vn2L        = vn2Loose[icent][iEP][iqn];
	double vn4L        = vn4Loose[icent][iEP][iqn];
	double vn6L        = vn6Loose[icent][iEP][iqn];
	double vn8L        = vn8Loose[icent][iEP][iqn];
	double gamma1expL  = gamma1expLoose[icent][iEP][iqn];
	double vn6vn4L     = vn6vn4Loose[icent][iEP][iqn];
	double vn8vn4L     = vn8vn4Loose[icent][iEP][iqn];
	double vn8vn6L     = vn8vn6Loose[icent][iEP][iqn];

	double rmsvnT      = rmsVnTight[icent][iEP][iqn];
        double meanvnT     = meanVnTight[icent][iEP][iqn];
        double stdevvnT    = stDevVnTight[icent][iEP][iqn];
	double relfluctvnT = relFluctVnTight[icent][iEP][iqn];
	double vn2T        = vn2Tight[icent][iEP][iqn];
	double vn4T        = vn4Tight[icent][iEP][iqn];
	double vn6T        = vn6Tight[icent][iEP][iqn];
	double vn8T        = vn8Tight[icent][iEP][iqn];
	double gamma1expT  = gamma1expTight[icent][iEP][iqn];
	double vn6vn4T     = vn6vn4Tight[icent][iEP][iqn];
	double vn8vn4T     = vn8vn4Tight[icent][iEP][iqn];
	double vn8vn6T     = vn8vn6Tight[icent][iEP][iqn];

	double rmsvnN      = rmsVnNominal[icent][iEP][iqn];
        double meanvnN     = meanVnNominal[icent][iEP][iqn];
        double stdevvnN    = stDevVnNominal[icent][iEP][iqn];
	double relfluctvnN = relFluctVnNominal[icent][iEP][iqn];
	double vn2N        = vn2Nominal[icent][iEP][iqn];
	double vn4N        = vn4Nominal[icent][iEP][iqn];
	double vn6N        = vn6Nominal[icent][iEP][iqn];
	double vn8N        = vn8Nominal[icent][iEP][iqn];
	double gamma1expN  = gamma1expNominal[icent][iEP][iqn];
	double vn6vn4N     = vn6vn4Nominal[icent][iEP][iqn];
	double vn8vn4N     = vn8vn4Nominal[icent][iEP][iqn];
	double vn8vn6N     = vn8vn6Nominal[icent][iEP][iqn];

	double rrmsL       = rmsvnL / rmsvnN;
	double rmeanL      = meanvnL / meanvnN;
	double rstdevL     = stdevvnL / stdevvnN;
	double rrelfluctL  = relfluctvnL / relfluctvnN;
	double r2L         = vn2L / vn2N;
	double r4L         = vn4L / vn4N;
	double r6L         = vn6L / vn6N;
	double r8L         = vn8L / vn8N;
	double rgamma1expL = gamma1expL / gamma1expN;
	double r64L        = vn6vn4L / vn6vn4N;
	double r84L        = vn8vn4L / vn8vn4N;
	double r86L        = vn8vn6L / vn8vn6N;

	double rrmsT       = rmsvnT / rmsvnN;
        double rmeanT      = meanvnT / meanvnN;
        double rstdevT     = stdevvnT / stdevvnN;
	double rrelfluctT  = relfluctvnT / relfluctvnN;
	double r2T         = vn2T / vn2N;
	double r4T         = vn4T / vn4N;
	double r6T         = vn6T / vn6N;
	double r8T         = vn8T / vn8N;
	double rgamma1expT = gamma1expT / gamma1expN;
	double r64T        = vn6vn4T / vn6vn4N;
	double r84T        = vn8vn4T / vn8vn4N;
	double r86T        = vn8vn6T / vn8vn6N;

	rmsLoose_RatiotoNominal[icent][iEP][iqn]       = rrmsL;
	meanLoose_RatiotoNominal[icent][iEP][iqn]      = rmeanL;
	stDevLoose_RatiotoNominal[icent][iEP][iqn]     = rstdevL;
	relFluctLoose_RatiotoNominal[icent][iEP][iqn]  = rrelfluctL;
	vn2Loose_RatiotoNominal[icent][iEP][iqn]       = r2L;
	vn4Loose_RatiotoNominal[icent][iEP][iqn]       = r4L;
	vn6Loose_RatiotoNominal[icent][iEP][iqn]       = r6L;
	vn8Loose_RatiotoNominal[icent][iEP][iqn]       = r8L;
	gamma1expLoose_RatiotoNominal[icent][iEP][iqn] = rgamma1expL;
	vn6vn4Loose_RatiotoNominal[icent][iEP][iqn]    = r64L;
	vn8vn4Loose_RatiotoNominal[icent][iEP][iqn]    = r84L;
	vn8vn6Loose_RatiotoNominal[icent][iEP][iqn]    = r86L;

	rmsTight_RatiotoNominal[icent][iEP][iqn]       = rrmsT;
        meanTight_RatiotoNominal[icent][iEP][iqn]      = rmeanT;
        stDevTight_RatiotoNominal[icent][iEP][iqn]     = rstdevT;
	relFluctTight_RatiotoNominal[icent][iEP][iqn]  = rrelfluctT;
	vn2Tight_RatiotoNominal[icent][iEP][iqn]       = r2T;
	vn4Tight_RatiotoNominal[icent][iEP][iqn]       = r4T;
	vn6Tight_RatiotoNominal[icent][iEP][iqn]       = r6T;
	vn8Tight_RatiotoNominal[icent][iEP][iqn]       = r8T;
	gamma1expTight_RatiotoNominal[icent][iEP][iqn] = rgamma1expT;
	vn6vn4Tight_RatiotoNominal[icent][iEP][iqn]    = r64T;
	vn8vn4Tight_RatiotoNominal[icent][iEP][iqn]    = r84T;
	vn8vn6Tight_RatiotoNominal[icent][iEP][iqn]    = r86T;

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	double pctRMSVn_L      = fabs(rmsVnLoose[icent][iEP][iqn] - rmsVnNominal[icent][iEP][iqn]) / fabs(rmsVnNominal[icent][iEP][iqn]);
	double pctMeanVn_L     = fabs(meanVnLoose[icent][iEP][iqn] - meanVnNominal[icent][iEP][iqn]) / fabs(meanVnNominal[icent][iEP][iqn]);
	double pctStDevVn_L    = fabs(stDevVnLoose[icent][iEP][iqn] - stDevVnNominal[icent][iEP][iqn]) / fabs(stDevVnNominal[icent][iEP][iqn]);
	double pctRelFluctVn_L = fabs(relFluctVnLoose[icent][iEP][iqn] - relFluctVnNominal[icent][iEP][iqn]) / fabs(relFluctVnNominal[icent][iEP][iqn]);
	double pctVn2_L        = fabs(vn2Loose[icent][iEP][iqn] - vn2Nominal[icent][iEP][iqn]) / fabs(vn2Nominal[icent][iEP][iqn]);
	double pctVn4_L        = fabs(vn4Loose[icent][iEP][iqn] - vn4Nominal[icent][iEP][iqn]) / fabs(vn4Nominal[icent][iEP][iqn]);
	double pctVn6_L        = fabs(vn6Loose[icent][iEP][iqn] - vn6Nominal[icent][iEP][iqn]) / fabs(vn6Nominal[icent][iEP][iqn]);
	double pctVn8_L        = fabs(vn8Loose[icent][iEP][iqn] - vn8Nominal[icent][iEP][iqn]) / fabs(vn8Nominal[icent][iEP][iqn]);
	double pctGamma1Exp_L  = fabs(gamma1expLoose[icent][iEP][iqn] - gamma1expNominal[icent][iEP][iqn]) / fabs(gamma1expNominal[icent][iEP][iqn]);
	double pctVn6Vn4_L     = fabs(vn6vn4Loose[icent][iEP][iqn] - vn6vn4Nominal[icent][iEP][iqn]) / fabs(vn6vn4Nominal[icent][iEP][iqn]);
	double pctVn8Vn4_L     = fabs(vn8vn4Loose[icent][iEP][iqn] - vn8vn4Nominal[icent][iEP][iqn]) / fabs(vn8vn4Nominal[icent][iEP][iqn]);
	double pctVn8Vn6_L     = fabs(vn8vn6Loose[icent][iEP][iqn] - vn8vn6Nominal[icent][iEP][iqn]) / fabs(vn8vn6Nominal[icent][iEP][iqn]);

	double pctRMSVn_T      = fabs(rmsVnTight[icent][iEP][iqn] - rmsVnNominal[icent][iEP][iqn]) / fabs(rmsVnNominal[icent][iEP][iqn]);
	double pctMeanVn_T     = fabs(meanVnTight[icent][iEP][iqn] - meanVnNominal[icent][iEP][iqn]) / fabs(meanVnNominal[icent][iEP][iqn]);
	double pctStDevVn_T    = fabs(stDevVnTight[icent][iEP][iqn] - stDevVnNominal[icent][iEP][iqn]) / fabs(stDevVnNominal[icent][iEP][iqn]);
	double pctRelFluctVn_T = fabs(relFluctVnTight[icent][iEP][iqn] - relFluctVnNominal[icent][iEP][iqn]) / fabs(relFluctVnNominal[icent][iEP][iqn]);
	double pctVn2_T        = fabs(vn2Tight[icent][iEP][iqn] - vn2Nominal[icent][iEP][iqn]) / fabs(vn2Nominal[icent][iEP][iqn]);
	double pctVn4_T        = fabs(vn4Tight[icent][iEP][iqn] - vn4Nominal[icent][iEP][iqn]) / fabs(vn4Nominal[icent][iEP][iqn]);
	double pctVn6_T        = fabs(vn6Tight[icent][iEP][iqn] - vn6Nominal[icent][iEP][iqn]) / fabs(vn6Nominal[icent][iEP][iqn]);
	double pctVn8_T        = fabs(vn8Tight[icent][iEP][iqn] - vn8Nominal[icent][iEP][iqn]) / fabs(vn8Nominal[icent][iEP][iqn]);
	double pctGamma1Exp_T  = fabs(gamma1expTight[icent][iEP][iqn] - gamma1expNominal[icent][iEP][iqn]) / fabs(gamma1expNominal[icent][iEP][iqn]);
	double pctVn6Vn4_T     = fabs(vn6vn4Tight[icent][iEP][iqn] - vn6vn4Nominal[icent][iEP][iqn]) / fabs(vn6vn4Nominal[icent][iEP][iqn]);
	double pctVn8Vn4_T     = fabs(vn8vn4Tight[icent][iEP][iqn] - vn8vn4Nominal[icent][iEP][iqn]) / fabs(vn8vn4Nominal[icent][iEP][iqn]);
	double pctVn8Vn6_T     = fabs(vn8vn6Tight[icent][iEP][iqn] - vn8vn6Nominal[icent][iEP][iqn]) / fabs(vn8vn6Nominal[icent][iEP][iqn]);

	double reRMSVn      = max(pctRMSVn_L, pctRMSVn_T);
	double reMeanVn     = max(pctMeanVn_L, pctMeanVn_T);
	double reStDevVn    = max(pctStDevVn_L, pctStDevVn_T);
	double reRelFluctVn = max(pctRelFluctVn_L, pctRelFluctVn_T);
	double reVn2        = max(pctVn2_L, pctVn2_T);
	double reVn4        = max(pctVn4_L, pctVn4_T);
	double reVn6        = max(pctVn6_L, pctVn6_T);
	double reVn8        = max(pctVn8_L, pctVn8_T);
	double reGamma1Exp  = max(pctGamma1Exp_L, pctGamma1Exp_T);
	double reVn6Vn4     = max(pctVn6Vn4_L, pctVn6Vn4_T);
	double reVn8Vn4     = max(pctVn8Vn4_L, pctVn8Vn4_T);
	double reVn8Vn6     = max(pctVn8Vn6_L, pctVn8Vn6_T);

	if( TMath::IsNaN(reVn2) || std::isinf(reVn2) )             reVn2 = 0;
	if( TMath::IsNaN(reVn4) || std::isinf(reVn4) )             reVn4 = 0;
	if( TMath::IsNaN(reVn6) || std::isinf(reVn6) )             reVn6 = 0;
	if( TMath::IsNaN(reVn8) || std::isinf(reVn8) )             reVn8 = 0;
	if( TMath::IsNaN(reGamma1Exp) || std::isinf(reGamma1Exp) ) reGamma1Exp = 0;
	if( TMath::IsNaN(reVn6Vn4) || std::isinf(reVn6Vn4) )       reVn6Vn4 = 0;
	if( TMath::IsNaN(reVn8Vn4) || std::isinf(reVn8Vn4) )       reVn8Vn4 = 0;
	if( TMath::IsNaN(reVn8Vn6) || std::isinf(reVn8Vn6) )       reVn8Vn6 = 0;

	relErrRMSVn[iEP][iqn]->SetBinContent(icent+1, reRMSVn);
	relErrMeanVn[iEP][iqn]->SetBinContent(icent+1, reMeanVn);
	relErrStDevVn[iEP][iqn]->SetBinContent(icent+1, reStDevVn);
	relErrRelFluctVn[iEP][iqn]->SetBinContent(icent+1, reRelFluctVn);
	relErrVn2[iEP][iqn]->SetBinContent(icent+1, reVn2);
	relErrVn4[iEP][iqn]->SetBinContent(icent+1, reVn4);
	relErrVn6[iEP][iqn]->SetBinContent(icent+1, reVn6);
	relErrVn8[iEP][iqn]->SetBinContent(icent+1, reVn8);
	relErrGamma1Exp[iEP][iqn]->SetBinContent(icent+1, reGamma1Exp);
	relErrVn6Vn4[iEP][iqn]->SetBinContent(icent+1, reVn6Vn4);
	relErrVn8Vn4[iEP][iqn]->SetBinContent(icent+1, reVn8Vn4);
	relErrVn8Vn6[iEP][iqn]->SetBinContent(icent+1, reVn8Vn6);

      } //-- End qn loop
    } //-- End EP loop
  } //-- End Cent loop

  fRelErr->Write();

}

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

void sysVtxCut(){

  int norder_ = 2;

  //-- Vertex < 3.0
  TFile * fAnaVtx3;
  TFile * fUnfVtx3;
  TH1D * hObsVtx3[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldVtx3[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldVtx3[NCENT][NEPSymm][NQN][NITER];

  double rmsVnVtx3[NCENT][NEPSymm][NQN];
  double meanVnVtx3[NCENT][NEPSymm][NQN];
  double stDevVnVtx3[NCENT][NEPSymm][NQN];
  double relFluctVnVtx3[NCENT][NEPSymm][NQN];
  double vn2Vtx3[NCENT][NEPSymm][NQN];
  double vn4Vtx3[NCENT][NEPSymm][NQN];
  double vn6Vtx3[NCENT][NEPSymm][NQN];
  double vn8Vtx3[NCENT][NEPSymm][NQN];
  double gamma1expVtx3[NCENT][NEPSymm][NQN];
  double vn6vn4Vtx3[NCENT][NEPSymm][NQN];
  double vn8vn4Vtx3[NCENT][NEPSymm][NQN];
  double vn8vn6Vtx3[NCENT][NEPSymm][NQN];

  double rmsVtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double meanVtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double stDevVtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double relFluctVtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn2Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn4Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn6Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double gamma1expVtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn6vn4Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8vn4Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8vn6Vtx3_RatiotoVtx15[NCENT][NEPSymm][NQN];

  bool iterCutoffVtx3[NCENT][NEPSymm][NQN];

  //-- 3.0 < Vertex < 15.0
  TFile * fAnaVtx3_15;
  TFile * fUnfVtx3_15;
  TH1D * hObsVtx3_15[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldVtx3_15[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldVtx3_15[NCENT][NEPSymm][NQN][NITER];

  double rmsVnVtx3_15[NCENT][NEPSymm][NQN];
  double meanVnVtx3_15[NCENT][NEPSymm][NQN];
  double stDevVnVtx3_15[NCENT][NEPSymm][NQN];
  double relFluctVnVtx3_15[NCENT][NEPSymm][NQN];
  double vn2Vtx3_15[NCENT][NEPSymm][NQN];
  double vn4Vtx3_15[NCENT][NEPSymm][NQN];
  double vn6Vtx3_15[NCENT][NEPSymm][NQN];
  double vn8Vtx3_15[NCENT][NEPSymm][NQN];
  double gamma1expVtx3_15[NCENT][NEPSymm][NQN];
  double vn6vn4Vtx3_15[NCENT][NEPSymm][NQN];
  double vn8vn4Vtx3_15[NCENT][NEPSymm][NQN];
  double vn8vn6Vtx3_15[NCENT][NEPSymm][NQN];

  double rmsVtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double meanVtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double stDevVtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double relFluctVtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn2Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn4Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn6Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double gamma1expVtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn6vn4Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8vn4Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];
  double vn8vn6Vtx3_15_RatiotoVtx15[NCENT][NEPSymm][NQN];

  bool iterCutoffVtx3_15[NCENT][NEPSymm][NQN];

  //-- Vertex < 15.0
  TFile * fAnaVtx15;
  TFile * fUnfVtx15;
  TH1D * hObsVtx15[NCENT][NEPSymm][NQN];
  TH1D * hUnfoldVtx15[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldVtx15[NCENT][NEPSymm][NQN][NITER];

  double rmsVnVtx15[NCENT][NEPSymm][NQN];
  double meanVnVtx15[NCENT][NEPSymm][NQN];
  double stDevVnVtx15[NCENT][NEPSymm][NQN];
  double relFluctVnVtx15[NCENT][NEPSymm][NQN];
  double vn2Vtx15[NCENT][NEPSymm][NQN];
  double vn4Vtx15[NCENT][NEPSymm][NQN];
  double vn6Vtx15[NCENT][NEPSymm][NQN];
  double vn8Vtx15[NCENT][NEPSymm][NQN];
  double gamma1expVtx15[NCENT][NEPSymm][NQN];
  double vn6vn4Vtx15[NCENT][NEPSymm][NQN];
  double vn8vn4Vtx15[NCENT][NEPSymm][NQN];
  double vn8vn6Vtx15[NCENT][NEPSymm][NQN];

  bool iterCutoffVtx15[NCENT][NEPSymm][NQN];

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
  fAnaVtx3    = new TFile( "vtx_leq_3/AnalyzerResults/CastleEbyE.root" );
  fUnfVtx3    = new TFile( Form("vtx_leq_3/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaVtx3_15 = new TFile( "vtx3_15/AnalyzerResults/CastleEbyE.root" );
  fUnfVtx3_15 = new TFile( Form("vtx3_15/UnfoldResults/dataResp/data%i.root", norder_) );

  fAnaVtx15   = new TFile( "../../AnalyzerResults/CastleEbyE.root"  );
  fUnfVtx15   = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorVtx.root", "recreate");
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

	//-- -------------------- Vertex < 3.0 --------------------
	hObsVtx3[icent][iEP][iqn] = (TH1D*) fAnaVtx3->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsVtx3[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Vtx3", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffVtx3[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldVtx3[icent][iEP][iqn][i] = (TH1D*) fUnfVtx3->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldVtx3[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Vtx3", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldVtx3[icent][iEP][iqn][i] = (TH1D*) fUnfVtx3->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldVtx3[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Vtx3", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldVtx3[icent][iEP][iqn][i]->Chi2Test(hObsVtx3[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffVtx3[icent][iEP][iqn] = i;
	
	    FixUnfold( hUnfoldVtx3[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx3[icent][iEP][iqn][i]);
	    vn2Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx3[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx3[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx3[icent][iEP][iqn] = 0;
	      vn8vn4Vtx3[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx3[icent][iEP][iqn] = vn6Vtx3[icent][iEP][iqn] / vn4Vtx3[icent][iEP][iqn];
	      vn8vn4Vtx3[icent][iEP][iqn] = vn8Vtx3[icent][iEP][iqn] / vn4Vtx3[icent][iEP][iqn];
	    }
	    if( vn6Vtx3[icent][iEP][iqn] == 0 ) vn8vn6Vtx3[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx3[icent][iEP][iqn] = vn8Vtx3[icent][iEP][iqn] / vn6Vtx3[icent][iEP][iqn];

	    double mean = hUnfoldVtx3[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldVtx3[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnVtx3[icent][iEP][iqn]   = rms;
            meanVnVtx3[icent][iEP][iqn]  = mean;
            stDevVnVtx3[icent][iEP][iqn] = stdev;
	    relFluctVnVtx3[icent][iEP][iqn] = relfluct;

	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffVtx3[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldVtx3[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx3[icent][iEP][iqn][i]);
	    vn2Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx3[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx3[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx3[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx3[icent][iEP][iqn] = 0;
	      vn8vn4Vtx3[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx3[icent][iEP][iqn] = vn6Vtx3[icent][iEP][iqn] / vn4Vtx3[icent][iEP][iqn];
	      vn8vn4Vtx3[icent][iEP][iqn] = vn8Vtx3[icent][iEP][iqn] / vn4Vtx3[icent][iEP][iqn];
	    }
	    if( vn6Vtx3[icent][iEP][iqn] == 0 ) vn8vn6Vtx3[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx3[icent][iEP][iqn] = vn8Vtx3[icent][iEP][iqn] / vn6Vtx3[icent][iEP][iqn];

	    double mean = hUnfoldVtx3[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldVtx3[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
	    rmsVnVtx3[icent][iEP][iqn]   = rms;
	    meanVnVtx3[icent][iEP][iqn]  = mean;
	    stDevVnVtx3[icent][iEP][iqn] = stdev;
	    relFluctVnVtx3[icent][iEP][iqn] = relfluct;
	    break;
	  }
 
	} //-- End iter loop

	//-- -------------------- 3.0 < Vertex < 15.0 --------------------
	hObsVtx3_15[icent][iEP][iqn] = (TH1D*) fAnaVtx3_15->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsVtx3_15[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Vtx3_15", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffVtx3_15[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldVtx3_15[icent][iEP][iqn][i] = (TH1D*) fUnfVtx3_15->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldVtx3_15[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Vtx3_15", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldVtx3_15[icent][iEP][iqn][i] = (TH1D*) fUnfVtx3_15->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldVtx3_15[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Vtx3_15", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldVtx3_15[icent][iEP][iqn][i]->Chi2Test(hObsVtx3_15[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffVtx3_15[icent][iEP][iqn] = i;
	    
	    FixUnfold( hUnfoldVtx3_15[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx3_15[icent][iEP][iqn][i]);
	    vn2Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx3_15[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx3_15[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx3_15[icent][iEP][iqn] = 0;
	      vn8vn4Vtx3_15[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx3_15[icent][iEP][iqn] = vn6Vtx3_15[icent][iEP][iqn] / vn4Vtx3_15[icent][iEP][iqn];
	      vn8vn4Vtx3_15[icent][iEP][iqn] = vn8Vtx3_15[icent][iEP][iqn] / vn4Vtx3_15[icent][iEP][iqn];
	    }
	    if( vn6Vtx3_15[icent][iEP][iqn] == 0 ) vn8vn6Vtx3_15[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx3_15[icent][iEP][iqn] = vn8Vtx3_15[icent][iEP][iqn] / vn6Vtx3_15[icent][iEP][iqn];

	    double mean = hUnfoldVtx3_15[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldVtx3_15[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnVtx3_15[icent][iEP][iqn]   = rms;
            meanVnVtx3_15[icent][iEP][iqn]  = mean;
            stDevVnVtx3_15[icent][iEP][iqn] = stdev;
	    relFluctVnVtx3_15[icent][iEP][iqn] = relfluct;

	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffVtx3_15[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldVtx3_15[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx3_15[icent][iEP][iqn][i]);
	    vn2Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx3_15[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx3_15[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx3_15[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx3_15[icent][iEP][iqn] = 0;
	      vn8vn4Vtx3_15[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx3_15[icent][iEP][iqn] = vn6Vtx3_15[icent][iEP][iqn] / vn4Vtx3_15[icent][iEP][iqn];
	      vn8vn4Vtx3_15[icent][iEP][iqn] = vn8Vtx3_15[icent][iEP][iqn] / vn4Vtx3_15[icent][iEP][iqn];
	    }
	    if( vn6Vtx3_15[icent][iEP][iqn] == 0 ) vn8vn6Vtx3_15[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx3_15[icent][iEP][iqn] = vn8Vtx3_15[icent][iEP][iqn] / vn6Vtx3_15[icent][iEP][iqn];

	    double mean = hUnfoldVtx3_15[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldVtx3_15[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
	    rmsVnVtx3_15[icent][iEP][iqn]   = rms;
	    meanVnVtx3_15[icent][iEP][iqn]  = mean;
	    stDevVnVtx3_15[icent][iEP][iqn] = stdev;
	    relFluctVnVtx3_15[icent][iEP][iqn] = relfluct;

	    break;
	  }
	} //-- End iter loop

	//-- -------------------- Vertex < 15.0 --------------------
	hObsVtx15[icent][iEP][iqn] = (TH1D*) fAnaVtx15->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObsVtx15[icent][iEP][iqn]->SetName( Form("hVnFull_%s_c%i_qbin%i_Vtx15", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Iter loop
	iterCutoffVtx15[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfoldVtx15[icent][iEP][iqn][i] = (TH1D*) fUnfVtx15->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldVtx15[icent][iEP][iqn][i]->SetName( Form("hreco%i_%s_c%i_qbin%i_Vtx15", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefoldVtx15[icent][iEP][iqn][i] = (TH1D*) fUnfVtx15->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldVtx15[icent][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i_Vtx15", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi squares
	  double chi2NDF_Refold = hRefoldVtx15[icent][iEP][iqn][i]->Chi2Test(hObsVtx15[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold < 1.2 ){
	    iterCutoffVtx15[icent][iEP][iqn] = i;
	    
	    FixUnfold( hUnfoldVtx15[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx15[icent][iEP][iqn][i]);
	    vn2Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx15[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx15[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx15[icent][iEP][iqn] = 0;
	      vn8vn4Vtx15[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx15[icent][iEP][iqn] = vn6Vtx15[icent][iEP][iqn] / vn4Vtx15[icent][iEP][iqn];
	      vn8vn4Vtx15[icent][iEP][iqn] = vn8Vtx15[icent][iEP][iqn] / vn4Vtx15[icent][iEP][iqn];
	    }
	    if( vn6Vtx15[icent][iEP][iqn] == 0 ) vn8vn6Vtx15[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx15[icent][iEP][iqn] = vn8Vtx15[icent][iEP][iqn] / vn6Vtx15[icent][iEP][iqn];

	    double mean = hUnfoldVtx15[icent][iEP][iqn][i]->GetMean();
            double stdev = hUnfoldVtx15[icent][iEP][iqn][i]->GetRMS();
            double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
            rmsVnVtx15[icent][iEP][iqn]   = rms;
            meanVnVtx15[icent][iEP][iqn]  = mean;
            stDevVnVtx15[icent][iEP][iqn] = stdev;
	    relFluctVnVtx15[icent][iEP][iqn] = relfluct;

	    break;
	  }
	  if( i == NITER - 1 ){
	    iterCutoffVtx15[icent][iEP][iqn] = i;

	    FixUnfold( hUnfoldVtx15[icent][iEP][iqn][i] );
	    EbyECumu cumu(hUnfoldVtx15[icent][iEP][iqn][i]);
	    vn2Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn2();
	    vn4Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn4();
	    vn6Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn6();
	    vn8Vtx15[icent][iEP][iqn]  = cumu.GetCumu_vn8();
	    gamma1expVtx15[icent][iEP][iqn] = cumu.GetGamma1Exp();

	    if( vn4Vtx15[icent][iEP][iqn] == 0 ){
	      vn6vn4Vtx15[icent][iEP][iqn] = 0;
	      vn8vn4Vtx15[icent][iEP][iqn] = 0;
	    }
	    else{
	      vn6vn4Vtx15[icent][iEP][iqn] = vn6Vtx15[icent][iEP][iqn] / vn4Vtx15[icent][iEP][iqn];
	      vn8vn4Vtx15[icent][iEP][iqn] = vn8Vtx15[icent][iEP][iqn] / vn4Vtx15[icent][iEP][iqn];
	    }
	    if( vn6Vtx15[icent][iEP][iqn] == 0 ) vn8vn6Vtx15[icent][iEP][iqn] = 0;
	    else                                vn8vn6Vtx15[icent][iEP][iqn] = vn8Vtx15[icent][iEP][iqn] / vn6Vtx15[icent][iEP][iqn];

	    double mean = hUnfoldVtx15[icent][iEP][iqn][i]->GetMean();
	    double stdev = hUnfoldVtx15[icent][iEP][iqn][i]->GetRMS();
	    double rms = sqrt( mean*mean + stdev*stdev );
	    double relfluct = stdev / mean;
	    rmsVnVtx15[icent][iEP][iqn]   = rms;
	    meanVnVtx15[icent][iEP][iqn]  = mean;
	    stDevVnVtx15[icent][iEP][iqn] = stdev;
	    relFluctVnVtx15[icent][iEP][iqn] = relfluct;
	    break;
	  }
 
	} //-- End iter loop

      }  //-- End qn loop
    } //-- End EP loop
  } //-- End Cent loop

  //-- Calculate ratios to the vtx < 15 case
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	double rmsvn3      = rmsVnVtx3[icent][iEP][iqn];
	double meanvn3     = meanVnVtx3[icent][iEP][iqn];
	double stdevvn3    = stDevVnVtx3[icent][iEP][iqn];
	double relfluctvn3 = relFluctVnVtx3[icent][iEP][iqn];
	double vn23        = vn2Vtx3[icent][iEP][iqn];
	double vn43        = vn4Vtx3[icent][iEP][iqn];
	double vn63        = vn6Vtx3[icent][iEP][iqn];
	double vn83        = vn8Vtx3[icent][iEP][iqn];
	double gamma1exp3  = gamma1expVtx3[icent][iEP][iqn];
	double vn6vn43     = vn6vn4Vtx3[icent][iEP][iqn];
	double vn8vn43     = vn8vn4Vtx3[icent][iEP][iqn];
	double vn8vn63     = vn8vn6Vtx3[icent][iEP][iqn];

	double rmsvn3_15      = rmsVnVtx3_15[icent][iEP][iqn];
        double meanvn3_15     = meanVnVtx3_15[icent][iEP][iqn];
        double stdevvn3_15    = stDevVnVtx3_15[icent][iEP][iqn];
	double relfluctvn3_15 = relFluctVnVtx3_15[icent][iEP][iqn];
	double vn23_15        = vn2Vtx3_15[icent][iEP][iqn];
	double vn43_15        = vn4Vtx3_15[icent][iEP][iqn];
	double vn63_15        = vn6Vtx3_15[icent][iEP][iqn];
	double vn83_15        = vn8Vtx3_15[icent][iEP][iqn];
	double gamma1exp3_15  = gamma1expVtx3_15[icent][iEP][iqn];
	double vn6vn43_15     = vn6vn4Vtx3_15[icent][iEP][iqn];
	double vn8vn43_15     = vn8vn4Vtx3_15[icent][iEP][iqn];
	double vn8vn63_15     = vn8vn6Vtx3_15[icent][iEP][iqn];

	double rmsvn15      = rmsVnVtx15[icent][iEP][iqn];
        double meanvn15     = meanVnVtx15[icent][iEP][iqn];
        double stdevvn15    = stDevVnVtx15[icent][iEP][iqn];
	double relfluctvn15 = relFluctVnVtx15[icent][iEP][iqn];
	double vn215        = vn2Vtx15[icent][iEP][iqn];
	double vn415        = vn4Vtx15[icent][iEP][iqn];
	double vn615        = vn6Vtx15[icent][iEP][iqn];
	double vn815        = vn8Vtx15[icent][iEP][iqn];
	double gamma1exp15  = gamma1expVtx15[icent][iEP][iqn];
	double vn6vn415     = vn6vn4Vtx15[icent][iEP][iqn];
	double vn8vn415     = vn8vn4Vtx15[icent][iEP][iqn];
	double vn8vn615     = vn8vn6Vtx15[icent][iEP][iqn];

	double rrms3       = rmsvn3 / rmsvn15;
	double rmean3      = meanvn3 / meanvn15;
	double rstdev3     = stdevvn3 / stdevvn15;
	double rrelfluct3  = relfluctvn3 / relfluctvn15;
	double r23         = vn23 / vn215;
	double r43         = vn43 / vn415;
	double r63         = vn63 / vn615;
	double r83         = vn83 / vn815;
	double rgamma1exp3 = gamma1exp3 / gamma1exp15;
	double r643        = vn6vn43 / vn6vn415;
	double r843        = vn8vn43 / vn8vn415;
	double r863        = vn8vn63 / vn8vn615;

	double rrms3_15       = rmsvn3_15 / rmsvn15;
        double rmean3_15      = meanvn3_15 / meanvn15;
        double rstdev3_15     = stdevvn3_15 / stdevvn15;
	double rrelfluct3_15  = relfluctvn3_15 / relfluctvn15;
	double r23_15         = vn23_15 / vn215;
	double r43_15         = vn43_15 / vn415;
	double r63_15         = vn63_15 / vn615;
	double r83_15         = vn83_15 / vn815;
	double rgamma1exp3_15 = gamma1exp3_15 / gamma1exp15;
	double r643_15        = vn6vn43_15 / vn6vn415;
	double r843_15        = vn8vn43_15 / vn8vn415;
	double r863_15        = vn8vn63_15 / vn8vn615;

	rmsVtx3_RatiotoVtx15[icent][iEP][iqn]       = rrms3;
	meanVtx3_RatiotoVtx15[icent][iEP][iqn]      = rmean3;
	stDevVtx3_RatiotoVtx15[icent][iEP][iqn]     = rstdev3;
	relFluctVtx3_RatiotoVtx15[icent][iEP][iqn]  = rrelfluct3;
	vn2Vtx3_RatiotoVtx15[icent][iEP][iqn]       = r23;
	vn4Vtx3_RatiotoVtx15[icent][iEP][iqn]       = r43;
	vn6Vtx3_RatiotoVtx15[icent][iEP][iqn]       = r63;
	vn8Vtx3_RatiotoVtx15[icent][iEP][iqn]       = r83;
	gamma1expVtx3_RatiotoVtx15[icent][iEP][iqn] = rgamma1exp3;
	vn6vn4Vtx3_RatiotoVtx15[icent][iEP][iqn]    = r643;
	vn8vn4Vtx3_RatiotoVtx15[icent][iEP][iqn]    = r843;
	vn8vn6Vtx3_RatiotoVtx15[icent][iEP][iqn]    = r863;

	rmsVtx3_15_RatiotoVtx15[icent][iEP][iqn]       = rrms3_15;
        meanVtx3_15_RatiotoVtx15[icent][iEP][iqn]      = rmean3_15;
        stDevVtx3_15_RatiotoVtx15[icent][iEP][iqn]     = rstdev3_15;
	relFluctVtx3_15_RatiotoVtx15[icent][iEP][iqn]  = rrelfluct3_15;
	vn2Vtx3_15_RatiotoVtx15[icent][iEP][iqn]       = r23_15;
	vn4Vtx3_15_RatiotoVtx15[icent][iEP][iqn]       = r43_15;
	vn6Vtx3_15_RatiotoVtx15[icent][iEP][iqn]       = r63_15;
	vn8Vtx3_15_RatiotoVtx15[icent][iEP][iqn]       = r83_15;
	gamma1expVtx3_15_RatiotoVtx15[icent][iEP][iqn] = rgamma1exp3_15;
	vn6vn4Vtx3_15_RatiotoVtx15[icent][iEP][iqn]    = r643_15;
	vn8vn4Vtx3_15_RatiotoVtx15[icent][iEP][iqn]    = r843_15;
	vn8vn6Vtx3_15_RatiotoVtx15[icent][iEP][iqn]    = r863_15;

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	double pctRMSVn_3      = fabs(rmsVnVtx3[icent][iEP][iqn] - rmsVnVtx15[icent][iEP][iqn]) / fabs(rmsVnVtx15[icent][iEP][iqn]);
	double pctMeanVn_3     = fabs(meanVnVtx3[icent][iEP][iqn] - meanVnVtx15[icent][iEP][iqn]) / fabs(meanVnVtx15[icent][iEP][iqn]);
	double pctStDevVn_3    = fabs(stDevVnVtx3[icent][iEP][iqn] - stDevVnVtx15[icent][iEP][iqn]) / fabs(stDevVnVtx15[icent][iEP][iqn]);
	double pctRelFluctVn_3 = fabs(relFluctVnVtx3[icent][iEP][iqn] - relFluctVnVtx15[icent][iEP][iqn]) / fabs(relFluctVnVtx15[icent][iEP][iqn]);
	double pctVn2_3        = fabs(vn2Vtx3[icent][iEP][iqn] - vn2Vtx15[icent][iEP][iqn]) / fabs(vn2Vtx15[icent][iEP][iqn]);
	double pctVn4_3        = fabs(vn4Vtx3[icent][iEP][iqn] - vn4Vtx15[icent][iEP][iqn]) / fabs(vn4Vtx15[icent][iEP][iqn]);
	double pctVn6_3        = fabs(vn6Vtx3[icent][iEP][iqn] - vn6Vtx15[icent][iEP][iqn]) / fabs(vn6Vtx15[icent][iEP][iqn]);
	double pctVn8_3        = fabs(vn8Vtx3[icent][iEP][iqn] - vn8Vtx15[icent][iEP][iqn]) / fabs(vn8Vtx15[icent][iEP][iqn]);
	double pctGamma1Exp_3  = fabs(gamma1expVtx3[icent][iEP][iqn] - gamma1expVtx15[icent][iEP][iqn]) / fabs(gamma1expVtx15[icent][iEP][iqn]);
	double pctVn6Vn4_3     = fabs(vn6vn4Vtx3[icent][iEP][iqn] - vn6vn4Vtx15[icent][iEP][iqn]) / fabs(vn6vn4Vtx15[icent][iEP][iqn]);
	double pctVn8Vn4_3     = fabs(vn8vn4Vtx3[icent][iEP][iqn] - vn8vn4Vtx15[icent][iEP][iqn]) / fabs(vn8vn4Vtx15[icent][iEP][iqn]);
	double pctVn8Vn6_3     = fabs(vn8vn6Vtx3[icent][iEP][iqn] - vn8vn6Vtx15[icent][iEP][iqn]) / fabs(vn8vn6Vtx15[icent][iEP][iqn]);

	double pctRMSVn_3_15      = fabs(rmsVnVtx3_15[icent][iEP][iqn] - rmsVnVtx15[icent][iEP][iqn]) / fabs(rmsVnVtx15[icent][iEP][iqn]);
	double pctMeanVn_3_15     = fabs(meanVnVtx3_15[icent][iEP][iqn] - meanVnVtx15[icent][iEP][iqn]) / fabs(meanVnVtx15[icent][iEP][iqn]);
	double pctStDevVn_3_15    = fabs(stDevVnVtx3_15[icent][iEP][iqn] - stDevVnVtx15[icent][iEP][iqn]) / fabs(stDevVnVtx15[icent][iEP][iqn]);
	double pctRelFluctVn_3_15 = fabs(relFluctVnVtx3_15[icent][iEP][iqn] - relFluctVnVtx15[icent][iEP][iqn]) / fabs(relFluctVnVtx15[icent][iEP][iqn]);
	double pctVn2_3_15        = fabs(vn2Vtx3_15[icent][iEP][iqn] - vn2Vtx15[icent][iEP][iqn]) / fabs(vn2Vtx15[icent][iEP][iqn]);
	double pctVn4_3_15        = fabs(vn4Vtx3_15[icent][iEP][iqn] - vn4Vtx15[icent][iEP][iqn]) / fabs(vn4Vtx15[icent][iEP][iqn]);
	double pctVn6_3_15        = fabs(vn6Vtx3_15[icent][iEP][iqn] - vn6Vtx15[icent][iEP][iqn]) / fabs(vn6Vtx15[icent][iEP][iqn]);
	double pctVn8_3_15        = fabs(vn8Vtx3_15[icent][iEP][iqn] - vn8Vtx15[icent][iEP][iqn]) / fabs(vn8Vtx15[icent][iEP][iqn]);
	double pctGamma1Exp_3_15  = fabs(gamma1expVtx3_15[icent][iEP][iqn] - gamma1expVtx15[icent][iEP][iqn]) / fabs(gamma1expVtx15[icent][iEP][iqn]);
	double pctVn6Vn4_3_15     = fabs(vn6vn4Vtx3_15[icent][iEP][iqn] - vn6vn4Vtx15[icent][iEP][iqn]) / fabs(vn6vn4Vtx15[icent][iEP][iqn]);
	double pctVn8Vn4_3_15     = fabs(vn8vn4Vtx3_15[icent][iEP][iqn] - vn8vn4Vtx15[icent][iEP][iqn]) / fabs(vn8vn4Vtx15[icent][iEP][iqn]);
	double pctVn8Vn6_3_15     = fabs(vn8vn6Vtx3_15[icent][iEP][iqn] - vn8vn6Vtx15[icent][iEP][iqn]) / fabs(vn8vn6Vtx15[icent][iEP][iqn]);

	double reRMSVn      = max(pctRMSVn_3, pctRMSVn_3_15);
	double reMeanVn     = max(pctMeanVn_3, pctMeanVn_3_15);
	double reStDevVn    = max(pctStDevVn_3, pctStDevVn_3_15);
	double reRelFluctVn = max(pctRelFluctVn_3, pctRelFluctVn_3_15);
	double reVn2        = max(pctVn2_3, pctVn2_3_15);
	double reVn4        = max(pctVn4_3, pctVn4_3_15);
	double reVn6        = max(pctVn6_3, pctVn6_3_15);
	double reVn8        = max(pctVn8_3, pctVn8_3_15);
	double reGamma1Exp  = max(pctGamma1Exp_3, pctGamma1Exp_3_15);
	double reVn6Vn4     = max(pctVn6Vn4_3, pctVn6Vn4_3_15);
	double reVn8Vn4     = max(pctVn8Vn4_3, pctVn8Vn4_3_15);
	double reVn8Vn6     = max(pctVn8Vn6_3, pctVn8Vn6_3_15);

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

#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TExec.h"
#include "TF1.h"
#include "TPaletteAxis.h"
#include "TDecompSVD.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void sysUnfoldDistns(int n = 2, double e = 1.0){

  int norder_  = n;
  double tkEta = e;

  bool smooth  = 1;
  bool drawObs = 1;

  TFile * fAna;
  TH1D * hObs[NCENT];
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  TFile * fStat;
  TH1D * hVarianceOfMean_Bin[NCENT];

  TFile * fUnfoldCov;
  TH2D * hUnfoldCov[NCENT][NITER];

  double refoldChi2[NCENT][NITER];
  TGraph * grRefoldChi2[NCENT];

  TFile * fUnfoldResp;
  TH1D * hUnfoldResp[NCENT][NITER];

  TFile * fAnaVtx3;
  TH1D * hObsVtx3[NCENT];
  TFile * fUnfoldVtx3;
  TH1D * hUnfoldVtx3[NCENT][NITER];
  TH1D * hRefoldVtx3[NCENT][NITER];

  TFile * fAnaVtx3_15;
  TH1D * hObsVtx3_15[NCENT];
  TFile * fUnfoldVtx3_15;
  TH1D * hUnfoldVtx3_15[NCENT][NITER];
  TH1D * hRefoldVtx3_15[NCENT][NITER];

  TFile * fAnaTkLoose;
  TH1D * hObsTkLoose[NCENT];
  TFile * fUnfoldTkLoose;
  TH1D * hUnfoldTkLoose[NCENT][NITER];
  TH1D * hRefoldTkLoose[NCENT][NITER];

  TFile * fAnaTkTight;
  TH1D * hObsTkTight[NCENT];
  TFile * fUnfoldTkTight;
  TH1D * hUnfoldTkTight[NCENT][NITER];
  TH1D * hRefoldTkTight[NCENT][NITER];

  TFile * fAnaNewCC;
  TH1D * hObsNewCC[NCENT];
  TFile * fUnfoldNewCC;
  TH1D * hUnfoldNewCC[NCENT][NITER];
  TH1D * hRefoldNewCC[NCENT][NITER];

  TFile * fUnfoldGaussResp;
  TH1D * hUnfoldGaussResp[NCENT][NITER];
  TH1D * hRefoldGaussResp[NCENT][NITER];

  TH1D * hFinalUnf[NCENT];
  TH1D * hFinalUnfSysReg[NCENT];
  TH1D * hFinalUnfSysResp[NCENT];
  TH1D * hFinalUnfSysVtx3[NCENT];
  TH1D * hFinalUnfSysVtx3_15[NCENT];
  TH1D * hFinalUnfSysTkLoose[NCENT];
  TH1D * hFinalUnfSysTkTight[NCENT];
  TH1D * hFinalUnfSysNewCC[NCENT];
  TH1D * hFinalUnfSysGaussResp[NCENT];

  TH1D * hSysReg_RatioToNominal[NCENT];
  TH1D * hSysVtx3_RatioToNominal[NCENT];
  TH1D * hSysVtx3_15_RatioToNominal[NCENT];
  TH1D * hSysTkLoose_RatioToNominal[NCENT];
  TH1D * hSysTkTight_RatioToNominal[NCENT];
  TH1D * hSysNewCC_RatioToNominal[NCENT];
  TH1D * hSysGaussResp_RatioToNominal[NCENT];

  TH1D * hSysErrPointToPoint_Reg[NCENT];
  TH1D * hSysErrPointToPoint_Vtx[NCENT];
  TH1D * hSysErrPointToPoint_TkQ[NCENT];
  TH1D * hSysErrPointToPoint_NCC[NCENT];
  TH1D * hSysErrPointToPoint_GRP[NCENT];
  TH1D * hSysErrPointToPoint_Tot[NCENT];
  TF1 * fSmooth[NCENT];
  double VNMAX[NCENT];

  TFile * fOut;
  TH1D * hFinalUnfoldSys[NCENT];
  TH1D * hFinalUnfoldStatAndSys[NCENT];
  TH1D * hNormFactor;

  TCanvas * cSysReg;
  TCanvas * cSysVtx;
  TCanvas * cSysTkQ;
  TCanvas * cSysNCC;
  TCanvas * cSysGRP;
  TCanvas * cSysResp;
  TCanvas * cFinalUnfold;
  TCanvas * cSysPointToPoint;

  TLine * upperCut_5Sigma[NCENT];
  TLine * upperCut_4Sigma[NCENT];
  TLine * upperCut_3Sigma[NCENT];

  //-- Save Response and Covariance Matrices:
  TH2D * hResponse[NCENT];
  TH2D * hCovMatrix[NCENT];
  TH2D * hCorrMatrix[NCENT];
  TCanvas * cCovRespUnf[NCENT];

  TLatex latex;
  TLatex latex2;
  TLatex latex3;
  TLatex latex4;

  //
  // MAIN
  // 
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  latex2.SetNDC();

  latex3.SetNDC();
  latex3.SetTextFont(43);
  latex3.SetTextSize(26);

  latex4.SetNDC();
  latex4.SetTextFont(43);
  latex4.SetTextSize(33);

  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");

  //-- Get the output files for the reported distributions
  fAna    = new TFile( "../AnalyzerResults/CastleEbyE.root" );
  fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the stat errors for the reported distributions
  fStat = new TFile( Form("../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertainties_v%i.root ", norder_, tkEta, norder_) );
  for(int icent = 0; icent < NCENT; icent++) hVarianceOfMean_Bin[icent] = (TH1D*) fStat->Get( Form("hVarianceOfMean_Bin_c%i", icent) );

  //-- Get the output files for SysResp
  fUnfoldResp = new TFile( Form("../UnfoldResults/dataResp/data%i_dosys.root", norder_) );

  //-- Get the output files for SysVtx3
  fAnaVtx3    = new TFile( "vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE.root" );
  fUnfoldVtx3 = new TFile( Form("vtxCut/vtx_leq_3/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the output files for SysVtx3_15
  fAnaVtx3_15    = new TFile( "vtxCut/vtx3_15/AnalyzerResults/CastleEbyE.root" );
  fUnfoldVtx3_15 = new TFile( Form("vtxCut/vtx3_15/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the output files for SysTkLoose
  fAnaTkLoose    = new TFile( "tkQuality/loose/AnalyzerResults/CastleEbyE.root" );
  fUnfoldTkLoose = new TFile( Form("tkQuality/loose/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the output files for SysTkTight
  fAnaTkTight    = new TFile( "tkQuality/tight/AnalyzerResults/CastleEbyE.root" );
  fUnfoldTkTight = new TFile( Form("tkQuality/tight/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the output files for sysNewCC
  fAnaNewCC    = new TFile( "clusCompatTune/newCCTune2pct/AnalyzerResults/CastleEbyE.root" );
  fUnfoldNewCC = new TFile( Form("clusCompatTune/newCCTune2pct/UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get the output files for sysGaussResp
  fUnfoldGaussResp = new TFile( Form("../UnfoldResults/dataResp/data%iGauss.root", norder_) );

  if(fAna->IsZombie() || fUnfold->IsZombie() ){
    std::cout << "WARNING! Unfolding procedure not run." << std::endl;
    std::cout << "Please run the unfolding first and then run this macro..." << std::endl;
    std::cout << "Exiting now..." << std::endl;
    exit(0);
  }
  if( fStat->IsZombie() ){
    std::cout << "WARNING! Statistic resampling not run." << std::endl;
    std::cout << "Please run the resampling first and then run this macro..." << std::endl;
    std::cout << "Exiting now..." << std::endl;
    exit(0);
  }
  if( fUnfoldResp->IsZombie()    ||
      fAnaVtx3->IsZombie()       ||
      fUnfoldVtx3->IsZombie()    ||
      fAnaVtx3_15->IsZombie()    || 
      fUnfoldVtx3_15->IsZombie() || 
      fAnaTkLoose->IsZombie()    || 
      fUnfoldTkLoose->IsZombie() || 
      fAnaTkTight->IsZombie()    || 
      fUnfoldTkTight->IsZombie() || 
      fAnaNewCC->IsZombie()      || 
      fUnfoldNewCC->IsZombie()   || 
      fUnfoldGaussResp->IsZombie() ){
    std::cout << "WARNING! Systematic studies not run." << std::endl;
    std::cout << "Please run the systematics first and then run this macro..." << std::endl;
    std::cout << "Exiting now..." << std::endl;
    exit(0);
  }

  //-- Set up the output file
  fOut = new TFile(Form("SysUnfoldDistns_v%i.root", norder_), "recreate");

  cSysReg = new TCanvas("cSysReg", "cSysReg", 2000, 1500);
  cSysReg->Divide(4,3);

  cSysVtx = new TCanvas("cSysVtx", "cSysVtx", 2000, 1500);
  cSysVtx->Divide(4,3);

  cSysTkQ = new TCanvas("cSysTkQ", "cSysTkQ", 2000, 1500);
  cSysTkQ->Divide(4,3);

  cSysNCC = new TCanvas("cSysNCC", "cSysNCC", 2000, 1500);
  cSysNCC->Divide(4,3);

  cSysGRP = new TCanvas("cSysGRP", "cSysGRP", 2000, 1500);
  cSysGRP->Divide(4,3);

  cFinalUnfold = new TCanvas("cFinalUnfold", "cFinalUnfold", 2000, 1500);
  cFinalUnfold->Divide(4,3);

  hNormFactor = new TH1D("hNormFactor","hNormFactor", NCENT, centbinsDefault);

  //-- Start fetching histograms....
  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the Response matrix
    hResponse[icent] = (TH2D*) fUnfold->Get( Form("hresp_c%i", icent) );

    //-- -------------- Real Dists --------------
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMarkerStyle(24);

    bool iterCut          = 0;
    bool iterCutReg       = 0;
    bool iterCutGaussResp = 0;
    for(int i = 0; i < NITER; i++){
      hUnfold[icent][i]          = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold[icent][i]->SetLineColor( col[i] );
      hUnfold[icent][i]->SetMarkerColor( col[i] );
      hUnfoldResp[icent][i]      = (TH1D*) fUnfoldResp->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldGaussResp[icent][i] = (TH1D*) fUnfoldGaussResp->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i]          = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold[icent][i]->SetLineColor(col[i]);
      hRefold[icent][i]->SetLineWidth(2);
      hRefoldGaussResp[icent][i] = (TH1D*) fUnfoldGaussResp->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2    = hRefold[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF");
      double chi2GRP = hRefoldGaussResp[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF");

      refoldChi2[icent][i] = chi2;

      if( chi2 < 1.2 && !iterCut ){
	hFinalUnf[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnf_c%i", icent) );
	hFinalUnfSysResp[icent] = (TH1D*) hUnfoldResp[icent][i]->Clone( Form("hFinalUnfSysResp_c%i", icent) );
	iterCut = 1;
	//-- Grab the cov matrix
	hCovMatrix[icent] = (TH2D*) fUnfold->Get( Form("hCovMat%i_c%i", iter[i], icent) );
      }
      if( i == NITER - 1 && !iterCut){
	hFinalUnf[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnf_c%i", icent) );
	hFinalUnfSysResp[icent] = (TH1D*) hUnfoldResp[icent][i]->Clone( Form("hFinalUnfSysResp_c%i", icent) );
	//-- Grab the cov matrix
	hCovMatrix[icent] = (TH2D*) fUnfold->Get( Form("hCovMat%i_c%i", iter[i], icent) );
      }

      //-- -------------- Reg --------------
      if( chi2 < 2. && !iterCutReg ){
        hFinalUnfSysReg[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfSysReg_c%i", icent) );
        iterCutReg = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysReg[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfSysReg_c%i", icent) );
      }

      //-- -------------- GaussResp --------------
      if( chi2GRP < 1.2 && !iterCutGaussResp ){
        hFinalUnfSysGaussResp[icent] = (TH1D*) hUnfoldGaussResp[icent][i]->Clone( Form("hFinalUnfSysGRP_c%i", icent) );
        iterCutGaussResp = 1;
      }
      if( i == NITER - 1 && !iterCutGaussResp){
	hFinalUnfSysGaussResp[icent] = (TH1D*) hUnfoldGaussResp[icent][i]->Clone( Form("hFinalUnfSysGRP_c%i", icent) );
      }

    }

    //-- -------------- Vtx3 --------------
    hObsVtx3[icent] = (TH1D*) fAnaVtx3->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldVtx3[icent][i] = (TH1D*) fUnfoldVtx3->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldVtx3[icent][i] = (TH1D*) fUnfoldVtx3->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldVtx3[icent][i]->Chi2Test(hObsVtx3[icent], "UWCHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysVtx3[icent] = (TH1D*) hUnfoldVtx3[icent][i]->Clone( Form("hFinalUnfSysVtx3_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysVtx3[icent] = (TH1D*) hUnfoldVtx3[icent][i]->Clone( Form("hFinalUnfSysVtx3_c%i", icent) );
      }
    }


    //-- -------------- Vtx3_15 --------------
    hObsVtx3_15[icent]= (TH1D*) fAnaVtx3_15->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldVtx3_15[icent][i] = (TH1D*) fUnfoldVtx3_15->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldVtx3_15[icent][i] = (TH1D*) fUnfoldVtx3_15->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldVtx3_15[icent][i]->Chi2Test(hObsVtx3_15[icent], "UWCHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysVtx3_15[icent] = (TH1D*) hUnfoldVtx3_15[icent][i]->Clone( Form("hFinalUnfSysVtx3_15_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysVtx3_15[icent] = (TH1D*) hUnfoldVtx3_15[icent][i]->Clone( Form("hFinalUnfSysVtx3_15_c%i", icent) );
      }
    }


    //-- -------------- TkLoose --------------
    hObsTkLoose[icent]= (TH1D*) fAnaTkLoose->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldTkLoose[icent][i] = (TH1D*) fUnfoldTkLoose->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldTkLoose[icent][i] = (TH1D*) fUnfoldTkLoose->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldTkLoose[icent][i]->Chi2Test(hObsTkLoose[icent], "UWCHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysTkLoose[icent] = (TH1D*) hUnfoldTkLoose[icent][i]->Clone( Form("hFinalUnfSysTkLoose_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysTkLoose[icent] = (TH1D*) hUnfoldTkLoose[icent][i]->Clone( Form("hFinalUnfSysTkLoose_c%i", icent) );
      }
    }

    //-- -------------- TkTight --------------
    hObsTkTight[icent]= (TH1D*) fAnaTkTight->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldTkTight[icent][i] = (TH1D*) fUnfoldTkTight->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldTkTight[icent][i] = (TH1D*) fUnfoldTkTight->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldTkTight[icent][i]->Chi2Test(hObsTkTight[icent], "UWCHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysTkTight[icent] = (TH1D*) hUnfoldTkTight[icent][i]->Clone( Form("hFinalUnfSysTkTight_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysTkTight[icent] = (TH1D*) hUnfoldTkTight[icent][i]->Clone( Form("hFinalUnfSysTkTight_c%i", icent) );
      }
    }

    //-- -------------- NewCC --------------
    hObsNewCC[icent]= (TH1D*) fAnaNewCC->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldNewCC[icent][i] = (TH1D*) fUnfoldNewCC->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldNewCC[icent][i] = (TH1D*) fUnfoldNewCC->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldNewCC[icent][i]->Chi2Test(hObsNewCC[icent], "UWCHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysNewCC[icent] = (TH1D*) hUnfoldNewCC[icent][i]->Clone( Form("hFinalUnfSysNewCC_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysNewCC[icent] = (TH1D*) hUnfoldNewCC[icent][i]->Clone( Form("hFinalUnfSysNewCC_c%i", icent) );
      }
    }

  } //-- End Cent loop

  //-- Done fetching files and determining the final unfolded distns, now let's calculate some ratios!
  for(int icent = 0; icent < NCENT; icent++){

    //-- Chi2 vs iter
    grRefoldChi2[icent] = new TGraph(NITER, diter, refoldChi2[icent]);
    formatGraph(grRefoldChi2[icent], "Iteration", 0.001, 100000., "#chi^{2}/NDF", 1, 20, Form("grRefoldChi2_c%i", icent));

    //-- Fix Unfolds
    FixUnfold( hFinalUnf[icent] );
    FixUnfold( hFinalUnfSysResp[icent] );
    FixUnfold( hFinalUnfSysReg[icent] );
    FixUnfold( hFinalUnfSysVtx3[icent] );
    FixUnfold( hFinalUnfSysVtx3_15[icent] );
    FixUnfold( hFinalUnfSysTkLoose[icent] );
    FixUnfold( hFinalUnfSysTkTight[icent] );
    FixUnfold( hFinalUnfSysNewCC[icent] );
    FixUnfold( hFinalUnfSysGaussResp[icent] );

    //-- Update Stat Uncertainty for hFinalUnf
    UpdateUnfoldBinErrors(hFinalUnf[icent], hVarianceOfMean_Bin[icent]);

    //-- Normalize
    hNormFactor->SetBinContent( icent+1, 1./hFinalUnf[icent]->Integral() );
    hFinalUnf[icent]->Scale( 1. / hFinalUnf[icent]->Integral() );
    hFinalUnfSysResp[icent]->Scale( 1. / hFinalUnfSysResp[icent]->Integral() );
    hFinalUnfSysReg[icent]->Scale( 1. / hFinalUnfSysReg[icent]->Integral() );
    hFinalUnfSysVtx3[icent]->Scale( 1. / hFinalUnfSysVtx3[icent]->Integral() );
    hFinalUnfSysVtx3_15[icent]->Scale( 1. / hFinalUnfSysVtx3_15[icent]->Integral() );
    hFinalUnfSysTkLoose[icent]->Scale( 1. / hFinalUnfSysTkLoose[icent]->Integral() );
    hFinalUnfSysTkTight[icent]->Scale( 1. / hFinalUnfSysTkTight[icent]->Integral() );
    hFinalUnfSysNewCC[icent]->Scale( 1. / hFinalUnfSysNewCC[icent]->Integral() );
    hFinalUnfSysGaussResp[icent]->Scale( 1. / hFinalUnfSysGaussResp[icent]->Integral() );

    //-- Clone originals
    hSysReg_RatioToNominal[icent]       = (TH1D*) hFinalUnfSysReg[icent]->Clone( Form("hSysReg_RatioToNominal_c%i", icent) );
    hSysVtx3_RatioToNominal[icent]      = (TH1D*) hFinalUnfSysVtx3[icent]->Clone( Form("hSysVtx3_RatioToNominal_c%i", icent) );
    hSysVtx3_15_RatioToNominal[icent]   = (TH1D*) hFinalUnfSysVtx3_15[icent]->Clone( Form("hSysVtx3_15_RatioToNominal_c%i", icent) );
    hSysTkLoose_RatioToNominal[icent]   = (TH1D*) hFinalUnfSysTkLoose[icent]->Clone( Form("hSysTkLoose_RatioToNominal_c%i", icent) );
    hSysTkTight_RatioToNominal[icent]   = (TH1D*) hFinalUnfSysTkTight[icent]->Clone( Form("hSysTkTight_RatioToNominal_c%i", icent) );
    hSysNewCC_RatioToNominal[icent]     = (TH1D*) hFinalUnfSysNewCC[icent]->Clone( Form("hSysNewCC_RatioToNominal_c%i", icent) );
    hSysGaussResp_RatioToNominal[icent] = (TH1D*) hFinalUnfSysGaussResp[icent]->Clone( Form("hSysGaussResp_RatioToNominal_c%i", icent) );


    //-- Point to point uncert for smoothing
    hSysErrPointToPoint_Reg[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Reg_c%i", icent) );
    hSysErrPointToPoint_Vtx[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Vtx_c%i", icent) );
    hSysErrPointToPoint_TkQ[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_TkQ_c%i", icent) );
    hSysErrPointToPoint_NCC[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_NCC_c%i", icent) );
    hSysErrPointToPoint_GRP[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_GRP_c%i", icent) );
    hSysErrPointToPoint_Tot[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Tot_c%i", icent) );

    hSysReg_RatioToNominal[icent]       -> Reset();
    hSysVtx3_RatioToNominal[icent]      -> Reset();
    hSysVtx3_15_RatioToNominal[icent]   -> Reset();
    hSysTkLoose_RatioToNominal[icent]   -> Reset();
    hSysTkTight_RatioToNominal[icent]   -> Reset();
    hSysNewCC_RatioToNominal[icent]     -> Reset();
    hSysGaussResp_RatioToNominal[icent] -> Reset();

    hSysErrPointToPoint_Reg[icent] -> Reset();
    hSysErrPointToPoint_Vtx[icent] -> Reset();
    hSysErrPointToPoint_TkQ[icent] -> Reset();
    hSysErrPointToPoint_NCC[icent] -> Reset();
    hSysErrPointToPoint_GRP[icent] -> Reset();
    hSysErrPointToPoint_Tot[icent] -> Reset();

    //-- Rescale to the same mean as nominal
    double meanvn_Nominal   = hFinalUnf[icent]->GetMean();
    double meanvn_Reg       = hFinalUnfSysReg[icent]->GetMean();
    double meanvn_Vtx3      = hFinalUnfSysVtx3[icent]->GetMean();
    double meanvn_Vtx3_15   = hFinalUnfSysVtx3_15[icent]->GetMean();
    double meanvn_TkLoose   = hFinalUnfSysTkLoose[icent]->GetMean();
    double meanvn_TkTight   = hFinalUnfSysTkTight[icent]->GetMean();
    double meanvn_NewCC     = hFinalUnfSysNewCC[icent]->GetMean();
    double meanvn_GaussResp = hFinalUnfSysGaussResp[icent]->GetMean();

    for(int i = 1; i <= NBins; i++){

      double vn = hFinalUnf[icent]->GetBinCenter(i);
      double pvn_Nominal   = hFinalUnf[icent]->GetBinContent(i);
      double pvn_Reg       = hFinalUnfSysReg[icent]->GetBinContent(i);
      double pvn_Vtx3      = hFinalUnfSysVtx3[icent]->GetBinContent(i);
      double pvn_Vtx3_15   = hFinalUnfSysVtx3_15[icent]->GetBinContent(i);
      double pvn_TkLoose   = hFinalUnfSysTkLoose[icent]->GetBinContent(i);
      double pvn_TkTight   = hFinalUnfSysTkTight[icent]->GetBinContent(i);
      double pvn_NewCC     = hFinalUnfSysNewCC[icent]->GetBinContent(i);
      double pvn_GaussResp = hFinalUnfSysGaussResp[icent]->GetBinContent(i);

      double pvne_Nominal   = hFinalUnf[icent]->GetBinError(i);
      double pvne_Reg       = hFinalUnfSysReg[icent]->GetBinError(i);
      double pvne_Vtx3      = hFinalUnfSysVtx3[icent]->GetBinError(i);
      double pvne_Vtx3_15   = hFinalUnfSysVtx3_15[icent]->GetBinError(i);
      double pvne_TkLoose   = hFinalUnfSysTkLoose[icent]->GetBinError(i);
      double pvne_TkTight   = hFinalUnfSysTkTight[icent]->GetBinError(i);
      double pvne_NewCC     = hFinalUnfSysNewCC[icent]->GetBinError(i);
      double pvne_GaussResp = hFinalUnfSysGaussResp[icent]->GetBinError(i);

      double shiftedvn_Reg       = vn + ( meanvn_Nominal - meanvn_Reg );
      double shiftedvn_Vtx3      = vn + ( meanvn_Nominal - meanvn_Vtx3 );
      double shiftedvn_Vtx3_15   = vn + ( meanvn_Nominal - meanvn_Vtx3_15 );
      double shiftedvn_TkLoose   = vn + ( meanvn_Nominal - meanvn_TkLoose );
      double shiftedvn_TkTight   = vn + ( meanvn_Nominal - meanvn_TkTight );
      double shiftedvn_NewCC     = vn + ( meanvn_Nominal - meanvn_NewCC );
      double shiftedvn_GaussResp = vn + ( meanvn_Nominal - meanvn_GaussResp );

      int shiftedvnBIN_Reg       = hSysReg_RatioToNominal[icent]->FindBin( shiftedvn_Reg );
      int shiftedvnBIN_Vtx3      = hSysVtx3_RatioToNominal[icent]->FindBin( shiftedvn_Vtx3 );
      int shiftedvnBIN_Vtx3_15   = hSysVtx3_15_RatioToNominal[icent]->FindBin( shiftedvn_Vtx3_15 );
      int shiftedvnBIN_TkLoose   = hSysTkLoose_RatioToNominal[icent]->FindBin( shiftedvn_TkLoose );
      int shiftedvnBIN_TkTight   = hSysTkTight_RatioToNominal[icent]->FindBin( shiftedvn_TkTight );
      int shiftedvnBIN_NewCC     = hSysNewCC_RatioToNominal[icent]->FindBin( shiftedvn_NewCC );
      int shiftedvnBIN_GaussResp = hSysGaussResp_RatioToNominal[icent]->FindBin( shiftedvn_GaussResp );

      hSysReg_RatioToNominal[icent]       -> SetBinContent(shiftedvnBIN_Reg, pvn_Reg);
      hSysVtx3_RatioToNominal[icent]      -> SetBinContent(shiftedvnBIN_Vtx3, pvn_Vtx3);
      hSysVtx3_15_RatioToNominal[icent]   -> SetBinContent(shiftedvnBIN_Vtx3_15, pvn_Vtx3_15);
      hSysTkLoose_RatioToNominal[icent]   -> SetBinContent(shiftedvnBIN_TkLoose, pvn_TkLoose);
      hSysTkTight_RatioToNominal[icent]   -> SetBinContent(shiftedvnBIN_TkTight, pvn_TkTight);
      hSysNewCC_RatioToNominal[icent]     -> SetBinContent(shiftedvnBIN_NewCC, pvn_NewCC);
      hSysGaussResp_RatioToNominal[icent] -> SetBinContent(shiftedvnBIN_GaussResp, pvn_GaussResp);

      hSysReg_RatioToNominal[icent]       -> SetBinError(shiftedvnBIN_Reg, pvne_Reg);
      hSysVtx3_RatioToNominal[icent]      -> SetBinError(shiftedvnBIN_Vtx3, pvne_Vtx3);
      hSysVtx3_15_RatioToNominal[icent]   -> SetBinError(shiftedvnBIN_Vtx3_15, pvne_Vtx3_15);
      hSysTkLoose_RatioToNominal[icent]   -> SetBinError(shiftedvnBIN_TkLoose, pvne_TkLoose);
      hSysTkTight_RatioToNominal[icent]   -> SetBinError(shiftedvnBIN_TkTight, pvne_TkTight);
      hSysNewCC_RatioToNominal[icent]     -> SetBinError(shiftedvnBIN_NewCC, pvne_NewCC);
      hSysGaussResp_RatioToNominal[icent] -> SetBinError(shiftedvnBIN_GaussResp, pvne_GaussResp);

    }


    //-- Calculate ratio to nominal
    hSysReg_RatioToNominal[icent]       -> Divide( hFinalUnf[icent] );
    hSysVtx3_RatioToNominal[icent]      -> Divide( hFinalUnf[icent] );
    hSysVtx3_15_RatioToNominal[icent]   -> Divide( hFinalUnf[icent] );
    hSysTkLoose_RatioToNominal[icent]   -> Divide( hFinalUnf[icent] );
    hSysTkTight_RatioToNominal[icent]   -> Divide( hFinalUnf[icent] );
    hSysNewCC_RatioToNominal[icent]     -> Divide( hFinalUnf[icent] );
    hSysGaussResp_RatioToNominal[icent] -> Divide( hFinalUnf[icent] );

    hSysReg_RatioToNominal[icent]       -> SetMinimum( 0.8 );
    hSysVtx3_RatioToNominal[icent]      -> SetMinimum( 0.8 );
    hSysVtx3_15_RatioToNominal[icent]   -> SetMinimum( 0.8 );
    hSysTkLoose_RatioToNominal[icent]   -> SetMinimum( 0.8 );
    hSysTkTight_RatioToNominal[icent]   -> SetMinimum( 0.8 );
    hSysNewCC_RatioToNominal[icent]     -> SetMinimum( 0.8 );
    hSysGaussResp_RatioToNominal[icent] -> SetMinimum( 0.8 );

    hSysReg_RatioToNominal[icent]       -> SetMaximum( 1.2 );
    hSysVtx3_RatioToNominal[icent]      -> SetMaximum( 1.2 );
    hSysVtx3_15_RatioToNominal[icent]   -> SetMaximum( 1.2 );
    hSysTkLoose_RatioToNominal[icent]   -> SetMaximum( 1.2 );
    hSysTkTight_RatioToNominal[icent]   -> SetMaximum( 1.2 );
    hSysNewCC_RatioToNominal[icent]     -> SetMaximum( 1.2 );
    hSysGaussResp_RatioToNominal[icent] -> SetMaximum( 1.2 );

    hSysErrPointToPoint_Reg[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_Vtx[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_TkQ[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_NCC[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_GRP[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_Tot[icent]->SetMaximum(0.01);

    hSysReg_RatioToNominal[icent]       -> SetLineColor(kOrange+2);
    hSysVtx3_RatioToNominal[icent]      -> SetLineColor(2);
    hSysVtx3_15_RatioToNominal[icent]   -> SetLineColor(4);
    hSysTkLoose_RatioToNominal[icent]   -> SetLineColor(8);
    hSysTkTight_RatioToNominal[icent]   -> SetLineColor(9);
    hSysNewCC_RatioToNominal[icent]     -> SetLineColor(kViolet-1);
    hSysGaussResp_RatioToNominal[icent] -> SetLineColor(kSpring+4);

    hSysReg_RatioToNominal[icent]       -> SetMarkerColor(kOrange+2);
    hSysVtx3_RatioToNominal[icent]      -> SetMarkerColor(2);
    hSysVtx3_15_RatioToNominal[icent]   -> SetMarkerColor(4);
    hSysTkLoose_RatioToNominal[icent]   -> SetMarkerColor(8);
    hSysTkTight_RatioToNominal[icent]   -> SetMarkerColor(9);
    hSysNewCC_RatioToNominal[icent]     -> SetMarkerColor(kViolet-1);
    hSysGaussResp_RatioToNominal[icent] -> SetMarkerColor(kSpring+4);

    hSysReg_RatioToNominal[icent]       -> SetMarkerStyle(34);
    hSysVtx3_RatioToNominal[icent]      -> SetMarkerStyle(20);
    hSysVtx3_15_RatioToNominal[icent]   -> SetMarkerStyle(24);
    hSysTkLoose_RatioToNominal[icent]   -> SetMarkerStyle(21);
    hSysTkTight_RatioToNominal[icent]   -> SetMarkerStyle(25);
    hSysNewCC_RatioToNominal[icent]     -> SetMarkerStyle(22);
    hSysGaussResp_RatioToNominal[icent] -> SetMarkerStyle(23);

    hSysReg_RatioToNominal[icent]       -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysVtx3_RatioToNominal[icent]      -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysVtx3_15_RatioToNominal[icent]   -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysTkLoose_RatioToNominal[icent]   -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysTkTight_RatioToNominal[icent]   -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysNewCC_RatioToNominal[icent]     -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysGaussResp_RatioToNominal[icent] -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );

    hSysReg_RatioToNominal[icent]       -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysVtx3_RatioToNominal[icent]      -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysVtx3_15_RatioToNominal[icent]   -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysTkLoose_RatioToNominal[icent]   -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysTkTight_RatioToNominal[icent]   -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysNewCC_RatioToNominal[icent]     -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysGaussResp_RatioToNominal[icent] -> GetYaxis() -> SetTitle( "Ratio to Default" );

    //-- Set up the sysErr final unfolded distn
    hFinalUnfoldSys[icent]        = (TH1D*) hFinalUnf[icent]->Clone( Form("hFinalUnfoldSys_c%i", icent) );
    hFinalUnfoldStatAndSys[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hFinalUnfoldStatAndSys_c%i", icent) );

    hFinalUnfoldSys[icent]->GetXaxis()->SetNdivisions(508);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitle( Form("p(v_{%i})", norder_) );

    hFinalUnfoldStatAndSys[icent]->GetXaxis()->SetNdivisions(508);
    hFinalUnfoldStatAndSys[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hFinalUnfoldStatAndSys[icent]->GetYaxis()->SetTitle( Form("p(v_{%i})", norder_) );

    double upC5Sig =  hFinalUnf[icent]->GetMean() + 5.*hFinalUnf[icent]->GetRMS();
    double upC4Sig =  hFinalUnf[icent]->GetMean() + 4.*hFinalUnf[icent]->GetRMS();
    double upC3Sig =  hFinalUnf[icent]->GetMean() + 3.*hFinalUnf[icent]->GetRMS();

    double UL = hFinalUnf[icent]->GetMaximum();

    upperCut_5Sigma[icent] = new TLine(upC5Sig, 0., upC5Sig, UL);
    upperCut_5Sigma[icent]->SetLineColor(2);
    upperCut_5Sigma[icent]->SetLineWidth(2);
    upperCut_4Sigma[icent] = new TLine(upC4Sig, 0., upC4Sig, UL);
    upperCut_4Sigma[icent]->SetLineColor(4);
    upperCut_4Sigma[icent]->SetLineWidth(2);
    upperCut_3Sigma[icent] = new TLine(upC3Sig, 0., upC3Sig, UL);
    upperCut_3Sigma[icent]->SetLineColor(8);
    upperCut_3Sigma[icent]->SetLineWidth(2);

    //-- Loop over bins and calculate total sysErr
    for(int i = 1; i <= NBins; i++){

      hFinalUnfoldSys[icent]->SetBinError(i, 0);
      hFinalUnfoldStatAndSys[icent]->SetBinError(i, 0);

      if( hFinalUnf[icent]->GetBinContent(i) == 0 ) continue;

      VNMAX[icent] = hFinalUnfoldSys[icent]->GetBinCenter(i);

      double relErr_Reg = hSysReg_RatioToNominal[icent]->GetBinContent(i);
      double relErr_Rege = hSysReg_RatioToNominal[icent]->GetBinError(i);
      relErr_Reg = fabs(1. - relErr_Reg);
      if( compatibleWithZero( relErr_Reg, relErr_Rege) ) relErr_Reg = 0.;

      double relErr_Vtx3 = hSysVtx3_RatioToNominal[icent]->GetBinContent(i);
      double relErr_Vtx3e = hSysVtx3_RatioToNominal[icent]->GetBinError(i);
      relErr_Vtx3 = fabs(1. - relErr_Vtx3);
      if( compatibleWithZero( relErr_Vtx3, relErr_Vtx3e) ) relErr_Vtx3 = 0.;

      double relErr_Vtx3_15 = hSysVtx3_15_RatioToNominal[icent]->GetBinContent(i);
      double relErr_Vtx3_15e = hSysVtx3_15_RatioToNominal[icent]->GetBinError(i);
      relErr_Vtx3_15 = fabs(1. - relErr_Vtx3_15);
      if( compatibleWithZero( relErr_Vtx3_15, relErr_Vtx3_15e) ) relErr_Vtx3_15 = 0.;

      double relErr_TkLoose = hSysTkLoose_RatioToNominal[icent]->GetBinContent(i);
      double relErr_TkLoosee = hSysTkLoose_RatioToNominal[icent]->GetBinError(i);
      relErr_TkLoose = fabs(1. - relErr_TkLoose);
      if( compatibleWithZero( relErr_TkLoose, relErr_TkLoosee) ) relErr_TkLoose = 0.;

      double relErr_TkTight = hSysTkTight_RatioToNominal[icent]->GetBinContent(i);
      double relErr_TkTighte = hSysTkTight_RatioToNominal[icent]->GetBinError(i);
      relErr_TkTight = fabs(1. - relErr_TkTight);
      if( compatibleWithZero( relErr_TkTight, relErr_TkTighte) ) relErr_TkTight = 0.;

      double relErr_NewCC = hSysNewCC_RatioToNominal[icent]->GetBinContent(i);
      double relErr_NewCCe = hSysNewCC_RatioToNominal[icent]->GetBinError(i);
      relErr_NewCC = fabs(1. - relErr_NewCC);
      if( compatibleWithZero( relErr_NewCC, relErr_NewCCe) ) relErr_NewCC = 0.;

      double relErr_GaussResp = hSysGaussResp_RatioToNominal[icent]->GetBinContent(i);
      double relErr_GaussRespe = hSysGaussResp_RatioToNominal[icent]->GetBinError(i);
      relErr_GaussResp = fabs(1. - relErr_GaussResp);
      if( compatibleWithZero( relErr_GaussResp, relErr_GaussRespe) ) relErr_GaussResp = 0.;

      double relErr_Vtx       = std::max(relErr_Vtx3, relErr_Vtx3_15);
      double relErr_TkQuality = std::max(relErr_TkLoose, relErr_TkLoose);

      double sysErrReg = relErr_Reg * hFinalUnf[icent]->GetBinContent(i);
      double sysErrVtx = relErr_Vtx * hFinalUnf[icent]->GetBinContent(i);
      double sysErrTkQ = relErr_TkQuality * hFinalUnf[icent]->GetBinContent(i);
      double sysErrNCC = relErr_NewCC * hFinalUnf[icent]->GetBinContent(i);
      double sysErrGRP = relErr_GaussResp * hFinalUnf[icent]->GetBinContent(i);

      double relErrTot = sqrt( pow(relErr_Reg, 2) + pow(relErr_Vtx, 2) + pow(relErr_TkQuality, 2) + pow(sysErrNCC, 2) + pow(sysErrGRP, 2) );
      double statErr    = hFinalUnf[icent]->GetBinError(i);
      double respErr    = hFinalUnfSysResp[icent]->GetBinError(i);
      double sysErrRVT  = relErrTot * hFinalUnf[icent]->GetBinContent(i);
      double sysErrResp = sqrt( pow(respErr,2) - pow(statErr,2) );
      double sysErrTot  = sqrt( sysErrRVT*sysErrRVT + sysErrResp*sysErrResp );
      double totErr     = sqrt(pow(sysErrTot, 2) + pow(statErr, 2));

      hSysErrPointToPoint_Reg[icent]->SetBinContent(i, sysErrReg );
      hSysErrPointToPoint_Vtx[icent]->SetBinContent(i, sysErrVtx );
      hSysErrPointToPoint_TkQ[icent]->SetBinContent(i, sysErrTkQ );
      hSysErrPointToPoint_NCC[icent]->SetBinContent(i, sysErrNCC );
      hSysErrPointToPoint_GRP[icent]->SetBinContent(i, sysErrGRP );
      hSysErrPointToPoint_Tot[icent]->SetBinContent(i, sysErrTot );

      hFinalUnfoldSys[icent]->SetBinError(i, sysErrTot);
      hFinalUnfoldStatAndSys[icent]->SetBinError(i, totErr);

    }

    if(smooth){

      hSysErrPointToPoint_Reg[icent]->Smooth();
      hSysErrPointToPoint_Vtx[icent]->Smooth();
      hSysErrPointToPoint_TkQ[icent]->Smooth();
      hSysErrPointToPoint_NCC[icent]->Smooth();
      hSysErrPointToPoint_GRP[icent]->Smooth();
      hSysErrPointToPoint_Tot[icent]->Smooth();

      for(int i = 1; i <= NBins; i++){
	if( hFinalUnf[icent]->GetBinContent(i) == 0 ) continue;
	double sysReg  = hSysErrPointToPoint_Reg[icent]->GetBinContent(i);
	double sysVtx  = hSysErrPointToPoint_Vtx[icent]->GetBinContent(i);
	double sysTkQ  = hSysErrPointToPoint_TkQ[icent]->GetBinContent(i);
	double sysNCC  = hSysErrPointToPoint_NCC[icent]->GetBinContent(i);
	double sysGRP  = hSysErrPointToPoint_GRP[icent]->GetBinContent(i);

	double statErr = hFinalUnf[icent]->GetBinError(i);
	double respErr = hFinalUnfSysResp[icent]->GetBinError(i);
	double sysResp = sqrt( pow(respErr,2) - pow(statErr,2) );

	double sys = sqrt( pow(sysResp,2) + pow(sysTkQ,2) + pow(sysVtx,2) + pow(sysReg,2) + pow(sysNCC,2) + pow(sysGRP,2) );
	double tot = sqrt( sys*sys + statErr*statErr );
	hFinalUnfoldSys[icent]->SetBinError(i, sys);
	hFinalUnfoldStatAndSys[icent]->SetBinError(i, tot);
      }
    }

    //fSmooth[icent] = new TF1(Form("fSmooth_c%i", icent), "pol2", 0.0, VNMAX[icent]);
    //hSysErrPointToPoint[icent]->Fit( Form("fSmooth_c%i", icent), "QL0", "", 0.0, VNMAX[icent]);

    fOut->cd();
    hFinalUnf[icent]->Write( Form( "hFinalUnfoldStat_c%i", icent) );
    hFinalUnfoldSys[icent]->Write();
    hFinalUnfoldStatAndSys[icent]->Write();
    hCovMatrix[icent]->Write( Form("finalCovMat_c%i", icent) );
    hResponse[icent]->Write();
    if(icent == NCENT-1) hNormFactor->Write();

    cSysReg->cd(icent+1);
    hSysReg_RatioToNominal[icent]     -> Draw();
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
    latex.DrawLatex(0.67, 0.8, "SysReg");

    cSysVtx->cd(icent+1);
    hSysVtx3_RatioToNominal[icent]    -> Draw();
    hSysVtx3_15_RatioToNominal[icent] -> Draw("same");
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
    latex.DrawLatex(0.67, 0.8, "SysVtx");

    cSysTkQ->cd(icent+1);
    hSysTkLoose_RatioToNominal[icent] -> Draw();
    hSysTkTight_RatioToNominal[icent] -> Draw("same");
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
    latex.DrawLatex(0.67, 0.8, "SysTkQ");

    cSysNCC->cd(icent+1);
    hSysNewCC_RatioToNominal[icent] -> Draw();
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
    latex.DrawLatex(0.67, 0.8, "SysNCC");

    cSysGRP->cd(icent+1);
    hSysGaussResp_RatioToNominal[icent] -> Draw();
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
    latex.DrawLatex(0.67, 0.8, "SysGRP");

  } // End Cent loop

  TLegend * legVtx = new TLegend(0.67, 0.6, 0.9, 0.79);
  legVtx->SetBorderSize(0);
  legVtx->SetFillStyle(0);
  legVtx->AddEntry(hSysVtx3_RatioToNominal[0], "|v_{z}| < 3 cm", "lp");
  legVtx->AddEntry(hSysVtx3_15_RatioToNominal[0], "3 < |v_{z}| < 15 cm", "lp");

  TLegend * legTkQ = new TLegend(0.67, 0.6, 0.9, 0.79);
  legTkQ->SetBorderSize(0);
  legTkQ->SetFillStyle(0);
  legTkQ->AddEntry(hSysTkLoose_RatioToNominal[0], "Loose Cuts", "lp");
  legTkQ->AddEntry(hSysTkTight_RatioToNominal[0], "Tight Cuts", "lp");

  cSysVtx->cd(1);
  legVtx->Draw("same");

  cSysTkQ->cd(1);
  legTkQ->Draw("same");

  cSysReg->SaveAs( Form("../plots/systematicStudies/sysReg_UnfoldDistn_v%i.pdf", norder_) );
  cSysVtx->SaveAs( Form("../plots/systematicStudies/sysVtx_UnfoldDistn_v%i.pdf", norder_) );
  cSysTkQ->SaveAs( Form("../plots/systematicStudies/sysTkQ_UnfoldDistn_v%i.pdf", norder_) );
  cSysNCC->SaveAs( Form("../plots/systematicStudies/sysNCC_UnfoldDistn_v%i.pdf", norder_) );
  cSysGRP->SaveAs( Form("../plots/systematicStudies/sysGRP_UnfoldDistn_v%i.pdf", norder_) );

  //-- Draw the resp and cov response matrices with the final unfolded distn
  TPaletteAxis * palette[NCENT];
  TPaletteAxis * paletteCov[NCENT];

  TLegend * legUnf[NCENT];

  TLegend * leg = new TLegend(0.1283, 0.2440, 0.8108, 0.8539);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hObs[0], Form("v_{%i}^{obs}", norder_), "lp");
  for(int i = 0; i < 6; i++) leg->AddEntry(hRefold[0][i], Form("D'Agostini iter = %i", iter[i]), "l");

  TLine * l1p2 = new TLine(grRefoldChi2[0]->GetXaxis()->GetXmin(), 1.2, grRefoldChi2[0]->GetXaxis()->GetXmax(), 1.2);
  l1p2->SetLineStyle(2);
  l1p2->SetLineWidth(2);

  for(int icent = 0; icent < NCENT; icent++){
    cCovRespUnf[icent] = new TCanvas(Form("cCovRespUnf_c%i", icent), Form("cCovRespUnf_c%i", icent), 1500, 1000);
    cCovRespUnf[icent]->Divide(3,2);

    cCovRespUnf[icent]->cd(1);
    cCovRespUnf[icent]->cd(1)->SetLogz();
    hResponse[icent]->GetZaxis()->SetNdivisions(508);
    //hResponse[icent]->Scale(1./hResponse[icent]->Integral());
    hResponse[icent]->Draw("colz");
    cCovRespUnf[icent]->Update();
    palette[icent] = (TPaletteAxis*) hResponse[icent]->GetListOfFunctions()->FindObject("palette");
    palette[icent]->SetX1NDC(0.75);
    palette[icent]->SetY1NDC(0.2);
    palette[icent]->SetX2NDC(0.8);
    palette[icent]->SetY2NDC(0.5);
    latex.DrawLatex(0.6, 0.89, "Response Matrix");
    latex.DrawLatex(0.18, 0.89, "A");

    cCovRespUnf[icent]->cd(2);
    cCovRespUnf[icent]->cd(2)->SetLogy();
    hObs[icent]->Scale(1./hObs[icent]->Integral());
    hObs[icent]->GetYaxis()->SetTitle( Form("p(v_{%i})", norder_) );
    hObs[icent]->Draw();
    hFinalUnfoldSys[icent]->Draw("e2same");
    setex2->Draw();
    hFinalUnf[icent]->Draw("same");
    legUnf[icent] = new TLegend(0.6481, 0.7578, 0.9933, 0.8940);
    legUnf[icent]->SetBorderSize(0);
    legUnf[icent]->SetFillStyle(0);
    legUnf[icent]->AddEntry(hObs[icent],       "Observed", "lp");
    legUnf[icent]->AddEntry(hFinalUnf[icent],  "Unfolded", "lp");
    legUnf[icent]->Draw("same");
    latex.DrawLatex(0.56, 0.89, "Unfold Distribution");
    latex.DrawLatex(0.18, 0.89, "B");

    //-- Convert cov matrix to corr matrix
    hCorrMatrix[icent] = new TH2D( Form("hCorrMatrix_c%i", icent), Form("hCorrMatrix_c%i", icent), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_] );
    hCorrMatrix[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hCorrMatrix[icent]->GetYaxis()->SetTitle( Form("v_{%i}", norder_) );
    hCorrMatrix[icent]->SetOption("colz");

    TMatrixD MCov = H2M( hCovMatrix[icent] );
    TMatrixD S(NBins, NBins);
    for(int i = 0; i<NBins; i++){
      S(i,i) = sqrt( MCov(i,i) );
    }
    TDecompSVD svd(S);
    TMatrixD Sinv = svd.Invert();
    TMatrixD MCov2 = MCov;
    MCov2 *= Sinv;
    TMatrixD MCorr = Sinv;
    MCorr *= MCov2;
    M2H(MCorr, hCorrMatrix[icent]);

    /*
    cCovRespUnf[icent]->cd(3);
    hCovMatrix[icent]->GetZaxis()->SetNdivisions(508);
    hCovMatrix[icent]->Draw("colz");
    cCovRespUnf[icent]->Update();
    paletteCov[icent] = (TPaletteAxis*) hCovMatrix[icent]->GetListOfFunctions()->FindObject("palette");
    paletteCov[icent]->SetX1NDC(0.75);
    paletteCov[icent]->SetY1NDC(0.2);
    paletteCov[icent]->SetX2NDC(0.8);
    paletteCov[icent]->SetY2NDC(0.5);
    latex.DrawLatex(0.58, 0.89, "Covariance Matrix");
    */
    cCovRespUnf[icent]->cd(3);
    hCorrMatrix[icent]->GetXaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->GetYaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->GetZaxis()->SetNdivisions(508);
    hCorrMatrix[icent]->Draw("colz");
    cCovRespUnf[icent]->Update();
    paletteCov[icent] = (TPaletteAxis*) hCorrMatrix[icent]->GetListOfFunctions()->FindObject("palette");
    paletteCov[icent]->SetX1NDC(0.75);
    paletteCov[icent]->SetY1NDC(0.2);
    paletteCov[icent]->SetX2NDC(0.8);
    paletteCov[icent]->SetY2NDC(0.5);
    latex.DrawLatex(0.58, 0.89, "Correlation Matrix");
    latex.DrawLatex(0.18, 0.89, "C");

    cCovRespUnf[icent]->cd(4);
    cCovRespUnf[icent]->cd(4)->SetLogx();
    grRefoldChi2[icent]->Draw("alp");
    l1p2->Draw("same");
    grRefoldChi2[icent]->GetYaxis()->SetRangeUser(0.1, 10000.);
    cCovRespUnf[icent]->cd(4)->SetLogy();
    latex.DrawLatex(0.65, 0.89, "Refold #chi^{2}/NDF");
    latex.DrawLatex(0.18, 0.89, "D");

    cCovRespUnf[icent]->cd(5);
    cCovRespUnf[icent]->cd(5)->SetLogy();
    hObs[icent]->Draw();
    for(int i = 0; i<6; i++){
      hRefold[icent][i]->Scale(1./hRefold[icent][i]->Integral());
      hRefold[icent][i]->Draw("same");
    }
    latex.DrawLatex(0.7, 0.89, "Refolding");
    latex.DrawLatex(0.18, 0.89, "E");

    cCovRespUnf[icent]->cd(6);
    leg->Draw();
    latex2.DrawLatex(0.1283, 0.87,  Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));

    cCovRespUnf[icent]->Update();
    cCovRespUnf[icent]->SaveAs( Form("../plots/unfolding/cCovRespUnf_cent%i_%i.pdf", cent_min[icent], cent_max[icent]) );
  }

  //-- Big unfold distns plots:
  for(int icent = 0; icent < NCENT; icent++){
    cFinalUnfold->cd(icent+1);
    cFinalUnfold->cd(icent+1)->SetLogy();
    hFinalUnfoldSys[icent]->SetLineColor(2);
    hFinalUnfoldSys[icent]->SetMarkerColor(2);
    hFinalUnfoldSys[icent]->SetFillColor(15);
    hFinalUnfoldSys[icent]->GetXaxis()->SetRange(1, hFinalUnfoldSys[icent]->FindBin(0.29));
    hFinalUnfoldSys[icent]->Draw("e2");
    setex2->Draw();
    hFinalUnf[icent]->SetLineColor(1);
    hFinalUnf[icent]->SetMarkerColor(1);
    hFinalUnf[icent]->Draw("same");
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
  }
  cFinalUnfold->SaveAs( Form("../plots/systematicStudies/sysFinalUnfoldDistns_v%i.pdf", norder_) );

  //-- 6 x 2
  /*
  //-- Draw Final Unfolded distns on one plot
  TLegend * legUnfObs = new TLegend(0.62, 0.06, 0.82, 0.26);
  legInit(legUnfObs);

  TCanvas * cFinalUnfoldMerged = new TCanvas("cFinalUnfoldMerged", "cFinalUnfoldMerged", 3000, 1000);
  cFinalUnfoldMerged->Divide(6,2,0,0);

  for(int icent = 0; icent < NCENT; icent++){

    cFinalUnfoldMerged->cd(icent+1);
    cFinalUnfoldMerged->cd(icent+1)->SetLogy();
    if(icent < 6) cFinalUnfoldMerged->cd(icent+1)->SetTopMargin(0.07);
    else          cFinalUnfoldMerged->cd(icent+1)->SetBottomMargin(0.22);
    if(icent == 0 || icent == 6)  cFinalUnfoldMerged->cd(icent+1)->SetLeftMargin(0.17);
    if(icent == 5 || icent == 11) cFinalUnfoldMerged->cd(icent+1)->SetRightMargin(0.01);

    //-- X axes 
    double m = 0.27;
    hObs[icent]->GetXaxis()->SetRange(0, hFinalUnfoldSys[icent]->FindBin(m));
    hFinalUnfoldSys[icent]->GetXaxis()->SetRange(0, hFinalUnfoldSys[icent]->FindBin(m));
    hFinalUnf[icent]->GetXaxis()->SetRange(0, hFinalUnfoldSys[icent]->FindBin(m));

    hObs[icent]->GetXaxis()->SetNdivisions(507);
    hFinalUnfoldSys[icent]->GetXaxis()->SetNdivisions(507);
    hFinalUnf[icent]->GetXaxis()->SetNdivisions(507);

    hObs[icent]->GetXaxis()->SetLabelFont(43);
    hFinalUnfoldSys[icent]->GetXaxis()->SetLabelFont(43);
    hFinalUnf[icent]->GetXaxis()->SetLabelFont(43);

    hObs[icent]->GetXaxis()->SetLabelSize(30);
    hFinalUnfoldSys[icent]->GetXaxis()->SetLabelSize(30);
    hFinalUnf[icent]->GetXaxis()->SetLabelSize(30);

    hObs[icent]->GetXaxis()->SetTitleFont(43);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitleFont(43);
    hFinalUnf[icent]->GetXaxis()->SetTitleFont(43);

    hObs[icent]->GetXaxis()->SetTitleSize(45);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitleSize(45);
    hFinalUnf[icent]->GetXaxis()->SetTitleSize(45);

    hObs[icent]->GetXaxis()->SetTitleOffset(1.3);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitleOffset(1.3);
    hFinalUnf[icent]->GetXaxis()->SetTitleOffset(1.3);

    //-- Y axes
    double mm = 0.5;
    hObs[icent]->SetMaximum(mm);
    hFinalUnfoldSys[icent]->SetMaximum(mm);
    hFinalUnf[icent]->SetMaximum(mm);

    double mmm = 2e-5;
    hObs[icent]->SetMinimum(mmm);
    hFinalUnfoldSys[icent]->SetMinimum(mmm);
    hFinalUnf[icent]->SetMinimum(mmm);

    hObs[icent]->GetYaxis()->SetLabelFont(43);
    hFinalUnfoldSys[icent]->GetYaxis()->SetLabelFont(43);
    hFinalUnf[icent]->GetYaxis()->SetLabelFont(43);

    hObs[icent]->GetYaxis()->SetLabelSize(30);
    hFinalUnfoldSys[icent]->GetYaxis()->SetLabelSize(30);
    hFinalUnf[icent]->GetYaxis()->SetLabelSize(30);

    hObs[icent]->GetYaxis()->SetTitleFont(43);
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitleFont(43);
    hFinalUnf[icent]->GetYaxis()->SetTitleFont(43);

    hObs[icent]->GetYaxis()->SetTitleSize(40);
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitleSize(40);
    hFinalUnf[icent]->GetYaxis()->SetTitleSize(40);

    hObs[icent]->GetYaxis()->SetTitleOffset(1.8);
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitleOffset(1.8);
    hFinalUnf[icent]->GetYaxis()->SetTitleOffset(1.8);

    //-- Colors/Cosmetics
    hObs[icent]->Scale(1./hObs[icent]->Integral());
    hObs[icent]->SetLineColor(1);
    hObs[icent]->SetMarkerColor(1);
    hObs[icent]->SetMarkerStyle(24);

    hFinalUnfoldSys[icent]->SetLineColor(4);
    hFinalUnfoldSys[icent]->SetMarkerColor(4);
    hFinalUnfoldSys[icent]->SetFillColorAlpha(kBlue-7, 0.6);

    hFinalUnf[icent]->SetLineColor(4);
    hFinalUnf[icent]->SetMarkerColor(4);

    if(drawObs){
      hFinalUnfoldSys[icent]->Draw("e2");
      hFinalUnf[icent]->Draw("same");
      hObs[icent]->Draw("same");

      //-- Add legend
      if(icent == 0){
	legUnfObs->AddEntry(hObs[icent],      "Observed p(v_{2})", "lp");
	legUnfObs->AddEntry(hFinalUnf[icent], "Unfold p(v_{2})",   "lp");
	legUnfObs->Draw("same");
	legUnfObs->SetTextFont(43);
	legUnfObs->SetTextSize(26);
      }
    }
    else{
      hFinalUnfoldSys[icent]->Draw("e2");
      hFinalUnf[icent]->Draw("same");
    }

    //-- Centrality tags
    if(icent == 0)      latex3.DrawLatex(0.21, 0.06, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent < 6)  latex3.DrawLatex(0.05, 0.06, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent == 6) latex3.DrawLatex(0.20, 0.27, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent < 12) latex3.DrawLatex(0.05, 0.27, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );

    //-- CMS, lumi pt and eta tags
    if(icent == 0){
      latex4.DrawLatex(0.19, 0.95, "#bf{CMS}");
      latex4.DrawLatex(0.65, 0.95, "PbPb 5.02 TeV");
      latex3.DrawLatex(0.57, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
      latex3.DrawLatex(0.80, 0.76, Form("|#eta| < %.1f", tkEta));
    }

  }

  cFinalUnfoldMerged->Update();
  cFinalUnfoldMerged->SaveAs("../plots/skew/cFinalUnfoldMerged.pdf");
  */

  //-- 4 x 3
  //-- Draw Final Unfolded distns on one plot
  TLegend * legUnfObs = new TLegend(0.55, 0.06, 0.75, 0.26);
  legInit(legUnfObs);

  TH1D * hDummy = new TH1D("hDummy", "hDummy", 152, 0-binw, vnMax[norder_]);

  TCanvas * cFinalUnfoldMerged = new TCanvas("cFinalUnfoldMerged", "cFinalUnfoldMerged", 2000, 1500);
  cFinalUnfoldMerged->SetLeftMargin(0.18);
  cFinalUnfoldMerged->SetRightMargin(0.01);
  cFinalUnfoldMerged->SetTopMargin(0.1); //-0.07
  cFinalUnfoldMerged->SetBottomMargin(0.19);  
  cFinalUnfoldMerged->Modified();
  cFinalUnfoldMerged->Update();
  cFinalUnfoldMerged->Divide(4,3,0,0);

  for(int icent = 0; icent < NCENT; icent++){

    cFinalUnfoldMerged->cd(icent+1);
    cFinalUnfoldMerged->cd(icent+1)->SetLogy();
    //if(icent < 4) cFinalUnfoldMerged->cd(icent+1)->SetTopMargin(0.07);
    //if(icent > 7) cFinalUnfoldMerged->cd(icent+1)->SetBottomMargin(0.22);
    //if(icent == 0 || icent == 4 || icent == 8)  cFinalUnfoldMerged->cd(icent+1)->SetLeftMargin(0.20);
    //if(icent == 3 || icent == 7 || icent == 11) cFinalUnfoldMerged->cd(icent+1)->SetRightMargin(0.0);

    //-- X axes 
    double m = 0.37;
    hDummy->GetXaxis()->SetTitle("v_{2}");
    hDummy->GetXaxis()->CenterTitle();
    hDummy->GetXaxis()->SetRange(0, hFinalUnfoldSys[icent]->FindBin(m));
    hDummy->GetXaxis()->SetNdivisions(507);
    hDummy->GetXaxis()->SetLabelFont(43);
    hDummy->GetXaxis()->SetLabelSize(38);
    hDummy->GetXaxis()->SetTitleFont(43);
    hDummy->GetXaxis()->SetTitleSize(47);
    hDummy->GetXaxis()->SetTitleOffset(2.5);


    //-- Y axes
    hDummy->GetYaxis()->SetTitle("p(v_{2})");
    hDummy->GetYaxis()->CenterTitle();
    hDummy->SetMaximum(0.5);
    hDummy->SetMinimum(2e-5);
    hDummy->GetYaxis()->SetLabelFont(43);
    hDummy->GetYaxis()->SetLabelSize(38);
    hDummy->GetYaxis()->SetTitleFont(43);
    hDummy->GetYaxis()->SetTitleSize(45);
    hDummy->GetYaxis()->SetTitleOffset(3.4);

    //-- Colors/Cosmetics
    hObs[icent]->Scale(1./hObs[icent]->Integral());
    hObs[icent]->SetLineColor(1);
    hObs[icent]->SetMarkerColor(1);
    hObs[icent]->SetMarkerStyle(24);

    hFinalUnfoldSys[icent]->SetLineColor(kRed+1);
    hFinalUnfoldSys[icent]->SetMarkerColor(kRed+1);
    hFinalUnfoldSys[icent]->SetFillColorAlpha(kRed-7, 0.6);

    hFinalUnf[icent]->SetLineColor(kRed+1);
    hFinalUnf[icent]->SetMarkerColor(kRed+1);

    if(drawObs){
      hDummy->Draw();
      hFinalUnfoldSys[icent]->Draw("e2same");
      hFinalUnf[icent]->Draw("same");
      hObs[icent]->Draw("same");

      //-- Add legend
      if(icent == 0){
	legUnfObs->AddEntry(hObs[icent],      "Observed p(v_{2})", "lp");
	legUnfObs->AddEntry(hFinalUnf[icent], "Unfold p(v_{2})",   "lp");
	legUnfObs->Draw("same");
	legUnfObs->SetTextFont(43);
	legUnfObs->SetTextSize(33);
      }
    }
    else{
      hDummy->Draw();
      hFinalUnfoldSys[icent]->Draw("e2same");
      hFinalUnf[icent]->Draw("same");
    }

    //-- Centrality tags
    if(icent == 0)                   latex4.DrawLatex(0.25, 0.06, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent < 4)               latex4.DrawLatex(0.05, 0.06, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent == 4)              latex4.DrawLatex(0.25, 0.06, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent > 4 && icent < 8)  latex4.DrawLatex(0.05, 0.06, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent == 8)              latex4.DrawLatex(0.25, 0.27, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    else if(icent > 8 && icent < 12) latex4.DrawLatex(0.05, 0.27, Form("#bf{%i - %i%s}", cent_min[icent], cent_max[icent], "%") );

    //-- pt, eta tags
    if(icent == 0){
      latex4.DrawLatex(0.455, 0.88, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
      latex4.DrawLatex(0.75, 0.78, Form("|#eta| < %.1f", tkEta));
    }

  }

  //-- CMS, PbPb, pt, eta tags
  cFinalUnfoldMerged->cd(0);
  latex4.DrawLatex(0.07, 0.973, "#bf{CMS}");
  latex4.DrawLatex(0.18, 0.973, "PbPb 5.02 TeV");

  cFinalUnfoldMerged->Update();
  cFinalUnfoldMerged->SaveAs("../plots/skew/cFinalUnfoldMerged.pdf");



}

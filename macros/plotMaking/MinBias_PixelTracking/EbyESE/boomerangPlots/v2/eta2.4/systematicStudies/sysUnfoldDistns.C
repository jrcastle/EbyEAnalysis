#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void sysUnfoldDistns(){

  int norder_     = 2;
  double tkEta    = 2.4;

  bool smooth = 1;

  TFile * fAna;
  TH1D * hObs[NCENT];
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

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

  TH1D * hFinalUnf[NCENT];
  TH1D * hFinalUnfSysReg[NCENT];
  TH1D * hFinalUnfSysResp[NCENT];
  TH1D * hFinalUnfSysVtx3[NCENT];
  TH1D * hFinalUnfSysVtx3_15[NCENT];
  TH1D * hFinalUnfSysTkLoose[NCENT];
  TH1D * hFinalUnfSysTkTight[NCENT];

  TH1D * hSysReg_RatioToNominal[NCENT];
  TH1D * hSysVtx3_RatioToNominal[NCENT];
  TH1D * hSysVtx3_15_RatioToNominal[NCENT];
  TH1D * hSysTkLoose_RatioToNominal[NCENT];
  TH1D * hSysTkTight_RatioToNominal[NCENT];

  TH1D * hSysErrPointToPoint_Reg[NCENT];
  TH1D * hSysErrPointToPoint_Vtx[NCENT];
  TH1D * hSysErrPointToPoint_TkQ[NCENT];
  TH1D * hSysErrPointToPoint_Tot[NCENT];
  TF1 * fSmooth[NCENT];
  double vnMax[NCENT];

  TFile * fOut;
  TH1D * hFinalUnfoldSys[NCENT];
  TH1D * hFinalUnfoldStatAndSys[NCENT];

  TCanvas * cSysReg;
  TCanvas * cSysVtx;
  TCanvas * cSysTkQ;
  TCanvas * cSysResp;
  TCanvas * cFinalUnfold;
  TCanvas * cSysPointToPoint;

  TLine * upperCut_5Sigma[NCENT];
  TLine * upperCut_4Sigma[NCENT];
  TLine * upperCut_3Sigma[NCENT];



  TLatex latex;

  //
  // MAIN
  // 
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");

  //-- Get the output files for the reported distributions
  fAna    = new TFile( "../AnalyzerResults/CastleEbyE.root" );
  fUnfold = new TFile( Form("../UnfoldResults/dataResp/data%i.root", norder_) );

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

  //-- Set up the output file
  fOut = new TFile(Form("SysUnfoldDistns_v%i.root", norder_), "recreate");

  cSysReg = new TCanvas("cSysReg", "cSysReg", 2000, 1500);
  cSysReg->Divide(4,3);

  cSysVtx = new TCanvas("cSysVtx", "cSysVtx", 2000, 1500);
  cSysVtx->Divide(4,3);

  cSysTkQ = new TCanvas("cSysTkQ", "cSysTkQ", 2000, 1500);
  cSysTkQ->Divide(4,3);

  cFinalUnfold = new TCanvas("cFinalUnfold", "cFinalUnfold", 2000, 1500);
  cFinalUnfold->Divide(4,3);

  //cSysPointToPoint = new TCanvas("cSysPointToPoint", "cSysPointToPoint", 2000, 1500);
  //cSysPointToPoint->Divide(4,3);

  //-- Start fetching histograms....
  for(int icent = 0; icent < NCENT; icent++){

    //-- -------------- Real Dists --------------
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    bool iterCut    = 0;
    bool iterCutReg = 0;
    for(int i = 0; i < NITER; i++){
      hUnfold[icent][i]     = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldResp[icent][i] = (TH1D*) fUnfoldResp->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i]     = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
	hFinalUnf[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnf_c%i", icent) );
	hFinalUnfSysResp[icent] = (TH1D*) hUnfoldResp[icent][i]->Clone( Form("hFinalUnfSysResp_c%i", icent) );
	iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
	hFinalUnf[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnf_c%i", icent) );
	hFinalUnfSysResp[icent] = (TH1D*) hUnfoldResp[icent][i]->Clone( Form("hFinalUnfSysResp_c%i", icent) );
      }

      //-- -------------- Reg --------------
      if( chi2 < 2. && !iterCutReg ){
        hFinalUnfSysReg[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfSysReg_c%i", icent) );
        iterCutReg = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysReg[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfSysReg_c%i", icent) );
      }

    }

    //-- -------------- Vtx3 --------------
    hObsVtx3[icent] = (TH1D*) fAnaVtx3->Get( Form("qwebye/hVnFull_c%i", icent) );
    iterCut= 0;
    for(int i = 0; i < NITER; i++){
      hUnfoldVtx3[icent][i] = (TH1D*) fUnfoldVtx3->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefoldVtx3[icent][i] = (TH1D*) fUnfoldVtx3->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefoldVtx3[icent][i]->Chi2Test(hObsVtx3[icent], "CHI2/NDF");

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
      double chi2 = hRefoldVtx3_15[icent][i]->Chi2Test(hObsVtx3_15[icent], "CHI2/NDF");

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
      double chi2 = hRefoldTkLoose[icent][i]->Chi2Test(hObsTkLoose[icent], "CHI2/NDF");

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
      double chi2 = hRefoldTkTight[icent][i]->Chi2Test(hObsTkTight[icent], "CHI2/NDF");

      if( chi2 < 1.2 && !iterCut ){
        hFinalUnfSysTkTight[icent] = (TH1D*) hUnfoldTkTight[icent][i]->Clone( Form("hFinalUnfSysTkTight_c%i", icent) );
        iterCut = 1;
      }
      if( i == NITER - 1 && !iterCut){
        hFinalUnfSysTkTight[icent] = (TH1D*) hUnfoldTkTight[icent][i]->Clone( Form("hFinalUnfSysTkTight_c%i", icent) );
      }
    }

  } //-- End Cent loop

  //-- Done fetching files and determining the final unfolded distns, now let's calculate some ratios!
  for(int icent = 0; icent < NCENT; icent++){


    //-- Fix Unfolds
    FixUnfold( hFinalUnf[icent] );
    FixUnfold( hFinalUnfSysResp[icent] );
    FixUnfold( hFinalUnfSysReg[icent] );
    FixUnfold( hFinalUnfSysVtx3[icent] );
    FixUnfold( hFinalUnfSysVtx3_15[icent] );
    FixUnfold( hFinalUnfSysTkLoose[icent] );
    FixUnfold( hFinalUnfSysTkTight[icent] );

    //-- Normalize
    hFinalUnf[icent]->Scale( 1. / hFinalUnf[icent]->Integral() );
    hFinalUnfSysResp[icent]->Scale( 1. / hFinalUnfSysResp[icent]->Integral() );
    hFinalUnfSysReg[icent]->Scale( 1. / hFinalUnfSysReg[icent]->Integral() );
    hFinalUnfSysVtx3[icent]->Scale( 1. / hFinalUnfSysVtx3[icent]->Integral() );
    hFinalUnfSysVtx3_15[icent]->Scale( 1. / hFinalUnfSysVtx3_15[icent]->Integral() );
    hFinalUnfSysTkLoose[icent]->Scale( 1. / hFinalUnfSysTkLoose[icent]->Integral() );
    hFinalUnfSysTkTight[icent]->Scale( 1. / hFinalUnfSysTkTight[icent]->Integral() );

    //-- Clone originals
    hSysReg_RatioToNominal[icent]     = (TH1D*) hFinalUnfSysReg[icent]->Clone( Form("hSysReg_RatioToNominal_c%i", icent) );
    hSysVtx3_RatioToNominal[icent]    = (TH1D*) hFinalUnfSysVtx3[icent]->Clone( Form("hSysVtx3_RatioToNominal_c%i", icent) );
    hSysVtx3_15_RatioToNominal[icent] = (TH1D*) hFinalUnfSysVtx3_15[icent]->Clone( Form("hSysVtx3_15_RatioToNominal_c%i", icent) );
    hSysTkLoose_RatioToNominal[icent] = (TH1D*) hFinalUnfSysTkLoose[icent]->Clone( Form("hSysTkLoose_RatioToNominal_c%i", icent) );
    hSysTkTight_RatioToNominal[icent] = (TH1D*) hFinalUnfSysTkTight[icent]->Clone( Form("hSysTkTight_RatioToNominal_c%i", icent) );

    //-- Point to point uncert for smoothing
    hSysErrPointToPoint_Reg[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Reg_c%i", icent) );
    hSysErrPointToPoint_Vtx[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Vtx_c%i", icent) );
    hSysErrPointToPoint_TkQ[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_TkQ_c%i", icent) );
    hSysErrPointToPoint_Tot[icent] = (TH1D*) hFinalUnf[icent]->Clone( Form("hSysErrPointToPoint_Tot_c%i", icent) );

    hSysReg_RatioToNominal[icent]     -> Reset();
    hSysVtx3_RatioToNominal[icent]    -> Reset();
    hSysVtx3_15_RatioToNominal[icent] -> Reset();
    hSysTkLoose_RatioToNominal[icent] -> Reset();
    hSysTkTight_RatioToNominal[icent] -> Reset();

    hSysErrPointToPoint_Reg[icent] -> Reset();
    hSysErrPointToPoint_Vtx[icent] -> Reset();
    hSysErrPointToPoint_TkQ[icent] -> Reset();
    hSysErrPointToPoint_Tot[icent] -> Reset();

    //-- Rescale to the same mean as nominal
    double meanvn_Nominal = hFinalUnf[icent]->GetMean();
    double meanvn_Reg     = hFinalUnfSysReg[icent]->GetMean();
    double meanvn_Vtx3    = hFinalUnfSysVtx3[icent]->GetMean();
    double meanvn_Vtx3_15 = hFinalUnfSysVtx3_15[icent]->GetMean();
    double meanvn_TkLoose = hFinalUnfSysTkLoose[icent]->GetMean();
    double meanvn_TkTight = hFinalUnfSysTkTight[icent]->GetMean();

    for(int i = 1; i <= NBins; i++){

      double vn = hFinalUnf[icent]->GetBinCenter(i);
      double pvn_Nominal = hFinalUnf[icent]->GetBinContent(i);
      double pvn_Reg     = hFinalUnfSysReg[icent]->GetBinContent(i);
      double pvn_Vtx3    = hFinalUnfSysVtx3[icent]->GetBinContent(i);
      double pvn_Vtx3_15 = hFinalUnfSysVtx3_15[icent]->GetBinContent(i);
      double pvn_TkLoose = hFinalUnfSysTkLoose[icent]->GetBinContent(i);
      double pvn_TkTight = hFinalUnfSysTkTight[icent]->GetBinContent(i);

      double pvne_Nominal = hFinalUnf[icent]->GetBinError(i);
      double pvne_Reg     = hFinalUnfSysReg[icent]->GetBinError(i);
      double pvne_Vtx3    = hFinalUnfSysVtx3[icent]->GetBinError(i);
      double pvne_Vtx3_15 = hFinalUnfSysVtx3_15[icent]->GetBinError(i);
      double pvne_TkLoose = hFinalUnfSysTkLoose[icent]->GetBinError(i);
      double pvne_TkTight = hFinalUnfSysTkTight[icent]->GetBinError(i);

      double shiftedvn_Reg     = vn + ( meanvn_Nominal - meanvn_Reg );
      double shiftedvn_Vtx3    = vn + ( meanvn_Nominal - meanvn_Vtx3 );
      double shiftedvn_Vtx3_15 = vn + ( meanvn_Nominal - meanvn_Vtx3_15 );
      double shiftedvn_TkLoose = vn + ( meanvn_Nominal - meanvn_TkLoose );
      double shiftedvn_TkTight = vn + ( meanvn_Nominal - meanvn_TkTight );

      int shiftedvnBIN_Reg     = hSysReg_RatioToNominal[icent]->FindBin( shiftedvn_Reg );
      int shiftedvnBIN_Vtx3    = hSysVtx3_RatioToNominal[icent]->FindBin( shiftedvn_Vtx3 );
      int shiftedvnBIN_Vtx3_15 = hSysVtx3_15_RatioToNominal[icent]->FindBin( shiftedvn_Vtx3_15 );
      int shiftedvnBIN_TkLoose = hSysTkLoose_RatioToNominal[icent]->FindBin( shiftedvn_TkLoose );
      int shiftedvnBIN_TkTight = hSysTkTight_RatioToNominal[icent]->FindBin( shiftedvn_TkTight );

      hSysReg_RatioToNominal[icent]     -> SetBinContent(shiftedvnBIN_Reg, pvn_Reg);
      hSysVtx3_RatioToNominal[icent]    -> SetBinContent(shiftedvnBIN_Vtx3, pvn_Vtx3);
      hSysVtx3_15_RatioToNominal[icent] -> SetBinContent(shiftedvnBIN_Vtx3_15, pvn_Vtx3_15);
      hSysTkLoose_RatioToNominal[icent] -> SetBinContent(shiftedvnBIN_TkLoose, pvn_TkLoose);
      hSysTkTight_RatioToNominal[icent] -> SetBinContent(shiftedvnBIN_TkTight, pvn_TkTight);

      hSysReg_RatioToNominal[icent]     -> SetBinError(shiftedvnBIN_Reg, pvne_Reg);
      hSysVtx3_RatioToNominal[icent]    -> SetBinError(shiftedvnBIN_Vtx3, pvne_Vtx3);
      hSysVtx3_15_RatioToNominal[icent] -> SetBinError(shiftedvnBIN_Vtx3_15, pvne_Vtx3_15);
      hSysTkLoose_RatioToNominal[icent] -> SetBinError(shiftedvnBIN_TkLoose, pvne_TkLoose);
      hSysTkTight_RatioToNominal[icent] -> SetBinError(shiftedvnBIN_TkTight, pvne_TkTight);

    }


    //-- Calculate ratio to nominal
    hSysReg_RatioToNominal[icent]     -> Divide( hFinalUnf[icent] );
    hSysVtx3_RatioToNominal[icent]    -> Divide( hFinalUnf[icent] );
    hSysVtx3_15_RatioToNominal[icent] -> Divide( hFinalUnf[icent] );
    hSysTkLoose_RatioToNominal[icent] -> Divide( hFinalUnf[icent] );
    hSysTkTight_RatioToNominal[icent] -> Divide( hFinalUnf[icent] );

    hSysReg_RatioToNominal[icent]     -> SetMinimum( 0.8 );
    hSysVtx3_RatioToNominal[icent]    -> SetMinimum( 0.8 );
    hSysVtx3_15_RatioToNominal[icent] -> SetMinimum( 0.8 );
    hSysTkLoose_RatioToNominal[icent] -> SetMinimum( 0.8 );
    hSysTkTight_RatioToNominal[icent] -> SetMinimum( 0.8 );

    hSysReg_RatioToNominal[icent]     -> SetMaximum( 1.2 );
    hSysVtx3_RatioToNominal[icent]    -> SetMaximum( 1.2 );
    hSysVtx3_15_RatioToNominal[icent] -> SetMaximum( 1.2 );
    hSysTkLoose_RatioToNominal[icent] -> SetMaximum( 1.2 );
    hSysTkTight_RatioToNominal[icent] -> SetMaximum( 1.2 );

    hSysErrPointToPoint_Reg[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_Vtx[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_TkQ[icent]->SetMaximum(0.01);
    hSysErrPointToPoint_Tot[icent]->SetMaximum(0.01);

    hSysReg_RatioToNominal[icent]     -> SetLineColor(kOrange+2);
    hSysVtx3_RatioToNominal[icent]    -> SetLineColor(2);
    hSysVtx3_15_RatioToNominal[icent] -> SetLineColor(4);
    hSysTkLoose_RatioToNominal[icent] -> SetLineColor(8);
    hSysTkTight_RatioToNominal[icent] -> SetLineColor(9);

    hSysReg_RatioToNominal[icent]     -> SetMarkerColor(kOrange+2);
    hSysVtx3_RatioToNominal[icent]    -> SetMarkerColor(2);
    hSysVtx3_15_RatioToNominal[icent] -> SetMarkerColor(4);
    hSysTkLoose_RatioToNominal[icent] -> SetMarkerColor(8);
    hSysTkTight_RatioToNominal[icent] -> SetMarkerColor(9);

    hSysReg_RatioToNominal[icent]     -> SetMarkerStyle(34);
    hSysVtx3_RatioToNominal[icent]    -> SetMarkerStyle(20);
    hSysVtx3_15_RatioToNominal[icent] -> SetMarkerStyle(24);
    hSysTkLoose_RatioToNominal[icent] -> SetMarkerStyle(21);
    hSysTkTight_RatioToNominal[icent] -> SetMarkerStyle(25);

    hSysReg_RatioToNominal[icent]     -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysVtx3_RatioToNominal[icent]    -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysVtx3_15_RatioToNominal[icent] -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysTkLoose_RatioToNominal[icent] -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );
    hSysTkTight_RatioToNominal[icent] -> GetXaxis() -> SetTitle( Form("v_{%i}", norder_ ) );

    hSysReg_RatioToNominal[icent]     -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysVtx3_RatioToNominal[icent]    -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysVtx3_15_RatioToNominal[icent] -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysTkLoose_RatioToNominal[icent] -> GetYaxis() -> SetTitle( "Ratio to Default" );
    hSysTkTight_RatioToNominal[icent] -> GetYaxis() -> SetTitle( "Ratio to Default" );


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

      vnMax[icent] = hFinalUnfoldSys[icent]->GetBinCenter(i);

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

      double relErr_Vtx       = std::max(relErr_Vtx3, relErr_Vtx3_15);
      double relErr_TkQuality = std::max(relErr_TkLoose, relErr_TkLoose);

      double sysErrReg = relErr_Reg * hFinalUnf[icent]->GetBinContent(i);
      double sysErrVtx = relErr_Vtx * hFinalUnf[icent]->GetBinContent(i);
      double sysErrTkQ = relErr_TkQuality * hFinalUnf[icent]->GetBinContent(i);

      double relErrTot = sqrt( pow(relErr_Reg, 2) + pow(relErr_Vtx, 2) + pow(relErr_TkQuality, 2) );
      double statErr    = hFinalUnf[icent]->GetBinError(i);
      double respErr    = hFinalUnfSysResp[icent]->GetBinError(i);
      double sysErrRVT  = relErrTot * hFinalUnf[icent]->GetBinContent(i);
      double sysErrResp = sqrt( pow(respErr,2) - pow(statErr,2) );
      double sysErrTot  = sqrt( sysErrRVT*sysErrRVT + sysErrResp*sysErrResp );
      double totErr     = sqrt(pow(sysErrTot, 2) + pow(statErr, 2));

      hSysErrPointToPoint_Reg[icent]->SetBinContent(i, sysErrReg );
      hSysErrPointToPoint_Vtx[icent]->SetBinContent(i, sysErrVtx );
      hSysErrPointToPoint_TkQ[icent]->SetBinContent(i, sysErrTkQ );
      hSysErrPointToPoint_Tot[icent]->SetBinContent(i, sysErrTot );

      hFinalUnfoldSys[icent]->SetBinError(i, sysErrTot);
      hFinalUnfoldStatAndSys[icent]->SetBinError(i, totErr);

    }

    if(smooth){

      hSysErrPointToPoint_Reg[icent]->Smooth();
      hSysErrPointToPoint_Vtx[icent]->Smooth();
      hSysErrPointToPoint_TkQ[icent]->Smooth();
      hSysErrPointToPoint_Tot[icent]->Smooth();

      for(int i = 1; i <= NBins; i++){
	if( hFinalUnf[icent]->GetBinContent(i) == 0 ) continue;
	double sysReg  = hSysErrPointToPoint_Reg[icent]->GetBinContent(i);
	double sysVtx  = hSysErrPointToPoint_Vtx[icent]->GetBinContent(i);
	double sysTkQ  = hSysErrPointToPoint_TkQ[icent]->GetBinContent(i);

	double statErr    = hFinalUnf[icent]->GetBinError(i);
	double respErr = hFinalUnfSysResp[icent]->GetBinError(i);
	double sysResp = sqrt( pow(respErr,2) - pow(statErr,2) );

	double sys = sqrt( pow(sysResp,2) + pow(sysTkQ,2) + pow(sysVtx,2) + pow(sysReg,2) );
	double tot = sqrt( sys*sys + statErr*statErr );
	hFinalUnfoldSys[icent]->SetBinError(i, sys);
	hFinalUnfoldStatAndSys[icent]->SetBinError(i, tot);
      }
    }

    //fSmooth[icent] = new TF1(Form("fSmooth_c%i", icent), "pol2", 0.0, vnMax[icent]);
    //hSysErrPointToPoint[icent]->Fit( Form("fSmooth_c%i", icent), "QL0", "", 0.0, vnMax[icent]);

    fOut->cd();
    hFinalUnf[icent]->Write( Form( "hFinalUnfoldStat_c%i", icent) );
    hFinalUnfoldSys[icent]->Write();
    hFinalUnfoldStatAndSys[icent]->Write();

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

    cFinalUnfold->cd(icent+1);
    cFinalUnfold->cd(icent+1)->SetLogy();
    hFinalUnfoldSys[icent]->SetLineColor(2);
    hFinalUnfoldSys[icent]->SetMarkerColor(2);
    hFinalUnfoldSys[icent]->SetFillColor(15);
    hFinalUnfoldSys[icent]->GetXaxis()->SetRange(1, hFinalUnfoldSys[icent]->FindBin(0.3));
    hFinalUnfoldSys[icent]->Draw("e2");
    setex2->Draw();
    hFinalUnf[icent]->Draw("same");
    upperCut_5Sigma[icent]->Draw("same");
    upperCut_4Sigma[icent]->Draw("same");
    upperCut_3Sigma[icent]->Draw("same");
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );

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
  cFinalUnfold->SaveAs( Form("../plots/systematicStudies/sysFinalUnfoldDistns_v%i.pdf", norder_) );

}

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void makePhysics(){

  int norder_     = 2;
  int QnBinOrder_ = 2;
  double tkEta    = 2.4;

  bool dosys_      = 1;
  double sysWidth_ = 0.01;

  double mar   = 0.2;
  double offsx = 1.2;
  double offsy = 1.6;

  bool looseChi2IterCut   = 0;
  bool nominalChi2IterCut = 1;
  bool tightChi2IterCut   = 0;

  double qnmin = 0;
  double qnmax = 0.26;

  double rmsVnMin      = 0.;
  double rmsVnMax      = 0.17;
  double meanVnMin     = 0.;
  double meanVnMax     = 0.17;
  double relFluctVnMin = 0.0;
  double relFluctVnMax = 0.79;
  double stdDevVnMin    = 0.0;
  double stdDevVnMax    = 0.041;

  double vn2Min = 0.;
  double vn2Max = 0.2;
  double vn4Min = 0.;
  double vn4Max = 0.2;
  double vn6Min = 0.;
  double vn6Max = 0.2;
  double vn8Min = 0.;
  double vn8Max = 0.2;

  double vn6vn4Min = 0.95;
  double vn6vn4Max = 1.05;
  double vn8vn4Min = 0.95;
  double vn8vn4Max = 1.05;
  double vn8vn6Min = 0.99;
  double vn8vn6Max = 1.01;

  double g1expMin = -1.;
  double g1expMax = 1.2;

  TLatex latex;

  //-- Qn binning
  TFile * fQn;
  TH1D * hqbins[NCENT][NEPSymm];
  TH1D * hqnHF_EP[NCENT][NEPSymm];
  double qnBinCenter[NCENT][NEPSymm][NQN];
  double qnBinCentere[NCENT][NEPSymm][NQN];

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT][NEPSymm][NQN];

  //-- Unfolding output
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NITER];

  //-- Statististical Errors
  TFile * fStatErr;
  TH1D * hVarianceOfMean_RMSVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_MeanVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_RelFluctVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_StdDevVn[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn2[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn6[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Gamma1Exp[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn6Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8Vn4[NEPSymm][NQN];
  TH1D * hVarianceOfMean_Vn8Vn6[NEPSymm][NQN];

  //-- Systematic Errors
  TFile * fRelErr_Regularization;
  TH1D * relErr_Regularization_RMSVn[NEPSymm][NQN];
  TH1D * relErr_Regularization_MeanVn[NEPSymm][NQN];
  TH1D * relErr_Regularization_StDevVn[NEPSymm][NQN];
  TH1D * relErr_Regularization_RelFluctVn[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn2[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn4[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn6[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn8[NEPSymm][NQN];
  TH1D * relErr_Regularization_Gamma1Exp[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn6Vn4[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn8Vn4[NEPSymm][NQN];
  TH1D * relErr_Regularization_Vn8Vn6[NEPSymm][NQN];

  TFile * fRelErr_Vtx;
  TH1D * relErr_Vtx_RMSVn[NEPSymm][NQN];
  TH1D * relErr_Vtx_MeanVn[NEPSymm][NQN];
  TH1D * relErr_Vtx_StDevVn[NEPSymm][NQN];
  TH1D * relErr_Vtx_RelFluctVn[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn2[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn4[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn6[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn8[NEPSymm][NQN];
  TH1D * relErr_Vtx_Gamma1Exp[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn6Vn4[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn8Vn4[NEPSymm][NQN];
  TH1D * relErr_Vtx_Vn8Vn6[NEPSymm][NQN];

  TFile * fRelErr_TkQuality;
  TH1D * relErr_TkQuality_RMSVn[NEPSymm][NQN];
  TH1D * relErr_TkQuality_MeanVn[NEPSymm][NQN];
  TH1D * relErr_TkQuality_StDevVn[NEPSymm][NQN];
  TH1D * relErr_TkQuality_RelFluctVn[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn2[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn4[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn6[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn8[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Gamma1Exp[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn6Vn4[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn8Vn4[NEPSymm][NQN];
  TH1D * relErr_TkQuality_Vn8Vn6[NEPSymm][NQN];

  //-- Final Unfolded iteration
  double chi2Cut;
  int iterCut[NCENT][NEPSymm][NQN];

  //-- Moment vs qn
  double rmsVn_vs_qn[NCENT][NEPSymm][NQN];
  double meanVn_vs_qn[NCENT][NEPSymm][NQN];
  double relFluctVn_vs_qn[NCENT][NEPSymm][NQN];
  double stdDevVn_vs_qn[NCENT][NEPSymm][NQN];

  double rmsVn_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double meanVn_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double relFluctVn_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double stdDevVn_vs_qn_statErr[NCENT][NEPSymm][NQN];

  double rmsVn_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double meanVn_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double relFluctVn_vs_qn_sysErr[NCENT][NEPSymm][NQN];

  TGraphErrors * grRMSVnVSQn[NCENT][NEPSymm];
  TGraphErrors * grMeanVnVSQn[NCENT][NEPSymm];
  TGraphErrors * grRelFluctVnVSQn[NCENT][NEPSymm];
  TGraphErrors * grStdDevVnVSQn[NCENT][NEPSymm];

  TGraphErrors * grRMSVnVSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grMeanVnVSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grRelFluctVnVSQn_DoSys[NCENT][NEPSymm];

  //-- Cumulant vs qn
  double vn2_vs_qn[NCENT][NEPSymm][NQN];
  double vn4_vs_qn[NCENT][NEPSymm][NQN];
  double vn6_vs_qn[NCENT][NEPSymm][NQN];
  double vn8_vs_qn[NCENT][NEPSymm][NQN];

  double vn2_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double vn4_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double vn6_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double vn8_vs_qn_statErr[NCENT][NEPSymm][NQN];

  double vn2_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double vn4_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double vn6_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double vn8_vs_qn_sysErr[NCENT][NEPSymm][NQN];

  TGraphErrors * grVn2VSQn[NCENT][NEPSymm];
  TGraphErrors * grVn4VSQn[NCENT][NEPSymm];
  TGraphErrors * grVn6VSQn[NCENT][NEPSymm];
  TGraphErrors * grVn8VSQn[NCENT][NEPSymm];

  TGraphErrors * grVn2VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grVn4VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grVn6VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grVn8VSQn_DoSys[NCENT][NEPSymm];

  //-- Cumu Ratio and Gamma1Exp vs qn
  double vn6vn4_vs_qn[NCENT][NEPSymm][NQN];
  double vn8vn4_vs_qn[NCENT][NEPSymm][NQN];
  double vn8vn6_vs_qn[NCENT][NEPSymm][NQN];
  double g1exp_vs_qn[NCENT][NEPSymm][NQN];

  double vn6vn4_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double vn8vn4_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double vn8vn6_vs_qn_statErr[NCENT][NEPSymm][NQN];
  double g1exp_vs_qn_statErr[NCENT][NEPSymm][NQN];

  double vn6vn4_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double vn8vn4_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double vn8vn6_vs_qn_sysErr[NCENT][NEPSymm][NQN];
  double g1exp_vs_qn_sysErr[NCENT][NEPSymm][NQN];

  TGraphErrors * grVn6Vn4VSQn[NCENT][NEPSymm];
  TGraphErrors * grVn8Vn4VSQn[NCENT][NEPSymm];
  TGraphErrors * grVn8Vn6VSQn[NCENT][NEPSymm];
  TGraphErrors * grG1ExpVSQn[NCENT][NEPSymm];

  TGraphErrors * grVn6Vn4VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grVn8Vn4VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grVn8Vn6VSQn_DoSys[NCENT][NEPSymm];
  TGraphErrors * grG1ExpVSQn_DoSys[NCENT][NEPSymm];

  //-- Canvases
  TCanvas * cMomentVSQn[NEPSymm];
  TCanvas * cRMSVSQn_Big[NEPSymm];
  TCanvas * cMeanVSQn_Big[NEPSymm];
  TCanvas * cRelFluctVSQn_Big[NEPSymm];

  TCanvas * cCumuVSQn[NEPSymm];
  TCanvas * cVn2VSQn_Big[NEPSymm];
  TCanvas * cVn4VSQn_Big[NEPSymm];
  TCanvas * cVn6VSQn_Big[NEPSymm];
  TCanvas * cVn8VSQn_Big[NEPSymm];

  TCanvas * cCumuRatioVSQn[NEPSymm];
  TCanvas * cVn6Vn4VSQn_Big[NEPSymm];
  TCanvas * cVn8Vn4VSQn_Big[NEPSymm];
  TCanvas * cVn8Vn6VSQn_Big[NEPSymm];
  TCanvas * cG1ExpVSQn_Big[NEPSymm];


  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();

  //-- Set the chi2 cutoff
  bool c2c1 = looseChi2IterCut   && nominalChi2IterCut;
  bool c2c2 = looseChi2IterCut   && tightChi2IterCut;
  bool c2c3 = nominalChi2IterCut && tightChi2IterCut;
  bool c2c4 = looseChi2IterCut   && nominalChi2IterCut && tightChi2IterCut;

  if( c2c1 || c2c2 || c2c3 || c2c4){
    std::cout<<"WARNING! More than one chi2 cutoff scenario defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  if( looseChi2IterCut )   chi2Cut = 1.5;
  if( nominalChi2IterCut ) chi2Cut = 1.2;
  if( tightChi2IterCut )   chi2Cut = 1.0;

  //-- Get the QN binning file
  fQn = new TFile( Form( "../../v%i/eta2.4/AnalyzerResults/q%iCuts.root", QnBinOrder_, QnBinOrder_) );

  //-- Get the Analyzer output file
  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );

  //-- Get Unfolding output
  fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get Statististical Errors
  fStatErr = new TFile( Form("../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", norder_, norder_) );
  for(int iEP = 0; iEP < NEPSymm; iEP++){
    if( iEP != EPSymmBin ) continue;
    for(int iqn = 0; iqn < NQN; iqn++){
      hVarianceOfMean_RMSVn[iEP][iqn]      = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_RMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_MeanVn[iEP][iqn]     = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_MeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_RelFluctVn[iEP][iqn] = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_RelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_StdDevVn[iEP][iqn]   = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_StDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn2[iEP][iqn]        = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn2_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn4[iEP][iqn]        = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn4_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn6[iEP][iqn]        = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn6_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn8[iEP][iqn]        = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn8_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Gamma1Exp[iEP][iqn]  = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Gamma1Exp_%s_qbin%i",  EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn6Vn4[iEP][iqn]     = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn6Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn8Vn4[iEP][iqn]     = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn8Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
      hVarianceOfMean_Vn8Vn6[iEP][iqn]     = (TH1D*) fStatErr->Get( Form("hVarianceOfMean_Vn8Vn6_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
    }
  }

  //-- Get Systematic Errors
  double sysWidth[NQN];
  if( dosys_ ){
    //-- Systematic Errors
    fRelErr_Regularization = new TFile("systematicStudies/chi2Cutoff/relErrorRegularization.root");
    fRelErr_Vtx            = new TFile("systematicStudies/vtxCut/relErrorVtx.root");
    fRelErr_TkQuality      = new TFile("systematicStudies/tkQuality/relErrorTkQuality.root");

    //-- Widths for systematic error bars
    for(int iqn = 0; iqn < NQN; iqn++) sysWidth[iqn] = sysWidth_;

    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	relErr_Regularization_RMSVn[iEP][iqn]      = (TH1D*) fRelErr_Regularization->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_MeanVn[iEP][iqn]     = (TH1D*) fRelErr_Regularization->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_StDevVn[iEP][iqn]    = (TH1D*) fRelErr_Regularization->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_RelFluctVn[iEP][iqn] = (TH1D*) fRelErr_Regularization->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn2[iEP][iqn]        = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn2_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn4[iEP][iqn]        = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn4_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn6[iEP][iqn]        = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn6_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn8[iEP][iqn]        = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn8_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Gamma1Exp[iEP][iqn]  = (TH1D*) fRelErr_Regularization->Get( Form("relErrGamma1Exp_%s_qbin%i",  EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn6Vn4[iEP][iqn]     = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn6Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn8Vn4[iEP][iqn]     = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn8Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Regularization_Vn8Vn6[iEP][iqn]     = (TH1D*) fRelErr_Regularization->Get( Form("relErrVn8Vn6_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );


	relErr_Vtx_RMSVn[iEP][iqn]      = (TH1D*) fRelErr_Vtx->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_MeanVn[iEP][iqn]     = (TH1D*) fRelErr_Vtx->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_StDevVn[iEP][iqn]    = (TH1D*) fRelErr_Vtx->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_RelFluctVn[iEP][iqn] = (TH1D*) fRelErr_Vtx->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn2[iEP][iqn]        = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn2_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn4[iEP][iqn]        = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn4_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn6[iEP][iqn]        = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn6_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn8[iEP][iqn]        = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn8_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Gamma1Exp[iEP][iqn]  = (TH1D*) fRelErr_Vtx->Get( Form("relErrGamma1Exp_%s_qbin%i",  EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn6Vn4[iEP][iqn]     = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn6Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn8Vn4[iEP][iqn]     = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn8Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_Vtx_Vn8Vn6[iEP][iqn]     = (TH1D*) fRelErr_Vtx->Get( Form("relErrVn8Vn6_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );


	relErr_TkQuality_RMSVn[iEP][iqn]      = (TH1D*) fRelErr_TkQuality->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_MeanVn[iEP][iqn]     = (TH1D*) fRelErr_TkQuality->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_StDevVn[iEP][iqn]    = (TH1D*) fRelErr_TkQuality->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_RelFluctVn[iEP][iqn] = (TH1D*) fRelErr_TkQuality->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn2[iEP][iqn]        = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn2_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn4[iEP][iqn]        = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn4_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn6[iEP][iqn]        = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn6_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn8[iEP][iqn]        = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn8_%s_qbin%i",        EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Gamma1Exp[iEP][iqn]  = (TH1D*) fRelErr_TkQuality->Get( Form("relErrGamma1Exp_%s_qbin%i",  EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn6Vn4[iEP][iqn]     = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn6Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn8Vn4[iEP][iqn]     = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn8Vn4_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	relErr_TkQuality_Vn8Vn6[iEP][iqn]     = (TH1D*) fRelErr_TkQuality->Get( Form("relErrVn8Vn6_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );

      }
    }

  }

  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){

    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;

      //-- Get the QN binning and the QN distributions
      hqbins[icent][iEP]   = (TH1D*) fQn->Get( Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent) );
      hqnHF_EP[icent][iEP] = (TH1D*) fQn->Get( Form("hqnHF_c%i_EP%i", icent, iEP) );

      for(int iqn = 0; iqn < NQN; iqn++){

	//-- Get the VN observed histogram
	hObs[icent][iEP][iqn] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Determine the SW qn values
	double qMin = hqbins[icent][iEP]->GetBinLowEdge(iqn+1);
	double qMax = hqbins[icent][iEP]->GetBinLowEdge(iqn+2);
	int binMin  = hqnHF_EP[icent][iEP]->FindBin(qMin);
	int binMax  = hqnHF_EP[icent][iEP]->FindBin(qMax);
	double qnSW = 0.;
	double sumw = 0;
	for(int ib = binMin; ib <= binMax; ib++){
	  double w  = hqnHF_EP[icent][iEP]->GetBinContent(ib);
	  double qn = hqnHF_EP[icent][iEP]->GetBinCenter(ib);
	  qnSW += w * qn;
	  sumw += w;
	}
	qnSW /= sumw;
	qnBinCenter[icent][iEP][iqn]  = qnSW;
	qnBinCentere[icent][iEP][iqn] = 0;

	iterCut[icent][iEP][iqn] = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the Refolded histograms
	  hRefold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Chi2 Cut
	  double chi2NDF_Refold = hRefold[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold <= chi2Cut ){
	    iterCut[icent][iEP][iqn] = i;
	    break;
	  }
	  if( i == NITER-1 ){
	    iterCut[icent][iEP][iqn] = i;
	    break;
	  }

	} //-- End unfold iteration loop


	//-- Fill Arrays for plots
	int it = iterCut[icent][iEP][iqn];
	FixUnfold( hUnfold[icent][iEP][iqn][it] );
	EbyECumu cumu( hUnfold[icent][iEP][iqn][it] );

	//-- Values
	double mean     = hUnfold[icent][iEP][iqn][it]->GetMean();
	if( icent == 6 && iqn == 9 ) std::cout<<"LOOK HERE!!!!!!!!\t"<< hUnfold[icent][iEP][iqn][it]->GetRMS() <<std::endl;

	double stdev    = hUnfold[icent][iEP][iqn][it]->GetRMS();
	double rms      = sqrt( mean*mean + stdev*stdev );
	double relFluct = stdev / mean;

	double vn2 = cumu.GetCumu_vn2();
	double vn4 = cumu.GetCumu_vn4();
	double vn6 = cumu.GetCumu_vn6();
	double vn8 = cumu.GetCumu_vn8();

	double vn6vn4;
	if( vn6 == 0 || vn4 == 0 ) vn6vn4 = -1;
	else                       vn6vn4 = vn6 / vn4;
	double vn8vn4;
	if( vn8 == 0 || vn4 == 0 ) vn8vn4 = -1;
	else                       vn8vn4 = vn8 / vn4;
	double vn8vn6;
	if( vn8 == 0 || vn6 == 0 ) vn8vn6 = -1;
	else                       vn8vn6 = vn8 / vn6;

	double g1exp = cumu.GetGamma1Exp();

	//-- Stat Errors
	double rmsStatErr      = sqrt( hVarianceOfMean_RMSVn[iEP][iqn]->GetBinContent(icent+1) );
	double meanStatErr     = sqrt( hVarianceOfMean_MeanVn[iEP][iqn]->GetBinContent(icent+1) );
	double relFluctStatErr = sqrt( hVarianceOfMean_RelFluctVn[iEP][iqn]->GetBinContent(icent+1) );
	double stdDevStatErr   = sqrt( hVarianceOfMean_StdDevVn[iEP][iqn]->GetBinContent(icent+1) );

	double vn2StatErr = sqrt( hVarianceOfMean_Vn2[iEP][iqn]->GetBinContent(icent+1) );
	double vn4StatErr = sqrt( hVarianceOfMean_Vn4[iEP][iqn]->GetBinContent(icent+1) );
	double vn6StatErr = sqrt( hVarianceOfMean_Vn6[iEP][iqn]->GetBinContent(icent+1) );
	double vn8StatErr = sqrt( hVarianceOfMean_Vn8[iEP][iqn]->GetBinContent(icent+1) );

	double vn6vn4StatErr = sqrt( hVarianceOfMean_Vn6Vn4[iEP][iqn]->GetBinContent(icent+1) );
	double vn8vn4StatErr = sqrt( hVarianceOfMean_Vn8Vn4[iEP][iqn]->GetBinContent(icent+1) );
	double vn8vn6StatErr = sqrt( hVarianceOfMean_Vn8Vn6[iEP][iqn]->GetBinContent(icent+1) );
	double g1expStatErr  = sqrt( hVarianceOfMean_Gamma1Exp[iEP][iqn]->GetBinContent(icent+1) );

	//-- Sys Errors
	if( dosys_ ){
	  //-- Individual studies
	  double re_Reg_RMSVn      = relErr_Regularization_RMSVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_MeanVn     = relErr_Regularization_MeanVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_StDevVn    = relErr_Regularization_StDevVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_RelFluctVn = relErr_Regularization_RelFluctVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn2        = relErr_Regularization_Vn2[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn4        = relErr_Regularization_Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn6        = relErr_Regularization_Vn6[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn8        = relErr_Regularization_Vn8[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Gamma1Exp  = relErr_Regularization_Gamma1Exp[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn6Vn4     = relErr_Regularization_Vn6Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn8Vn4     = relErr_Regularization_Vn8Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Reg_Vn8Vn6     = relErr_Regularization_Vn8Vn6[iEP][iqn]->GetBinContent(icent+1);

	  double re_Vtx_RMSVn      = relErr_Vtx_RMSVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_MeanVn     = relErr_Vtx_MeanVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_StDevVn    = relErr_Vtx_StDevVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_RelFluctVn = relErr_Vtx_RelFluctVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn2        = relErr_Vtx_Vn2[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn4        = relErr_Vtx_Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn6        = relErr_Vtx_Vn6[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn8        = relErr_Vtx_Vn8[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Gamma1Exp  = relErr_Vtx_Gamma1Exp[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn6Vn4     = relErr_Vtx_Vn6Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn8Vn4     = relErr_Vtx_Vn8Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_Vtx_Vn8Vn6     = relErr_Vtx_Vn8Vn6[iEP][iqn]->GetBinContent(icent+1);

	  double re_TkQ_RMSVn      = relErr_TkQuality_RMSVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_MeanVn     = relErr_TkQuality_MeanVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_StDevVn    = relErr_TkQuality_StDevVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_RelFluctVn = relErr_TkQuality_RelFluctVn[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn2        = relErr_TkQuality_Vn2[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn4        = relErr_TkQuality_Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn6        = relErr_TkQuality_Vn6[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn8        = relErr_TkQuality_Vn8[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Gamma1Exp  = relErr_TkQuality_Gamma1Exp[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn6Vn4     = relErr_TkQuality_Vn6Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn8Vn4     = relErr_TkQuality_Vn8Vn4[iEP][iqn]->GetBinContent(icent+1);
	  double re_TkQ_Vn8Vn6     = relErr_TkQuality_Vn8Vn6[iEP][iqn]->GetBinContent(icent+1);

	  //-- Merge
	  double sys_RMSVn      = rms * sqrt( pow(re_Reg_RMSVn,2) + pow(re_Vtx_RMSVn,2) + pow(re_TkQ_RMSVn,2) );
          double sys_MeanVn     = mean * sqrt( pow(re_Reg_MeanVn,2) + pow(re_Vtx_MeanVn,2) + pow(re_TkQ_MeanVn,2) );
          double sys_StDevVn    = stdev * sqrt( pow(re_Reg_StDevVn,2) + pow(re_Vtx_StDevVn,2) + pow(re_TkQ_StDevVn,2) );
          double sys_RelFluctVn = relFluct * sqrt( pow(re_Reg_RelFluctVn,2) + pow(re_Vtx_RelFluctVn,2) + pow(re_TkQ_RelFluctVn,2) );
          double sys_Vn2        = vn2 * sqrt( pow(re_Reg_Vn2,2) + pow(re_Vtx_Vn2,2) + pow(re_TkQ_Vn2,2) );
          double sys_Vn4        = vn4 * sqrt( pow(re_Reg_Vn4,2) + pow(re_Vtx_Vn4,2) + pow(re_TkQ_Vn4,2) );
          double sys_Vn6        = vn6 * sqrt( pow(re_Reg_Vn6,2) + pow(re_Vtx_Vn6,2) + pow(re_TkQ_Vn6,2) );
          double sys_Vn8        = vn8 * sqrt( pow(re_Reg_Vn8,2) + pow(re_Vtx_Vn8,2) + pow(re_TkQ_Vn8,2) );
          double sys_Gamma1Exp  = g1exp * sqrt( pow(re_Reg_Gamma1Exp,2) + pow(re_Vtx_Gamma1Exp,2) + pow(re_TkQ_Gamma1Exp,2) );
          double sys_Vn6Vn4     = vn6vn4 * sqrt( pow(re_Reg_Vn6Vn4,2) + pow(re_Vtx_Vn6Vn4,2) + pow(re_TkQ_Vn6Vn4,2) );
          double sys_Vn8Vn4     = vn8vn4 * sqrt( pow(re_Reg_Vn8Vn4,2) + pow(re_Vtx_Vn8Vn4,2) + pow(re_TkQ_Vn8Vn4,2) );
          double sys_Vn8Vn6     = vn8vn6 * sqrt( pow(re_Reg_Vn8Vn6,2) + pow(re_Vtx_Vn8Vn6,2) + pow(re_TkQ_Vn8Vn6,2) );

	  //-- Set
	  rmsVn_vs_qn_sysErr[icent][iEP][iqn]      = sys_RMSVn;
	  meanVn_vs_qn_sysErr[icent][iEP][iqn]     = sys_MeanVn;
	  relFluctVn_vs_qn_sysErr[icent][iEP][iqn] = sys_RelFluctVn;

	  vn2_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn2;
	  vn4_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn4;
	  vn6_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn6;
	  vn8_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn8;

	  vn6vn4_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn6Vn4;
	  vn8vn4_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn8Vn4;
	  vn8vn6_vs_qn_sysErr[icent][iEP][iqn] = sys_Vn8Vn6;
	  g1exp_vs_qn_sysErr[icent][iEP][iqn]  = sys_Gamma1Exp;


	} //End if( dosys_ )

	rmsVn_vs_qn[icent][iEP][iqn]      = rms;
	meanVn_vs_qn[icent][iEP][iqn]     = mean;
	relFluctVn_vs_qn[icent][iEP][iqn] = relFluct;
	stdDevVn_vs_qn[icent][iEP][iqn]   = stdev;

	rmsVn_vs_qn_statErr[icent][iEP][iqn]      = rmsStatErr;
	meanVn_vs_qn_statErr[icent][iEP][iqn]     = meanStatErr;
	relFluctVn_vs_qn_statErr[icent][iEP][iqn] = relFluctStatErr;
	stdDevVn_vs_qn_statErr[icent][iEP][iqn]   = stdDevStatErr;

	vn2_vs_qn[icent][iEP][iqn] = vn2;
	vn4_vs_qn[icent][iEP][iqn] = vn4;
	vn6_vs_qn[icent][iEP][iqn] = vn6;
	vn8_vs_qn[icent][iEP][iqn] = vn8;

	vn2_vs_qn_statErr[icent][iEP][iqn] = vn2StatErr;
	vn4_vs_qn_statErr[icent][iEP][iqn] = vn4StatErr;
	vn6_vs_qn_statErr[icent][iEP][iqn] = vn6StatErr;
	vn8_vs_qn_statErr[icent][iEP][iqn] = vn8StatErr;

	vn6vn4_vs_qn[icent][iEP][iqn] = vn6vn4;
	vn8vn4_vs_qn[icent][iEP][iqn] = vn8vn4;
	vn8vn6_vs_qn[icent][iEP][iqn] = vn8vn6;
	g1exp_vs_qn[icent][iEP][iqn]  = g1exp;

	vn6vn4_vs_qn_statErr[icent][iEP][iqn] = vn6vn4StatErr;
	vn8vn4_vs_qn_statErr[icent][iEP][iqn] = vn8vn4StatErr;
	vn8vn6_vs_qn_statErr[icent][iEP][iqn] = vn8vn6StatErr;
	g1exp_vs_qn_statErr[icent][iEP][iqn]  = g1expStatErr;

      } //-- End QN loop

      //-- Initialize TGraphErrors
      grRMSVnVSQn[icent][iEP]      = new TGraphErrors(NQN, qnBinCenter[icent][iEP], rmsVn_vs_qn[icent][iEP],      qnBinCentere[icent][iEP], rmsVn_vs_qn_statErr[icent][iEP]);
      grMeanVnVSQn[icent][iEP]     = new TGraphErrors(NQN, qnBinCenter[icent][iEP], meanVn_vs_qn[icent][iEP],     qnBinCentere[icent][iEP], meanVn_vs_qn_statErr[icent][iEP]);
      grRelFluctVnVSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], relFluctVn_vs_qn[icent][iEP], qnBinCentere[icent][iEP], relFluctVn_vs_qn_statErr[icent][iEP]);
      grStdDevVnVSQn[icent][iEP]   = new TGraphErrors(NQN, qnBinCenter[icent][iEP], stdDevVn_vs_qn[icent][iEP],   qnBinCentere[icent][iEP], stdDevVn_vs_qn_statErr[icent][iEP]);

      grVn2VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn2_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn2_vs_qn_statErr[icent][iEP]);
      grVn4VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn4_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn4_vs_qn_statErr[icent][iEP]);
      grVn6VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn6_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn6_vs_qn_statErr[icent][iEP]);
      grVn8VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn8_vs_qn_statErr[icent][iEP]);

      grVn6Vn4VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn6vn4_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn6vn4_vs_qn_statErr[icent][iEP]);
      grVn8Vn4VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8vn4_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn8vn4_vs_qn_statErr[icent][iEP]);
      grVn8Vn6VSQn[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8vn6_vs_qn[icent][iEP], qnBinCentere[icent][iEP], vn8vn6_vs_qn_statErr[icent][iEP]);
      grG1ExpVSQn[icent][iEP]  = new TGraphErrors(NQN, qnBinCenter[icent][iEP], g1exp_vs_qn[icent][iEP],  qnBinCentere[icent][iEP], g1exp_vs_qn_statErr[icent][iEP]);

      //-- Format Graphs
      formatGraph(grRMSVnVSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), rmsVnMin, rmsVnMax, Form("#sqrt{#LTv_{%i}^{2}#GT}", norder_), centCol[icent], centMark[icent], Form("grRMSVnVSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grMeanVnVSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), meanVnMin, meanVnMax, Form("#LTv_{%i}#GT", norder_), centCol[icent], centMark[icent], Form("grMeanVnVSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grRelFluctVnVSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), relFluctVnMin, relFluctVnMax, Form("#sigma_{v_{%i}}/#LTv_{%i}#GT", norder_, norder_), centCol[icent], centMark[icent], Form("grRelFluctVnVSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grStdDevVnVSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), stdDevVnMin, stdDevVnMax, Form("#sigma_{v_{%i}}", norder_), centCol[icent], centMark[icent], Form("grStdDevVnVSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn2VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn2Min, vn2Max, Form("v_{%i}{2}", norder_), centCol[icent], centMark[icent], Form("grVn2VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn4VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn4Min, vn4Max, Form("v_{%i}{4}", norder_), centCol[icent], centMark[icent], Form("grVn4VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn6VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn6Min, vn6Max, Form("v_{%i}{6}", norder_), centCol[icent], centMark[icent], Form("grVn6VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn8VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8Min, vn8Max, Form("v_{%i}{8}", norder_), centCol[icent], centMark[icent], Form("grVn8VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));

      formatGraph(grVn6Vn4VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn6vn4Min, vn6vn4Max, Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_), centCol[icent], centMark[icent], Form("grVn6Vn4VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn8Vn4VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8vn4Min, vn8vn4Max, Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_), centCol[icent], centMark[icent], Form("grVn8Vn4VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grVn8Vn6VSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8vn6Min, vn8vn6Max, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_), centCol[icent], centMark[icent], Form("grVn8Vn6VSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));
      formatGraph(grG1ExpVSQn[icent][iEP], Form("q_{%i}", QnBinOrder_), g1expMin, g1expMax, "#gamma_{1}^{exp}", centCol[icent], centMark[icent], Form("grG1ExpVSQn_%s_c%i", EPSymmNames[iEP].data(), icent ));


      grRMSVnVSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grMeanVnVSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grRelFluctVnVSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grStdDevVnVSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);

      grVn2VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grVn4VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grVn6VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grVn8VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);

      grVn6Vn4VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grVn8Vn4VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grVn8Vn6VSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
      grG1ExpVSQn[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);

      grRMSVnVSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grMeanVnVSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grRelFluctVnVSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grStdDevVnVSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);

      grVn2VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grVn4VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grVn6VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grVn8VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);

      grVn6Vn4VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grVn8Vn4VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grVn8Vn6VSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);
      grG1ExpVSQn[icent][iEP]->GetXaxis()->SetNdivisions(509);

      grRMSVnVSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grMeanVnVSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grRelFluctVnVSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grStdDevVnVSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);

      grVn2VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grVn4VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grVn6VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grVn8VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);

      grVn6Vn4VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grVn8Vn4VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grVn8Vn6VSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);
      grG1ExpVSQn[icent][iEP]->GetYaxis()->SetNdivisions(509);

      if( dosys_ ){

	//-- Initialize TGraphErrors
	grRMSVnVSQn_DoSys[icent][iEP]      = new TGraphErrors(NQN, qnBinCenter[icent][iEP], rmsVn_vs_qn[icent][iEP],      sysWidth, rmsVn_vs_qn_sysErr[icent][iEP]);
	grMeanVnVSQn_DoSys[icent][iEP]     = new TGraphErrors(NQN, qnBinCenter[icent][iEP], meanVn_vs_qn[icent][iEP],     sysWidth, meanVn_vs_qn_sysErr[icent][iEP]);
	grRelFluctVnVSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], relFluctVn_vs_qn[icent][iEP], sysWidth, relFluctVn_vs_qn_sysErr[icent][iEP]);

	grVn2VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn2_vs_qn[icent][iEP], sysWidth, vn2_vs_qn_sysErr[icent][iEP]);
	grVn4VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn4_vs_qn[icent][iEP], sysWidth, vn4_vs_qn_sysErr[icent][iEP]);
	grVn6VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn6_vs_qn[icent][iEP], sysWidth, vn6_vs_qn_sysErr[icent][iEP]);
	grVn8VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8_vs_qn[icent][iEP], sysWidth, vn8_vs_qn_sysErr[icent][iEP]);

	grVn6Vn4VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn6vn4_vs_qn[icent][iEP], sysWidth, vn6vn4_vs_qn_sysErr[icent][iEP]);
	grVn8Vn4VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8vn4_vs_qn[icent][iEP], sysWidth, vn8vn4_vs_qn_sysErr[icent][iEP]);
	grVn8Vn6VSQn_DoSys[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], vn8vn6_vs_qn[icent][iEP], sysWidth, vn8vn6_vs_qn_sysErr[icent][iEP]);
	grG1ExpVSQn_DoSys[icent][iEP]  = new TGraphErrors(NQN, qnBinCenter[icent][iEP], g1exp_vs_qn[icent][iEP],  sysWidth, g1exp_vs_qn_sysErr[icent][iEP]);

	//-- Format Graphs
	formatGraph(grRMSVnVSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), rmsVnMin, rmsVnMax, Form("#sqrt{#LTv_{%i}^{2}#GT}", norder_), 17, centMark[icent], Form("grRMSVnVSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grMeanVnVSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), meanVnMin, meanVnMax, Form("#LTv_{%i}#GT", norder_), 17, centMark[icent], Form("grMeanVnVSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grRelFluctVnVSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), relFluctVnMin, relFluctVnMax, Form("#sigma_{v_{%i}}/#LTv_{%i}#GT", norder_, norder_), 17, centMark[icent], Form("grRelFluctVnVSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));

	formatGraph(grVn2VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn2Min, vn2Max, Form("v_{%i}{2}", norder_), 17, centMark[icent], Form("grVn2VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grVn4VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn4Min, vn4Max, Form("v_{%i}{4}", norder_), 17, centMark[icent], Form("grVn4VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grVn6VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn6Min, vn6Max, Form("v_{%i}{6}", norder_), 17, centMark[icent], Form("grVn6VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grVn8VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8Min, vn8Max, Form("v_{%i}{8}", norder_), 17, centMark[icent], Form("grVn8VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));

	formatGraph(grVn6Vn4VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn6vn4Min, vn6vn4Max, Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_), 17, centMark[icent], Form("grVn6Vn4VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grVn8Vn4VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8vn4Min, vn8vn4Max, Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_), 17, centMark[icent], Form("grVn8Vn4VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grVn8Vn6VSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), vn8vn6Min, vn8vn6Max, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_), 17, centMark[icent], Form("grVn8Vn6VSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));
	formatGraph(grG1ExpVSQn_DoSys[icent][iEP], Form("q_{%i}", QnBinOrder_), g1expMin, g1expMax, "#gamma_{1}^{exp}", 17, centMark[icent], Form("grG1ExpVSQn_DoSys_%s_c%i", EPSymmNames[iEP].data(), icent ));


	grRMSVnVSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grMeanVnVSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grRelFluctVnVSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);

	grVn2VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grVn4VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grVn6VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grVn8VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);

	grVn6Vn4VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grVn8Vn4VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	grVn8Vn6VSQn_DoSys[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);


	grRMSVnVSQn_DoSys[icent][iEP]->SetFillColor(17);
	grMeanVnVSQn_DoSys[icent][iEP]->SetFillColor(17);
	grRelFluctVnVSQn_DoSys[icent][iEP]->SetFillColor(17);
	
	grVn2VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grVn4VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grVn6VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grVn8VSQn_DoSys[icent][iEP]->SetFillColor(17);
	
	grVn6Vn4VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grVn8Vn4VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grVn8Vn6VSQn_DoSys[icent][iEP]->SetFillColor(17);
	grG1ExpVSQn_DoSys[icent][iEP]->SetFillColor(17);
 
      } //-- End if( dosys_ )

    } //-- End EP loop
  } //-- End cent loop

  //-- DRAW!!!!!
  TLegend * legCent = new TLegend(0.46, 0.61, 0.95, 0.90);
  legCent->SetFillStyle(0);
  legCent->SetBorderSize(0);
  legCent->SetNColumns(2);
  for(int icent = 0; icent < NCENT; icent++) legCent->AddEntry(grRMSVnVSQn[icent][EPSymmBin], Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"), "lp");

  double BGrf = sqrt( 4./TMath::Pi() - 1);
  TLine * rfLine = new TLine(grRelFluctVnVSQn[0][EPSymmBin]->GetXaxis()->GetXmin(), BGrf, grRelFluctVnVSQn[0][EPSymmBin]->GetXaxis()->GetXmax(), BGrf);
  rfLine->SetLineColor(1);
  rfLine->SetLineStyle(2);
  rfLine->SetLineWidth(2);;

  //-- Initialize Canvases
  cMomentVSQn[EPSymmBin]       = new TCanvas(Form("cMomentVSQn_%s",       EPSymmNames[EPSymmBin].data()), Form("cMomentVSQn_%s",       EPSymmNames[EPSymmBin].data()), 1000, 1000);
  cRMSVSQn_Big[EPSymmBin]      = new TCanvas(Form("cRMSVSQn_Big_%s",      EPSymmNames[EPSymmBin].data()), Form("cRMSVSQn_Big_%s",      EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cMeanVSQn_Big[EPSymmBin]     = new TCanvas(Form("cMeanVSQn_Big_%s",     EPSymmNames[EPSymmBin].data()), Form("cMeanVSQn_Big_%s",     EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cRelFluctVSQn_Big[EPSymmBin] = new TCanvas(Form("cRelFluctVSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cRelFluctVSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);

  cCumuVSQn[EPSymmBin]    = new TCanvas(Form("cCumuVSQn_%s",    EPSymmNames[EPSymmBin].data()), Form("cCumuVSQn_%s",    EPSymmNames[EPSymmBin].data()), 1500, 1000);
  cVn2VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn2VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn2VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cVn4VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cVn6VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn6VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn6VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cVn8VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn8VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn8VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);

  cCumuRatioVSQn[EPSymmBin]  = new TCanvas(Form("cCumuRatioVSQn_%s", EPSymmNames[EPSymmBin].data()), Form("cCumuRatioVSQn_%s",   EPSymmNames[EPSymmBin].data()), 1500, 1000);
  cVn6Vn4VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn6Vn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn6Vn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cVn8Vn4VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn8Vn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn8Vn4VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cVn8Vn6VSQn_Big[EPSymmBin] = new TCanvas(Form("cVn8Vn6VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cVn8Vn6VSQn_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cG1ExpVSQn_Big[EPSymmBin]  = new TCanvas(Form("cG1ExpVSQn_Big_%s",  EPSymmNames[EPSymmBin].data()), Form("cG1ExpVSQn_Big_%s",  EPSymmNames[EPSymmBin].data()), 1500, 2000);

  cMomentVSQn[EPSymmBin]       -> Divide(2,2);
  cRMSVSQn_Big[EPSymmBin]      -> Divide(3,4);
  cMeanVSQn_Big[EPSymmBin]     -> Divide(3,4);
  cRelFluctVSQn_Big[EPSymmBin] -> Divide(3,4);

  cCumuVSQn[EPSymmBin]    -> Divide(3,2);
  cVn2VSQn_Big[EPSymmBin] -> Divide(3,4);
  cVn4VSQn_Big[EPSymmBin] -> Divide(3,4);
  cVn6VSQn_Big[EPSymmBin] -> Divide(3,4);
  cVn8VSQn_Big[EPSymmBin] -> Divide(3,4);

  cCumuRatioVSQn[EPSymmBin]  -> Divide(3,2);
  cVn6Vn4VSQn_Big[EPSymmBin] -> Divide(3,4);
  cVn8Vn4VSQn_Big[EPSymmBin] -> Divide(3,4);
  cVn8Vn6VSQn_Big[EPSymmBin] -> Divide(3,4);
  cG1ExpVSQn_Big[EPSymmBin]  -> Divide(3,4);


  //-- ========================== MomentVSQn plot ==========================
  for(int icent = 0; icent < NCENT; icent++){
  
    grMeanVnVSQn[icent][EPSymmBin]->GetXaxis()->SetTitleOffset(offsx);
    grRMSVnVSQn[icent][EPSymmBin]->GetXaxis()->SetTitleOffset(offsx);
    grRelFluctVnVSQn[icent][EPSymmBin]->GetXaxis()->SetTitleOffset(offsx);
    grStdDevVnVSQn[icent][EPSymmBin]->GetXaxis()->SetTitleOffset(offsx);

    grMeanVnVSQn[icent][EPSymmBin]->GetYaxis()->SetTitleOffset(offsy);
    grRMSVnVSQn[icent][EPSymmBin]->GetYaxis()->SetTitleOffset(offsy);
    grRelFluctVnVSQn[icent][EPSymmBin]->GetYaxis()->SetTitleOffset(offsy);
    grStdDevVnVSQn[icent][EPSymmBin]->GetYaxis()->SetTitleOffset(offsy);

    cMomentVSQn[EPSymmBin]->cd(1)->SetLeftMargin(0.2);
    if( icent == 0 ) grMeanVnVSQn[icent][EPSymmBin]->Draw("alp");
    else             grMeanVnVSQn[icent][EPSymmBin]->Draw("lpsame");

    cMomentVSQn[EPSymmBin]->cd(2)->SetLeftMargin(0.2);
    if( icent == 0 ) grStdDevVnVSQn[icent][EPSymmBin]->Draw("alp");
    else             grStdDevVnVSQn[icent][EPSymmBin]->Draw("lpsame");

    cMomentVSQn[EPSymmBin]->cd(3)->SetLeftMargin(0.2);
    if( icent == 0 ) grRMSVnVSQn[icent][EPSymmBin]->Draw("alp");
    else             grRMSVnVSQn[icent][EPSymmBin]->Draw("lpsame");

    cMomentVSQn[EPSymmBin]->cd(4)->SetLeftMargin(0.2);
    if( icent == 0 ) grRelFluctVnVSQn[icent][EPSymmBin]->Draw("alp");
    else             grRelFluctVnVSQn[icent][EPSymmBin]->Draw("lpsame");
  }

  cMomentVSQn[EPSymmBin]->cd(4);
  legCent->Draw();

  cMomentVSQn[EPSymmBin]->cd(1);
  latex.DrawLatex(0.25, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.25, 0.82, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.25, 0.76, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  cMomentVSQn[EPSymmBin]->SaveAs( Form("plots/cMomentVSQn_%s.pdf", EPSymmNames[EPSymmBin].data()) );

  //-- ========================== CumuVSQn plot ==========================
  for(int icent = 0; icent < NCENT; icent++){
    cCumuVSQn[EPSymmBin]->cd(1);
    if( icent == 0 ) grVn2VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn2VSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuVSQn[EPSymmBin]->cd(2);
    if( icent == 0 ) grVn4VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn4VSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuVSQn[EPSymmBin]->cd(3);
    if( icent == 0 ) grVn6VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn6VSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuVSQn[EPSymmBin]->cd(4);
    if( icent == 0 ) grVn8VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn8VSQn[icent][EPSymmBin]->Draw("lpsame");
  }

  cCumuVSQn[EPSymmBin]->cd(5);
  legCent->Draw();
  latex.DrawLatex(0.0775, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.0775, 0.8, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.0775, 0.72, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  cCumuVSQn[EPSymmBin]->SaveAs( Form("plots/physics/cCumuVSQn_%s.pdf", EPSymmNames[EPSymmBin].data()) );

  //-- ========================== CumuRatioVSQn plot ==========================
  for(int icent = 0; icent < NCENT; icent++){
    cCumuRatioVSQn[EPSymmBin]->cd(1);
    if( icent == 0 ) grG1ExpVSQn[icent][EPSymmBin]->Draw("alp");
    else             grG1ExpVSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuRatioVSQn[EPSymmBin]->cd(2);
    if( icent == 0 ) grVn6Vn4VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn6Vn4VSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuRatioVSQn[EPSymmBin]->cd(3);
    if( icent == 0 ) grVn8Vn4VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn8Vn4VSQn[icent][EPSymmBin]->Draw("lpsame");

    cCumuRatioVSQn[EPSymmBin]->cd(4);
    if( icent == 0 ) grVn8Vn6VSQn[icent][EPSymmBin]->Draw("alp");
    else             grVn8Vn6VSQn[icent][EPSymmBin]->Draw("lpsame");
  }

  cCumuRatioVSQn[EPSymmBin]->cd(5);
  legCent->Draw();
  latex.DrawLatex(0.0775, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.0775, 0.8, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.0775, 0.72, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  cCumuRatioVSQn[EPSymmBin]->SaveAs( Form("plots/physics/cCumuRatioVSQn_%s.pdf", EPSymmNames[EPSymmBin].data()) );

  //-- ========================== *_Big plots ==========================
  for(int icent = 0; icent < NCENT; icent++){

    //-- cRMSVSQn_Big
    cRMSVSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grRMSVnVSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grRMSVnVSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grRMSVnVSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cMeanVSQn_Big
    cMeanVSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grMeanVnVSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grMeanVnVSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grMeanVnVSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //--c RelFluctVSQn_Big
    cRelFluctVSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grRelFluctVnVSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grRelFluctVnVSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grRelFluctVnVSQn[icent][EPSymmBin]->Draw("alp");
    rfLine->Draw("same");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn2VSQn_Big
    cVn2VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn2VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn2VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn2VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn4VSQn_Big
    cVn4VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn4VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn4VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn4VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn6VSQn_Big
    cVn6VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn6VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn6VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn6VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn8VSQn_Big
    cVn8VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn8VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn8VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn8VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn6Vn4VSQn_Big
    cVn6Vn4VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn6Vn4VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn6Vn4VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn6Vn4VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn8Vn4VSQn_Big
    cVn8Vn4VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn8Vn4VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn8Vn4VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn8Vn4VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cVn8Vn6VSQn_Big
    cVn8Vn6VSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grVn8Vn4VSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grVn8Vn4VSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grVn8Vn4VSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cG1ExpVSQn_Big
    cG1ExpVSQn_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ){
      grG1ExpVSQn_DoSys[icent][EPSymmBin]->Draw("apE2");
      grG1ExpVSQn[icent][EPSymmBin]->Draw("lpsame");
    }
    else         grG1ExpVSQn[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

  }

  if( dosys_ ){
    cRMSVSQn_Big[EPSymmBin]      -> SaveAs( Form("plots/physics/cRMSVSQn_Big_%s_DOSYS.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanVSQn_Big[EPSymmBin]     -> SaveAs( Form("plots/physics/cMeanVSQn_Big_%s_DOSYS.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctVSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cRelFluctVSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    
    cVn2VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn2VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn4VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn6VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn6VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    
    cVn6Vn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn6Vn4VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8Vn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8Vn4VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8Vn6VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8Vn6VSQn_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
    cG1ExpVSQn_Big[EPSymmBin]  -> SaveAs( Form("plots/physics/cG1ExpVSQn_Big_%s_DOSYS.pdf",  EPSymmNames[EPSymmBin].data()) );
  }
  else{
    cRMSVSQn_Big[EPSymmBin]      -> SaveAs( Form("plots/physics/cRMSVSQn_Big_%s.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanVSQn_Big[EPSymmBin]     -> SaveAs( Form("plots/physics/cMeanVSQn_Big_%s.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctVSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cRelFluctVSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );

    cVn2VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn2VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn4VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn6VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn6VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );

    cVn6Vn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn6Vn4VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8Vn4VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8Vn4VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cVn8Vn6VSQn_Big[EPSymmBin] -> SaveAs( Form("plots/physics/cVn8Vn6VSQn_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
    cG1ExpVSQn_Big[EPSymmBin]  -> SaveAs( Form("plots/physics/cG1ExpVSQn_Big_%s.pdf",  EPSymmNames[EPSymmBin].data()) );
  }

}

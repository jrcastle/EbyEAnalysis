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


void crossHarmonicCorr(){

  int norder_     = 2;
  double tkEta    = 2.4;

  bool dosys_      = 0;
  double sysWidth_ = 0.1;

  bool looseChi2IterCut   = 0;
  bool nominalChi2IterCut = 1;
  bool tightChi2IterCut   = 0;

  double rmsV2Min      = 0.;
  double rmsV2Max      = 0.2;
  double meanV2Min     = 0.;
  double meanV2Max     = 0.2;
  double relFluctV2Min = 0.0;
  double relFluctV2Max = 0.7;

  double rmsV3Min      = 0.;
  double rmsV3Max      = 0.05;
  double meanV3Min     = 0.;
  double meanV3Max     = 0.05;
  double relFluctV3Min = 0.4;
  double relFluctV3Max = 0.7;

  double rmsV4Min      = 0.;
  double rmsV4Max      = 0.05;
  double meanV4Min     = 0.;
  double meanV4Max     = 0.05;
  double relFluctV4Min = 0.4;
  double relFluctV4Max = 0.7;

  TLatex latex;

  //-- Analyzer Output
  TFile * fAna[NVN];
  TH1D * hObs[NVN][NCENT][NEPSymm][NQN];

  //-- Unfolding output
  TFile * fUnfold[NVN];
  TH1D * hUnfold[NVN][NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NVN][NCENT][NEPSymm][NQN][NITER];

  //-- Statististical Errors
  TFile * fStatErr[NVN];
  TH1D * hVarianceOfMean_RMSVn[NVN][NEPSymm][NQN];
  TH1D * hVarianceOfMean_MeanVn[NVN][NEPSymm][NQN];
  TH1D * hVarianceOfMean_RelFluctVn[NVN][NEPSymm][NQN];

  //-- Systematic Errors
  TFile * fRelErr_Regularization[NVN];
  TH1D * relErr_Regularization_RMSVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Regularization_MeanVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Regularization_StDevVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Regularization_RelFluctVn[NVN][NEPSymm][NQN];

  TFile * fRelErr_RespEl[NVN];
  TH1D * relErr_RespEl_RMSVn[NVN][NEPSymm][NQN];
  TH1D * relErr_RespEl_MeanVn[NVN][NEPSymm][NQN];
  TH1D * relErr_RespEl_StDevVn[NVN][NEPSymm][NQN];
  TH1D * relErr_RespEl_RelFluctVn[NVN][NEPSymm][NQN];

  TFile * fRelErr_Vtx[NVN];
  TH1D * relErr_Vtx_RMSVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Vtx_MeanVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Vtx_StDevVn[NVN][NEPSymm][NQN];
  TH1D * relErr_Vtx_RelFluctVn[NVN][NEPSymm][NQN];

  TFile * fRelErr_TkEff[NVN];
  TH1D * relErr_TkEff_RMSVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkEff_MeanVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkEff_StDevVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkEff_RelFluctVn[NVN][NEPSymm][NQN];

  TFile * fRelErr_TkQuality[NVN];
  TH1D * relErr_TkQuality_RMSVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkQuality_MeanVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkQuality_StDevVn[NVN][NEPSymm][NQN];
  TH1D * relErr_TkQuality_RelFluctVn[NVN][NEPSymm][NQN];

  //-- Final Unfolded iteration
  int iterCut[NVN][NCENT][NEPSymm][NQN];

  //-- Moment vs vm 
  double rmsVn_vs_vm[NVN][NCENT][NEPSymm][NQN];
  double meanVn_vs_vm[NVN][NCENT][NEPSymm][NQN];
  double relFluctVn_vs_vm[NVN][NCENT][NEPSymm][NQN];

  double rmsVn_vs_vm_statErr[NVN][NCENT][NEPSymm][NQN];
  double meanVn_vs_vm_statErr[NVN][NCENT][NEPSymm][NQN];
  double relFluctVn_vs_vm_statErr[NVN][NCENT][NEPSymm][NQN];

  double rmsVn_vs_vm_sysErr[NVN][NCENT][NEPSymm][NQN];
  double meanVn_vs_vm_sysErr[NVN][NCENT][NEPSymm][NQN];
  double relFluctVn_vs_vm_sysErr[NVN][NCENT][NEPSymm][NQN];

  TGraphErrors * grRMSV3VSV2[NCENT][NEPSymm];
  TGraphErrors * grMeanV3VSV2[NCENT][NEPSymm];
  TGraphErrors * grRelFluctV3VSV2[NCENT][NEPSymm];

  TGraphErrors * grRMSV3VSV2_DoSys[NCENT][NEPSymm];
  TGraphErrors * grMeanV3VSV2_DoSys[NCENT][NEPSymm];
  TGraphErrors * grRelFluctV3VSV2_DoSys[NCENT][NEPSymm];

  TGraphErrors * grRMSV4VSV2[NCENT][NEPSymm];
  TGraphErrors * grMeanV4VSV2[NCENT][NEPSymm];
  TGraphErrors * grRelFluctV4VSV2[NCENT][NEPSymm];

  TGraphErrors * grRMSV4VSV2_DoSys[NCENT][NEPSymm];
  TGraphErrors * grMeanV4VSV2_DoSys[NCENT][NEPSymm];
  TGraphErrors * grRelFluctV4VSV2_DoSys[NCENT][NEPSymm];


  //-- Canvases
  TCanvas * cMoment_V3vsV2[NEPSymm];
  TCanvas * cMoment_V4vsV2[NEPSymm];

  TCanvas * cRMSV3vsV2_Big[NEPSymm];
  TCanvas * cMeanV3vsV2_Big[NEPSymm];
  TCanvas * cRelFluctV3vsV2_Big[NEPSymm];

  TCanvas * cRMSV4vsV2_Big[NEPSymm];
  TCanvas * cMeanV4vsV2_Big[NEPSymm];
  TCanvas * cRelFluctV4vsV2_Big[NEPSymm];
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

  double chi2Cut;
  if( looseChi2IterCut )   chi2Cut = 1.5;
  if( nominalChi2IterCut ) chi2Cut = 1.2;
  if( tightChi2IterCut )   chi2Cut = 1.0;

  //double sysWidth[NQN];
  //-- Widths for systematic error bars
  //for(intiqn = 0; iqn < NQN; iqn++) sysWidth[iqn] = sysWidth_;

  for(int ivn = 0; ivn < NVN; ivn++){

    //-- Get the Analyzer output file
    fAna[ivn] = new TFile( Form("../../v%i/eta2.4/AnalyzerResults/CastleEbyE.root", vn_[ivn]) );

    //-- Get Unfolding output
    fUnfold[ivn] = new TFile( Form("../../v%i/eta2.4/UnfoldResults/dataResp/data%i.root", vn_[ivn], vn_[ivn]) );

    //-- Get Statististical Errors
    fStatErr[ivn] = new TFile( Form("../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", vn_[ivn], vn_[ivn]) );
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){
	hVarianceOfMean_RMSVn[ivn][iEP][iqn]      = (TH1D*) fStatErr[ivn]->Get( Form("hVarianceOfMean_RMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	hVarianceOfMean_MeanVn[ivn][iEP][iqn]     = (TH1D*) fStatErr[ivn]->Get( Form("hVarianceOfMean_MeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	hVarianceOfMean_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fStatErr[ivn]->Get( Form("hVarianceOfMean_RelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
      }
    }

    //-- Get Systematic Errors
    if( dosys_ ){
      //-- Systematic Errors
      fRelErr_Regularization[ivn] = new TFile( Form("../../v%i/eta2.4/systematicStudies/chi2Cutoff/relErrorRegularization.root", vn_[ivn]) );
      fRelErr_RespEl[ivn]         = new TFile( Form("../../v%i/eta2.4/systematicStudies/responseElements/relErrorResponse.root", vn_[ivn]) );
      fRelErr_Vtx[ivn]            = new TFile( Form("../../v%i/eta2.4/systematicStudies/vtxCut/relErrorVtx.root", vn_[ivn]) );
      fRelErr_TkEff[ivn]          = new TFile( Form("../../v%i/eta2.4/systematicStudies/tkEff/relErrorTkEff.root", vn_[ivn]) );
      fRelErr_TkQuality[ivn]      = new TFile( Form("../../v%i/eta2.4/systematicStudies/tkQuality/relErrorTkQuality.root", vn_[ivn]) );

      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;
	for(int iqn = 0; iqn < NQN; iqn++){

	  relErr_Regularization_RMSVn[ivn][iEP][iqn]      = (TH1D*) fRelErr_Regularization[ivn]->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	  relErr_Regularization_MeanVn[ivn][iEP][iqn]     = (TH1D*) fRelErr_Regularization[ivn]->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	  relErr_Regularization_StDevVn[ivn][iEP][iqn]    = (TH1D*) fRelErr_Regularization[ivn]->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	  relErr_Regularization_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fRelErr_Regularization[ivn]->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );

	  relErr_RespEl_RMSVn[ivn][iEP][iqn]      = (TH1D*) fRelErr_RespEl[ivn]->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	  relErr_RespEl_MeanVn[ivn][iEP][iqn]     = (TH1D*) fRelErr_RespEl[ivn]->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	  relErr_RespEl_StDevVn[ivn][iEP][iqn]    = (TH1D*) fRelErr_RespEl[ivn]->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	  relErr_RespEl_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fRelErr_RespEl[ivn]->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );

	  relErr_Vtx_RMSVn[ivn][iEP][iqn]      = (TH1D*) fRelErr_Vtx[ivn]->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	  relErr_Vtx_MeanVn[ivn][iEP][iqn]     = (TH1D*) fRelErr_Vtx[ivn]->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	  relErr_Vtx_StDevVn[ivn][iEP][iqn]    = (TH1D*) fRelErr_Vtx[ivn]->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	  relErr_Vtx_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fRelErr_Vtx[ivn]->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );

	  relErr_TkEff_RMSVn[ivn][iEP][iqn]      = (TH1D*) fRelErr_TkEff[ivn]->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	  relErr_TkEff_MeanVn[ivn][iEP][iqn]     = (TH1D*) fRelErr_TkEff[ivn]->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	  relErr_TkEff_StDevVn[ivn][iEP][iqn]    = (TH1D*) fRelErr_TkEff[ivn]->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	  relErr_TkEff_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fRelErr_TkEff[ivn]->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );

	  relErr_TkQuality_RMSVn[ivn][iEP][iqn]      = (TH1D*) fRelErr_TkQuality[ivn]->Get( Form("relErrRMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	  relErr_TkQuality_MeanVn[ivn][iEP][iqn]     = (TH1D*) fRelErr_TkQuality[ivn]->Get( Form("relErrMeanVn_%s_qbin%i",     EPSymmNames[iEP].data(), iqn) );
	  relErr_TkQuality_StDevVn[ivn][iEP][iqn]    = (TH1D*) fRelErr_TkQuality[ivn]->Get( Form("relErrStDevVn_%s_qbin%i",    EPSymmNames[iEP].data(), iqn) );
	  relErr_TkQuality_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fRelErr_TkQuality[ivn]->Get( Form("relErrRelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );

	}
      }

    } // End if( dosys_ )

    //-- Start looping over the data...
    for(int icent = 0; icent < NCENT; icent++){

      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;

	for(int iqn = 0; iqn < NQN; iqn++){

	  //-- Get the VN observed histogram
	  hObs[ivn][icent][iEP][iqn] = (TH1D*) fAna[ivn]->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	  iterCut[ivn][icent][iEP][iqn] = 0;
	  for(int i = 0; i < NITER; i++){

	    //-- Get the unfolded histograms
	    hUnfold[ivn][icent][iEP][iqn][i] = (TH1D*) fUnfold[ivn]->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	    //-- Get the Refolded histograms
	    hRefold[ivn][icent][iEP][iqn][i] = (TH1D*) fUnfold[ivn]->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	    //-- Chi2 Cut
	    double chi2NDF_Refold = hRefold[ivn][icent][iEP][iqn][i]->Chi2Test(hObs[ivn][icent][iEP][iqn], "CHI2/NDF");

	    if( chi2NDF_Refold <= chi2Cut ){
	      iterCut[ivn][icent][iEP][iqn] = i;
	      break;
	    }
	    if( i == NITER-1 ){
	      iterCut[ivn][icent][iEP][iqn] = i;
	      break;
	    }

	  } //-- End unfold iteration loop


	  //-- Fill Arrays for plots
	  int it = iterCut[ivn][icent][iEP][iqn];
	  FixUnfold( hUnfold[ivn][icent][iEP][iqn][it] );
	  EbyECumu cumu( hUnfold[ivn][icent][iEP][iqn][it] );

	  //-- Values
	  double mean     = hUnfold[ivn][icent][iEP][iqn][it]->GetMean();
	  double stdev    = hUnfold[ivn][icent][iEP][iqn][it]->GetRMS();
	  double rms      = sqrt( mean*mean + stdev*stdev );
	  double relFluct = stdev / mean;

	  //-- Stat Errors
	  double rmsStatErr      = sqrt( hVarianceOfMean_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1) );
	  double meanStatErr     = sqrt( hVarianceOfMean_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1) );
	  double relFluctStatErr = sqrt( hVarianceOfMean_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1) );

	  //-- Sys Errors
	  if( dosys_ ){

	    //-- Individual studies
	    double re_Reg_RMSVn      = relErr_Regularization_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Reg_MeanVn     = relErr_Regularization_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Reg_StDevVn    = relErr_Regularization_StDevVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Reg_RelFluctVn = relErr_Regularization_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1);

	    double re_Resp_RMSVn      = relErr_RespEl_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Resp_MeanVn     = relErr_RespEl_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Resp_StDevVn    = relErr_RespEl_StDevVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Resp_RelFluctVn = relErr_RespEl_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1);

	    double re_Vtx_RMSVn      = relErr_Vtx_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Vtx_MeanVn     = relErr_Vtx_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Vtx_StDevVn    = relErr_Vtx_StDevVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_Vtx_RelFluctVn = relErr_Vtx_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1);

	    double re_TkEff_RMSVn      = relErr_TkEff_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkEff_MeanVn     = relErr_TkEff_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkEff_StDevVn    = relErr_TkEff_StDevVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkEff_RelFluctVn = relErr_TkEff_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1);

	    double re_TkQ_RMSVn      = relErr_TkQuality_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkQ_MeanVn     = relErr_TkQuality_MeanVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkQ_StDevVn    = relErr_TkQuality_StDevVn[ivn][iEP][iqn]->GetBinContent(icent+1);
	    double re_TkQ_RelFluctVn = relErr_TkQuality_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1);

	    //-- Merge
	    double sys_RMSVn      = rms * sqrt( pow(re_Reg_RMSVn,2) + pow(re_Resp_RMSVn,2) + pow(re_Vtx_RMSVn,2) + pow(re_TkEff_RMSVn,2) + pow(re_TkQ_RMSVn,2) );
	    double sys_MeanVn     = mean * sqrt( pow(re_Reg_MeanVn,2) + pow(re_Resp_MeanVn,2) + pow(re_Vtx_MeanVn,2) + pow(re_TkEff_MeanVn,2) + pow(re_TkQ_MeanVn,2) );
	    double sys_StDevVn    = stdev * sqrt( pow(re_Reg_StDevVn,2) + pow(re_Resp_StDevVn,2) + pow(re_Vtx_StDevVn,2) + pow(re_TkEff_StDevVn,2) + pow(re_TkQ_StDevVn,2) );
	    double sys_RelFluctVn = relFluct * sqrt( pow(re_Reg_RelFluctVn,2) + pow(re_Resp_RelFluctVn,2) + pow(re_Vtx_RelFluctVn,2) + pow(re_TkEff_RelFluctVn,2) + pow(re_TkQ_RelFluctVn,2) );

	    //-- Set
	    rmsVn_vs_vm_sysErr[ivn][icent][iEP][iqn]      = sys_RMSVn;
	    meanVn_vs_vm_sysErr[ivn][icent][iEP][iqn]     = sys_MeanVn;
	    relFluctVn_vs_vm_sysErr[ivn][icent][iEP][iqn] = sys_RelFluctVn;

	  } //End if( dosys_ )

	  rmsVn_vs_vm[ivn][icent][iEP][iqn]      = rms;
	  meanVn_vs_vm[ivn][icent][iEP][iqn]     = mean;
	  relFluctVn_vs_vm[ivn][icent][iEP][iqn] = relFluct;

	  rmsVn_vs_vm_statErr[ivn][icent][iEP][iqn]      = rmsStatErr;
	  meanVn_vs_vm_statErr[ivn][icent][iEP][iqn]     = meanStatErr;
	  relFluctVn_vs_vm_statErr[ivn][icent][iEP][iqn] = relFluctStatErr;

	} //-- End QN loop

	//-- Initialize TGraphErrors
	if( ivn == NVN-1 ){

	  grRMSV3VSV2[icent][iEP]      = new TGraphErrors(NQN, rmsVn_vs_vm[0][icent][iEP],      rmsVn_vs_vm[1][icent][iEP],      rmsVn_vs_vm_statErr[0][icent][iEP],      rmsVn_vs_vm_statErr[1][icent][iEP]);
	  grMeanV3VSV2[icent][iEP]     = new TGraphErrors(NQN, meanVn_vs_vm[0][icent][iEP],     meanVn_vs_vm[1][icent][iEP],     meanVn_vs_vm_statErr[0][icent][iEP],     meanVn_vs_vm_statErr[1][icent][iEP]);
	  grRelFluctV3VSV2[icent][iEP] = new TGraphErrors(NQN, relFluctVn_vs_vm[0][icent][iEP], relFluctVn_vs_vm[1][icent][iEP], relFluctVn_vs_vm_statErr[0][icent][iEP], relFluctVn_vs_vm_statErr[1][icent][iEP]);

	  grRMSV4VSV2[icent][iEP]      = new TGraphErrors(NQN, rmsVn_vs_vm[0][icent][iEP],      rmsVn_vs_vm[2][icent][iEP],      rmsVn_vs_vm_statErr[0][icent][iEP],      rmsVn_vs_vm_statErr[2][icent][iEP]);
          grMeanV4VSV2[icent][iEP]     = new TGraphErrors(NQN, meanVn_vs_vm[0][icent][iEP],     meanVn_vs_vm[2][icent][iEP],     meanVn_vs_vm_statErr[0][icent][iEP],     meanVn_vs_vm_statErr[2][icent][iEP]);
          grRelFluctV4VSV2[icent][iEP] = new TGraphErrors(NQN, relFluctVn_vs_vm[0][icent][iEP], relFluctVn_vs_vm[2][icent][iEP], relFluctVn_vs_vm_statErr[0][icent][iEP], relFluctVn_vs_vm_statErr[2][icent][iEP]);

	  //-- Format Graphs

	  //-- V3 vs V2
	  grRMSV3VSV2[icent][iEP]->GetXaxis()->SetTitle( "#sqrt{#LTv_{2}^{2}#GT}" );
	  grRMSV3VSV2[icent][iEP]->GetXaxis()->SetLimits(rmsV2Min, rmsV2Max);
	  grRMSV3VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
	  grRMSV3VSV2[icent][iEP]->GetYaxis()->SetTitle( "#sqrt{#LTv_{3}^{2}#GT}" );
	  grRMSV3VSV2[icent][iEP]->GetYaxis()->SetRangeUser(rmsV3Min, rmsV3Max);
	  grRMSV3VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
	  grRMSV3VSV2[icent][iEP]->SetLineColor( centCol[icent] );
	  grRMSV3VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
	  grRMSV3VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );

	  grMeanV3VSV2[icent][iEP]->GetXaxis()->SetTitle( "#LTv_{2}#GT" );
          grMeanV3VSV2[icent][iEP]->GetXaxis()->SetLimits(meanV2Min, meanV2Max);
          grMeanV3VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
          grMeanV3VSV2[icent][iEP]->GetYaxis()->SetTitle( "#LTv_{3}#GT" );
          grMeanV3VSV2[icent][iEP]->GetYaxis()->SetRangeUser(meanV3Min, meanV3Max);
	  grMeanV3VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grMeanV3VSV2[icent][iEP]->SetLineColor( centCol[icent] );
          grMeanV3VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
          grMeanV3VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );

	  grRelFluctV3VSV2[icent][iEP]->GetXaxis()->SetTitle( "#sigma_{v_{2}} / #LTv_{2}#GT" );
          grRelFluctV3VSV2[icent][iEP]->GetXaxis()->SetLimits(relFluctV2Min, relFluctV2Max);
          grRelFluctV3VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
          grRelFluctV3VSV2[icent][iEP]->GetYaxis()->SetTitle( "#sigma_{v_{3}} / #LTv_{3}#GT" );
	  grRelFluctV3VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grRelFluctV3VSV2[icent][iEP]->GetYaxis()->SetRangeUser(relFluctV3Min, relFluctV3Max);
          grRelFluctV3VSV2[icent][iEP]->SetLineColor( centCol[icent] );
          grRelFluctV3VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
          grRelFluctV3VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );

	  //-- V4 vs V2
	  grRMSV4VSV2[icent][iEP]->GetXaxis()->SetTitle( "#sqrt{#LTv_{2}^{2}#GT}" );
          grRMSV4VSV2[icent][iEP]->GetXaxis()->SetLimits(rmsV2Min, rmsV2Max);
          grRMSV4VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
          grRMSV4VSV2[icent][iEP]->GetYaxis()->SetTitle( "#sqrt{#LTv_{4}^{2}#GT}" );
	  grRMSV4VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grRMSV4VSV2[icent][iEP]->GetYaxis()->SetRangeUser(rmsV4Min, rmsV4Max);
          grRMSV4VSV2[icent][iEP]->SetLineColor( centCol[icent] );
          grRMSV4VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
          grRMSV4VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );

          grMeanV4VSV2[icent][iEP]->GetXaxis()->SetTitle( "#LTv_{2}#GT" );
          grMeanV4VSV2[icent][iEP]->GetXaxis()->SetLimits(meanV2Min, meanV2Max);
          grMeanV4VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
          grMeanV4VSV2[icent][iEP]->GetYaxis()->SetTitle( "#LTv_{4}#GT" );
          grMeanV4VSV2[icent][iEP]->GetYaxis()->SetRangeUser(meanV4Min, meanV4Max);
	  grMeanV4VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grMeanV4VSV2[icent][iEP]->SetLineColor( centCol[icent] );
          grMeanV4VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
          grMeanV4VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );

          grRelFluctV4VSV2[icent][iEP]->GetXaxis()->SetTitle( "#sigma_{v_{2}} / #LTv_{2}#GT" );
          grRelFluctV4VSV2[icent][iEP]->GetXaxis()->SetLimits(relFluctV2Min, relFluctV2Max);
          grRelFluctV4VSV2[icent][iEP]->GetXaxis()->SetNdivisions(509);
          grRelFluctV4VSV2[icent][iEP]->GetYaxis()->SetTitle( "#sigma_{v_{4}} / #LTv_{4}#GT" );
          grRelFluctV4VSV2[icent][iEP]->GetYaxis()->SetRangeUser(relFluctV4Min, relFluctV4Max);
	  grRelFluctV4VSV2[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grRelFluctV4VSV2[icent][iEP]->SetLineColor( centCol[icent] );
          grRelFluctV4VSV2[icent][iEP]->SetMarkerColor( centCol[icent] );
          grRelFluctV4VSV2[icent][iEP]->SetMarkerStyle( centMark[icent] );


	  if( dosys_ ){

	    grRMSV3VSV2_DoSys[icent][iEP]      = new TGraphErrors(NQN, rmsVn_vs_vm[0][icent][iEP],      rmsVn_vs_vm[1][icent][iEP],      rmsVn_vs_vm_sysErr[0][icent][iEP],      rmsVn_vs_vm_sysErr[1][icent][iEP]);
	    grMeanV3VSV2_DoSys[icent][iEP]     = new TGraphErrors(NQN, meanVn_vs_vm[0][icent][iEP],     meanVn_vs_vm[1][icent][iEP],     meanVn_vs_vm_sysErr[0][icent][iEP],     meanVn_vs_vm_sysErr[1][icent][iEP]);
	    grRelFluctV3VSV2_DoSys[icent][iEP] = new TGraphErrors(NQN, relFluctVn_vs_vm[0][icent][iEP], relFluctVn_vs_vm[1][icent][iEP], relFluctVn_vs_vm_sysErr[0][icent][iEP], relFluctVn_vs_vm_sysErr[1][icent][iEP]);

	    grRMSV4VSV2_DoSys[icent][iEP]      = new TGraphErrors(NQN, rmsVn_vs_vm[0][icent][iEP],      rmsVn_vs_vm[2][icent][iEP],      rmsVn_vs_vm_sysErr[0][icent][iEP],      rmsVn_vs_vm_sysErr[2][icent][iEP]);
	    grMeanV4VSV2_DoSys[icent][iEP]     = new TGraphErrors(NQN, meanVn_vs_vm[0][icent][iEP],     meanVn_vs_vm[2][icent][iEP],     meanVn_vs_vm_sysErr[0][icent][iEP],     meanVn_vs_vm_sysErr[2][icent][iEP]);
	    grRelFluctV4VSV2_DoSys[icent][iEP] = new TGraphErrors(NQN, relFluctVn_vs_vm[0][icent][iEP], relFluctVn_vs_vm[2][icent][iEP], relFluctVn_vs_vm_sysErr[0][icent][iEP], relFluctVn_vs_vm_sysErr[2][icent][iEP]);

	    //-- Format Graphs

	    //-- V3 vs V2
	    grRMSV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#sqrt{#LTv_{2}^{2}#GT}" );
	    grRMSV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(rmsV2Min, rmsV2Max);
	    grRMSV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grRMSV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#sqrt{#LTv_{3}^{2}#GT}" );
	    grRMSV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(rmsV3Min, rmsV3Max);
	    grRMSV3VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grRMSV3VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grRMSV3VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grRMSV3VSV2_DoSys[icent][iEP]->SetFillColor(17);

	    grMeanV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#LTv_{2}#GT" );
	    grMeanV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(meanV2Min, meanV2Max);
	    grMeanV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grMeanV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#LTv_{3}#GT" );
	    grMeanV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(meanV3Min, meanV3Max);
	    grMeanV3VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grMeanV3VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grMeanV3VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grMeanV3VSV2_DoSys[icent][iEP]->SetFillColor(17);

	    grRelFluctV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#sigma_{v_{2}} / #LTv_{2}#GT" );
	    grRelFluctV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(relFluctV2Min, relFluctV2Max);
	    grRelFluctV3VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grRelFluctV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#sigma_{v_{3}} / #LTv_{3}#GT" );
	    grRelFluctV3VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(relFluctV3Min, relFluctV3Max);
	    grRelFluctV3VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grRelFluctV3VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grRelFluctV3VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grRelFluctV3VSV2_DoSys[icent][iEP]->SetFillColor(17);

	    //-- V4 vs V2
	    grRMSV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#sqrt{#LTv_{2}^{2}#GT}" );
	    grRMSV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(rmsV2Min, rmsV2Max);
	    grRMSV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grRMSV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#sqrt{#LTv_{4}^{2}#GT}" );
	    grRMSV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(rmsV4Min, rmsV4Max);
	    grRMSV4VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grRMSV4VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grRMSV4VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grRMSV4VSV2_DoSys[icent][iEP]->SetFillColor(17);

	    grMeanV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#LTv_{2}#GT" );
	    grMeanV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(meanV2Min, meanV2Max);
	    grMeanV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grMeanV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#LTv_{4}#GT" );
	    grMeanV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(meanV4Min, meanV4Max);
	    grMeanV4VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grMeanV4VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grMeanV4VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grMeanV4VSV2_DoSys[icent][iEP]->SetFillColor(17);

	    grRelFluctV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetTitle( "#sigma_{v_{2}} / #LTv_{2}#GT" );
	    grRelFluctV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetLimits(relFluctV2Min, relFluctV2Max);
	    grRelFluctV4VSV2_DoSys[icent][iEP]->GetXaxis()->SetNdivisions(509);
	    grRelFluctV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetTitle( "#sigma_{v_{4}} / #LTv_{4}#GT" );
	    grRelFluctV4VSV2_DoSys[icent][iEP]->GetYaxis()->SetRangeUser(relFluctV4Min, relFluctV4Max);
	    grRelFluctV4VSV2_DoSys[icent][iEP]->SetLineColor( centCol[icent] );
	    grRelFluctV4VSV2_DoSys[icent][iEP]->SetMarkerColor( centCol[icent] );
	    grRelFluctV4VSV2_DoSys[icent][iEP]->SetMarkerStyle( centMark[icent] );
	    grRelFluctV4VSV2_DoSys[icent][iEP]->SetFillColor(17);

	  } //-- End if( dosys_ )

	} //-- End if( ivn == NVN-1 )

      } //-- End EP loop
    } //-- End cent loop
  } //-- End vn loop

  //-- DRAW!!!!!
  TLegend * legCent = new TLegend(0.0775, 0.1812, 0.9936, 0.6799);
  legCent->SetFillStyle(0);
  legCent->SetBorderSize(0);
  legCent->SetNColumns(2);
  for(int icent = 0; icent < NCENT; icent++) legCent->AddEntry(grRMSV3VSV2[icent][EPSymmBin], Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"), "lp");

  double BGrf = sqrt( 4./TMath::Pi() - 1);
  TLine * rfLine = new TLine(grRelFluctV3VSV2[0][EPSymmBin]->GetXaxis()->GetXmin(), BGrf, grRelFluctV3VSV2[0][EPSymmBin]->GetXaxis()->GetXmax(), BGrf);
  rfLine->SetLineColor(1);
  rfLine->SetLineStyle(2);
  rfLine->SetLineWidth(2);;

  //-- Initialize Canvases


  cMoment_V3vsV2[EPSymmBin] = new TCanvas(Form("cMoment V3vsV2_%s", EPSymmNames[EPSymmBin].data()), Form("cMoment V3vsV2_%s", EPSymmNames[EPSymmBin].data()), 1000, 1000);
  cMoment_V4vsV2[EPSymmBin] = new TCanvas(Form("cMoment V4vsV2_%s", EPSymmNames[EPSymmBin].data()), Form("cMoment V4vsV2_%s", EPSymmNames[EPSymmBin].data()), 1000, 1000);

  cRMSV3vsV2_Big[EPSymmBin]      = new TCanvas(Form("cRMSV3vsV2_Big_%s",      EPSymmNames[EPSymmBin].data()), Form("cRMSV3vsV2_Big_%s",      EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cMeanV3vsV2_Big[EPSymmBin]     = new TCanvas(Form("cMeanV3vsV2_Big_%s",     EPSymmNames[EPSymmBin].data()), Form("cMeanV3vsV2_Big_%s",     EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cRelFluctV3vsV2_Big[EPSymmBin] = new TCanvas(Form("cRelFluctV3vsV2_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cRelFluctV3vsV2_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);

  cRMSV4vsV2_Big[EPSymmBin]      = new TCanvas(Form("cRMSV4vsV2_Big_%s",      EPSymmNames[EPSymmBin].data()), Form("cRMSV4vsV2_Big_%s",      EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cMeanV4vsV2_Big[EPSymmBin]     = new TCanvas(Form("cMeanV4vsV2_Big_%s",     EPSymmNames[EPSymmBin].data()), Form("cMeanV4vsV2_Big_%s",     EPSymmNames[EPSymmBin].data()), 1500, 2000);
  cRelFluctV4vsV2_Big[EPSymmBin] = new TCanvas(Form("cRelFluctV4vsV2_Big_%s", EPSymmNames[EPSymmBin].data()), Form("cRelFluctV4vsV2_Big_%s", EPSymmNames[EPSymmBin].data()), 1500, 2000);

  cMoment_V3vsV2[EPSymmBin] -> Divide(2,2);
  cMoment_V4vsV2[EPSymmBin] -> Divide(2,2); 

  cRMSV3vsV2_Big[EPSymmBin]      -> Divide(3,4);
  cMeanV3vsV2_Big[EPSymmBin]     -> Divide(3,4);
  cRelFluctV3vsV2_Big[EPSymmBin] -> Divide(3,4);

  cRMSV4vsV2_Big[EPSymmBin]      -> Divide(3,4);
  cMeanV4vsV2_Big[EPSymmBin]     -> Divide(3,4);
  cRelFluctV4vsV2_Big[EPSymmBin] -> Divide(3,4);

  //-- ========================== Moment_V3vsV2 plot ==========================
  for(int icent = 0; icent < NCENT; icent++){
    cMoment_V3vsV2[EPSymmBin]->cd(1);
    if( icent == 0 ) grMeanV3VSV2[icent][EPSymmBin]->Draw("alp");
    else             grMeanV3VSV2[icent][EPSymmBin]->Draw("lpsame");

    cMoment_V3vsV2[EPSymmBin]->cd(2);
    if( icent == 0 ) grRMSV3VSV2[icent][EPSymmBin]->Draw("alp");
    else             grRMSV3VSV2[icent][EPSymmBin]->Draw("lpsame");

    cMoment_V3vsV2[EPSymmBin]->cd(3);
    if( icent == 0 ) grRelFluctV3VSV2[icent][EPSymmBin]->Draw("alp");
    else             grRelFluctV3VSV2[icent][EPSymmBin]->Draw("lpsame");
  }

  cMoment_V3vsV2[EPSymmBin]->cd(3);
  rfLine->Draw("same");

  cMoment_V3vsV2[EPSymmBin]->cd(4);
  legCent->Draw();
  latex.DrawLatex(0.0775, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.0775, 0.8, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.0775, 0.72, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  cMoment_V3vsV2[EPSymmBin]->SaveAs( Form("plots/cMoment_V3vsV2_%s.pdf", EPSymmNames[EPSymmBin].data()) );

  //-- ========================== Moment_V4vsV2 plot ==========================
  for(int icent = 0; icent < NCENT; icent++){
    cMoment_V4vsV2[EPSymmBin]->cd(1);
    if( icent == 0 ) grMeanV4VSV2[icent][EPSymmBin]->Draw("alp");
    else             grMeanV4VSV2[icent][EPSymmBin]->Draw("lpsame");

    cMoment_V4vsV2[EPSymmBin]->cd(2);
    if( icent == 0 ) grRMSV4VSV2[icent][EPSymmBin]->Draw("alp");
    else             grRMSV4VSV2[icent][EPSymmBin]->Draw("lpsame");

    cMoment_V4vsV2[EPSymmBin]->cd(3);
    if( icent == 0 ) grRelFluctV4VSV2[icent][EPSymmBin]->Draw("alp");
    else             grRelFluctV4VSV2[icent][EPSymmBin]->Draw("lpsame");
  }

  cMoment_V4vsV2[EPSymmBin]->cd(3);
  rfLine->Draw("same");

  cMoment_V4vsV2[EPSymmBin]->cd(4);
  legCent->Draw();
  latex.DrawLatex(0.0775, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.0775, 0.8, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.0775, 0.72, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  cMoment_V4vsV2[EPSymmBin]->SaveAs( Form("plots/cMoment_V4vsV2_%s.pdf", EPSymmNames[EPSymmBin].data()) );

  //-- ========================== *_Big plots ==========================
  for(int icent = 0; icent < NCENT; icent++){

    //-- cRMSV3vsV2_Big
    cRMSV3vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grRMSV3VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grRMSV3VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cMeanV3vsV2_Big 
    cMeanV3vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grMeanV3VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grMeanV3VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cRelFluctV3vsV2_Big 
    cRelFluctV3vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grRelFluctV3VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grRelFluctV3VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cRMSV4vsV2_Big
    cRMSV4vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grRMSV4VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grRMSV4VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cMeanV4vsV2_Big
    cMeanV4vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grMeanV4VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grMeanV4VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

    //-- cRelFluctV4vsV2_Big
    cRelFluctV4vsV2_Big[EPSymmBin]->cd(icent+1);
    if( dosys_ ) grRelFluctV4VSV2_DoSys[icent][EPSymmBin]->Draw("alp");
    else         grRelFluctV4VSV2[icent][EPSymmBin]->Draw("alp");
    latex.DrawLatex(0.54, 0.88, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%") );

  }

  if( dosys_ ){
    cRMSV3vsV2_Big[EPSymmBin]      -> SaveAs( Form("plots/cRMSV3vsV2_Big_%s_DOSYS.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanV3vsV2_Big[EPSymmBin]     -> SaveAs( Form("plots/cMeanV3vsV2_Big_%s_DOSYS.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctV3vsV2_Big[EPSymmBin] -> SaveAs( Form("plots/cRelFluctV3vsV2_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );

    cRMSV4vsV2_Big[EPSymmBin]      -> SaveAs( Form("plots/cRMSV4vsV2_Big_%s_DOSYS.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanV4vsV2_Big[EPSymmBin]     -> SaveAs( Form("plots/cMeanV4vsV2_Big_%s_DOSYS.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctV4vsV2_Big[EPSymmBin] -> SaveAs( Form("plots/cRelFluctV4vsV2_Big_%s_DOSYS.pdf", EPSymmNames[EPSymmBin].data()) );
  }
  else{
    cRMSV3vsV2_Big[EPSymmBin]      -> SaveAs( Form("plots/cRMSV3vsV2_Big_%s.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanV3vsV2_Big[EPSymmBin]     -> SaveAs( Form("plots/cMeanV3vsV2_Big_%s.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctV3vsV2_Big[EPSymmBin] -> SaveAs( Form("plots/cRelFluctV3vsV2_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );

    cRMSV4vsV2_Big[EPSymmBin]      -> SaveAs( Form("plots/cRMSV4vsV2_Big_%s.pdf",      EPSymmNames[EPSymmBin].data()) );
    cMeanV4vsV2_Big[EPSymmBin]     -> SaveAs( Form("plots/cMeanV4vsV2_Big_%s.pdf",     EPSymmNames[EPSymmBin].data()) );
    cRelFluctV4vsV2_Big[EPSymmBin] -> SaveAs( Form("plots/cRelFluctV4vsV2_Big_%s.pdf", EPSymmNames[EPSymmBin].data()) );
  }

}

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
#include "TError.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;
using namespace hi;

void sysChi2Cutoff(){

  int norder_     = 3;
  double tkEta    = 2.4;

  bool dosys     = 0;
  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  TLatex latex;

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT][NEPSymm][NQN];

  //-- Unfolding output
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NITER];

  //-- Chi2 Refold vs iter
  TGraphErrors * grChi2NDFvsIter[NCENT][NEPSymm][NQN];
  double chi2NDF[NCENT][NEPSymm][NQN][NITER];

  //-- Cumu vs chi2 Cutoff
  double unfold_RMSVnVSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_MeanVnVSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_StDevVnVSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_RelFluctVnVSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn2VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn4VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn6VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn8VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Gamma1ExpVSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn6Vn4VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn8Vn4VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];
  double unfold_Vn8Vn6VSChi2Cut[NCENT][NEPSymm][NQN][NCHI2];

  //-- Ratio to the Chi2=1 cutoff scenario
  double unfoldRMSVn_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldMeanVn_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldStDevVn_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldRelFluctVn_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn2_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn4_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn6_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn8_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldGamma1Exp_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn6Vn4_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn8Vn4_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];
  double unfoldVn8Vn6_RatioToChi21[NCENT][NEPSymm][NQN][NCHI2];

  bool iterCutoff[NCHI2];

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

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  gErrorIgnoreLevel = kWarning;

  //-- Get the Analyzer output file
  fAna = new TFile( "../../AnalyzerResults/CastleEbyE.root" );

  //-- Get the unfolding output file
  bool R1 = dataResp && studTResp && gaussResp;
  bool R2 = dataResp && studTResp;
  bool R3 = dataResp && gaussResp;
  bool R4 = studTResp && gaussResp;

  if( R1 || R2 || R3 || R4){
    std::cout<<"WARNING! More than one response function defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  fUnfold = 0;
  if( gaussResp ){
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/gaussResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/gaussResp/data%i_dosys.root", norder_) );
  }
  if( studTResp ){
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/studTResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/studTResp/data%i_dosys.root", norder_) );
  }
  if( dataResp ){
    if( !dosys ) fUnfold = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("../../UnfoldResults/dataResp/data%i_dosys.root", norder_) );
  }

  fRelErr   = new TFile("relErrorRegularization.root", "recreate");
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

  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	//-- Get the VN observed histogram
	hObs[icent][iEP][iqn] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObs[icent][iEP][iqn]->SetMaximum( 10.*hObs[icent][iEP][iqn]->GetMaximum() );

	//-- Reset the flags that say whether or not the chi2 cutoff has been reached
	for(int ichi = 0; ichi < NCHI2; ichi++) iterCutoff[ichi] = 0;

	//-- Loop over iterations
	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
          hUnfold[icent][iEP][iqn][i]->SetLineColor(col[i]);
          hUnfold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the Refolded histograms
	  hRefold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
          hRefold[icent][iEP][iqn][i]->SetLineWidth(2);
          hRefold[icent][iEP][iqn][i]->SetLineColor(col[i]);
          hRefold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Calculate the cumulants of each unfolded histogram
	  FixUnfold( hUnfold[icent][iEP][iqn][i] );
	  double meanvn     = hUnfold[icent][iEP][iqn][i]->GetMean();
	  double stdevvn    = hUnfold[icent][iEP][iqn][i]->GetRMS();
	  double rmsvn      = sqrt( meanvn*meanvn + stdevvn*stdevvn );
	  double relfluctvn = stdevvn / meanvn;
	  EbyECumu cumu(hUnfold[icent][iEP][iqn][i]);
	  double vn2       = cumu.GetCumu_vn2();
	  double vn4       = cumu.GetCumu_vn4();
	  double vn6       = cumu.GetCumu_vn6();
	  double vn8       = cumu.GetCumu_vn8();
	  double gamma1exp = cumu.GetGamma1Exp();
	  double vn6vn4;
	  if(vn4 == 0) vn6vn4 = 0;
	  else         vn6vn4 = vn6 / vn4;
	  double vn8vn4;
	  if(vn4 == 0) vn8vn4 = 0;
	  else         vn8vn4 = vn8 / vn4;
	  double vn8vn6;
	  if(vn6 == 0) vn8vn6 = 0;
	  else         vn8vn6 = vn8 / vn6;

	  //-- Calculate the refold chi2
	  double chi2NDF_Refold = hRefold[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");
	  chi2NDF[icent][iEP][iqn][i] = chi2NDF_Refold;

	  //-- loop over the scenarios and choose the iteration for each
	  for(int ichi = 0; ichi < NCHI2; ichi++){

	    if( chi2NDF_Refold <= chi2Cutoff[ichi] && !iterCutoff[ichi] ){
	      unfold_RMSVnVSChi2Cut[icent][iEP][iqn][ichi]      = rmsvn;
	      unfold_MeanVnVSChi2Cut[icent][iEP][iqn][ichi]     = meanvn;
	      unfold_StDevVnVSChi2Cut[icent][iEP][iqn][ichi]    = stdevvn;
	      unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][ichi] = relfluctvn;
	      unfold_Vn2VSChi2Cut[icent][iEP][iqn][ichi]        = vn2;
	      unfold_Vn4VSChi2Cut[icent][iEP][iqn][ichi]        = vn4;
	      unfold_Vn6VSChi2Cut[icent][iEP][iqn][ichi]        = vn6;
	      unfold_Vn8VSChi2Cut[icent][iEP][iqn][ichi]        = vn8;
	      unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][ichi]  = gamma1exp;
	      unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][ichi]     = vn6vn4;
	      unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][ichi]     = vn8vn4;
	      unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][ichi]     = vn8vn6;

	      iterCutoff[ichi] = 1;
	    }
	    if( i == NITER-1 && !iterCutoff[ichi] ){
	      unfold_RMSVnVSChi2Cut[icent][iEP][iqn][ichi]      = rmsvn;
              unfold_MeanVnVSChi2Cut[icent][iEP][iqn][ichi]     = meanvn;
              unfold_StDevVnVSChi2Cut[icent][iEP][iqn][ichi]    = stdevvn;
	      unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][ichi] = relfluctvn;
              unfold_Vn2VSChi2Cut[icent][iEP][iqn][ichi]        = vn2;
              unfold_Vn4VSChi2Cut[icent][iEP][iqn][ichi]        = vn4;
              unfold_Vn6VSChi2Cut[icent][iEP][iqn][ichi]        = vn6;
              unfold_Vn8VSChi2Cut[icent][iEP][iqn][ichi]        = vn8;
              unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][ichi]  = gamma1exp;
              unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][ichi]     = vn6vn4;
              unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][ichi]     = vn8vn4;
              unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][ichi]     = vn8vn6;

              iterCutoff[ichi] = 1;
	    }

	  } //-- End chi2 scenario loop

	} //-- End unfold iteration loop

	//-- Ratio to the Chi2=1 cutoff scenario
	for(int ichi = 0; ichi < NCHI2; ichi++){

	  double rmsvn_i      = unfold_RMSVnVSChi2Cut[icent][iEP][iqn][ichi];
	  double rmsvn_0      = unfold_RMSVnVSChi2Cut[icent][iEP][iqn][0];
	  double meanvn_i     = unfold_MeanVnVSChi2Cut[icent][iEP][iqn][ichi];
          double meanvn_0     = unfold_MeanVnVSChi2Cut[icent][iEP][iqn][0];
	  double stdevvn_i    = unfold_StDevVnVSChi2Cut[icent][iEP][iqn][ichi];
          double stdevvn_0    = unfold_StDevVnVSChi2Cut[icent][iEP][iqn][0];
	  double relfluctvn_i = unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][ichi];
          double relfluctvn_0 = unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][0];

	  double vn2_i = unfold_Vn2VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn2_0 = unfold_Vn2VSChi2Cut[icent][iEP][iqn][0];
	  double vn4_i = unfold_Vn4VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn4_0 = unfold_Vn4VSChi2Cut[icent][iEP][iqn][0];
	  double vn6_i = unfold_Vn6VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn6_0 = unfold_Vn6VSChi2Cut[icent][iEP][iqn][0];
	  double vn8_i = unfold_Vn8VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn8_0 = unfold_Vn8VSChi2Cut[icent][iEP][iqn][0];

	  double gamma1exp_i = unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][ichi];
	  double gamma1exp_0 = unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][0];
	  double vn6vn4_i    = unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn6vn4_0    = unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][0];
	  double vn8vn4_i    = unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn8vn4_0    = unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][0];
	  double vn8vn6_i    = unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][ichi];
	  double vn8vn6_0    = unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][0];

	  double rmsvnRat      = rmsvn_i / rmsvn_0;
	  double meanvnRat     = meanvn_i / meanvn_0;
	  double stdevvnRat    = stdevvn_i / stdevvn_0;
	  double relfluctvnRat = relfluctvn_i / relfluctvn_0;
	  double vn2Rat        = vn2_i / vn2_0;
	  double vn4Rat        = vn4_i / vn4_0;
	  double vn6Rat        = vn6_i / vn6_0;
	  double vn8Rat        = vn8_i / vn8_0;
	  double g1expRat      = gamma1exp_i / gamma1exp_0;
	  double vn6vn4Rat     = vn6vn4_i / vn6vn4_0;
	  double vn8vn4Rat     = vn8vn4_i / vn8vn4_0;
	  double vn8vn6Rat     = vn8vn6_i / vn8vn6_0;

	  unfoldRMSVn_RatioToChi21[icent][iEP][iqn][ichi]      = rmsvnRat;
	  unfoldMeanVn_RatioToChi21[icent][iEP][iqn][ichi]     = meanvnRat;
	  unfoldStDevVn_RatioToChi21[icent][iEP][iqn][ichi]    = stdevvnRat;
	  unfoldRelFluctVn_RatioToChi21[icent][iEP][iqn][ichi] = relfluctvnRat;
	  unfoldVn2_RatioToChi21[icent][iEP][iqn][ichi]        = vn2Rat;
	  unfoldVn4_RatioToChi21[icent][iEP][iqn][ichi]        = vn4Rat;
	  unfoldVn6_RatioToChi21[icent][iEP][iqn][ichi]        = vn6Rat;
	  unfoldVn8_RatioToChi21[icent][iEP][iqn][ichi]        = vn8Rat;
	  unfoldGamma1Exp_RatioToChi21[icent][iEP][iqn][ichi]  = g1expRat;
	  unfoldVn6Vn4_RatioToChi21[icent][iEP][iqn][ichi]     = vn6vn4Rat;
	  unfoldVn8Vn4_RatioToChi21[icent][iEP][iqn][ichi]     = vn8vn4Rat;
	  unfoldVn8Vn6_RatioToChi21[icent][iEP][iqn][ichi]     = vn8vn6Rat;

	} //-- End chi2 loop for ratios

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //---------------------------------------------------------------------------------------------------- 
  //-- Save Relative Error to Root File
  vector<double> rrmsvn;
  vector<double> rmeanvn;
  vector<double> rstdevvn;
  vector<double> rrelfluctvn;
  vector<double> rvn2;
  vector<double> rvn4;
  vector<double> rvn6;
  vector<double> rvn8;
  vector<double> rg1exp;
  vector<double> rvn6vn4;
  vector<double> rvn8vn4;
  vector<double> rvn8vn6;

  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	rrmsvn.clear();
	rmeanvn.clear();
	rstdevvn.clear();
	rrelfluctvn.clear();
	rvn2.clear();
	rvn4.clear();
	rvn6.clear();
	rvn8.clear();
	rg1exp.clear();
	rvn6vn4.clear();
	rvn8vn4.clear();
	rvn8vn6.clear();

	for(int i = 0; i < 4; i++){
      
	  double rmsvnPct      = fabs(unfold_RMSVnVSChi2Cut[icent][iEP][iqn][i] - unfold_RMSVnVSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_RMSVnVSChi2Cut[icent][iEP][iqn][0]);
	  double meanvnPct     = fabs(unfold_MeanVnVSChi2Cut[icent][iEP][iqn][i] - unfold_MeanVnVSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_MeanVnVSChi2Cut[icent][iEP][iqn][0]);
	  double stdevvnPct    = fabs(unfold_StDevVnVSChi2Cut[icent][iEP][iqn][i] - unfold_StDevVnVSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_StDevVnVSChi2Cut[icent][iEP][iqn][0]);
	  double relfluctvnPct = fabs(unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][i] - unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_RelFluctVnVSChi2Cut[icent][iEP][iqn][0]);
	  double vn2Pct        = fabs(unfold_Vn2VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn2VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn2VSChi2Cut[icent][iEP][iqn][0]);
	  double vn4Pct        = fabs(unfold_Vn4VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn4VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn4VSChi2Cut[icent][iEP][iqn][0]);
	  double vn6Pct        = fabs(unfold_Vn6VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn6VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn6VSChi2Cut[icent][iEP][iqn][0]);
	  double vn8Pct        = fabs(unfold_Vn8VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn8VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn8VSChi2Cut[icent][iEP][iqn][0]);
	  double g1expPct      = fabs(unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][i] - unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Gamma1ExpVSChi2Cut[icent][iEP][iqn][0]);
	  double vn6vn4Pct     = fabs(unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn6Vn4VSChi2Cut[icent][iEP][iqn][0]);
	  double vn8vn4Pct     = fabs(unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn8Vn4VSChi2Cut[icent][iEP][iqn][0]);
	  double vn8vn6Pct     = fabs(unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][i] - unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][0]) / fabs(unfold_Vn8Vn6VSChi2Cut[icent][iEP][iqn][0]);

	  if( TMath::IsNaN( vn2Pct ) )    vn2Pct    = 0;
	  if( TMath::IsNaN( vn4Pct ) )    vn4Pct    = 0;
	  if( TMath::IsNaN( vn6Pct ) )    vn6Pct    = 0;
	  if( TMath::IsNaN( vn8Pct ) )    vn8Pct    = 0;
	  if( TMath::IsNaN( g1expPct ) )  g1expPct  = 0;
	  if( TMath::IsNaN( vn6vn4Pct ) ) vn6vn4Pct = 0;
	  if( TMath::IsNaN( vn8vn4Pct ) ) vn8vn4Pct = 0;
	  if( TMath::IsNaN( vn8vn6Pct ) ) vn8vn6Pct = 0;

	  if( std::isinf( vn2Pct ) )    vn2Pct    = 0;
	  if( std::isinf( vn4Pct ) )    vn4Pct    = 0;
	  if( std::isinf( vn6Pct ) )    vn6Pct    = 0;
	  if( std::isinf( vn8Pct ) )    vn8Pct    = 0;
	  if( std::isinf( g1expPct ) )  g1expPct  = 0;
	  if( std::isinf( vn6vn4Pct ) ) vn6vn4Pct = 0;
	  if( std::isinf( vn8vn4Pct ) ) vn8vn4Pct = 0;
	  if( std::isinf( vn8vn6Pct ) ) vn8vn6Pct = 0;

	  rrmsvn.push_back(rmsvnPct);
	  rmeanvn.push_back(meanvnPct);
	  rstdevvn.push_back(stdevvnPct);
	  rrelfluctvn.push_back(relfluctvnPct);
	  rvn2.push_back(vn2Pct);
	  rvn4.push_back(vn4Pct);
	  rvn6.push_back(vn6Pct);
	  rvn8.push_back(vn8Pct);
	  rg1exp.push_back(g1expPct);
	  rvn6vn4.push_back(vn6vn4Pct);
	  rvn8vn4.push_back(vn8vn4Pct);
	  rvn8vn6.push_back(vn8vn6Pct);
	  
	}

	double reRMSVn      = *max_element(std::begin(rrmsvn), std::end(rrmsvn));
	double reMeanVn     = *max_element(std::begin(rmeanvn), std::end(rmeanvn));
	double reStDevVn    = *max_element(std::begin(rstdevvn), std::end(rstdevvn));
	double reRelFluctVn = *max_element(std::begin(rrelfluctvn), std::end(rrelfluctvn));
	double reVn2        = *max_element(std::begin(rvn2), std::end(rvn2));
	double reVn4        = *max_element(std::begin(rvn4), std::end(rvn4));
	double reVn6        = *max_element(std::begin(rvn6), std::end(rvn6));
	double reVn8        = *max_element(std::begin(rvn8), std::end(rvn8));
	double reG1Exp      = *max_element(std::begin(rg1exp), std::end(rg1exp));    
	double reVn6Vn4     = *max_element(std::begin(rvn6vn4), std::end(rvn6vn4));
	double reVn8Vn4     = *max_element(std::begin(rvn8vn4), std::end(rvn8vn4));
	double reVn8Vn6     = *max_element(std::begin(rvn8vn6), std::end(rvn8vn6));

	relErrRMSVn[iEP][iqn]->SetBinContent(icent+1, reRMSVn);
	relErrMeanVn[iEP][iqn]->SetBinContent(icent+1, reMeanVn);
	relErrStDevVn[iEP][iqn]->SetBinContent(icent+1, reStDevVn);
	relErrRelFluctVn[iEP][iqn]->SetBinContent(icent+1, reRelFluctVn);
	relErrVn2[iEP][iqn]->SetBinContent(icent+1, reVn2);
	relErrVn4[iEP][iqn]->SetBinContent(icent+1, reVn4);
	relErrVn6[iEP][iqn]->SetBinContent(icent+1, reVn6);
	relErrVn8[iEP][iqn]->SetBinContent(icent+1, reVn8);
	relErrGamma1Exp[iEP][iqn]->SetBinContent(icent+1, reG1Exp);
	relErrVn6Vn4[iEP][iqn]->SetBinContent(icent+1, reVn6Vn4);
	relErrVn8Vn4[iEP][iqn]->SetBinContent(icent+1, reVn8Vn4);
	relErrVn8Vn6[iEP][iqn]->SetBinContent(icent+1, reVn8Vn6);

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

  fRelErr->Write();
  
}

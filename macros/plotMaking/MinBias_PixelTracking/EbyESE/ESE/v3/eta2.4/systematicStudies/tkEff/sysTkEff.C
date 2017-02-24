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

void sysTkEff(){

  int norder_ = 3;

  TFile * fAna;
  TH1D * hObs[NCENT][NEPSymm][NQN];

  //-- Standard Unfolding
  TFile * fUnf;
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NITER];

  //-- NoTeff Unfolding
  TFile * fUnfNoTeff;
  TH1D * hUnfoldNoTeff[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldNoTeff[NCENT][NEPSymm][NQN][NITER];

  //-- Iter Cutoffs
  bool iterStopFound[NCENT][NEPSymm][NQN];
  int iterCutoff[NCENT][NEPSymm][NQN];

  bool iterStopFoundNoTeff[NCENT][NEPSymm][NQN];
  int iterCutoffNoTeff[NCENT][NEPSymm][NQN];

  TLatex latex;

  //-- Systematics
  double rmsVnNoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double meanVnNoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double stDevVnNoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double relFluctVnNoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn2NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn4NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn6NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double gamma1expNoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn6vn4NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8vn4NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8vn6NoTeff_RatioToNominal[NCENT][NEPSymm][NQN];

  double rmsVnNoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double meanVnNoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double stDevVnNoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double relFluctVnNoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn2NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn4NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn6NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double gamma1expNoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn6vn4NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8vn4NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8vn6NoTeff_PctDiffToNominal[NCENT][NEPSymm][NQN];

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
  //-- MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  gErrorIgnoreLevel = kWarning;

  fAna       = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fUnf       = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnfNoTeff = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );

  fRelErr   = new TFile("relErrorTkEff.root", "recreate");
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

	//-- Get the VN observed histogram
	hObs[icent][iEP][iqn] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hObs[icent][iEP][iqn]->SetMaximum( 10.*hObs[icent][iEP][iqn]->GetMaximum() );

	//-- Initialize iteration cutoff values/booleans
	iterStopFound[icent][iEP][iqn] = 0;
	iterCutoff[icent][iEP][iqn]    = 0;
	
	iterStopFoundNoTeff[icent][iEP][iqn] = 0;
	iterCutoffNoTeff[icent][iEP][iqn]    = 0;

	for(int i = 0; i < NITER; i++){

	  //-- Get the unfolded histograms
	  hUnfold[icent][iEP][iqn][i] = (TH1D*) fUnf->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfold[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hUnfold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the Refolded histograms
	  hRefold[icent][iEP][iqn][i] = (TH1D*) fUnf->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefold[icent][iEP][iqn][i]->SetLineWidth(2);
	  hRefold[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hRefold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the NoTeff unfolded histograms
	  hUnfoldNoTeff[icent][iEP][iqn][i] = (TH1D*) fUnfNoTeff->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldNoTeff[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hUnfoldNoTeff[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the NoTeff refolded histograms
	  hRefoldNoTeff[icent][iEP][iqn][i] = (TH1D*) fUnfNoTeff->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldNoTeff[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hRefoldNoTeff[icent][iEP][iqn][i]->SetMarkerColor(col[i]);


	  //-- Chi squares
	  double chi2NDF_Refold        = hRefold[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");
	  double chi2NDF_Refold_NoTeff = hRefoldNoTeff[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");


	  //-- Normal unfolding chi2 check
	  if( chi2NDF_Refold < 1.2 && !iterStopFound[icent][iEP][iqn] ){
	    iterStopFound[icent][iEP][iqn] = 1;
	    iterCutoff[icent][iEP][iqn]    = i;
	  }
	  if( i == (NITER - 1) && !iterStopFound[icent][iEP][iqn] ){
	    iterStopFound[icent][iEP][iqn] = 1;
	    iterCutoff[icent][iEP][iqn]    = i;
	  }

	  //-- NoTeff unfolding chi2 check
	  if( chi2NDF_Refold_NoTeff < 1.2 && !iterStopFoundNoTeff[icent][iEP][iqn] ){
	    iterStopFoundNoTeff[icent][iEP][iqn] = 1;
	    iterCutoffNoTeff[icent][iEP][iqn]    = i;
	  }
	  if( i == (NITER - 1) && !iterStopFoundNoTeff[icent][iEP][iqn] ){
	    iterStopFoundNoTeff[icent][iEP][iqn] = 1;
	    iterCutoffNoTeff[icent][iEP][iqn]    = i;
	  }

	} //-- End iter loop

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

  //-- Figure out the response matrix uncertainty component
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	std::cout << Form("=================== Cent %i %s Qbin %i ===================", icent, EPSymmNames[iEP].data(), iqn) << std::endl;

	int i       = iterCutoff[icent][iEP][iqn];
	int iNoTeff = iterCutoffNoTeff[icent][iEP][iqn];

	FixUnfold( hUnfold[icent][iEP][iqn][i] );
	FixUnfold( hUnfoldNoTeff[icent][iEP][iqn][iNoTeff] );

	EbyECumu cumu( hUnfold[icent][iEP][iqn][i] );
	EbyECumu cumuNoTeff( hUnfoldNoTeff[icent][iEP][iqn][iNoTeff] );

	double meanVn           = hUnfold[icent][iEP][iqn][i]->GetMean();
	double stDevVn          = hUnfold[icent][iEP][iqn][i]->GetRMS();
	double rmsVn            = sqrt( meanVn*meanVn + stDevVn*stDevVn );
	double relFluctVn       = stDevVn / meanVn;
	double meanVnNoTeff     = hUnfoldNoTeff[icent][iEP][iqn][i]->GetMean();
        double stDevVnNoTeff    = hUnfoldNoTeff[icent][iEP][iqn][i]->GetRMS();
        double rmsVnNoTeff      = sqrt( meanVnNoTeff*meanVnNoTeff + stDevVnNoTeff*stDevVnNoTeff );
	double relFluctVnNoTeff = stDevVnNoTeff / meanVnNoTeff;


	double vn2         = cumu.GetCumu_vn2();
	double vn4         = cumu.GetCumu_vn4();
	double vn6         = cumu.GetCumu_vn6();
	double vn8         = cumu.GetCumu_vn8();
	double vn2NoTeff   = cumuNoTeff.GetCumu_vn2();
	double vn4NoTeff   = cumuNoTeff.GetCumu_vn4();
	double vn6NoTeff   = cumuNoTeff.GetCumu_vn6();
	double vn8NoTeff   = cumuNoTeff.GetCumu_vn8();

	double gamma1exp       = cumu.GetGamma1Exp();
	double gamma1expNoTeff = cumuNoTeff.GetGamma1Exp();

	double vn6vn4;
	double vn6vn4NoTeff;
	if( vn6 == 0 || vn4 == 0 ) vn6vn4 = 0;
	else                       vn6vn4 = vn6 / vn4;
	if( vn6NoTeff == 0 || vn4NoTeff == 0) vn6vn4NoTeff = 0;
	else                                  vn6vn4NoTeff = vn6NoTeff / vn4NoTeff;

	double vn8vn4;
	double vn8vn4NoTeff;
	if( vn8 == 0 || vn4 == 0 ) vn8vn4 = 0;
	else                       vn8vn4 = vn8 / vn4;
	if( vn8NoTeff == 0 || vn4NoTeff == 0) vn8vn4NoTeff = 0;
	else                                  vn8vn4NoTeff = vn8NoTeff / vn4NoTeff;

	double vn8vn6;
	double vn8vn6NoTeff;
	if( vn8 == 0 || vn6 == 0 ) vn8vn6 = 0;
	else                       vn8vn6 = vn8 / vn6;
	if( vn8NoTeff == 0 || vn6NoTeff == 0) vn8vn6NoTeff = 0;
	else                                  vn8vn6NoTeff = vn8NoTeff / vn6NoTeff;

	//-- Calculate ratios
	rmsVnNoTeff_RatioToNominal[icent][iEP][iqn]   = rmsVnNoTeff / rmsVn;
	meanVnNoTeff_RatioToNominal[icent][iEP][iqn]  = meanVnNoTeff / meanVn;
	stDevVnNoTeff_RatioToNominal[icent][iEP][iqn] = stDevVnNoTeff / stDevVn;
	relFluctVnNoTeff_RatioToNominal[icent][iEP][iqn] = relFluctVnNoTeff / relFluctVn;

	if( vn2 == 0 || vn2NoTeff == 0 ) vn2NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                             vn2NoTeff_RatioToNominal[icent][iEP][iqn] = vn2NoTeff / vn2;

	if( vn4 == 0 || vn4NoTeff == 0 ) vn4NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                             vn4NoTeff_RatioToNominal[icent][iEP][iqn] = vn4NoTeff / vn4;

	if( vn6 == 0 || vn6NoTeff == 0 ) vn6NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                             vn6NoTeff_RatioToNominal[icent][iEP][iqn] = vn6NoTeff / vn6;

	if( vn8 == 0 || vn8NoTeff == 0 ) vn8NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                             vn8NoTeff_RatioToNominal[icent][iEP][iqn] = vn8NoTeff / vn8;

	if( gamma1exp == 0 || gamma1expNoTeff == 0 ) gamma1expNoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                                         gamma1expNoTeff_RatioToNominal[icent][iEP][iqn] = gamma1expNoTeff / gamma1exp;

	if( vn6vn4 == 0 || vn6vn4NoTeff == 0 ) vn6vn4NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                                   vn6vn4NoTeff_RatioToNominal[icent][iEP][iqn] = vn6vn4NoTeff / vn6vn4;

	if( vn8vn4 == 0 || vn8vn4NoTeff == 0 ) vn8vn4NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                                   vn8vn4NoTeff_RatioToNominal[icent][iEP][iqn] = vn8vn4NoTeff / vn8vn4;

	if( vn8vn6 == 0 || vn8vn6NoTeff == 0 ) vn8vn6NoTeff_RatioToNominal[icent][iEP][iqn] = 0;
	else                                   vn8vn6NoTeff_RatioToNominal[icent][iEP][iqn] = vn8vn6NoTeff / vn8vn6;

	//-- Calculate pct difference relative to nominal
	rmsVnNoTeff_PctDiffToNominal[icent][iEP][iqn]   = fabs(rmsVnNoTeff - rmsVn) / fabs(rmsVn);
	meanVnNoTeff_PctDiffToNominal[icent][iEP][iqn]  = fabs(meanVnNoTeff - meanVn) / fabs(meanVn);
	stDevVnNoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs(stDevVnNoTeff - stDevVn) / fabs(stDevVn);
	relFluctVnNoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs(relFluctVnNoTeff - relFluctVn) / fabs(relFluctVn);

	if( vn2 == 0 || vn2NoTeff == 0 ) vn2NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                             vn2NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn2NoTeff - vn2 ) / fabs( vn2 );

	if( vn4 == 0 || vn4NoTeff == 0 ) vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                             vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn4NoTeff - vn4 ) / fabs( vn4 );

	if( vn6 == 0 || vn6NoTeff == 0 ) vn6NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                             vn6NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn6NoTeff - vn6 ) / fabs( vn6 );

	if( vn8 == 0 || vn8NoTeff == 0 ) vn8NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                             vn8NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8NoTeff - vn8 ) / fabs( vn8 );

	if( gamma1exp == 0 || gamma1expNoTeff == 0 ) gamma1expNoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                         gamma1expNoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( gamma1expNoTeff - gamma1exp ) / fabs( gamma1exp );

	if( vn6vn4 == 0 || vn6vn4NoTeff == 0 ) vn6vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                   vn6vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn6vn4NoTeff - vn6vn4 ) / fabs( vn6vn4 );

	if( vn8vn4 == 0 || vn8vn4NoTeff == 0 ) vn8vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                   vn8vn4NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8vn4NoTeff - vn8vn4 ) / fabs( vn8vn4 );

	if( vn8vn6 == 0 || vn8vn6NoTeff == 0 ) vn8vn6NoTeff_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                   vn8vn6NoTeff_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8vn6NoTeff - vn8vn6 ) / fabs( vn8vn6 );

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){
	relErrRMSVn[iEP][iqn]->SetBinContent(icent+1, rmsVnNoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrMeanVn[iEP][iqn]->SetBinContent(icent+1, meanVnNoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrStDevVn[iEP][iqn]->SetBinContent(icent+1, stDevVnNoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrRelFluctVn[iEP][iqn]->SetBinContent(icent+1, relFluctVnNoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn2[iEP][iqn]->SetBinContent(icent+1, vn2NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn4[iEP][iqn]->SetBinContent(icent+1, vn4NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn6[iEP][iqn]->SetBinContent(icent+1, vn6NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8[iEP][iqn]->SetBinContent(icent+1, vn6NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrGamma1Exp[iEP][iqn]->SetBinContent(icent+1, gamma1expNoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn6Vn4[iEP][iqn]->SetBinContent(icent+1, vn6vn4NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8Vn4[iEP][iqn]->SetBinContent(icent+1, vn8vn4NoTeff_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8Vn6[iEP][iqn]->SetBinContent(icent+1, vn8vn6NoTeff_PctDiffToNominal[icent][iEP][iqn]);
      }
    }
  }

  fRelErr->Write();

}

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

void sysRespEl(){

  int norder_ = 3;

  TFile * fAna;
  TH1D * hObs[NCENT][NEPSymm][NQN];

  //-- Standard Unfolding
  TFile * fUnf;
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NITER];

  //-- DoSys Unfolding
  TFile * fUnfDoSys;
  TH1D * hUnfoldDoSys[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldDoSys[NCENT][NEPSymm][NQN][NITER];

  //-- Iter Cutoffs
  bool iterStopFound[NCENT][NEPSymm][NQN];
  int iterCutoff[NCENT][NEPSymm][NQN];

  bool iterStopFoundDoSys[NCENT][NEPSymm][NQN];
  int iterCutoffDoSys[NCENT][NEPSymm][NQN];

  TLatex latex;

  //-- Systematics
  double rmsVnDoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double meanVnDoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double stDevVnDoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double relFluctVnDoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn2DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn4DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn6DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double gamma1expDoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn6vn4DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8vn4DoSys_RatioToNominal[NCENT][NEPSymm][NQN];
  double vn8vn6DoSys_RatioToNominal[NCENT][NEPSymm][NQN];

  double rmsVnDoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double meanVnDoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double stDevVnDoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double relFluctVnDoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn2DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn4DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn6DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double gamma1expDoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn6vn4DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8vn4DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];
  double vn8vn6DoSys_PctDiffToNominal[NCENT][NEPSymm][NQN];

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

  fAna      = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fUnf      = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fUnfDoSys = new TFile( Form("../../UnfoldResults/dataResp/data%i_dosys.root", norder_) );

  fRelErr   = new TFile("relErrorResponse.root", "recreate");
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
	
	iterStopFoundDoSys[icent][iEP][iqn] = 0;
	iterCutoffDoSys[icent][iEP][iqn]    = 0;

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

	  //-- Get the DoSys unfolded histograms
	  hUnfoldDoSys[icent][iEP][iqn][i] = (TH1D*) fUnfDoSys->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfoldDoSys[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hUnfoldDoSys[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the DoSys refolded histograms
	  hRefoldDoSys[icent][iEP][iqn][i] = (TH1D*) fUnfDoSys->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldDoSys[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hRefoldDoSys[icent][iEP][iqn][i]->SetMarkerColor(col[i]);


	  //-- Chi squares
	  double chi2NDF_Refold       = hRefold[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");
	  double chi2NDF_Refold_DoSys = hRefoldDoSys[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");


	  //-- Normal unfolding chi2 check
	  if( chi2NDF_Refold < 1.2 && !iterStopFound[icent][iEP][iqn] ){
	    iterStopFound[icent][iEP][iqn] = 1;
	    iterCutoff[icent][iEP][iqn]    = i;
	  }
	  if( i == (NITER - 1) && !iterStopFound[icent][iEP][iqn] ){
	    iterStopFound[icent][iEP][iqn] = 1;
	    iterCutoff[icent][iEP][iqn]    = i;
	  }

	  //-- DoSys unfolding chi2 check
	  if( chi2NDF_Refold_DoSys < 1.2 && !iterStopFoundDoSys[icent][iEP][iqn] ){
	    iterStopFoundDoSys[icent][iEP][iqn] = 1;
	    iterCutoffDoSys[icent][iEP][iqn]    = i;
	  }
	  if( i == (NITER - 1) && !iterStopFoundDoSys[icent][iEP][iqn] ){
	    iterStopFoundDoSys[icent][iEP][iqn] = 1;
	    iterCutoffDoSys[icent][iEP][iqn]    = i;
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

	int i      = iterCutoff[icent][iEP][iqn];
	int iDoSys = iterCutoffDoSys[icent][iEP][iqn];

	FixUnfold( hUnfold[icent][iEP][iqn][i] );
	FixUnfold( hUnfoldDoSys[icent][iEP][iqn][iDoSys] );

	EbyECumu cumu( hUnfold[icent][iEP][iqn][i] );
	EbyECumu cumuDoSys( hUnfoldDoSys[icent][iEP][iqn][iDoSys] );

	double meanVn          = hUnfold[icent][iEP][iqn][i]->GetMean();
	double stDevVn         = hUnfold[icent][iEP][iqn][i]->GetRMS();
	double rmsVn           = sqrt( meanVn*meanVn + stDevVn*stDevVn );
	double relFluctVn      = stDevVn / meanVn;
	double meanVnDoSys     = hUnfoldDoSys[icent][iEP][iqn][i]->GetMean();
        double stDevVnDoSys    = hUnfoldDoSys[icent][iEP][iqn][i]->GetRMS();
        double rmsVnDoSys      = sqrt( meanVnDoSys*meanVnDoSys + stDevVnDoSys*stDevVnDoSys );
	double relFluctVnDoSys = stDevVnDoSys / meanVnDoSys;

	double vn2        = cumu.GetCumu_vn2();
	double vn4        = cumu.GetCumu_vn4();
	double vn6        = cumu.GetCumu_vn6();
	double vn8        = cumu.GetCumu_vn8();
	double vn2DoSys   = cumuDoSys.GetCumu_vn2();
	double vn4DoSys   = cumuDoSys.GetCumu_vn4();
	double vn6DoSys   = cumuDoSys.GetCumu_vn6();
	double vn8DoSys   = cumuDoSys.GetCumu_vn8();

	double gamma1exp      = cumu.GetGamma1Exp();
	double gamma1expDoSys = cumuDoSys.GetGamma1Exp();

	double vn6vn4;
	double vn6vn4DoSys;
	if( vn6 == 0 || vn4 == 0 ) vn6vn4 = 0;
	else                       vn6vn4 = vn6 / vn4;
	if( vn6DoSys == 0 || vn4DoSys == 0) vn6vn4DoSys = 0;
	else                                vn6vn4DoSys = vn6DoSys / vn4DoSys;

	double vn8vn4;
	double vn8vn4DoSys;
	if( vn8 == 0 || vn4 == 0 ) vn8vn4 = 0;
	else                       vn8vn4 = vn8 / vn4;
	if( vn8DoSys == 0 || vn4DoSys == 0) vn8vn4DoSys = 0;
	else                                vn8vn4DoSys = vn8DoSys / vn4DoSys;

	double vn8vn6;
	double vn8vn6DoSys;
	if( vn8 == 0 || vn6 == 0 ) vn8vn6 = 0;
	else                       vn8vn6 = vn8 / vn6;
	if( vn8DoSys == 0 || vn6DoSys == 0) vn8vn6DoSys = 0;
	else                                vn8vn6DoSys = vn8DoSys / vn6DoSys;

	//-- Calculate ratios
	rmsVnDoSys_RatioToNominal[icent][iEP][iqn]      = rmsVnDoSys / rmsVn;
	meanVnDoSys_RatioToNominal[icent][iEP][iqn]     = meanVnDoSys / meanVn;
	stDevVnDoSys_RatioToNominal[icent][iEP][iqn]    = stDevVnDoSys / stDevVn;
	relFluctVnDoSys_RatioToNominal[icent][iEP][iqn] = relFluctVnDoSys / relFluctVn;

	if( vn2 == 0 || vn2DoSys == 0 ) vn2DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                            vn2DoSys_RatioToNominal[icent][iEP][iqn] = vn2DoSys / vn2;

	if( vn4 == 0 || vn4DoSys == 0 ) vn4DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                            vn4DoSys_RatioToNominal[icent][iEP][iqn] = vn4DoSys / vn4;

	if( vn6 == 0 || vn6DoSys == 0 ) vn6DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                            vn6DoSys_RatioToNominal[icent][iEP][iqn] = vn6DoSys / vn6;

	if( vn8 == 0 || vn8DoSys == 0 ) vn8DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                            vn8DoSys_RatioToNominal[icent][iEP][iqn] = vn8DoSys / vn8;

	if( gamma1exp == 0 || gamma1expDoSys == 0 ) gamma1expDoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                                        gamma1expDoSys_RatioToNominal[icent][iEP][iqn] = gamma1expDoSys / gamma1exp;

	if( vn6vn4 == 0 || vn6vn4DoSys == 0 ) vn6vn4DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                                  vn6vn4DoSys_RatioToNominal[icent][iEP][iqn] = vn6vn4DoSys / vn6vn4;

	if( vn8vn4 == 0 || vn8vn4DoSys == 0 ) vn8vn4DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                                  vn8vn4DoSys_RatioToNominal[icent][iEP][iqn] = vn8vn4DoSys / vn8vn4;

	if( vn8vn6 == 0 || vn8vn6DoSys == 0 ) vn8vn6DoSys_RatioToNominal[icent][iEP][iqn] = 0;
	else                                  vn8vn6DoSys_RatioToNominal[icent][iEP][iqn] = vn8vn6DoSys / vn8vn6;

	//-- Calculate pct difference relative to nominal
	rmsVnDoSys_PctDiffToNominal[icent][iEP][iqn]      = fabs(rmsVnDoSys - rmsVn) / fabs(rmsVn);
	meanVnDoSys_PctDiffToNominal[icent][iEP][iqn]     = fabs(meanVnDoSys - meanVn) / fabs(meanVn);
	stDevVnDoSys_PctDiffToNominal[icent][iEP][iqn]    = fabs(stDevVnDoSys - stDevVn) / fabs(stDevVn);
	relFluctVnDoSys_PctDiffToNominal[icent][iEP][iqn] = fabs(relFluctVnDoSys - relFluctVn) / fabs(relFluctVn);

	if( vn2 == 0 || vn2DoSys == 0 ) vn2DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                            vn2DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn2DoSys - vn2 ) / fabs( vn2 );

	if( vn4 == 0 || vn4DoSys == 0 ) vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                            vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn4DoSys - vn4 ) / fabs( vn4 );

	if( vn6 == 0 || vn6DoSys == 0 ) vn6DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                            vn6DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn6DoSys - vn6 ) / fabs( vn6 );

	if( vn8 == 0 || vn8DoSys == 0 ) vn8DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                            vn8DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8DoSys - vn8 ) / fabs( vn8 );

	if( gamma1exp == 0 || gamma1expDoSys == 0 ) gamma1expDoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                        gamma1expDoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( gamma1expDoSys - gamma1exp ) / fabs( gamma1exp );

	if( vn6vn4 == 0 || vn6vn4DoSys == 0 ) vn6vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                  vn6vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn6vn4DoSys - vn6vn4 ) / fabs( vn6vn4 );

	if( vn8vn4 == 0 || vn8vn4DoSys == 0 ) vn8vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                  vn8vn4DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8vn4DoSys - vn8vn4 ) / fabs( vn8vn4 );

	if( vn8vn6 == 0 || vn8vn6DoSys == 0 ) vn8vn6DoSys_PctDiffToNominal[icent][iEP][iqn] = 0;
	else                                  vn8vn6DoSys_PctDiffToNominal[icent][iEP][iqn] = fabs( vn8vn6DoSys - vn8vn6 ) / fabs( vn8vn6 );

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop
  gErrorIgnoreLevel = kError;

  //-- Save Relative Errors
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){
	relErrRMSVn[iEP][iqn]->SetBinContent(icent+1, rmsVnDoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrMeanVn[iEP][iqn]->SetBinContent(icent+1, meanVnDoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrStDevVn[iEP][iqn]->SetBinContent(icent+1, stDevVnDoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrRelFluctVn[iEP][iqn]->SetBinContent(icent+1, relFluctVnDoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn2[iEP][iqn]->SetBinContent(icent+1, vn2DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn4[iEP][iqn]->SetBinContent(icent+1, vn4DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn6[iEP][iqn]->SetBinContent(icent+1, vn6DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8[iEP][iqn]->SetBinContent(icent+1, vn6DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrGamma1Exp[iEP][iqn]->SetBinContent(icent+1, gamma1expDoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn6Vn4[iEP][iqn]->SetBinContent(icent+1, vn6vn4DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8Vn4[iEP][iqn]->SetBinContent(icent+1, vn8vn4DoSys_PctDiffToNominal[icent][iEP][iqn]);
	relErrVn8Vn6[iEP][iqn]->SetBinContent(icent+1, vn8vn6DoSys_PctDiffToNominal[icent][iEP][iqn]);
      }
    }
  }

  fRelErr->Write();

}

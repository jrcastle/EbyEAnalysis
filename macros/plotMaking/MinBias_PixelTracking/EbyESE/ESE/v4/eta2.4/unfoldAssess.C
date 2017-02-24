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
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void unfoldAssess(){

  int norder_     = 4;
  int QnBinOrder_ = 2;
  double tkEta    = 2.4;

  bool RFplot     = 1;
  bool SaveSumPl  = 1;

  bool dosys     = 0;

  bool looseChi2IterCut   = 0;
  bool nominalChi2IterCut = 1;
  bool tightChi2IterCut   = 0;

  double momentMin = 0.0;
  double momentMax = 0.39;

  double momentTitleSize = 0.05;
  double momentTitleOffS = 1.5;
  double momentAxisLabS  = 0.04;
  double qnmin           = 0;
  double qnmax           = 0.26;
  double chi2Min         = 0.1;
  double chi2Max         = 10000.;

  TLatex latex;
  TLatex latex2;

  //-- Qn binning
  TFile * fQn;
  TH1D * hqbins[NCENT][NEPSymm];
  TH1D * hqnHF_EP[NCENT][NEPSymm];

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT][NEPSymm][NQN];

  //-- Unfolding output
  TFile * fUnfold;
  TH2D * hResp[NCENT][NEPSymm][NQN];
  TH1D * hUnfold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefoldRatio[NCENT][NEPSymm][NQN][NITER];

  //-- Unfolding moments
  double unfold_mean[NCENT][NEPSymm][NQN][NITER];
  double unfold_stDev[NCENT][NEPSymm][NQN][NITER];
  double unfold_RMS[NCENT][NEPSymm][NQN][NITER];
  double unfold_RelFluct[NCENT][NEPSymm][NQN][NITER];

  double unfold_meane[NCENT][NEPSymm][NQN][NITER];
  double unfold_stDeve[NCENT][NEPSymm][NQN][NITER];
  double unfold_RMSe[NCENT][NEPSymm][NQN][NITER];
  double unfold_RelFlucte[NCENT][NEPSymm][NQN][NITER];

  //-- Moment vs iter graphs
  TMultiGraph * grMoments[NCENT][NEPSymm][NQN];
  TGraphErrors * grUnfoldMean[NCENT][NEPSymm][NQN];
  TGraphErrors * grUnfoldStDev[NCENT][NEPSymm][NQN];
  TGraphErrors * grUnfoldRMS[NCENT][NEPSymm][NQN];

  //-- Chi Squares for iteration cutoffs
  double chi2NDF_Refold[NCENT][NEPSymm][NQN][NITER];

  //-- Chi square plots
  double chi2Cut;
  TGraphErrors * grChi2NDF_Refold[NCENT][NEPSymm][NQN];

  //-- Canvases
  TCanvas * cSummary[NCENT][NEPSymm][NQN];

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(15);

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

  //-- Get the unfolding output file
  fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );


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
	hObs[icent][iEP][iqn]->SetMaximum( 10.*hObs[icent][iEP][iqn]->GetMaximum() );
	hObs[icent][iEP][iqn]->SetMarkerStyle(24);

	bool cutoff = 0;
	for(int i = 0; i < NITER; i++){

	  //-- Get the Response histogram
	  hResp[icent][iEP][iqn] = (TH2D*) fUnfold->Get( Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Get the unfolded histograms
	  hUnfold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hUnfold[icent][iEP][iqn][i]->SetLineColor(col[i]);
	  hUnfold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Get the Refolded histograms
	  hRefold[icent][iEP][iqn][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefold[icent][iEP][iqn][i]->SetLineWidth(2);
	  hRefold[icent][iEP][iqn][i]->SetLineColor(col[i]);
          hRefold[icent][iEP][iqn][i]->SetMarkerColor(col[i]);

	  //-- Make the refold/observed histograms
	  hRefoldRatio[icent][iEP][iqn][i] = (TH1D*) hRefold[icent][iEP][iqn][i]->Clone( Form("hRefoldRatio%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );
	  hRefoldRatio[icent][iEP][iqn][i]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
	  hRefoldRatio[icent][iEP][iqn][i]->GetYaxis()->SetTitle( Form("v_{%i}^{Refold} / v_{%i}^{Obs}", norder_, norder_) );
	  hRefoldRatio[icent][iEP][iqn][i]->Divide( hObs[icent][iEP][iqn] );

	  //-- Calculate the moments of each unfolded histogram
	  double uMu   = hUnfold[icent][iEP][iqn][i]->GetMean();
	  double uMue  = hUnfold[icent][iEP][iqn][i]->GetMeanError();
	  double uSig  = hUnfold[icent][iEP][iqn][i]->GetRMS();
	  double uSige = hUnfold[icent][iEP][iqn][i]->GetRMSError();
	  double uRMS  = TMath::Sqrt( pow(uMu, 2) + pow(uSig, 2) );
	  double uRMSe = TMath::Sqrt( ( pow(uMu * uMue, 2) + pow(uSig * uSige, 2) ) / ( pow(uMu,2) + pow(uSig,2) ) );
	  double uRF   = uSig / uMu;
	  double uRFe  = TMath::Sqrt( pow(uSig*uMue/uMu/uMu, 2) + pow(uSige/uMu, 2) );

	  unfold_mean[icent][iEP][iqn][i]      = uMu;
	  unfold_RMS[icent][iEP][iqn][i]       = uRMS;
	  unfold_stDev[icent][iEP][iqn][i]     = uSig;
	  unfold_RelFluct[icent][iEP][iqn][i]  = uRF;

	  unfold_meane[icent][iEP][iqn][i]     = uMue;
	  unfold_RMSe[icent][iEP][iqn][i]      = uRMSe;
	  unfold_stDeve[icent][iEP][iqn][i]    = uSige;
	  unfold_RelFlucte[icent][iEP][iqn][i] = uRFe;

	  //-- Chi squares
	  chi2NDF_Refold[icent][iEP][iqn][i] = hRefold[icent][iEP][iqn][i]->Chi2Test(hObs[icent][iEP][iqn], "CHI2/NDF");

	  if( chi2NDF_Refold[icent][iEP][iqn][i] <= chi2Cut && !cutoff ){
	    cutoff = 1;
	  }
	  if( i == NITER-1 && !cutoff ){
	    cutoff = 1;
	  }

	} //-- End unfold iteration loop

	//-- Set up all the wonderful graphs
	grUnfoldMean[icent][iEP][iqn]  = new TGraphErrors(NITER, diter, unfold_mean[icent][iEP][iqn], iterErr, unfold_meane[icent][iEP][iqn] );
	grUnfoldMean[icent][iEP][iqn]->GetXaxis()->SetTitle("Iteration");
	grUnfoldMean[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	grUnfoldMean[icent][iEP][iqn]->GetYaxis()->SetTitle("Moment");
	grUnfoldMean[icent][iEP][iqn]->SetLineColor(8);
	grUnfoldMean[icent][iEP][iqn]->SetMarkerColor(8);
	grUnfoldMean[icent][iEP][iqn]->SetMarkerStyle(20);

	grUnfoldStDev[icent][iEP][iqn] = new TGraphErrors(NITER, diter, unfold_stDev[icent][iEP][iqn], iterErr, unfold_stDeve[icent][iEP][iqn] ); //46
	grUnfoldStDev[icent][iEP][iqn]->GetXaxis()->SetTitle("Iteration");
	grUnfoldStDev[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	grUnfoldStDev[icent][iEP][iqn]->GetYaxis()->SetTitle("Moment");
        grUnfoldStDev[icent][iEP][iqn]->SetLineColor(46);
	grUnfoldStDev[icent][iEP][iqn]->SetMarkerColor(46);
	grUnfoldStDev[icent][iEP][iqn]->SetMarkerStyle(20);

	grUnfoldRMS[icent][iEP][iqn]   = new TGraphErrors(NITER, diter, unfold_RMS[icent][iEP][iqn], iterErr, unfold_RMSe[icent][iEP][iqn] ); //9
	grUnfoldRMS[icent][iEP][iqn]->GetXaxis()->SetTitle("Iteration");
	grUnfoldRMS[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	grUnfoldRMS[icent][iEP][iqn]->GetYaxis()->SetTitle("Moment");
        grUnfoldRMS[icent][iEP][iqn]->SetLineColor(9);
	grUnfoldRMS[icent][iEP][iqn]->SetMarkerColor(9);
	grUnfoldRMS[icent][iEP][iqn]->SetMarkerStyle(20);

	grChi2NDF_Refold[icent][iEP][iqn]    = new TGraphErrors(NITER, diter, chi2NDF_Refold[icent][iEP][iqn],    iterErr, iterErr);
	grChi2NDF_Refold[icent][iEP][iqn]->GetXaxis()->SetTitle("Iteration");
	grChi2NDF_Refold[icent][iEP][iqn]->GetYaxis()->SetTitle("#chi^{2}/NDF");
	grChi2NDF_Refold[icent][iEP][iqn]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
	grChi2NDF_Refold[icent][iEP][iqn]->SetLineColor(4);
	grChi2NDF_Refold[icent][iEP][iqn]->SetMarkerColor(4);
	grChi2NDF_Refold[icent][iEP][iqn]->SetMarkerStyle(20);

      } //-- End QN loop
    } //-- End EP loop
  } //-- End cent loop

  TLegend * leg = new TLegend(0.1283, 0.1440, 0.8108, 0.7539);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hObs[0][0][0], Form("v_{%i}^{obs}", norder_), "lp");
  for(int i = 0; i < NITER; i++) leg->AddEntry(hUnfold[0][0][0][i], Form("D'Agostini iter = %i", iter[i]), "lp");

  TLegend * legMoment = new TLegend(0.6025, 0.7404, 0.9742, 0.8824);
  legMoment->SetFillStyle(0);
  legMoment->SetBorderSize(0);
  legMoment->SetNColumns(2);
  legMoment->AddEntry(grUnfoldRMS[0][0][0],   "#sqrt{#LTv_{2}^{2}#GT}", "lp");
  legMoment->AddEntry(grUnfoldMean[0][0][0],  "#LTv_{2}#GT", "lp");
  legMoment->AddEntry(grUnfoldStDev[0][0][0], "#sigma_{v_{2}}", "lp");

  TLine * line = new TLine(iter[0], 1., iter[NITER-1], 1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  //-- DRAW!!!!
  for(int icent = 0; icent < NCENT; icent++){

    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	cSummary[icent][iEP][iqn] = new TCanvas( Form("summary_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("summary_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 1000, 600);
	cSummary[icent][iEP][iqn]->Divide(3,2);

	cSummary[icent][iEP][iqn]->cd(1);
	cSummary[icent][iEP][iqn]->cd(1)->SetLogy();
	hObs[icent][iEP][iqn]->Draw();
	for(int i = 0; i < NITER; i++) hUnfold[icent][iEP][iqn][i]->Draw("same");
	latex.DrawLatex(0.18, 0.88, "Unfolding");

	cSummary[icent][iEP][iqn]->cd(2);
        cSummary[icent][iEP][iqn]->cd(2)->SetLogy();
        hObs[icent][iEP][iqn]->Draw();
        for(int i = 0; i < NITER; i++) hRefold[icent][iEP][iqn][i]->Draw("same");
	hObs[icent][iEP][iqn]->Draw("same");
	latex.DrawLatex(0.18, 0.88, "Refolding");

	cSummary[icent][iEP][iqn]->cd(3);
        cSummary[icent][iEP][iqn]->cd(3)->SetLogy();
        for(int i = 0; i < NITER; i++){
          if(i==0) hRefoldRatio[icent][iEP][iqn][i]->Draw();
          else     hRefoldRatio[icent][iEP][iqn][i]->Draw("same");
        }

	cSummary[icent][iEP][iqn]->cd(4);
	cSummary[icent][iEP][iqn]->cd(4)->SetLogx();
	grUnfoldMean[icent][iEP][iqn]->Draw("alp");
	grUnfoldStDev[icent][iEP][iqn]->Draw("lpsame");
	grUnfoldRMS[icent][iEP][iqn]->Draw("lpsame");

	grUnfoldMean[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	grUnfoldRMS[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	grUnfoldStDev[icent][iEP][iqn]->GetYaxis()->SetRangeUser(momentMin, momentMax);
	legMoment->Draw("same");

	cSummary[icent][iEP][iqn]->cd(5);
	cSummary[icent][iEP][iqn]->cd(5)->SetLogx();
	grChi2NDF_Refold[icent][iEP][iqn]->Draw("alp");

	grChi2NDF_Refold[icent][iEP][iqn]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
	cSummary[icent][iEP][iqn]->cd(5)->SetLogy();
	line->Draw("same");

	cSummary[icent][iEP][iqn]->cd(6);
	cSummary[icent][iEP][iqn]->cd(6)->SetLogy();
	leg->Draw();
	int qmin = 100*(iqn)/NQN;
	int qmax = 100*(iqn+1)/NQN;
	latex.DrawLatex(0.18, 0.92, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%"));
	latex.DrawLatex(0.18, 0.85, Form("|#eta| < %.1f", tkEta) );
	latex.DrawLatex(0.18, 0.78,  Form("q_{%i} %i - %i%s", norder_, qmin, qmax, "%"));

	cSummary[icent][iEP][iqn]->Update();
	if( SaveSumPl ) cSummary[icent][iEP][iqn]->SaveAs( Form("plots/unfolding/summary_%s_c%i_qbin%i.pdf", EPSymmNames[iEP].data(), icent, iqn) );

      } //-- end Qn loop
    } //-- End EP loop
  } //-- End Cent loop

}






#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void unfoldAssess(){

  int centbin     = 4;
  int EPSymmBin   = 0;
  int norder_     = 4;
  int QnBinOrder_ = 2;

  bool dosys     = 0;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double momentMin = 0.00;
  double momentMax = 0.29;
  double chi2Min   = 0.1;
  double chi2Max   = 1000000.;

  TLatex latex;
  TLatex latex2;

  //-- Analyzer Output
  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Unfolding output
  TFile * fUnfold;
  TH2D * hResp[NCENT];
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];
  TH1D * hRefoldRatio[NCENT][NITER];

  //-- Save final unfolded distributions (cleaned up) to a ROOT file
  TFile * fSave;
  TH1D * hUnfoldFinal[NCENT];

  //-- Unfolding moments
  double unfold_mean[NCENT][NITER];
  double unfold_stDev[NCENT][NITER];
  double unfold_RMS[NCENT][NITER];

  double unfold_meane[NCENT][NITER];
  double unfold_stDeve[NCENT][NITER];
  double unfold_RMSe[NCENT][NITER];

  //-- Moment vs iter graphs
  TMultiGraph * grMoments[NCENT];
  TGraphErrors * grUnfoldMean[NCENT];
  TGraphErrors * grUnfoldStDev[NCENT];
  TGraphErrors * grUnfoldRMS[NCENT];

  //-- Chi Squares for iteration cutoffs
  double chi2NDF_Iteration[NCENT][NITER];
  double chi2NDF_Refold[NCENT][NITER];

  //-- Chi square plots
  TGraphErrors * grChi2NDF_Iteration[NCENT];
  TGraphErrors * grChi2NDF_Refold[NCENT];

  bool iterCutoff[NCENT];
  int iterCut[NCENT];

  //-- Canvases
  TCanvas * cSummary[NCENT];
  TCanvas * cUnfoldDistsBig;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();

  //-- Get the Analyzer output file
  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );

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
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/gaussResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/gaussResp/data%i_dosys.root", norder_) );
  }
  if( studTResp ){
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/studTResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/studTResp/data%i_dosys.root", norder_) );
  }
  if( dataResp ){
    if( !dosys ) fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );
    else         fUnfold = new TFile( Form("UnfoldResults/dataResp/data%i_dosys.root", norder_) );
  }

  //-- Initialize the file that will hold the final unfolded distns
  fSave = new TFile("finalUnfoldDistns.root", "recreate");

  //-- Start looping over the data...
  for(int icent = 0; icent < NCENT; icent++){

    std::cout<< "--------------- CentBin = "<< icent << " ---------------" <<std::endl;

    //-- Get the VN observed histogram
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs[icent]->SetMaximum( 10.*hObs[icent]->GetMaximum() );

    iterCutoff[icent] = 0;
    for(int i = 0; i < NITER; i++){

      //-- Get the Response histogram
      hResp[icent] = (TH2D*) fUnfold->Get( Form("hresp_c%i", icent) );

      //-- Get the unfolded histograms
      hUnfold[icent][i] = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold[icent][i]->SetLineColor(col[i]);
      hUnfold[icent][i]->SetMarkerColor(col[i]);

      //-- Get the Refolded histograms
      hRefold[icent][i] = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      hRefold[icent][i]->SetLineWidth(2);
      hRefold[icent][i]->SetLineColor(col[i]);
      hRefold[icent][i]->SetMarkerColor(col[i]);

      //-- Make the refold/observed histograms
      hRefoldRatio[icent][i] = (TH1D*) hRefold[icent][i]->Clone( Form("hRefoldRatio%i_c%i", iter[i], icent) );
      hRefoldRatio[icent][i]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
      hRefoldRatio[icent][i]->GetYaxis()->SetTitle( Form("v_{%i}^{Refold} / v_{%i}^{Obs}", norder_, norder_) );
      hRefoldRatio[icent][i]->Divide( hObs[icent] );

      //-- Calculate the moments of each unfolded histogram
      double uMu   = hUnfold[icent][i]->GetMean();
      double uMue  = hUnfold[icent][i]->GetMeanError();
      double uSig  = hUnfold[icent][i]->GetRMS();
      double uSige = hUnfold[icent][i]->GetRMSError();
      double uRMS  = TMath::Sqrt( pow(uMu, 2) + pow(uSig, 2) );
      double uRMSe = TMath::Sqrt( ( pow(uMu * uMue, 2) + pow(uSig * uSige, 2) ) / ( pow(uMu,2) + pow(uSig,2) ) );

      unfold_mean[icent][i]   = uMu;
      unfold_RMS[icent][i]    = uRMS;
      unfold_stDev[icent][i]  = uSig;
      
      unfold_meane[icent][i]  = uMue;
      unfold_RMSe[icent][i]   = uRMSe;
      unfold_stDeve[icent][i] = uSige;

      //-- Chi squares
      if( i == 0 ) chi2NDF_Iteration[icent][i] = hUnfold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      else         chi2NDF_Iteration[icent][i] = hUnfold[icent][i]->Chi2Test(hUnfold[icent][i-1], "CHI2/NDF");;
      chi2NDF_Refold[icent][i] = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF"); 

      if( chi2NDF_Refold[icent][i] <= 1.2 && !iterCutoff[icent] ){
	iterCutoff[icent] = true;
	iterCut[icent]    = i;
      }
      if( i == NITER-1 && !iterCutoff[icent] ){
	iterCutoff[icent] = true;
	iterCut[icent]    = i;
      }

    } //-- End unfold iteration loop

    //-- Set up all the wonderful graphs
    grUnfoldMean[icent]  = new TGraphErrors(NITER, diter, unfold_mean[icent], iterErr, unfold_meane[icent] );
    grUnfoldMean[icent]->GetXaxis()->SetTitle("Iteration");
    grUnfoldMean[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    grUnfoldMean[icent]->GetYaxis()->SetTitle("Moment");
    grUnfoldMean[icent]->SetLineColor(8);
    grUnfoldMean[icent]->SetMarkerColor(8);
    grUnfoldMean[icent]->SetMarkerStyle(20);

    grUnfoldStDev[icent] = new TGraphErrors(NITER, diter, unfold_stDev[icent], iterErr, unfold_stDeve[icent] ); //46
    grUnfoldStDev[icent]->GetXaxis()->SetTitle("Iteration");
    grUnfoldStDev[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    grUnfoldStDev[icent]->GetYaxis()->SetTitle("Moment");
    grUnfoldStDev[icent]->SetLineColor(46);
    grUnfoldStDev[icent]->SetMarkerColor(46);
    grUnfoldStDev[icent]->SetMarkerStyle(20);

    grUnfoldRMS[icent]   = new TGraphErrors(NITER, diter, unfold_RMS[icent], iterErr, unfold_RMSe[icent] ); //9
    grUnfoldRMS[icent]->GetXaxis()->SetTitle("Iteration");
    grUnfoldRMS[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    grUnfoldRMS[icent]->GetYaxis()->SetTitle("Moment");
    grUnfoldRMS[icent]->SetLineColor(9);
    grUnfoldRMS[icent]->SetMarkerColor(9);
    grUnfoldRMS[icent]->SetMarkerStyle(20);

    grChi2NDF_Iteration[icent] = new TGraphErrors(NITER, diter, chi2NDF_Iteration[icent], iterErr, iterErr);
    grChi2NDF_Iteration[icent]->GetXaxis()->SetTitle("Iteration");
    grChi2NDF_Iteration[icent]->GetYaxis()->SetTitle("#chi^{2}/NDF");
    grChi2NDF_Iteration[icent]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
    grChi2NDF_Iteration[icent]->SetLineColor(2);
    grChi2NDF_Iteration[icent]->SetMarkerColor(2);
    grChi2NDF_Iteration[icent]->SetMarkerStyle(20);

    grChi2NDF_Refold[icent]    = new TGraphErrors(NITER, diter, chi2NDF_Refold[icent],    iterErr, iterErr);
    grChi2NDF_Refold[icent]->GetXaxis()->SetTitle("Iteration");
    grChi2NDF_Refold[icent]->GetYaxis()->SetTitle("#chi^{2}/NDF");
    grChi2NDF_Refold[icent]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
    grChi2NDF_Refold[icent]->SetLineColor(4);
    grChi2NDF_Refold[icent]->SetMarkerColor(4);
    grChi2NDF_Refold[icent]->SetMarkerStyle(20);

  } //-- End cent loop

  //-- Doctor the final unfolded distns by throwing out bin contents less than one
  for(int icent = 0; icent < NCENT; icent++){

    int i = iterCut[icent];
    hUnfoldFinal[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hUnfoldFinal_c%i", icent) );
    FixUnfold( hUnfoldFinal[icent] );

    fSave->cd();
    hUnfoldFinal[icent]->Write();

  }

  //-- DRAW!!!! 
  TLegend * leg = new TLegend(0.1283, 0.1440, 0.8108, 0.7539);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hObs[0], Form("v_{%i}^{obs}", norder_), "lp");
  for(int i = 0; i < NITER; i++) leg->AddEntry(hUnfold[0][i], Form("D'Agostini iter = %i", iter[i]), "lp");

  TLegend * legMoment = new TLegend(0.6025, 0.7404, 0.9742, 0.8824);
  legMoment->SetFillStyle(0);
  legMoment->SetBorderSize(0);
  legMoment->SetNColumns(2);
  legMoment->AddEntry(grUnfoldRMS[0],   "#sqrt{#LTv_{2}^{2}#GT}", "lp");
  legMoment->AddEntry(grUnfoldMean[0],  "#LTv_{2}#GT", "lp");
  legMoment->AddEntry(grUnfoldStDev[0], "#sigma_{v_{2}}", "lp");

  TLegend * legChi2 = new TLegend(0.2179, 0.8205, 0.9164, 0.9698);
  legChi2->SetFillStyle(0);
  legChi2->SetBorderSize(0);
  legChi2->SetNColumns(2);
  legChi2->AddEntry(grChi2NDF_Iteration[0], "Iteration #chi^{2}/NDF", "lp");
  legChi2->AddEntry(grChi2NDF_Refold[0],    "Refold #chi^{2}/NDF","lp");

  TLine * line = new TLine(iter[0], 1., iter[NITER-1], 1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  for(int icent = 0; icent < NCENT; icent++){

    cSummary[icent] = new TCanvas( Form("summary_c%i", icent), Form("summary_c%i", icent), 1000, 600);
    cSummary[icent]->Divide(3,2);

    cSummary[icent]->cd(1);
    cSummary[icent]->cd(1)->SetLogy();
    hObs[icent]->Draw();
    for(int i = 0; i < NITER; i++) hUnfold[icent][i]->Draw("same");
    latex.DrawLatex(0.18, 0.88, "Unfolding");

    cSummary[icent]->cd(2);
    cSummary[icent]->cd(2)->SetLogy();
    hObs[icent]->Draw();
    for(int i = 0; i < NITER; i++) hRefold[icent][i]->Draw("same");
    latex.DrawLatex(0.18, 0.88, "Refolding");

    cSummary[icent]->cd(3);
    cSummary[icent]->cd(3)->SetLogy();
    for(int i = 0; i < NITER; i++){
      if(i==0) hRefoldRatio[icent][i]->Draw();
      else     hRefoldRatio[icent][i]->Draw("same");
    }

    cSummary[icent]->cd(4);
    cSummary[icent]->cd(4)->SetLogx();
    grUnfoldMean[icent]->Draw("alp");
    grUnfoldStDev[icent]->Draw("lpsame");
    grUnfoldRMS[icent]->Draw("lpsame");

    grUnfoldMean[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    grUnfoldRMS[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    grUnfoldStDev[icent]->GetYaxis()->SetRangeUser(momentMin, momentMax);
    legMoment->Draw("same");

    cSummary[icent]->cd(5);
    cSummary[icent]->cd(5)->SetLogx();
    grChi2NDF_Iteration[icent]->Draw("alp");
    grChi2NDF_Refold[icent]->Draw("lpsame");

    grChi2NDF_Iteration[icent]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
    grChi2NDF_Refold[icent]->GetYaxis()->SetRangeUser(chi2Min, chi2Max);
    cSummary[icent]->cd(5)->SetLogy();
    legChi2->Draw("same");
    line->Draw("same");

    cSummary[icent]->cd(6);
    cSummary[icent]->cd(6)->SetLogy();
    latex.DrawLatex(0.18, 0.8, Form("Centrality %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    leg->Draw();

    cSummary[icent]->Update();
    cSummary[icent]->SaveAs( Form("plots/unfolding/summary_c%i.pdf", icent) );

  } //-- End Cent loop

  cUnfoldDistsBig = new TCanvas("cUnfoldDistsBig", "cUnfoldDistsBig", 2000, 1500);
  cUnfoldDistsBig->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    cUnfoldDistsBig->cd(icent+1);
    cUnfoldDistsBig->cd(icent+1)->SetLogy();
    int i = iterCut[icent];
    hUnfold[icent][i]->Draw();
    latex.DrawLatex(0.65, 0.88, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
  }

  cUnfoldDistsBig->SaveAs("plots/unfolding/finalUnfDist_big.pdf");


}


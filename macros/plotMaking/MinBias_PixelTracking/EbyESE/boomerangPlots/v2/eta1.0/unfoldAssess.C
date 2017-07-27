#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"

using namespace hi;
using namespace ebyese;

void unfoldAssess(){

  int centbin     = 4;
  int EPSymmBin   = 0;
  int norder_     = 2;
  double tkEta    = 1.0;
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
  TH2D * h2Vn2D0v1[NCENT];
  TH1D * h2Vn2D0v1_X[NCENT];
  TH1D * h2Vn2D0v1_Y[NCENT];
  TH1D * hMult[NCENT];
  TF1 * fGausX[NCENT];
  TF1 * fGausY[NCENT];
  TF1 * fStudTX[NCENT];
  TF1 * fStudTY[NCENT];

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
  latex.SetTextFont(43);
  latex.SetTextSize(20);
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(23);

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

      //-- 2SE Hists
      h2Vn2D0v1[icent]   = (TH2D*) fAna->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );
      h2Vn2D0v1_X[icent] = (TH1D*) h2Vn2D0v1[icent]->ProjectionX( Form("h2Vn2D0v1_X_c%i", icent) );
      h2Vn2D0v1_Y[icent] = (TH1D*) h2Vn2D0v1[icent]->ProjectionY( Form("h2Vn2D0v1_Y_c%i", icent) );
      hMult[icent]       = (TH1D*) fAna->Get( Form("qwebye/Mult_c%i", icent) );

      fGausX[icent]  = new TF1(Form("fGausX_c%i", icent), "gaus", -vnMax[norder_], vnMax[norder_] );
      fGausX[icent]->SetLineColor(2);
      fGausX[icent]->SetLineWidth(2);

      fGausY[icent]  = new TF1(Form("fGausY_c%i", icent), "gaus", -vnMax[norder_], vnMax[norder_] );;
      fGausY[icent]->SetLineColor(2);
      fGausY[icent]->SetLineWidth(2);

      double ndf = hMult[icent]->GetMean() - 1;
      fStudTX[icent] = new TF1(Form("fStudTX_c%i", icent), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
      fStudTX[icent]->SetParName(0,"Norm");
      fStudTX[icent]->SetParName(1,"Mean");
      fStudTX[icent]->SetParName(2,"Sigma");
      fStudTX[icent]->SetParName(3,"nu");
      fStudTX[icent]->SetParameters(h2Vn2D0v1_X[icent]->GetMaximum(), 0, ((ndf-2.)/ndf)*h2Vn2D0v1_X[icent]->GetRMS(),ndf);
      fStudTX[icent]->FixParameter(3,ndf);
      fStudTX[icent]->SetLineColor(kCyan);
      fStudTX[icent]->SetLineWidth(2);

      fStudTY[icent] = new TF1(Form("fStudTY_c%i", icent), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
      fStudTY[icent]->SetParName(0,"Norm");
      fStudTY[icent]->SetParName(1,"Mean");
      fStudTY[icent]->SetParName(2,"Sigma");
      fStudTY[icent]->SetParName(3,"nu");
      fStudTY[icent]->SetParameters(h2Vn2D0v1_Y[icent]->GetMaximum(), 0, ((ndf-2.)/ndf)*h2Vn2D0v1_Y[icent]->GetRMS(),ndf);
      fStudTY[icent]->FixParameter(3,ndf);
      fStudTY[icent]->SetLineColor(kCyan);
      fStudTY[icent]->SetLineWidth(2);

      double sigmaFit = 2.;
      h2Vn2D0v1_X[icent]->Fit(Form("fGausX_c%i", icent), "0", "", h2Vn2D0v1_X[icent]->GetMean() - sigmaFit * h2Vn2D0v1_X[icent]->GetRMS(), h2Vn2D0v1_X[icent]->GetMean() + sigmaFit * h2Vn2D0v1_X[icent]->GetRMS());
      h2Vn2D0v1_Y[icent]->Fit(Form("fGausY_c%i", icent), "0", "", h2Vn2D0v1_Y[icent]->GetMean() - sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS(), h2Vn2D0v1_Y[icent]->GetMean() + sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS());

      if(icent == 11){
	h2Vn2D0v1_X[icent]->Fit(Form("fStudTX_c%i", icent), "0", "", h2Vn2D0v1_X[icent]->GetMean() - sigmaFit * h2Vn2D0v1_X[icent]->GetRMS(), h2Vn2D0v1_X[icent]->GetMean() + sigmaFit * h2Vn2D0v1_X[icent]->GetRMS());
	h2Vn2D0v1_Y[icent]->Fit(Form("fStudTY_c%i", icent), "0", "", h2Vn2D0v1_Y[icent]->GetMean() - sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS(), h2Vn2D0v1_Y[icent]->GetMean() + sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS());
      }
      else{
	h2Vn2D0v1_X[icent]->Fit(Form("fStudTX_c%i", icent), "L0", "", h2Vn2D0v1_X[icent]->GetMean() - sigmaFit * h2Vn2D0v1_X[icent]->GetRMS(), h2Vn2D0v1_X[icent]->GetMean() + sigmaFit * h2Vn2D0v1_X[icent]->GetRMS());
        h2Vn2D0v1_Y[icent]->Fit(Form("fStudTY_c%i", icent), "L0", "", h2Vn2D0v1_Y[icent]->GetMean() - sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS(), h2Vn2D0v1_Y[icent]->GetMean() + sigmaFit * h2Vn2D0v1_Y[icent]->GetRMS());
      }
      h2Vn2D0v1_X[icent]->SetStats(0);
      h2Vn2D0v1_X[icent]->GetXaxis()->SetTitleOffset(1.2);
      h2Vn2D0v1_X[icent]->GetXaxis()->SetLabelSize(0.06);
      h2Vn2D0v1_X[icent]->GetXaxis()->SetRange(h2Vn2D0v1_X[icent]->FindBin(-0.59), h2Vn2D0v1_X[icent]->FindBin(0.59));
      h2Vn2D0v1_X[icent]->GetYaxis()->SetTitle("Events");
      h2Vn2D0v1_X[icent]->GetYaxis()->SetTitleOffset(1.5);
      h2Vn2D0v1_X[icent]->GetYaxis()->SetLabelSize(0.06);

      h2Vn2D0v1_Y[icent]->SetStats(0);
      h2Vn2D0v1_Y[icent]->GetXaxis()->SetTitleOffset(1.2);
      h2Vn2D0v1_Y[icent]->GetXaxis()->SetLabelSize(0.06);
      h2Vn2D0v1_Y[icent]->GetXaxis()->SetRange(h2Vn2D0v1_Y[icent]->FindBin(-0.59), h2Vn2D0v1_Y[icent]->FindBin(0.59));
      h2Vn2D0v1_Y[icent]->GetYaxis()->SetTitle("Events");
      h2Vn2D0v1_Y[icent]->GetYaxis()->SetTitleOffset(1.5);
      h2Vn2D0v1_Y[icent]->GetYaxis()->SetLabelSize(0.06);

      h2Vn2D0v1[icent]->SetStats(0);
      h2Vn2D0v1[icent]->GetXaxis()->SetTitleOffset(1.2);
      h2Vn2D0v1[icent]->GetXaxis()->SetLabelSize(0.06);
      h2Vn2D0v1[icent]->GetXaxis()->SetRange(h2Vn2D0v1[icent]->GetXaxis()->FindBin(-0.59), h2Vn2D0v1[icent]->GetXaxis()->FindBin(0.59));
      h2Vn2D0v1[icent]->GetYaxis()->SetTitleOffset(1.5);
      h2Vn2D0v1[icent]->GetYaxis()->SetLabelSize(0.06);
      h2Vn2D0v1[icent]->GetYaxis()->SetRange(h2Vn2D0v1[icent]->GetYaxis()->FindBin(-0.59), h2Vn2D0v1[icent]->GetYaxis()->FindBin(0.59));

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
      if( i == 0 ) chi2NDF_Iteration[icent][i] = hUnfold[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF");
      else         chi2NDF_Iteration[icent][i] = hUnfold[icent][i]->Chi2Test(hUnfold[icent][i-1], "UWCHI2/NDF");;
      chi2NDF_Refold[icent][i] = hRefold[icent][i]->Chi2Test(hObs[icent], "UWCHI2/NDF"); 

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

  int c;
  TLegend * leg2SE = new TLegend(0.2, 0.75, 0.6, 0.88);
  legInit( leg2SE );
  leg2SE->AddEntry(fGausX[0], "Gaussian", "l");
  leg2SE->AddEntry(fStudTX[0], "Student's T", "l");


  //-- PAS Fig. 1
  TCanvas * c2SE = new TCanvas("c2SE", "c2SE", 1500, 1000);
  c2SE->Divide(3,2);
  c2SE->SetLeftMargin(0.4);

  c = 1;
  c2SE->cd(2);
  c2SE->cd(2)->SetLogy();
  c2SE->cd(2)->SetLeftMargin(0.19);
  c2SE->cd(2)->SetRightMargin(0.07);
  c2SE->cd(2)->SetTopMargin(0.1);
  h2Vn2D0v1_X[c]->SetMaximum( 20.*h2Vn2D0v1_X[c]->GetMaximum() );
  h2Vn2D0v1_X[c]->Draw();
  fGausX[c]->Draw("Same");
  fStudTX[c]->Draw("Same");
  leg2SE->Draw("same");
  latex.DrawLatex(0.63, 0.83, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  c2SE->cd(3);
  c2SE->cd(3)->SetLogy();
  c2SE->cd(3)->SetLeftMargin(0.19);
  c2SE->cd(3)->SetRightMargin(0.07);
  c2SE->cd(3)->SetTopMargin(0.1);
  h2Vn2D0v1_Y[c]->SetMaximum( 20.*h2Vn2D0v1_Y[c]->GetMaximum() );
  h2Vn2D0v1_Y[c]->Draw();
  fGausY[c]->Draw("Same");
  fStudTY[c]->Draw("Same");
  latex.DrawLatex(0.63, 0.83, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  c = 11;
  c2SE->cd(5);
  c2SE->cd(5)->SetLogy();
  c2SE->cd(5)->SetLeftMargin(0.19);
  c2SE->cd(5)->SetRightMargin(0.07);
  c2SE->cd(5)->SetTopMargin(0.1);
  h2Vn2D0v1_X[c]->SetMaximum( 20.*h2Vn2D0v1_X[c]->GetMaximum() );
  h2Vn2D0v1_X[c]->Draw();
  fGausX[c]->Draw("Same");
  fStudTX[c]->Draw("Same");
  leg2SE->Draw("same");
  latex.DrawLatex(0.63, 0.83, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  c2SE->cd(6);
  c2SE->cd(6)->SetLogy();
  c2SE->cd(6)->SetLeftMargin(0.19);
  c2SE->cd(6)->SetRightMargin(0.07);
  c2SE->cd(6)->SetTopMargin(0.1);
  h2Vn2D0v1_Y[c]->SetMaximum( 20.*h2Vn2D0v1_Y[c]->GetMaximum() );
  h2Vn2D0v1_Y[c]->Draw();
  fGausY[c]->Draw("Same");
  fStudTY[c]->Draw("Same");
  latex.DrawLatex(0.63, 0.83, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  leg2SE->SetTextFont(43);
  leg2SE->SetTextSize(23);

  c = 1;
  c2SE->cd(1);
  c2SE->cd(1)->SetLogz();
  c2SE->cd(1)->SetLeftMargin(0.19);
  c2SE->cd(1)->SetRightMargin(0.07);
  c2SE->cd(1)->SetTopMargin(0.1);
  h2Vn2D0v1[c]->Draw();
  c2SE->Update();
  TPaletteAxis * p1 = (TPaletteAxis*) h2Vn2D0v1[c]->GetListOfFunctions()->FindObject("palette");
  p1->SetX1NDC(0.75);
  p1->SetY1NDC(0.25);
  p1->SetX2NDC(0.8);
  p1->SetY2NDC(0.55);
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.23, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.23, 0.76, Form("|#eta| < %.1f", tkEta));
  latex2.DrawLatex(0.23, 0.21, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  c = 11;
  c2SE->cd(4);
  c2SE->cd(4)->SetLogz();
  c2SE->cd(4)->SetLeftMargin(0.19);
  c2SE->cd(4)->SetRightMargin(0.07);
  c2SE->cd(4)->SetTopMargin(0.1);
  h2Vn2D0v1[c]->Draw();
  c2SE->Update();
  TPaletteAxis * p2 = (TPaletteAxis*) h2Vn2D0v1[c]->GetListOfFunctions()->FindObject("palette");
  p2->SetX1NDC(0.75);
  p2->SetY1NDC(0.25);
  p2->SetX2NDC(0.8);
  p2->SetY2NDC(0.55);
  latex.DrawLatex(0.20, 0.915, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.66, 0.915, "PbPb 5.02 TeV");
  latex2.DrawLatex(0.23, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex2.DrawLatex(0.23, 0.76, Form("|#eta| < %.1f", tkEta));
  latex2.DrawLatex(0.23,0.21, Form("Cent. %i - %i%s", cent_min[c], cent_max[c], "%") );

  c2SE->Update();
  c2SE->SaveAs( "plots/unfolding/c2SE.pdf" );

  for(int icent = 0; icent < NCENT; icent++){

    double chi2GX = fGausX[icent]->GetChisquare() / (double)fGausX[icent]->GetNDF();
    double chi2GY = fGausY[icent]->GetChisquare() / (double)fGausY[icent]->GetNDF();
    double chi2SX = fStudTX[icent]->GetChisquare() / (double)fStudTX[icent]->GetNDF();
    double chi2SY = fStudTY[icent]->GetChisquare() / (double)fStudTY[icent]->GetNDF();

    std::cout << Form("c%i\t", icent) 
	      << Form("Chi2/NDF GX = %.1f\t", chi2GX) 
	      << Form("Chi2/NDF GY = %.1f\t", chi2GY) 
	      << Form("Chi2/NDF SX = %.1f\t", chi2SX) 
	      << Form("Chi2/NDF SY = %.1f", chi2SY) 
	      << std::endl;

  }

}


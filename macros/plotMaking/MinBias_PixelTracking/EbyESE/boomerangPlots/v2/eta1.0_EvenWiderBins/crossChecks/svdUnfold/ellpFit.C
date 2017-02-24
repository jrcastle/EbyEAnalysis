#include "TExec.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

double pEllP(double * x, double * par){

  int nbin = 180;
  //int nbin = 360;

  //-- [0] = e0
  //-- [1] = alpha
  //-- [2] = kn
  //-- [3] = Scale 
  //-- xx  = vn

  double e0    = par[0];
  double alpha = par[1];
  double kn    = par[2];
  double scale = par[3];
  double eccn  = x[0] / kn;


  double pi = TMath::Pi();
  double p1 = (scale * 2. * alpha * eccn / pi / kn) * pow( 1 - e0*e0, alpha + 0.5);
  //double p1 = (2. * alpha * eccn / pi / kn) * pow( 1 - e0*e0, alpha + 0.5);

  //-- [0] = xx
  //-- [1] = e0
  //-- [2] = alpha
  //-- [3] = kn
  //TF1 f("f", "pow(1 - [0]*[0]/[3]/[3], [2] - 1) / pow(1 - [1]*[0]*TMath::Cos(x)/[3]/[3], 2*[2]+1)", 0., 2*pi);
  //f.SetParameters(xx, e0, alpha, kn);
  //double integ = f.Integral(0, pi);
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 0; i <= nbin; i++){
    double phi = 0. + dphi*(double)(i);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  //integ += 0.5 * dphi * ( pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(0), -2.*alpha-1 ) + pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(pi), -2.*alpha-1 ) );
  double ellp = p1 * integ;

  return ellp;

}

void ellpFit(){

  const int norder_ = 2;
  bool TRUNC_ = 0;

  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.43,  0.41, 0.40, 0.39, 0.37, 0.35, 0.32, 0.28, 0.23};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 102.0, 70.0, 59.0, 47.0, 34.0, 25.0, 18.0, 11.0, 7.0};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.15,  0.18, 0.20, 0.22, 0.24, 0.25, 0.27, 0.30, 0.35};
  double vnmax[NCENT]   = {0.3,   0.3,   0.3,   0.3,   0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.25,  0.27};

  //-- Ana
  TFile * fAna;
  TH1D * hObs[NCENT];

  //-- Reg
  TFile * fDAG;
  TH1D * hUDAG[NCENT][NITER];
  TH1D * hRDAG[NCENT][NITER];
  TH1D * hUnfoldDAG[NCENT];

  TF1 * fEllP_DAG[NCENT];
  double fitE0_DAG[NCENT];
  double fitAlpha_DAG[NCENT];
  double fitKn_DAG[NCENT];
  double fitE0_DAG_err[NCENT];
  double fitAlpha_DAG_err[NCENT];
  double fitKn_DAG_err[NCENT];

  TGaxis * axEcc_DAG[NCENT];
  TLine * asmpt_DAG[NCENT];

  //-- Resp
  TFile * fSVD;
  TH1D * hUSVD[NCENT][NKREG];
  TH1D * hUnfoldSVD[NCENT];
  TH1D * hKreg;

  TF1 * fEllP_SVD[NCENT];
  double fitE0_SVD[NCENT];
  double fitAlpha_SVD[NCENT];
  double fitKn_SVD[NCENT];
  double fitE0_SVD_err[NCENT];
  double fitAlpha_SVD_err[NCENT];
  double fitKn_SVD_err[NCENT];

  TGaxis * axEcc_SVD[NCENT];
  TLine * asmpt_SVD[NCENT];

  TLatex latex;

  //
  // MAIN
  // 
  setTDRStyle();
  latex.SetNDC();

  //-- Get Files 
  fAna = new TFile("../../AnalyzerResults/CastleEbyE.root");
  fDAG = new TFile("../../UnfoldResults/dataResp/data2.root");
  fSVD = new TFile("../../UnfoldResults/dataResp/data2_svd.root");

  hKreg = (TH1D*) fSVD->Get( "hKreg" );
  hKreg->SetBinContent(6, hKreg->GetBinContent(6)+1);
  hKreg->SetBinContent(7, hKreg->GetBinContent(7)+1);

  //-- Canvases
  TCanvas * cFits_DAG = new TCanvas("cFits_DAG", "cFits_DAG", 1500, 1500);
  cFits_DAG->Divide(3,3);

  TCanvas * cFits_SVD = new TCanvas("cFits_SVD", "cFits_SVD", 1500, 1500);
  cFits_SVD->Divide(3,3);

  //-- Cent loop
  for(int icent = 0; icent < NCENT; icent++){
    if(icent < 3) continue;

    double e0, alpha, kn;
    double xmin = 0;
    double xmax = 0.3;
    double ymax;

    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    //-- ------------ DAG ------------
    std::cout << Form("============ DAG c%i ============", icent) << std::endl;

    int finalIter = 0;
    for(int i = 0; i < NITER; i++){
      hUDAG[icent][i] = (TH1D*) fDAG->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRDAG[icent][i] = (TH1D*) fDAG->Get( Form("hrefold%i_c%i", iter[i], icent) );

      double chi2 = hRDAG[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      if( chi2 < 1.2 ){
	finalIter = i;
	break;
      }
      if(i == NITER){
	finalIter = i;
	break;
      }

    }
    hUnfoldDAG[icent] = (TH1D*) hUDAG[icent][finalIter]->Clone( Form("hFinalUnfold_DAG_c%i", icent) );
    hUnfoldDAG[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldDAG[icent]->GetXaxis()->SetRange(1, hUnfoldDAG[icent]->FindBin(0.6));
    hUnfoldDAG[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    hUnfoldDAG[icent]->Scale(1./hUnfoldDAG[icent]->Integral());
    hUnfoldDAG[icent]->SetStats(0);

    fEllP_DAG[icent] = new TF1(Form("fEllP_DAG_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP_DAG[icent]->SetLineColor(2);
    fEllP_DAG[icent]->SetLineWidth(2);
    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP_DAG[icent]->SetParLimits(0, 0., 1.);
    fEllP_DAG[icent]->SetParLimits(2, 0., 1.);
    fEllP_DAG[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], 1.0);
    fEllP_DAG[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");

    hUnfoldDAG[icent]->Fit( Form("fEllP_DAG_c%i", icent), "R0", "", 0.0, vnmax[icent]);
    fitE0_DAG[icent]        = fEllP_DAG[icent]->GetParameter(0);
    fitAlpha_DAG[icent]     = fEllP_DAG[icent]->GetParameter(1);
    fitKn_DAG[icent]        = fEllP_DAG[icent]->GetParameter(2);
    fitE0_DAG_err[icent]    = fEllP_DAG[icent]->GetParError(0);
    fitAlpha_DAG_err[icent] = fEllP_DAG[icent]->GetParError(1);
    fitKn_DAG_err[icent]    = fEllP_DAG[icent]->GetParError(2);

    e0    = fEllP_DAG[icent]->GetParameter(0);
    alpha = fEllP_DAG[icent]->GetParameter(1);
    kn    = fEllP_DAG[icent]->GetParameter(2);
    ymax  = 1.9*hUnfoldDAG[icent]->GetMaximum();
    axEcc_DAG[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_DAG[icent]->SetLabelSize(0.055);
    axEcc_DAG[icent]->SetTitleSize(0.055);
    axEcc_DAG[icent]->SetTitle("#epsilon_{2}");

    asmpt_DAG[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_DAG[icent]->SetLineStyle(2);

    cFits_DAG->cd(icent-2);
    cFits_DAG->cd(icent-2)->SetLogy();
    cFits_DAG->cd(icent-2)->SetTopMargin(0.15);
    cFits_DAG->cd(icent-2)->SetTickx(0);
    hUnfoldDAG[icent]->Draw();
    fEllP_DAG[icent]->Draw("same");
    axEcc_DAG[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_DAG[icent]->Draw("same");

    //-- ------------ SVD ------------
    std::cout << Form("============ SVD c%i ============", icent) << std::endl;

    int kreg = hKreg->GetBinContent(icent+1);

    hUnfoldSVD[icent] = (TH1D*) fSVD->Get( Form("hrecokreg%i_c%i", kreg, icent) );
    hUnfoldSVD[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldSVD[icent]->GetXaxis()->SetRange(1, hUnfoldSVD[icent]->FindBin(0.6));
    hUnfoldSVD[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    hUnfoldSVD[icent]->Scale(1./hUnfoldSVD[icent]->Integral());
    hUnfoldSVD[icent]->SetStats(0);

    fEllP_SVD[icent] = new TF1(Form("fEllP_SVD_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP_SVD[icent]->SetLineColor(2);
    fEllP_SVD[icent]->SetLineWidth(2);
    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP_SVD[icent]->SetParLimits(0, 0., 1.);
    fEllP_SVD[icent]->SetParLimits(2, 0., 1.);
    fEllP_SVD[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], 1.0);
    fEllP_SVD[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");

    hUnfoldSVD[icent]->Fit( Form("fEllP_SVD_c%i", icent), "R0", "", 0.0, vnmax[icent]);
    fitE0_SVD[icent]        = fEllP_SVD[icent]->GetParameter(0);
    fitAlpha_SVD[icent]     = fEllP_SVD[icent]->GetParameter(1);
    fitKn_SVD[icent]        = fEllP_SVD[icent]->GetParameter(2);
    fitE0_SVD_err[icent]    = fEllP_SVD[icent]->GetParError(0);
    fitAlpha_SVD_err[icent] = fEllP_SVD[icent]->GetParError(1);
    fitKn_SVD_err[icent]    = fEllP_SVD[icent]->GetParError(2);

    e0    = fEllP_SVD[icent]->GetParameter(0);
    alpha = fEllP_SVD[icent]->GetParameter(1);
    kn    = fEllP_SVD[icent]->GetParameter(2);
    ymax  = 1.9*hUnfoldSVD[icent]->GetMaximum();
    axEcc_SVD[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_SVD[icent]->SetLabelSize(0.055);
    axEcc_SVD[icent]->SetTitleSize(0.055);
    axEcc_SVD[icent]->SetTitle("#epsilon_{2}");

    asmpt_SVD[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_SVD[icent]->SetLineStyle(2);

    cFits_SVD->cd(icent-2);
    cFits_SVD->cd(icent-2)->SetLogy();
    cFits_SVD->cd(icent-2)->SetTopMargin(0.15);
    cFits_SVD->cd(icent-2)->SetTickx(0);
    hUnfoldSVD[icent]->Draw();
    fEllP_SVD[icent]->Draw("same");
    axEcc_SVD[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_SVD[icent]->Draw("same");

  } //-- End Cent loop

  //cFits_DAG->SaveAs("../plots/systematicStudies/cFits_DAG.pdf");
  //cFits_SVD->SaveAs("../plots/systematicStudies/cFits_SVD.pdf");

}

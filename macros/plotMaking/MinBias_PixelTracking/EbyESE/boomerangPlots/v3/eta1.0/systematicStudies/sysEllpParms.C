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

  int nbin = 50;
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
  for(int i = 1; i <= nbin; i++){
    double phi = (2*(double)i-1.)*pi/(2.*(double)nbin);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  //integ += 0.5 * dphi * ( pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(0), -2.*alpha-1 ) + pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(pi), -2.*alpha-1 ) );
  double ellp = p1 * integ;

  return ellp;

}

void sysEllpParms(){

  const int norder_ = 3;

  bool moscowFits = 0;
  bool fixKn       = 0;

  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.40,    0.38, 0.37, 0.37, 0.34, 0.32, 0.29, 0.28, 0.24};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 71.1,    63.0, 47.0, 41.0, 30.0, 20.0, 12.0, 11.0, 6.8};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.17,    0.19, 0.22, 0.23, 0.26, 0.29, 0.30, 0.31, 0.31};
  double vnmax[NCENT]   = {0.25,   0.25,   0.25,   0.25,   0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.23};

  if( fixKn ){
    alGuess[0]  = 0.0;    alGuess[1]  = 0.0;   alGuess[2] = 0.0;   alGuess[3] = 120.0;  alGuess[4] = 98.0;
    alGuess[5]  = 81.0;   alGuess[6]  = 70.0;  alGuess[7] = 58.0;  alGuess[8] = 51.0;   alGuess[9] = 46.0;
    alGuess[10] = 43.0;   alGuess[11] = 43.0;

    e0Guess[0]  = 0.0;    e0Guess[1]  = 0.0;   e0Guess[2] = 0.0;   e0Guess[3] = 0.14;  e0Guess[4] = 0.16;
    e0Guess[5]  = 0.17;   e0Guess[6]  = 0.18;  e0Guess[7] = 0.19;  e0Guess[8] = 0.19;  e0Guess[9] = 0.19;
    e0Guess[10] = 0.19;   e0Guess[11] = 0.19;
  }
  double sysMin = 0.8;
  double sysMax = 1.2;

  //-- Nominal
  TFile * fNominal;
  TH1D * hUnfoldNominal[NCENT];
  TH1D * hNormFactor;

  TF1 * fEllP_Nominal[NCENT];
  double fitE0_Nominal[NCENT];
  double fitAlpha_Nominal[NCENT];
  double fitKn_Nominal[NCENT];
  double fitE0_Nominal_err[NCENT];
  double fitAlpha_Nominal_err[NCENT];
  double fitKn_Nominal_err[NCENT];

  TGaxis * axEcc_Nominal[NCENT];
  TLine * asmpt_Nominal[NCENT];


  //-- Reg
  TFile * fSysReg;
  TH1D * hUnfoldSysReg[NCENT];

  TF1 * fEllP_SysReg[NCENT];
  double fitE0_SysReg[NCENT];
  double fitAlpha_SysReg[NCENT];
  double fitKn_SysReg[NCENT];
  double fitE0_SysReg_err[NCENT];
  double fitAlpha_SysReg_err[NCENT];
  double fitKn_SysReg_err[NCENT];

  TGaxis * axEcc_SysReg[NCENT];
  TLine * asmpt_SysReg[NCENT];

  double fitE0_RatioToNom_SysReg[NCENT];
  double fitAlpha_RatioToNom_SysReg[NCENT];
  double fitKn_RatioToNom_SysReg[NCENT];
  double fitE0_RatioToNom_SysReg_err[NCENT];
  double fitAlpha_RatioToNom_SysReg_err[NCENT];
  double fitKn_RatioToNom_SysReg_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_SysReg;
  TGraphErrors * grFitAlpha_RatioToNom_SysReg;
  TGraphErrors * grFitKn_RatioToNom_SysReg;

  //-- Resp
  TFile * fSysResp;
  TH1D * hUnfoldSysResp[NCENT];

  TF1 * fEllP_SysResp[NCENT];
  double fitE0_SysResp[NCENT];
  double fitAlpha_SysResp[NCENT];
  double fitKn_SysResp[NCENT];
  double fitE0_SysResp_err[NCENT];
  double fitAlpha_SysResp_err[NCENT];
  double fitKn_SysResp_err[NCENT];

  TGaxis * axEcc_SysResp[NCENT];
  TLine * asmpt_SysResp[NCENT];

  double fitE0_RatioToNom_SysResp[NCENT];
  double fitAlpha_RatioToNom_SysResp[NCENT];
  double fitKn_RatioToNom_SysResp[NCENT];
  double fitE0_RatioToNom_SysResp_err[NCENT];
  double fitAlpha_RatioToNom_SysResp_err[NCENT];
  double fitKn_RatioToNom_SysResp_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_SysResp;
  TGraphErrors * grFitAlpha_RatioToNom_SysResp;
  TGraphErrors * grFitKn_RatioToNom_SysResp;

  //-- NewCC
  TFile * fSysNewCC;
  TH1D * hUnfoldSysNewCC[NCENT];

  TF1 * fEllP_SysNewCC[NCENT];
  double fitE0_SysNewCC[NCENT];
  double fitAlpha_SysNewCC[NCENT];
  double fitKn_SysNewCC[NCENT];
  double fitE0_SysNewCC_err[NCENT];
  double fitAlpha_SysNewCC_err[NCENT];
  double fitKn_SysNewCC_err[NCENT];

  TGaxis * axEcc_SysNewCC[NCENT];
  TLine * asmpt_SysNewCC[NCENT];

  double fitE0_RatioToNom_SysNewCC[NCENT];
  double fitAlpha_RatioToNom_SysNewCC[NCENT];
  double fitKn_RatioToNom_SysNewCC[NCENT];
  double fitE0_RatioToNom_SysNewCC_err[NCENT];
  double fitAlpha_RatioToNom_SysNewCC_err[NCENT];
  double fitKn_RatioToNom_SysNewCC_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_SysNewCC;
  TGraphErrors * grFitAlpha_RatioToNom_SysNewCC;
  TGraphErrors * grFitKn_RatioToNom_SysNewCC;

  //-- TkQ
  TFile * fSysTkQ;

  TH1D * hUnfoldLooseSysTkQ[NCENT];
  TF1 * fEllP_LooseSysTkQ[NCENT];
  double fitE0_LooseSysTkQ[NCENT];
  double fitAlpha_LooseSysTkQ[NCENT];
  double fitKn_LooseSysTkQ[NCENT];
  double fitE0_LooseSysTkQ_err[NCENT];
  double fitAlpha_LooseSysTkQ_err[NCENT];
  double fitKn_LooseSysTkQ_err[NCENT];

  TGaxis * axEcc_LooseSysTkQ[NCENT];
  TLine * asmpt_LooseSysTkQ[NCENT];

  double fitE0_RatioToNom_LooseSysTkQ[NCENT];
  double fitAlpha_RatioToNom_LooseSysTkQ[NCENT];
  double fitKn_RatioToNom_LooseSysTkQ[NCENT];
  double fitE0_RatioToNom_LooseSysTkQ_err[NCENT];
  double fitAlpha_RatioToNom_LooseSysTkQ_err[NCENT];
  double fitKn_RatioToNom_LooseSysTkQ_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_LooseSysTkQ;
  TGraphErrors * grFitAlpha_RatioToNom_LooseSysTkQ;
  TGraphErrors * grFitKn_RatioToNom_LooseSysTkQ;

  TH1D * hUnfoldTightSysTkQ[NCENT];
  TF1 * fEllP_TightSysTkQ[NCENT];
  double fitE0_TightSysTkQ[NCENT];
  double fitAlpha_TightSysTkQ[NCENT];
  double fitKn_TightSysTkQ[NCENT];
  double fitE0_TightSysTkQ_err[NCENT];
  double fitAlpha_TightSysTkQ_err[NCENT];
  double fitKn_TightSysTkQ_err[NCENT];

  TGaxis * axEcc_TightSysTkQ[NCENT];
  TLine * asmpt_TightSysTkQ[NCENT];

  double fitE0_RatioToNom_TightSysTkQ[NCENT];
  double fitAlpha_RatioToNom_TightSysTkQ[NCENT];
  double fitKn_RatioToNom_TightSysTkQ[NCENT];
  double fitE0_RatioToNom_TightSysTkQ_err[NCENT];
  double fitAlpha_RatioToNom_TightSysTkQ_err[NCENT];
  double fitKn_RatioToNom_TightSysTkQ_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_TightSysTkQ;
  TGraphErrors * grFitAlpha_RatioToNom_TightSysTkQ;
  TGraphErrors * grFitKn_RatioToNom_TightSysTkQ;

  //-- Vtx
  TFile * fSysVtx;

  TH1D * hUnfoldVtx3SysVtx[NCENT];
  TF1 * fEllP_Vtx3SysVtx[NCENT];
  double fitE0_Vtx3SysVtx[NCENT];
  double fitAlpha_Vtx3SysVtx[NCENT];
  double fitKn_Vtx3SysVtx[NCENT];
  double fitE0_Vtx3SysVtx_err[NCENT];
  double fitAlpha_Vtx3SysVtx_err[NCENT];
  double fitKn_Vtx3SysVtx_err[NCENT];

  TGaxis * axEcc_Vtx3SysVtx[NCENT];
  TLine * asmpt_Vtx3SysVtx[NCENT];

  double fitE0_RatioToNom_Vtx3SysVtx[NCENT];
  double fitAlpha_RatioToNom_Vtx3SysVtx[NCENT];
  double fitKn_RatioToNom_Vtx3SysVtx[NCENT];
  double fitE0_RatioToNom_Vtx3SysVtx_err[NCENT];
  double fitAlpha_RatioToNom_Vtx3SysVtx_err[NCENT];
  double fitKn_RatioToNom_Vtx3SysVtx_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_Vtx3SysVtx;
  TGraphErrors * grFitAlpha_RatioToNom_Vtx3SysVtx;
  TGraphErrors * grFitKn_RatioToNom_Vtx3SysVtx;

  TH1D * hUnfoldVtx3_15SysVtx[NCENT];
  TF1 * fEllP_Vtx3_15SysVtx[NCENT];
  double fitE0_Vtx3_15SysVtx[NCENT];
  double fitAlpha_Vtx3_15SysVtx[NCENT];
  double fitKn_Vtx3_15SysVtx[NCENT];
  double fitE0_Vtx3_15SysVtx_err[NCENT];
  double fitAlpha_Vtx3_15SysVtx_err[NCENT];
  double fitKn_Vtx3_15SysVtx_err[NCENT];

  TGaxis * axEcc_Vtx3_15SysVtx[NCENT];
  TLine * asmpt_Vtx3_15SysVtx[NCENT];

  double fitE0_RatioToNom_Vtx3_15SysVtx[NCENT];
  double fitAlpha_RatioToNom_Vtx3_15SysVtx[NCENT];
  double fitKn_RatioToNom_Vtx3_15SysVtx[NCENT];
  double fitE0_RatioToNom_Vtx3_15SysVtx_err[NCENT];
  double fitAlpha_RatioToNom_Vtx3_15SysVtx_err[NCENT];
  double fitKn_RatioToNom_Vtx3_15SysVtx_err[NCENT];

  TGraphErrors * grFitE0_RatioToNom_Vtx3_15SysVtx;
  TGraphErrors * grFitAlpha_RatioToNom_Vtx3_15SysVtx;
  TGraphErrors * grFitKn_RatioToNom_Vtx3_15SysVtx;

  //-- Moscow Fits
  TFile * fRussiaSysReg;
  TFile * fRussiaSysResp;
  TFile * fRussiaSysNewCC;
  TFile * fRussiaSysTkQ1;
  TFile * fRussiaSysTkQ2;
  TFile * fRussiaSysVtx1;
  TFile * fRussiaSysVtx2;

  TLatex latex;

  //
  // MAIN
  // 
  setTDRStyle();
  latex.SetNDC();

  //-- Get Files 
  fNominal  = new TFile("SysUnfoldDistns_v2.root");
  fSysReg   = new TFile("chi2Cutoff/SysReg.root");
  fSysResp  = new TFile("responseElements/SysResp.root");
  fSysNewCC = new TFile("clusCompatTune/SysNewCC.root");
  fSysTkQ   = new TFile("tkQuality/SysTkQuality.root");
  fSysVtx   = new TFile("vtxCut/SysVtx.root");

  //-- Moscow Fits
  fRussiaSysReg   = new TFile("chi2Cutoff/MoscowFits_SysReg.root");
  fRussiaSysResp  = new TFile("responseElements/MoscowFits_SysResp.root");
  fRussiaSysNewCC = new TFile("clusCompatTune/MoscowFits_SysNewCC.root");
  fRussiaSysTkQ1  = new TFile("tkQuality/MoscowFits_SysTrk1.root");
  fRussiaSysTkQ2  = new TFile("tkQuality/MoscowFits_SysTrk2.root");
  fRussiaSysVtx1  = new TFile("vtxCut/SysFits_SysVtx1.root");
  fRussiaSysVtx2  = new TFile("vtxCut/SysFits_SysVtx2.root");


  //-- Normalization factor for nominal
  hNormFactor = (TH1D*) fNominal->Get( "hNormFactor");

  //-- Canvases
  TCanvas * cFits_Nominal = new TCanvas("cFits_Nominal", "cFits_Nominal", 1500, 1500);
  cFits_Nominal->Divide(3,3);

  TCanvas * cFits_SysReg = new TCanvas("cFits_SysReg", "cFits_SysReg", 1500, 1500);
  cFits_SysReg->Divide(3,3);

  TCanvas * cFits_SysResp = new TCanvas("cFits_SysResp", "cFits_SysResp", 1500, 1500);
  cFits_SysResp->Divide(3,3);

  TCanvas * cFits_SysNewCC = new TCanvas("cFits_SysNewCC", "cFits_SysNewCC", 1500, 1500);
  cFits_SysNewCC->Divide(3,3);

  TCanvas * cFits_LooseSysTkQ = new TCanvas("cFits_LooseSysTkQ", "cFits_LooseSysTkQ", 1500, 1500);
  cFits_LooseSysTkQ->Divide(3,3);

  TCanvas * cFits_TightSysTkQ = new TCanvas("cFits_TightSysTkQ", "cFits_TightSysTkQ", 1500, 1500);
  cFits_TightSysTkQ->Divide(3,3);

  TCanvas * cFits_Vtx3SysVtx = new TCanvas("cFits_Vtx3SysVtx", "cFits_Vtx3SysVtx", 1500, 1500);
  cFits_Vtx3SysVtx->Divide(3,3);

  TCanvas * cFits_Vtx3_15SysVtx = new TCanvas("cFits_Vtx3_15SysVtx", "cFits_Vtx3_15SysVtx", 1500, 1500);
  cFits_Vtx3_15SysVtx->Divide(3,3);

  //-- Cent loop
  for(int icent = 0; icent < NCENT; icent++){
    if(icent < 3) continue;

    std::cout << Form("============ SysReg c%i ============", icent) << std::endl;
    double e0, alpha, kn;
    double xmin = 0;
    double xmax = vnmax[icent];
    double ymax;

    //-- ------------ Nominal ------------
    hUnfoldNominal[icent] = (TH1D*) fNominal->Get( Form("hFinalUnfoldStat_c%i", icent) );
    hUnfoldNominal[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldNominal[icent]->GetXaxis()->SetRange(1, hUnfoldNominal[icent]->FindBin(vnmax[icent]));
    hUnfoldNominal[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    hUnfoldNominal[icent]->Scale(1./hNormFactor->GetBinContent(icent+1));
    hUnfoldNominal[icent]->SetStats(0);

    fEllP_Nominal[icent] = new TF1(Form("fEllP_Nominal_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP_Nominal[icent]->SetLineColor(2);
    fEllP_Nominal[icent]->SetLineWidth(2);
    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP_Nominal[icent]->SetParLimits(0, 0., 1.);
    fEllP_Nominal[icent]->SetParLimits(2, 0., 1.);
    fEllP_Nominal[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldNominal[icent]->GetMaximum());
    fEllP_Nominal[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
    if( fixKn ) fEllP_Nominal[icent]->FixParameter(2, 0.37);
    hUnfoldNominal[icent]->Fit( Form("fEllP_Nominal_c%i", icent), "L0", "", 0.0, vnmax[icent]);

    fitE0_Nominal[icent]        = fEllP_Nominal[icent]->GetParameter(0);
    fitAlpha_Nominal[icent]     = fEllP_Nominal[icent]->GetParameter(1);
    fitKn_Nominal[icent]        = fEllP_Nominal[icent]->GetParameter(2);
    fitE0_Nominal_err[icent]    = fEllP_Nominal[icent]->GetParError(0);
    fitAlpha_Nominal_err[icent] = fEllP_Nominal[icent]->GetParError(1);
    fitKn_Nominal_err[icent]    = fEllP_Nominal[icent]->GetParError(2);

    e0    = fEllP_Nominal[icent]->GetParameter(0);
    alpha = fEllP_Nominal[icent]->GetParameter(1);
    kn    = fEllP_Nominal[icent]->GetParameter(2);
    ymax  = 1.9*hUnfoldNominal[icent]->GetMaximum();
    axEcc_Nominal[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_Nominal[icent]->SetLabelSize(0.055);
    axEcc_Nominal[icent]->SetTitleSize(0.055);
    axEcc_Nominal[icent]->SetTitle("#epsilon_{2}");

    asmpt_Nominal[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_Nominal[icent]->SetLineStyle(2);

    cFits_Nominal->cd(icent-2);
    cFits_Nominal->cd(icent-2)->SetLogy();
    cFits_Nominal->cd(icent-2)->SetTopMargin(0.15);
    cFits_Nominal->cd(icent-2)->SetTickx(0);
    hUnfoldNominal[icent]->Draw();
    fEllP_Nominal[icent]->Draw("same");
    axEcc_Nominal[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_Nominal[icent]->Draw("same");

    //-- ------------ SysReg ------------ 
    hUnfoldSysReg[icent] = (TH1D*) fSysReg->Get( Form("hFinalUnfold_SysReg_c%i", icent) );
    hUnfoldSysReg[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldSysReg[icent]->GetXaxis()->SetRange(1, hUnfoldSysReg[icent]->FindBin(vnmax[icent]));
    hUnfoldSysReg[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldSysReg[icent]->Scale(1./hUnfoldSysReg[icent]->Integral());
    hUnfoldSysReg[icent]->SetStats(0);

    if(moscowFits){
      fEllP_SysReg[icent] = (TF1*) fRussiaSysReg->Get( Form("elliptic%i", icent) );
      fEllP_SysReg[icent]->SetLineColor(2);
      fEllP_SysReg[icent]->SetLineWidth(2);

      fitE0_SysReg[icent]        = fEllP_SysReg[icent]->GetParameter(2);
      fitAlpha_SysReg[icent]     = fEllP_SysReg[icent]->GetParameter(0);
      fitKn_SysReg[icent]        = fEllP_SysReg[icent]->GetParameter(1);
      fitE0_SysReg_err[icent]    = fEllP_SysReg[icent]->GetParError(2);
      fitAlpha_SysReg_err[icent] = fEllP_SysReg[icent]->GetParError(0);
      fitKn_SysReg_err[icent]    = fEllP_SysReg[icent]->GetParError(1);

      e0    = fEllP_SysReg[icent]->GetParameter(2);
      alpha = fEllP_SysReg[icent]->GetParameter(0);
      kn    = fEllP_SysReg[icent]->GetParameter(1);
    }
    else{
      fEllP_SysReg[icent] = new TF1(Form("fEllP_SysReg_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_SysReg[icent]->SetLineColor(2);
      fEllP_SysReg[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_SysReg[icent]->SetParLimits(0, 0., 1.);
      fEllP_SysReg[icent]->SetParLimits(2, 0., 1.);
      fEllP_SysReg[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldSysReg[icent]->GetMaximum());
      fEllP_SysReg[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_SysReg[icent]->FixParameter(2,0.37);
      hUnfoldSysReg[icent]->Fit( Form("fEllP_SysReg_c%i", icent), "L0", "", 0.0, vnmax[icent]);

      fitE0_SysReg[icent]        = fEllP_SysReg[icent]->GetParameter(0);
      fitAlpha_SysReg[icent]     = fEllP_SysReg[icent]->GetParameter(1);
      fitKn_SysReg[icent]        = fEllP_SysReg[icent]->GetParameter(2);
      fitE0_SysReg_err[icent]    = fEllP_SysReg[icent]->GetParError(0);
      fitAlpha_SysReg_err[icent] = fEllP_SysReg[icent]->GetParError(1);
      fitKn_SysReg_err[icent]    = fEllP_SysReg[icent]->GetParError(2);

      e0    = fEllP_SysReg[icent]->GetParameter(0);
      alpha = fEllP_SysReg[icent]->GetParameter(1);
      kn    = fEllP_SysReg[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldSysReg[icent]->GetMaximum();
    axEcc_SysReg[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_SysReg[icent]->SetLabelSize(0.055);
    axEcc_SysReg[icent]->SetTitleSize(0.055);
    axEcc_SysReg[icent]->SetTitle("#epsilon_{2}");

    asmpt_SysReg[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_SysReg[icent]->SetLineStyle(2);

    cFits_SysReg->cd(icent-2);
    cFits_SysReg->cd(icent-2)->SetLogy();
    cFits_SysReg->cd(icent-2)->SetTopMargin(0.15);
    cFits_SysReg->cd(icent-2)->SetTickx(0);
    hUnfoldSysReg[icent]->Draw();
    fEllP_SysReg[icent]->Draw("same");
    axEcc_SysReg[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_SysReg[icent]->Draw("same");

    //-- ------------ SysResp ------------
    std::cout << Form("============ SysResp c%i ============", icent) << std::endl;

    hUnfoldSysResp[icent] = (TH1D*) fSysResp->Get( Form("hFinalUnfold_SysResp_c%i", icent) );
    hUnfoldSysResp[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldSysResp[icent]->GetXaxis()->SetRange(1, hUnfoldSysResp[icent]->FindBin(vnmax[icent]));
    hUnfoldSysResp[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldSysResp[icent]->Scale(1./hUnfoldSysResp[icent]->Integral());
    hUnfoldSysResp[icent]->SetStats(0);

    if(moscowFits){
      fEllP_SysResp[icent] = (TF1*) fRussiaSysResp->Get( Form("elliptic%i", icent) );
      fEllP_SysResp[icent]->SetLineColor(2);
      fEllP_SysResp[icent]->SetLineWidth(2);

      fitE0_SysResp[icent]        = fEllP_SysResp[icent]->GetParameter(2);
      fitAlpha_SysResp[icent]     = fEllP_SysResp[icent]->GetParameter(0);
      fitKn_SysResp[icent]        = fEllP_SysResp[icent]->GetParameter(1);
      fitE0_SysResp_err[icent]    = fEllP_SysResp[icent]->GetParError(2);
      fitAlpha_SysResp_err[icent] = fEllP_SysResp[icent]->GetParError(0);
      fitKn_SysResp_err[icent]    = fEllP_SysResp[icent]->GetParError(1);

      e0    = fEllP_SysResp[icent]->GetParameter(2);
      alpha = fEllP_SysResp[icent]->GetParameter(0);
      kn    = fEllP_SysResp[icent]->GetParameter(1);
    }
    else{
      fEllP_SysResp[icent] = new TF1(Form("fEllP_SysResp_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_SysResp[icent]->SetLineColor(2);
      fEllP_SysResp[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_SysResp[icent]->SetParLimits(0, 0., 1.);
      fEllP_SysResp[icent]->SetParLimits(2, 0., 1.);
      fEllP_SysResp[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldSysResp[icent]->GetMaximum());
      fEllP_SysResp[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_SysResp[icent]->FixParameter(2,0.37);
      hUnfoldSysResp[icent]->Fit( Form("fEllP_SysResp_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_SysResp[icent]        = fEllP_SysResp[icent]->GetParameter(0);
      fitAlpha_SysResp[icent]     = fEllP_SysResp[icent]->GetParameter(1);
      fitKn_SysResp[icent]        = fEllP_SysResp[icent]->GetParameter(2);
      fitE0_SysResp_err[icent]    = fEllP_SysResp[icent]->GetParError(0);
      fitAlpha_SysResp_err[icent] = fEllP_SysResp[icent]->GetParError(1);
      fitKn_SysResp_err[icent]    = fEllP_SysResp[icent]->GetParError(2);
      
      e0    = fEllP_SysResp[icent]->GetParameter(0);
      alpha = fEllP_SysResp[icent]->GetParameter(1);
      kn    = fEllP_SysResp[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldSysResp[icent]->GetMaximum();
    axEcc_SysResp[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_SysResp[icent]->SetLabelSize(0.055);
    axEcc_SysResp[icent]->SetTitleSize(0.055);
    axEcc_SysResp[icent]->SetTitle("#epsilon_{2}");

    asmpt_SysResp[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_SysResp[icent]->SetLineStyle(2);

    cFits_SysResp->cd(icent-2);
    cFits_SysResp->cd(icent-2)->SetLogy();
    cFits_SysResp->cd(icent-2)->SetTopMargin(0.15);
    cFits_SysResp->cd(icent-2)->SetTickx(0);
    hUnfoldSysResp[icent]->Draw();
    fEllP_SysResp[icent]->Draw("same");
    axEcc_SysResp[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_SysResp[icent]->Draw("same");

    //-- ------------ SysNewCC ------------
    std::cout << Form("============ SysNewCC c%i ============", icent) << std::endl;

    hUnfoldSysNewCC[icent] = (TH1D*) fSysNewCC->Get( Form("hFinalUnfold_SysNewCC_c%i", icent) );
    hUnfoldSysNewCC[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldSysNewCC[icent]->GetXaxis()->SetRange(1, hUnfoldSysNewCC[icent]->FindBin(vnmax[icent]));
    hUnfoldSysNewCC[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldSysNewCC[icent]->Scale(1./hUnfoldSysNewCC[icent]->Integral());
    hUnfoldSysNewCC[icent]->SetStats(0);

    if(moscowFits){
      fEllP_SysNewCC[icent] = (TF1*) fRussiaSysNewCC->Get( Form("elliptic%i", icent) );
      fEllP_SysNewCC[icent]->SetLineColor(2);
      fEllP_SysNewCC[icent]->SetLineWidth(2);

      fitE0_SysNewCC[icent]        = fEllP_SysNewCC[icent]->GetParameter(2);
      fitAlpha_SysNewCC[icent]     = fEllP_SysNewCC[icent]->GetParameter(0);
      fitKn_SysNewCC[icent]        = fEllP_SysNewCC[icent]->GetParameter(1);
      fitE0_SysNewCC_err[icent]    = fEllP_SysNewCC[icent]->GetParError(2);
      fitAlpha_SysNewCC_err[icent] = fEllP_SysNewCC[icent]->GetParError(0);
      fitKn_SysNewCC_err[icent]    = fEllP_SysNewCC[icent]->GetParError(1);

      e0    = fEllP_SysNewCC[icent]->GetParameter(2);
      alpha = fEllP_SysNewCC[icent]->GetParameter(0);
      kn    = fEllP_SysNewCC[icent]->GetParameter(1);
    }
    else{
      fEllP_SysNewCC[icent] = new TF1(Form("fEllP_SysNewCC_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_SysNewCC[icent]->SetLineColor(2);
      fEllP_SysNewCC[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_SysNewCC[icent]->SetParLimits(0, 0., 1.);
      fEllP_SysNewCC[icent]->SetParLimits(2, 0., 1.);
      fEllP_SysNewCC[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldSysNewCC[icent]->GetMaximum());
      fEllP_SysNewCC[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_SysNewCC[icent]->FixParameter(2,0.37);
      hUnfoldSysNewCC[icent]->Fit( Form("fEllP_SysNewCC_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_SysNewCC[icent]        = fEllP_SysNewCC[icent]->GetParameter(0);
      fitAlpha_SysNewCC[icent]     = fEllP_SysNewCC[icent]->GetParameter(1);
      fitKn_SysNewCC[icent]        = fEllP_SysNewCC[icent]->GetParameter(2);
      fitE0_SysNewCC_err[icent]    = fEllP_SysNewCC[icent]->GetParError(0);
      fitAlpha_SysNewCC_err[icent] = fEllP_SysNewCC[icent]->GetParError(1);
      fitKn_SysNewCC_err[icent]    = fEllP_SysNewCC[icent]->GetParError(2);
      
      e0    = fEllP_SysNewCC[icent]->GetParameter(0);
      alpha = fEllP_SysNewCC[icent]->GetParameter(1);
      kn    = fEllP_SysNewCC[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldSysNewCC[icent]->GetMaximum();
    axEcc_SysNewCC[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_SysNewCC[icent]->SetLabelSize(0.055);
    axEcc_SysNewCC[icent]->SetTitleSize(0.055);
    axEcc_SysNewCC[icent]->SetTitle("#epsilon_{2}");

    asmpt_SysNewCC[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_SysNewCC[icent]->SetLineStyle(2);

    cFits_SysNewCC->cd(icent-2);
    cFits_SysNewCC->cd(icent-2)->SetLogy();
    cFits_SysNewCC->cd(icent-2)->SetTopMargin(0.15);
    cFits_SysNewCC->cd(icent-2)->SetTickx(0);
    hUnfoldSysNewCC[icent]->Draw();
    fEllP_SysNewCC[icent]->Draw("same");
    axEcc_SysNewCC[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_SysNewCC[icent]->Draw("same");

    //-- ------------ LooseSysTkQ ------------
    std::cout << Form("============ LooseSysTkQ c%i ============", icent) << std::endl;

    hUnfoldLooseSysTkQ[icent] = (TH1D*) fSysTkQ->Get( Form("hFinalUnfold_LooseSysTkQ_c%i", icent) );
    hUnfoldLooseSysTkQ[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldLooseSysTkQ[icent]->GetXaxis()->SetRange(1, hUnfoldLooseSysTkQ[icent]->FindBin(vnmax[icent]));
    hUnfoldLooseSysTkQ[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldLooseSysTkQ[icent]->Scale(1./hUnfoldLooseSysTkQ[icent]->Integral());
    hUnfoldLooseSysTkQ[icent]->SetStats(0);

    if(moscowFits){
      fEllP_LooseSysTkQ[icent] = (TF1*) fRussiaSysTkQ1->Get( Form("elliptic%i", icent) );
      fEllP_LooseSysTkQ[icent]->SetLineColor(2);
      fEllP_LooseSysTkQ[icent]->SetLineWidth(2);

      fitE0_LooseSysTkQ[icent]        = fEllP_LooseSysTkQ[icent]->GetParameter(2);
      fitAlpha_LooseSysTkQ[icent]     = fEllP_LooseSysTkQ[icent]->GetParameter(0);
      fitKn_LooseSysTkQ[icent]        = fEllP_LooseSysTkQ[icent]->GetParameter(1);
      fitE0_LooseSysTkQ_err[icent]    = fEllP_LooseSysTkQ[icent]->GetParError(2);
      fitAlpha_LooseSysTkQ_err[icent] = fEllP_LooseSysTkQ[icent]->GetParError(0);
      fitKn_LooseSysTkQ_err[icent]    = fEllP_LooseSysTkQ[icent]->GetParError(1);

      e0    = fEllP_LooseSysTkQ[icent]->GetParameter(2);
      alpha = fEllP_LooseSysTkQ[icent]->GetParameter(0);
      kn    = fEllP_LooseSysTkQ[icent]->GetParameter(1);
    }
    else{
      fEllP_LooseSysTkQ[icent] = new TF1(Form("fEllP_LooseSysTkQ_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_LooseSysTkQ[icent]->SetLineColor(2);
      fEllP_LooseSysTkQ[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_LooseSysTkQ[icent]->SetParLimits(0, 0., 1.);
      fEllP_LooseSysTkQ[icent]->SetParLimits(2, 0., 1.);
      fEllP_LooseSysTkQ[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldLooseSysTkQ[icent]->GetMaximum());
      fEllP_LooseSysTkQ[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_LooseSysTkQ[icent]->FixParameter(2,0.37);
      hUnfoldLooseSysTkQ[icent]->Fit( Form("fEllP_LooseSysTkQ_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_LooseSysTkQ[icent]        = fEllP_LooseSysTkQ[icent]->GetParameter(0);
      fitAlpha_LooseSysTkQ[icent]     = fEllP_LooseSysTkQ[icent]->GetParameter(1);
      fitKn_LooseSysTkQ[icent]        = fEllP_LooseSysTkQ[icent]->GetParameter(2);
      fitE0_LooseSysTkQ_err[icent]    = fEllP_LooseSysTkQ[icent]->GetParError(0);
      fitAlpha_LooseSysTkQ_err[icent] = fEllP_LooseSysTkQ[icent]->GetParError(1);
      fitKn_LooseSysTkQ_err[icent]    = fEllP_LooseSysTkQ[icent]->GetParError(2);
      
      e0    = fEllP_LooseSysTkQ[icent]->GetParameter(0);
      alpha = fEllP_LooseSysTkQ[icent]->GetParameter(1);
      kn    = fEllP_LooseSysTkQ[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldLooseSysTkQ[icent]->GetMaximum();
    axEcc_LooseSysTkQ[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_LooseSysTkQ[icent]->SetLabelSize(0.055);
    axEcc_LooseSysTkQ[icent]->SetTitleSize(0.055);
    axEcc_LooseSysTkQ[icent]->SetTitle("#epsilon_{2}");

    asmpt_LooseSysTkQ[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_LooseSysTkQ[icent]->SetLineStyle(2);

    cFits_LooseSysTkQ->cd(icent-2);
    cFits_LooseSysTkQ->cd(icent-2)->SetLogy();
    cFits_LooseSysTkQ->cd(icent-2)->SetTopMargin(0.15);
    cFits_LooseSysTkQ->cd(icent-2)->SetTickx(0);
    hUnfoldLooseSysTkQ[icent]->Draw();
    fEllP_LooseSysTkQ[icent]->Draw("same");
    axEcc_LooseSysTkQ[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_LooseSysTkQ[icent]->Draw("same");

    //-- ------------ TightSysTkQ ------------
    std::cout << Form("============ TightSysTkQ c%i ============", icent) << std::endl;

    hUnfoldTightSysTkQ[icent] = (TH1D*) fSysTkQ->Get( Form("hFinalUnfold_TightSysTkQ_c%i", icent) );
    hUnfoldTightSysTkQ[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldTightSysTkQ[icent]->GetXaxis()->SetRange(1, hUnfoldTightSysTkQ[icent]->FindBin(vnmax[icent]));
    hUnfoldTightSysTkQ[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldTightSysTkQ[icent]->Scale(1./hUnfoldTightSysTkQ[icent]->Integral());
    hUnfoldTightSysTkQ[icent]->SetStats(0);

    if(moscowFits){
      fEllP_TightSysTkQ[icent] = (TF1*) fRussiaSysTkQ2->Get( Form("elliptic%i", icent) );
      fEllP_TightSysTkQ[icent]->SetLineColor(2);
      fEllP_TightSysTkQ[icent]->SetLineWidth(2);

      fitE0_TightSysTkQ[icent]        = fEllP_TightSysTkQ[icent]->GetParameter(2);
      fitAlpha_TightSysTkQ[icent]     = fEllP_TightSysTkQ[icent]->GetParameter(0);
      fitKn_TightSysTkQ[icent]        = fEllP_TightSysTkQ[icent]->GetParameter(1);
      fitE0_TightSysTkQ_err[icent]    = fEllP_TightSysTkQ[icent]->GetParError(2);
      fitAlpha_TightSysTkQ_err[icent] = fEllP_TightSysTkQ[icent]->GetParError(0);
      fitKn_TightSysTkQ_err[icent]    = fEllP_TightSysTkQ[icent]->GetParError(1);

      e0    = fEllP_TightSysTkQ[icent]->GetParameter(2);
      alpha = fEllP_TightSysTkQ[icent]->GetParameter(0);
      kn    = fEllP_TightSysTkQ[icent]->GetParameter(1);
    }
    else{
      fEllP_TightSysTkQ[icent] = new TF1(Form("fEllP_TightSysTkQ_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_TightSysTkQ[icent]->SetLineColor(2);
      fEllP_TightSysTkQ[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_TightSysTkQ[icent]->SetParLimits(0, 0., 1.);
      fEllP_TightSysTkQ[icent]->SetParLimits(2, 0., 1.);
      fEllP_TightSysTkQ[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldTightSysTkQ[icent]->GetMaximum());
      fEllP_TightSysTkQ[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_TightSysTkQ[icent]->FixParameter(2,0.37);
      hUnfoldTightSysTkQ[icent]->Fit( Form("fEllP_TightSysTkQ_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_TightSysTkQ[icent]        = fEllP_TightSysTkQ[icent]->GetParameter(0);
      fitAlpha_TightSysTkQ[icent]     = fEllP_TightSysTkQ[icent]->GetParameter(1);
      fitKn_TightSysTkQ[icent]        = fEllP_TightSysTkQ[icent]->GetParameter(2);
      fitE0_TightSysTkQ_err[icent]    = fEllP_TightSysTkQ[icent]->GetParError(0);
      fitAlpha_TightSysTkQ_err[icent] = fEllP_TightSysTkQ[icent]->GetParError(1);
      fitKn_TightSysTkQ_err[icent]    = fEllP_TightSysTkQ[icent]->GetParError(2);
      
      e0    = fEllP_TightSysTkQ[icent]->GetParameter(0);
      alpha = fEllP_TightSysTkQ[icent]->GetParameter(1);
      kn    = fEllP_TightSysTkQ[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldTightSysTkQ[icent]->GetMaximum();
    axEcc_TightSysTkQ[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_TightSysTkQ[icent]->SetLabelSize(0.055);
    axEcc_TightSysTkQ[icent]->SetTitleSize(0.055);
    axEcc_TightSysTkQ[icent]->SetTitle("#epsilon_{2}");

    asmpt_TightSysTkQ[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_TightSysTkQ[icent]->SetLineStyle(2);

    cFits_TightSysTkQ->cd(icent-2);
    cFits_TightSysTkQ->cd(icent-2)->SetLogy();
    cFits_TightSysTkQ->cd(icent-2)->SetTopMargin(0.15);
    cFits_TightSysTkQ->cd(icent-2)->SetTickx(0);
    hUnfoldTightSysTkQ[icent]->Draw();
    fEllP_TightSysTkQ[icent]->Draw("same");
    axEcc_TightSysTkQ[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_TightSysTkQ[icent]->Draw("same");

    //-- ------------ Vtx3SysVtx ------------
    std::cout << Form("============ Vtx3SysVtx c%i ============", icent) << std::endl;

    hUnfoldVtx3SysVtx[icent] = (TH1D*) fSysVtx->Get( Form("hFinalUnfold_Vtx3SysVtx_c%i", icent) );
    hUnfoldVtx3SysVtx[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldVtx3SysVtx[icent]->GetXaxis()->SetRange(1, hUnfoldVtx3SysVtx[icent]->FindBin(vnmax[icent]));
    hUnfoldVtx3SysVtx[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldVtx3SysVtx[icent]->Scale(1./hUnfoldVtx3SysVtx[icent]->Integral());
    hUnfoldVtx3SysVtx[icent]->SetStats(0);

    if(moscowFits){
      fEllP_Vtx3SysVtx[icent] = (TF1*) fRussiaSysVtx1->Get( Form("elliptic%i", icent) );
      fEllP_Vtx3SysVtx[icent]->SetLineColor(2);
      fEllP_Vtx3SysVtx[icent]->SetLineWidth(2);

      fitE0_Vtx3SysVtx[icent]        = fEllP_Vtx3SysVtx[icent]->GetParameter(2);
      fitAlpha_Vtx3SysVtx[icent]     = fEllP_Vtx3SysVtx[icent]->GetParameter(0);
      fitKn_Vtx3SysVtx[icent]        = fEllP_Vtx3SysVtx[icent]->GetParameter(1);
      fitE0_Vtx3SysVtx_err[icent]    = fEllP_Vtx3SysVtx[icent]->GetParError(2);
      fitAlpha_Vtx3SysVtx_err[icent] = fEllP_Vtx3SysVtx[icent]->GetParError(0);
      fitKn_Vtx3SysVtx_err[icent]    = fEllP_Vtx3SysVtx[icent]->GetParError(1);

      e0    = fEllP_Vtx3SysVtx[icent]->GetParameter(2);
      alpha = fEllP_Vtx3SysVtx[icent]->GetParameter(0);
      kn    = fEllP_Vtx3SysVtx[icent]->GetParameter(1);
    }
    else{
      fEllP_Vtx3SysVtx[icent] = new TF1(Form("fEllP_Vtx3SysVtx_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_Vtx3SysVtx[icent]->SetLineColor(2);
      fEllP_Vtx3SysVtx[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_Vtx3SysVtx[icent]->SetParLimits(0, 0., 1.);
      fEllP_Vtx3SysVtx[icent]->SetParLimits(2, 0., 1.);
      fEllP_Vtx3SysVtx[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldVtx3SysVtx[icent]->GetMaximum());
      fEllP_Vtx3SysVtx[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_Vtx3SysVtx[icent]->FixParameter(2,0.37);
      hUnfoldVtx3SysVtx[icent]->Fit( Form("fEllP_Vtx3SysVtx_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_Vtx3SysVtx[icent]        = fEllP_Vtx3SysVtx[icent]->GetParameter(0);
      fitAlpha_Vtx3SysVtx[icent]     = fEllP_Vtx3SysVtx[icent]->GetParameter(1);
      fitKn_Vtx3SysVtx[icent]        = fEllP_Vtx3SysVtx[icent]->GetParameter(2);
      fitE0_Vtx3SysVtx_err[icent]    = fEllP_Vtx3SysVtx[icent]->GetParError(0);
      fitAlpha_Vtx3SysVtx_err[icent] = fEllP_Vtx3SysVtx[icent]->GetParError(1);
      fitKn_Vtx3SysVtx_err[icent]    = fEllP_Vtx3SysVtx[icent]->GetParError(2);
      
      e0    = fEllP_Vtx3SysVtx[icent]->GetParameter(0);
      alpha = fEllP_Vtx3SysVtx[icent]->GetParameter(1);
      kn    = fEllP_Vtx3SysVtx[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldVtx3SysVtx[icent]->GetMaximum();
    axEcc_Vtx3SysVtx[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_Vtx3SysVtx[icent]->SetLabelSize(0.055);
    axEcc_Vtx3SysVtx[icent]->SetTitleSize(0.055);
    axEcc_Vtx3SysVtx[icent]->SetTitle("#epsilon_{2}");

    asmpt_Vtx3SysVtx[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_Vtx3SysVtx[icent]->SetLineStyle(2);

    cFits_Vtx3SysVtx->cd(icent-2);
    cFits_Vtx3SysVtx->cd(icent-2)->SetLogy();
    cFits_Vtx3SysVtx->cd(icent-2)->SetTopMargin(0.15);
    cFits_Vtx3SysVtx->cd(icent-2)->SetTickx(0);
    hUnfoldVtx3SysVtx[icent]->Draw();
    fEllP_Vtx3SysVtx[icent]->Draw("same");
    axEcc_Vtx3SysVtx[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_Vtx3SysVtx[icent]->Draw("same");

    //-- ------------ Vtx3_15SysVtx ------------
    std::cout << Form("============ Vtx3_15SysVtx c%i ============", icent) << std::endl;

    hUnfoldVtx3_15SysVtx[icent] = (TH1D*) fSysVtx->Get( Form("hFinalUnfold_Vtx3_15SysVtx_c%i", icent) );
    hUnfoldVtx3_15SysVtx[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hUnfoldVtx3_15SysVtx[icent]->GetXaxis()->SetRange(1, hUnfoldVtx3_15SysVtx[icent]->FindBin(vnmax[icent]));
    hUnfoldVtx3_15SysVtx[icent]->GetYaxis()->SetTitle(Form("p(v_{%i}), p(#epsilon_{%i})", norder_, norder_) );
    //hUnfoldVtx3_15SysVtx[icent]->Scale(1./hUnfoldVtx3_15SysVtx[icent]->Integral());
    hUnfoldVtx3_15SysVtx[icent]->SetStats(0);

    if(moscowFits){
      fEllP_Vtx3_15SysVtx[icent] = (TF1*) fRussiaSysVtx2->Get( Form("elliptic%i", icent) );
      fEllP_Vtx3_15SysVtx[icent]->SetLineColor(2);
      fEllP_Vtx3_15SysVtx[icent]->SetLineWidth(2);

      fitE0_Vtx3_15SysVtx[icent]        = fEllP_Vtx3_15SysVtx[icent]->GetParameter(2);
      fitAlpha_Vtx3_15SysVtx[icent]     = fEllP_Vtx3_15SysVtx[icent]->GetParameter(0);
      fitKn_Vtx3_15SysVtx[icent]        = fEllP_Vtx3_15SysVtx[icent]->GetParameter(1);
      fitE0_Vtx3_15SysVtx_err[icent]    = fEllP_Vtx3_15SysVtx[icent]->GetParError(2);
      fitAlpha_Vtx3_15SysVtx_err[icent] = fEllP_Vtx3_15SysVtx[icent]->GetParError(0);
      fitKn_Vtx3_15SysVtx_err[icent]    = fEllP_Vtx3_15SysVtx[icent]->GetParError(1);

      e0    = fEllP_Vtx3_15SysVtx[icent]->GetParameter(2);
      alpha = fEllP_Vtx3_15SysVtx[icent]->GetParameter(0);
      kn    = fEllP_Vtx3_15SysVtx[icent]->GetParameter(1);
    }
    else{
      fEllP_Vtx3_15SysVtx[icent] = new TF1(Form("fEllP_Vtx3_15SysVtx_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
      fEllP_Vtx3_15SysVtx[icent]->SetLineColor(2);
      fEllP_Vtx3_15SysVtx[icent]->SetLineWidth(2);
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP_Vtx3_15SysVtx[icent]->SetParLimits(0, 0., 1.);
      fEllP_Vtx3_15SysVtx[icent]->SetParLimits(2, 0., 1.);
      fEllP_Vtx3_15SysVtx[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hUnfoldVtx3_15SysVtx[icent]->GetMaximum());
      fEllP_Vtx3_15SysVtx[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
      if(fixKn )fEllP_Vtx3_15SysVtx[icent]->FixParameter(2,0.37);
      hUnfoldVtx3_15SysVtx[icent]->Fit( Form("fEllP_Vtx3_15SysVtx_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      fitE0_Vtx3_15SysVtx[icent]        = fEllP_Vtx3_15SysVtx[icent]->GetParameter(0);
      fitAlpha_Vtx3_15SysVtx[icent]     = fEllP_Vtx3_15SysVtx[icent]->GetParameter(1);
      fitKn_Vtx3_15SysVtx[icent]        = fEllP_Vtx3_15SysVtx[icent]->GetParameter(2);
      fitE0_Vtx3_15SysVtx_err[icent]    = fEllP_Vtx3_15SysVtx[icent]->GetParError(0);
      fitAlpha_Vtx3_15SysVtx_err[icent] = fEllP_Vtx3_15SysVtx[icent]->GetParError(1);
      fitKn_Vtx3_15SysVtx_err[icent]    = fEllP_Vtx3_15SysVtx[icent]->GetParError(2);
      
      e0    = fEllP_Vtx3_15SysVtx[icent]->GetParameter(0);
      alpha = fEllP_Vtx3_15SysVtx[icent]->GetParameter(1);
      kn    = fEllP_Vtx3_15SysVtx[icent]->GetParameter(2);
    }

    ymax  = 1.9*hUnfoldVtx3_15SysVtx[icent]->GetMaximum();
    axEcc_Vtx3_15SysVtx[icent] = new TGaxis(xmin, ymax, xmax, ymax, xmin/kn, xmax/kn, 507, "-");
    axEcc_Vtx3_15SysVtx[icent]->SetLabelSize(0.055);
    axEcc_Vtx3_15SysVtx[icent]->SetTitleSize(0.055);
    axEcc_Vtx3_15SysVtx[icent]->SetTitle("#epsilon_{2}");

    asmpt_Vtx3_15SysVtx[icent] = new TLine(kn, 0, kn, ymax);
    asmpt_Vtx3_15SysVtx[icent]->SetLineStyle(2);

    cFits_Vtx3_15SysVtx->cd(icent-2);
    cFits_Vtx3_15SysVtx->cd(icent-2)->SetLogy();
    cFits_Vtx3_15SysVtx->cd(icent-2)->SetTopMargin(0.15);
    cFits_Vtx3_15SysVtx->cd(icent-2)->SetTickx(0);
    hUnfoldVtx3_15SysVtx[icent]->Draw();
    fEllP_Vtx3_15SysVtx[icent]->Draw("same");
    axEcc_Vtx3_15SysVtx[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    asmpt_Vtx3_15SysVtx[icent]->Draw("same");

  } //-- End Cent loop

  cFits_SysReg->SaveAs("../plots/systematicStudies/cFits_SysReg.pdf");
  cFits_SysResp->SaveAs("../plots/systematicStudies/cFits_SysResp.pdf");
  cFits_SysNewCC->SaveAs("../plots/systematicStudies/cFits_SysNewCC.pdf");
  cFits_LooseSysTkQ->SaveAs("../plots/systematicStudies/cFits_LooseSysTkQ.pdf");
  cFits_TightSysTkQ->SaveAs("../plots/systematicStudies/cFits_TightSysTkQ.pdf");
  cFits_Vtx3SysVtx->SaveAs("../plots/systematicStudies/cFits_Vtx3SysVtx.pdf");
  cFits_Vtx3_15SysVtx->SaveAs("../plots/systematicStudies/cFits_Vtx3_15SysVtx.pdf");

  //-- Let's do some systematics!
  for(int icent = 3; icent < NCENT; icent++){

    //-- Nominal
    double E0 = fitE0_Nominal[icent];
    double AL = fitAlpha_Nominal[icent];
    double KN = fitKn_Nominal[icent];

    double E0e = fitE0_Nominal_err[icent];
    double ALe = fitAlpha_Nominal_err[icent];
    double KNe = fitKn_Nominal_err[icent];

    //-- SysReg
    double e0 = fitE0_SysReg[icent];
    double al = fitAlpha_SysReg[icent];
    double kn = fitKn_SysReg[icent];

    double e0e = fitE0_SysReg_err[icent];
    double ale = fitAlpha_SysReg_err[icent];
    double kne = fitKn_SysReg_err[icent];

    double r_e0 = e0 / E0;
    double r_al = al / AL; 
    double r_kn = kn / KN;

    double r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    double r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    double r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_SysReg[icent]        = r_e0;
    fitAlpha_RatioToNom_SysReg[icent]     = r_al;
    fitKn_RatioToNom_SysReg[icent]        = r_kn;
    fitE0_RatioToNom_SysReg_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_SysReg_err[icent] = r_ale;
    fitKn_RatioToNom_SysReg_err[icent]    = r_kne;

    //-- SysResp
    e0 = fitE0_SysResp[icent];
    al = fitAlpha_SysResp[icent];
    kn = fitKn_SysResp[icent];

    e0e = fitE0_SysResp_err[icent];
    ale = fitAlpha_SysResp_err[icent];
    kne = fitKn_SysResp_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_SysResp[icent]        = r_e0;
    fitAlpha_RatioToNom_SysResp[icent]     = r_al;
    fitKn_RatioToNom_SysResp[icent]        = r_kn;
    fitE0_RatioToNom_SysResp_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_SysResp_err[icent] = r_ale;
    fitKn_RatioToNom_SysResp_err[icent]    = r_kne;

    //-- SysNewCC
    e0 = fitE0_SysNewCC[icent];
    al = fitAlpha_SysNewCC[icent];
    kn = fitKn_SysNewCC[icent];

    e0e = fitE0_SysNewCC_err[icent];
    ale = fitAlpha_SysNewCC_err[icent];
    kne = fitKn_SysNewCC_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_SysNewCC[icent]        = r_e0;
    fitAlpha_RatioToNom_SysNewCC[icent]     = r_al;
    fitKn_RatioToNom_SysNewCC[icent]        = r_kn;
    fitE0_RatioToNom_SysNewCC_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_SysNewCC_err[icent] = r_ale;
    fitKn_RatioToNom_SysNewCC_err[icent]    = r_kne;

    //-- LooseSysTkQ
    e0 = fitE0_LooseSysTkQ[icent];
    al = fitAlpha_LooseSysTkQ[icent];
    kn = fitKn_LooseSysTkQ[icent];

    e0e = fitE0_LooseSysTkQ_err[icent];
    ale = fitAlpha_LooseSysTkQ_err[icent];
    kne = fitKn_LooseSysTkQ_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_LooseSysTkQ[icent]        = r_e0;
    fitAlpha_RatioToNom_LooseSysTkQ[icent]     = r_al;
    fitKn_RatioToNom_LooseSysTkQ[icent]        = r_kn;
    fitE0_RatioToNom_LooseSysTkQ_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_LooseSysTkQ_err[icent] = r_ale;
    fitKn_RatioToNom_LooseSysTkQ_err[icent]    = r_kne;

    //-- TightSysTkQ
    e0 = fitE0_TightSysTkQ[icent];
    al = fitAlpha_TightSysTkQ[icent];
    kn = fitKn_TightSysTkQ[icent];

    e0e = fitE0_TightSysTkQ_err[icent];
    ale = fitAlpha_TightSysTkQ_err[icent];
    kne = fitKn_TightSysTkQ_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_TightSysTkQ[icent]        = r_e0;
    fitAlpha_RatioToNom_TightSysTkQ[icent]     = r_al;
    fitKn_RatioToNom_TightSysTkQ[icent]        = r_kn;
    fitE0_RatioToNom_TightSysTkQ_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_TightSysTkQ_err[icent] = r_ale;
    fitKn_RatioToNom_TightSysTkQ_err[icent]    = r_kne;

    //-- Vtx3SysVtx
    e0 = fitE0_Vtx3SysVtx[icent];
    al = fitAlpha_Vtx3SysVtx[icent];
    kn = fitKn_Vtx3SysVtx[icent];

    e0e = fitE0_Vtx3SysVtx_err[icent];
    ale = fitAlpha_Vtx3SysVtx_err[icent];
    kne = fitKn_Vtx3SysVtx_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_Vtx3SysVtx[icent]        = r_e0;
    fitAlpha_RatioToNom_Vtx3SysVtx[icent]     = r_al;
    fitKn_RatioToNom_Vtx3SysVtx[icent]        = r_kn;
    fitE0_RatioToNom_Vtx3SysVtx_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_Vtx3SysVtx_err[icent] = r_ale;
    fitKn_RatioToNom_Vtx3SysVtx_err[icent]    = r_kne;

    //-- Vtx3_15SysVtx
    e0 = fitE0_Vtx3_15SysVtx[icent];
    al = fitAlpha_Vtx3_15SysVtx[icent];
    kn = fitKn_Vtx3_15SysVtx[icent];

    e0e = fitE0_Vtx3_15SysVtx_err[icent];
    ale = fitAlpha_Vtx3_15SysVtx_err[icent];
    kne = fitKn_Vtx3_15SysVtx_err[icent];

    r_e0 = e0 / E0;
    r_al = al / AL; 
    r_kn = kn / KN;

    r_e0e = sqrt( pow(e0e/E0, 2) + pow(e0*E0e/E0/E0, 2) );
    r_ale = sqrt( pow(ale/AL, 2) + pow(al*ALe/AL/AL, 2) );
    r_kne = sqrt( pow(kne/KN, 2) + pow(kn*KNe/KN/KN, 2) );

    fitE0_RatioToNom_Vtx3_15SysVtx[icent]        = r_e0;
    fitAlpha_RatioToNom_Vtx3_15SysVtx[icent]     = r_al;
    fitKn_RatioToNom_Vtx3_15SysVtx[icent]        = r_kn;
    fitE0_RatioToNom_Vtx3_15SysVtx_err[icent]    = r_e0e;
    fitAlpha_RatioToNom_Vtx3_15SysVtx_err[icent] = r_ale;
    fitKn_RatioToNom_Vtx3_15SysVtx_err[icent]    = r_kne;

  }

  //-- Make graphs

  //-- SysReg/Nominal
  grFitE0_RatioToNom_SysReg    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_SysReg, CERR,    fitE0_RatioToNom_SysReg_err);
  grFitAlpha_RatioToNom_SysReg = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_SysReg, CERR, fitAlpha_RatioToNom_SysReg_err);
  grFitKn_RatioToNom_SysReg    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_SysReg, CERR,    fitKn_RatioToNom_SysReg_err);

  formatGraph(grFitE0_RatioToNom_SysReg,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 34, "grFitE0_RatioToNom_SysReg");
  formatGraph(grFitAlpha_RatioToNom_SysReg, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 21, "grFitAlpha_RatioToNom_SysReg");
  formatGraph(grFitKn_RatioToNom_SysReg,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 20, "grFitKn_RatioToNom_SysReg");

  //-- SysResp/Nominal
  grFitE0_RatioToNom_SysResp    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_SysResp, CERR,    fitE0_RatioToNom_SysResp_err);
  grFitAlpha_RatioToNom_SysResp = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_SysResp, CERR, fitAlpha_RatioToNom_SysResp_err);
  grFitKn_RatioToNom_SysResp    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_SysResp, CERR,    fitKn_RatioToNom_SysResp_err);

  formatGraph(grFitE0_RatioToNom_SysResp,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 34, "grFitE0_RatioToNom_SysResp");
  formatGraph(grFitAlpha_RatioToNom_SysResp, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 21, "grFitAlpha_RatioToNom_SysResp");
  formatGraph(grFitKn_RatioToNom_SysResp,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 20, "grFitKn_RatioToNom_SysResp");

  //-- SysNewCC/Nominal
  grFitE0_RatioToNom_SysNewCC    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_SysNewCC, CERR,    fitE0_RatioToNom_SysNewCC_err);
  grFitAlpha_RatioToNom_SysNewCC = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_SysNewCC, CERR, fitAlpha_RatioToNom_SysNewCC_err);
  grFitKn_RatioToNom_SysNewCC    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_SysNewCC, CERR,    fitKn_RatioToNom_SysNewCC_err);

  formatGraph(grFitE0_RatioToNom_SysNewCC,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 34, "grFitE0_RatioToNom_SysNewCC");
  formatGraph(grFitAlpha_RatioToNom_SysNewCC, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 21, "grFitAlpha_RatioToNom_SysNewCC");
  formatGraph(grFitKn_RatioToNom_SysNewCC,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 20, "grFitKn_RatioToNom_SysNewCC");

  //-- LooseSysTkQ/Nominal
  grFitE0_RatioToNom_LooseSysTkQ    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_LooseSysTkQ, CERR,    fitE0_RatioToNom_LooseSysTkQ_err);
  grFitAlpha_RatioToNom_LooseSysTkQ = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_LooseSysTkQ, CERR, fitAlpha_RatioToNom_LooseSysTkQ_err);
  grFitKn_RatioToNom_LooseSysTkQ    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_LooseSysTkQ, CERR,    fitKn_RatioToNom_LooseSysTkQ_err);

  formatGraph(grFitE0_RatioToNom_LooseSysTkQ,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 34, "grFitE0_RatioToNom_LooseSysTkQ");
  formatGraph(grFitAlpha_RatioToNom_LooseSysTkQ, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 21, "grFitAlpha_RatioToNom_LooseSysTkQ");
  formatGraph(grFitKn_RatioToNom_LooseSysTkQ,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 20, "grFitKn_RatioToNom_LooseSysTkQ");

  //-- TightSysTkQ/Nominal
  grFitE0_RatioToNom_TightSysTkQ    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_TightSysTkQ, CERR,    fitE0_RatioToNom_TightSysTkQ_err);
  grFitAlpha_RatioToNom_TightSysTkQ = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_TightSysTkQ, CERR, fitAlpha_RatioToNom_TightSysTkQ_err);
  grFitKn_RatioToNom_TightSysTkQ    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_TightSysTkQ, CERR,    fitKn_RatioToNom_TightSysTkQ_err);

  formatGraph(grFitE0_RatioToNom_TightSysTkQ,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 28, "grFitE0_RatioToNom_TightSysTkQ");
  formatGraph(grFitAlpha_RatioToNom_TightSysTkQ, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 25, "grFitAlpha_RatioToNom_TightSysTkQ");
  formatGraph(grFitKn_RatioToNom_TightSysTkQ,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 24, "grFitKn_RatioToNom_TightSysTkQ");

  //-- Vtx3SysVtx/Nominal
  grFitE0_RatioToNom_Vtx3SysVtx    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_Vtx3SysVtx, CERR,    fitE0_RatioToNom_Vtx3SysVtx_err);
  grFitAlpha_RatioToNom_Vtx3SysVtx = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_Vtx3SysVtx, CERR, fitAlpha_RatioToNom_Vtx3SysVtx_err);
  grFitKn_RatioToNom_Vtx3SysVtx    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_Vtx3SysVtx, CERR,    fitKn_RatioToNom_Vtx3SysVtx_err);

  formatGraph(grFitE0_RatioToNom_Vtx3SysVtx,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 34, "grFitE0_RatioToNom_Vtx3SysVtx");
  formatGraph(grFitAlpha_RatioToNom_Vtx3SysVtx, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 21, "grFitAlpha_RatioToNom_Vtx3SysVtx");
  formatGraph(grFitKn_RatioToNom_Vtx3SysVtx,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 20, "grFitKn_RatioToNom_Vtx3SysVtx");

  //-- Vtx3_15SysVtx/Nominal
  grFitE0_RatioToNom_Vtx3_15SysVtx    = new TGraphErrors(NCENT, centBinCenter, fitE0_RatioToNom_Vtx3_15SysVtx, CERR,    fitE0_RatioToNom_Vtx3_15SysVtx_err);
  grFitAlpha_RatioToNom_Vtx3_15SysVtx = new TGraphErrors(NCENT, centBinCenter, fitAlpha_RatioToNom_Vtx3_15SysVtx, CERR, fitAlpha_RatioToNom_Vtx3_15SysVtx_err);
  grFitKn_RatioToNom_Vtx3_15SysVtx    = new TGraphErrors(NCENT, centBinCenter, fitKn_RatioToNom_Vtx3_15SysVtx, CERR,    fitKn_RatioToNom_Vtx3_15SysVtx_err);

  formatGraph(grFitE0_RatioToNom_Vtx3_15SysVtx,    "Centrality %", sysMin, sysMax, "#epsilon_{0}: Ratio to Nominal", 4, 28, "grFitE0_RatioToNom_Vtx3_15SysVtx");
  formatGraph(grFitAlpha_RatioToNom_Vtx3_15SysVtx, "Centrality %", sysMin, sysMax, "#alpha: Ratio to Nominal",       2, 25, "grFitAlpha_RatioToNom_Vtx3_15SysVtx");
  formatGraph(grFitKn_RatioToNom_Vtx3_15SysVtx,    "Centrality %", sysMin, sysMax, "k_{n}: Ratio to Nominal",        1, 24, "grFitKn_RatioToNom_Vtx3_15SysVtx");

  //-- Write Graphs to a file:
  TFile * fOut = new TFile("SysEllpParms.root", "recreate");
  fOut->cd();
  grFitE0_RatioToNom_SysReg    -> Write();
  grFitAlpha_RatioToNom_SysReg -> Write();
  grFitKn_RatioToNom_SysReg    -> Write();

  grFitE0_RatioToNom_SysResp    -> Write();
  grFitAlpha_RatioToNom_SysResp -> Write();
  grFitKn_RatioToNom_SysResp    -> Write();

  grFitE0_RatioToNom_SysNewCC    -> Write();
  grFitAlpha_RatioToNom_SysNewCC -> Write();
  grFitKn_RatioToNom_SysNewCC    -> Write();

  grFitE0_RatioToNom_LooseSysTkQ    -> Write();
  grFitAlpha_RatioToNom_LooseSysTkQ -> Write();
  grFitKn_RatioToNom_LooseSysTkQ    -> Write();

  grFitE0_RatioToNom_TightSysTkQ    -> Write();
  grFitAlpha_RatioToNom_TightSysTkQ -> Write();
  grFitKn_RatioToNom_TightSysTkQ    -> Write();

  grFitE0_RatioToNom_Vtx3SysVtx    -> Write();
  grFitAlpha_RatioToNom_Vtx3SysVtx -> Write();
  grFitKn_RatioToNom_Vtx3SysVtx    -> Write();

  grFitE0_RatioToNom_Vtx3_15SysVtx    -> Write();
  grFitAlpha_RatioToNom_Vtx3_15SysVtx -> Write();
  grFitKn_RatioToNom_Vtx3_15SysVtx    -> Write();

  //-- Draw
  TLine * lone = new TLine(0., 1., grFitKn_RatioToNom_SysReg->GetXaxis()->GetXmax(), 1.);
  lone->SetLineStyle(2);
  lone->SetLineWidth(2);

  //-- SysReg
  TCanvas * cSysReg = new TCanvas("cSysReg", "cSysReg", 1500, 500);
  cSysReg->Divide(3,1);

  cSysReg->cd(1);
  grFitKn_RatioToNom_SysReg->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysReg");

  cSysReg->cd(2);
  grFitE0_RatioToNom_SysReg->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysReg");

  cSysReg->cd(3);
  grFitAlpha_RatioToNom_SysReg->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysReg");

  cSysReg->SaveAs("../plots/systematicStudies/cSysFitParms_SysReg.pdf");

  //-- SysResp
  TCanvas * cSysResp = new TCanvas("cSysResp", "cSysResp", 1500, 500);
  cSysResp->Divide(3,1);

  cSysResp->cd(1);
  grFitKn_RatioToNom_SysResp->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysResp");

  cSysResp->cd(2);
  grFitE0_RatioToNom_SysResp->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysResp");

  cSysResp->cd(3);
  grFitAlpha_RatioToNom_SysResp->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysResp");

  cSysReg->SaveAs("../plots/systematicStudies/cSysFitParms_SysResp.pdf");

  //-- SysNewCC
  TCanvas * cSysNewCC = new TCanvas("cSysNewCC", "cSysNewCC", 1500, 500);
  cSysNewCC->Divide(3,1);

  cSysNewCC->cd(1);
  grFitKn_RatioToNom_SysNewCC->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysNewCC");

  cSysNewCC->cd(2);
  grFitE0_RatioToNom_SysNewCC->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysNewCC");

  cSysNewCC->cd(3);
  grFitAlpha_RatioToNom_SysNewCC->Draw("ap");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysNewCC");

  cSysReg->SaveAs("../plots/systematicStudies/cSysFitParms_SysNewCC.pdf");

  //-- SysTkQ
  TLegend * legKnTkQ = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legKnTkQ);
  legKnTkQ->AddEntry(grFitKn_RatioToNom_LooseSysTkQ, "Loose/Nominal", "lp");
  legKnTkQ->AddEntry(grFitKn_RatioToNom_TightSysTkQ, "Tight/Nominal", "lp");

  TLegend * legE0TkQ = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legE0TkQ);
  legE0TkQ->AddEntry(grFitE0_RatioToNom_LooseSysTkQ, "Loose/Nominal", "lp");
  legE0TkQ->AddEntry(grFitE0_RatioToNom_TightSysTkQ, "Tight/Nominal", "lp");

  TLegend * legAlphaTkQ = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legAlphaTkQ);
  legAlphaTkQ->AddEntry(grFitAlpha_RatioToNom_LooseSysTkQ, "Loose/Nominal", "lp");
  legAlphaTkQ->AddEntry(grFitAlpha_RatioToNom_TightSysTkQ, "Tight/Nominal", "lp");

  TCanvas * cSysTkQ = new TCanvas("cSysTkQ", "cSysTkQ", 1500, 500);
  cSysTkQ->Divide(3,1);

  cSysTkQ->cd(1);
  grFitKn_RatioToNom_LooseSysTkQ->Draw("ap");
  grFitKn_RatioToNom_TightSysTkQ->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysTkQ");
  legKnTkQ->Draw("same");

  cSysTkQ->cd(2);
  grFitE0_RatioToNom_LooseSysTkQ->Draw("ap");
  grFitE0_RatioToNom_TightSysTkQ->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysTkQ");
  legE0TkQ->Draw("same");

  cSysTkQ->cd(3);
  grFitAlpha_RatioToNom_LooseSysTkQ->Draw("ap");
  grFitAlpha_RatioToNom_TightSysTkQ->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysTkQ");
  legAlphaTkQ->Draw("same");

  cSysReg->SaveAs("../plots/systematicStudies/cSysFitParms_SysTkQ.pdf");

  //-- SysVtx
  TLegend * legKnVtx = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legKnVtx);
  legKnVtx->AddEntry(grFitKn_RatioToNom_Vtx3SysVtx, "Vtx3/Nominal", "lp");
  legKnVtx->AddEntry(grFitKn_RatioToNom_Vtx3_15SysVtx, "Vtx3_15/Nominal", "lp");

  TLegend * legE0Vtx = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legE0Vtx);
  legE0Vtx->AddEntry(grFitE0_RatioToNom_Vtx3SysVtx, "Vtx3/Nominal", "lp");
  legE0Vtx->AddEntry(grFitE0_RatioToNom_Vtx3_15SysVtx, "Vtx3_15/Nominal", "lp");

  TLegend * legAlphaVtx = new TLegend(0.20, 0.77, 0.64, 0.94);
  legInit(legAlphaVtx);
  legAlphaVtx->AddEntry(grFitAlpha_RatioToNom_Vtx3SysVtx, "Vtx3/Nominal", "lp");
  legAlphaVtx->AddEntry(grFitAlpha_RatioToNom_Vtx3_15SysVtx, "Vtx3_15/Nominal", "lp");

  TCanvas * cSysVtx = new TCanvas("cSysVtx", "cSysVtx", 1500, 500);
  cSysVtx->Divide(3,1);

  cSysVtx->cd(1);
  grFitKn_RatioToNom_Vtx3SysVtx->Draw("ap");
  grFitKn_RatioToNom_Vtx3_15SysVtx->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysVtx");
  legKnVtx->Draw("same");

  cSysVtx->cd(2);
  grFitE0_RatioToNom_Vtx3SysVtx->Draw("ap");
  grFitE0_RatioToNom_Vtx3_15SysVtx->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysVtx");
  legE0Vtx->Draw("same");

  cSysVtx->cd(3);
  grFitAlpha_RatioToNom_Vtx3SysVtx->Draw("ap");
  grFitAlpha_RatioToNom_Vtx3_15SysVtx->Draw("psame");
  lone->Draw("same");
  latex.DrawLatex(0.2, 0.2, "SysVtx");
  legAlphaVtx->Draw("same");

  cSysReg->SaveAs("../plots/systematicStudies/cSysFitParms_SysVtx.pdf");

}

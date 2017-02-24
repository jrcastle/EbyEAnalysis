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
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

int c = 11;

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

void FitPvn(){

  //double knGuess[NCENT] = {16.2, 16.2, 16.2, 0.8308, 0.5066, 0.3792, 0.3383, 0.3235, 0.2939, 0.28, 0.2423, 0.21};
  //double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 511.2, 147.5, 63.9, 40.7, 30.7, 21.2, 20., 10.0, 6.5};
  //double e0Guess[NCENT] = {0.0, 0.0, 0.0, 0.076, 0.1441, 0.2109, 0.2482, 0.2673, 0.2954, 0.3, 0.3303, 0.33};
  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.43,  0.41, 0.40, 0.39, 0.37, 0.35, 0.32, 0.28, 0.23};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 102.0, 70.0, 59.0, 47.0, 34.0, 25.0, 18.0, 11.0, 7.0};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.15,  0.18, 0.20, 0.22, 0.24, 0.25, 0.27, 0.30, 0.35};
  double vnmax[NCENT]   = {0.3,   0.3,   0.3,   0.3,   0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.2,  0.24};

  int norder_  = 2;
  double tkEta = 1.0;


  TLatex latex;
  TLatex latex2;

  bool russianFits = 1;
  TFile * fRussianFits;

  //-- Unfold histos
  TFile * fFinalUnf;
  TH1D * hFinalUnfold[NCENT];
  TH1D * hFinalUnfoldStat[NCENT];
  TH1D * hFinalUnfoldSys[NCENT];

  TF1 * fBG[NCENT];
  TF1 * fEllP[NCENT];

  //-- Parms vs npart
  double fitKn[NCENT];
  double fitAlpha[NCENT];
  double fitE0[NCENT];

  double fitKn_err[NCENT];
  double fitAlpha_err[NCENT];
  double fitE0_err[NCENT];

  TGraphErrors * grFitKn;
  TGraphErrors * grFitAlpha;
  TGraphErrors * grFitE0;

  //-- ATLAS
  TGraphErrors * grATLASKn;
  TGraphErrors * grATLASAlpha;
  TGraphErrors * grATLASEcc0;

  //
  // MAIN
  //

  latex.SetNDC();
  latex2.SetNDC();
  //latex2.SetTextFont(43);
  latex2.SetTextSize(0.08);

  setTDRStyle();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");

  fFinalUnf    = new TFile( Form("systematicStudies/SysUnfoldDistns_v%i.root", norder_) );
  fRussianFits = new TFile("RussianEllpFits.root");

  for(int icent = 0; icent < NCENT; icent++){

    hFinalUnfold[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldStatAndSys_c%i", icent) );
    hFinalUnfold[icent]->SetLineColor(1);
    hFinalUnfold[icent]->SetMarkerColor(1);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->SetFillColor(15);

    hFinalUnfoldStat[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldStat_c%i", icent) );
    hFinalUnfoldStat[icent]->SetLineColor(1);
    hFinalUnfoldStat[icent]->SetMarkerColor(1);
    hFinalUnfoldStat[icent]->SetMarkerStyle(20);

    hFinalUnfoldSys[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldSys_c%i", icent) );
    hFinalUnfoldSys[icent]->SetMarkerColor(1);
    hFinalUnfoldSys[icent]->SetLineColor(1);
    hFinalUnfoldSys[icent]->SetMarkerStyle(20);
    hFinalUnfoldSys[icent]->SetFillColor(15);
    hFinalUnfoldSys[icent]->SetStats(0);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitleSize(0.06);
    hFinalUnfoldSys[icent]->GetXaxis()->SetLabelSize(0.06);
    hFinalUnfoldSys[icent]->GetXaxis()->SetRange(1, 85);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitle("v_{2}");
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitle("p(#epsilon_{2}), p(v_{2})");


    //if( icent != 3 && icent != 5 && icent !=7 ) continue;
    if( icent < 3 ) continue;
    std::cout<<Form("-------- Centbin %i --------", icent)<<std::endl;

    fBG[icent] = new TF1(Form("fBG_c%i", icent), "[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))", 0.0, 0.4);
    fBG[icent]->SetLineColor(7);
    fBG[icent]->SetLineWidth(2);

    EbyECumu cumu(hFinalUnfold[icent]);
    double vn2 = cumu.GetCumu_vn2();
    double vn4 = cumu.GetCumu_vn4();
    double BGmu    = vn4;
    double BGDelta = sqrt( (vn2*vn2 - vn4*vn4)/2. ) ;
    if(icent == NCENT-1){
      BGmu = hFinalUnfold[icent]->GetMean();
      BGDelta = hFinalUnfold[icent]->GetRMS();
    }

    fBG[icent]->SetParameters(1., BGmu, BGDelta);

    if(russianFits) fEllP[icent] = (TF1*) fRussianFits->Get( Form("elliptic%i", icent) );
    else            fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP[icent]->SetLineColor(2);
    fEllP[icent]->SetLineWidth(2);

    if(!russianFits){
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP[icent]->SetParLimits(0, 0., 1.);
      //fEllP[icent]->SetParLimits(2, 0., 10.);
      fEllP[icent]->SetParLimits(2, 0., 1.);
      double ag = 1./2./BGDelta/BGDelta;
      //fEllP[icent]->SetParameters(1., BGmu, ag, knGuess[icent]);
      fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], 1.0);
      fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
    }
    if(russianFits){
      hFinalUnfold[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfold[icent]->GetMaximum());
      hFinalUnfoldStat[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfoldStat[icent]->GetMaximum());
      hFinalUnfoldSys[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfoldSys[icent]->GetMaximum());
    }
    else{
      hFinalUnfold[icent]->Scale(1./hFinalUnfold[icent]->Integral());
      hFinalUnfoldStat[icent]->Scale(1./hFinalUnfoldStat[icent]->Integral());
      hFinalUnfoldSys[icent]->Scale(1./hFinalUnfoldSys[icent]->Integral());
    }
    hFinalUnfold[icent]->Fit( Form("fBG_c%i", icent), "0", "", 0.0, 0.3);
    if(!russianFits) hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "R0", "", 0.0, vnmax[icent]);

    if(russianFits){
      fitKn[icent]    = fEllP[icent]->GetParameter(1);
      fitAlpha[icent] = fEllP[icent]->GetParameter(0);
      fitE0[icent]    = fEllP[icent]->GetParameter(2);

      fitKn_err[icent]    = fEllP[icent]->GetParError(1);
      fitAlpha_err[icent] = fEllP[icent]->GetParError(0);
      fitE0_err[icent]    = fEllP[icent]->GetParError(2);
    }
    else{
      fitKn[icent]    = fEllP[icent]->GetParameter(2);
      fitAlpha[icent] = fEllP[icent]->GetParameter(1);
      fitE0[icent]    = fEllP[icent]->GetParameter(0);

      fitKn_err[icent]    = fEllP[icent]->GetParError(2);
      fitAlpha_err[icent] = fEllP[icent]->GetParError(1);
      fitE0_err[icent]    = fEllP[icent]->GetParError(0);
    }

  } // End cent loop

  //-- Make Graphs
  double c_err[NCENT];
  for(int icent = 0; icent < NCENT; icent++) c_err[icent] = 0;
  /*
  grFitKn    = new TGraphErrors(NCENT, Npart, fitKn,    c_err, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, Npart, fitAlpha, c_err, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, Npart, fitE0,    c_err, fitE0_err);

  formatGraph(grFitKn,    "N_{Part}", 0.1, 1.0,  "Fit k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "N_{Part}", 0.1, 120,  "Fit #alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "N_{Part}", 0.1, 0.5,  "Fit #epsilon_{0}", 4, 34, "grFitE0");
  */
  grFitKn    = new TGraphErrors(NCENT, centBinCenter, fitKn,    c_err, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, centBinCenter, fitAlpha, c_err, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, centBinCenter, fitE0,    c_err, fitE0_err);

  formatGraph(grFitKn,    "Centrality %", 0.1, 1.0,  "k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "Centrality %", 0.1, 120,  "#alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "Centrality %", 0.1, 0.5,  "#epsilon_{0}", 4, 34, "grFitE0");

  //-- ATLAS
  grATLASKn    = new TGraphErrors("k_n_ATLAS_green.txt", "%lg %lg %lg");
  grATLASAlpha = new TGraphErrors("alpha_ATLAS_green.txt", "%lg %lg %lg");
  grATLASEcc0  = new TGraphErrors("ecc0_ATLAS_green.txt", "%lg %lg %lg");

  formatGraph(grATLASKn,    "Centrality %", 0.1, 1.0,  "k_{n}",        1, 24, "grATLASKn");
  formatGraph(grATLASAlpha, "Centrality %", 0.1, 120,  "#alpha",       2, 25, "grATLASAlpha");
  formatGraph(grATLASEcc0,  "Centrality %", 0.1, 0.5,  "#epsilon_{0}", 4, 26, "grATLASE0");

  TLegend * legKn = new TLegend(0.5, 0.75, 0.95, 0.9);
  legKn->SetBorderSize(0);
  legKn->SetFillStyle(0);
  legKn->AddEntry(grFitKn,   "CMS",   "lp");
  legKn->AddEntry(grATLASKn, "ATLAS 2.76 TeV", "lp");

  TLegend * legAlpha = new TLegend(0.5, 0.75, 0.95, 0.9);
  legAlpha->SetBorderSize(0);
  legAlpha->SetFillStyle(0);
  legAlpha->AddEntry(grFitAlpha,   "CMS",   "lp");
  legAlpha->AddEntry(grATLASAlpha, "ATLAS 2.76 TeV", "lp");

  TLegend * legEcc0 = new TLegend(0.2, 0.54, 0.65, 0.69);
  legEcc0->SetBorderSize(0);
  legEcc0->SetFillStyle(0);
  legEcc0->AddEntry(grFitE0,     "CMS",   "lp");
  legEcc0->AddEntry(grATLASEcc0, "ATLAS 2.76 TeV", "lp");

  //-- Fit parms vs npart
  TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 1500, 500);
  cParmSummary->Divide(3,1);

  cParmSummary->cd(1);
  grFitKn->Draw("ap");
  grATLASKn->Draw("psame");
  legKn->Draw("same");

  cParmSummary->cd(2);
  grFitE0->Draw("ap");
  grATLASEcc0->Draw("psame");
  legEcc0->Draw("same");
  latex.DrawLatex(0.2, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.82, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.76, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

  cParmSummary->cd(3);
  grFitAlpha->Draw("ap");
  grATLASAlpha->Draw("psame");
  legAlpha->Draw("same");

  cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_v%i.pdf",norder_));

  //-- Fit Summaries
  TLegend * leg2 = new TLegend(0.1730, 0.2718, 0.8615, 0.4666);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(fEllP[3], "Elliptic Power", "l");
  leg2->AddEntry(fBG[3], "Bessel Gaussian", "l");

  TCanvas * cUnfoldDistsBig = new TCanvas("cUnfoldDistsBig", "cUnfoldDistsBig", 1500, 1000);
  cUnfoldDistsBig->Divide(3,2);

  cUnfoldDistsBig->cd(1);
  cUnfoldDistsBig->cd(1)->SetLogy();
  hFinalUnfoldSys[3]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[3]->Draw("same");
  fBG[3]->Draw("same");
  fEllP[3]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[3], cent_max[3], "%"));
  double alpha3;
  double alpha3e;
  double eps03;
  double eps03e;
  double kn3;
  double kn3e;
  if(russianFits){
    alpha3  = fEllP[3]->GetParameter(0);
    alpha3e = fEllP[3]->GetParError(0);
    eps03   = fEllP[3]->GetParameter(2);
    eps03e  = fEllP[3]->GetParError(2);
    kn3     = fEllP[3]->GetParameter(1);
    kn3e    = fEllP[3]->GetParError(1);
  }
  else{
    alpha3  = fEllP[3]->GetParameter(1);
    alpha3e = fEllP[3]->GetParError(1);
    eps03   = fEllP[3]->GetParameter(0);
    eps03e  = fEllP[3]->GetParError(0);
    kn3     = fEllP[3]->GetParameter(2);
    kn3e    = fEllP[3]->GetParError(2);
  }
  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f", eps03,  eps03e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha3, alpha3e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",       kn3,    kn3e));

  cUnfoldDistsBig->cd(1)->SetTopMargin(0.15);
  cUnfoldDistsBig->cd(1)->SetTickx(0);
  double xmin3 = 0;
  double xmax3 = 0.34;
  double ymax3 = 1.9*hFinalUnfold[3]->GetMaximum();

  //-- POOP!

  TGaxis * axEcc3 = new TGaxis(xmin3, ymax3, xmax3, ymax3, xmin3/kn3, xmax3/kn3, 509, "-");
  axEcc3->SetLabelSize(0.055);
  axEcc3->SetTitleSize(0.055);
  axEcc3->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(1);
  axEcc3->Draw("same");

  cUnfoldDistsBig->cd(2);
  cUnfoldDistsBig->cd(2)->SetLogy();
  hFinalUnfoldSys[5]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[5]->Draw("same");
  fBG[5]->Draw("same");
  fEllP[5]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[5], cent_max[5], "%"));
  double alpha5;
  double alpha5e;
  double eps05;
  double eps05e;
  double kn5;
  double kn5e;
  if(russianFits){
    alpha5  = fEllP[5]->GetParameter(0);
    alpha5e = fEllP[5]->GetParError(0);
    eps05   = fEllP[5]->GetParameter(2);
    eps05e  = fEllP[5]->GetParError(2);
    kn5     = fEllP[5]->GetParameter(1);
    kn5e    = fEllP[5]->GetParError(1);
  }
  else{
    alpha5  = fEllP[5]->GetParameter(1);
    alpha5e = fEllP[5]->GetParError(1);
    eps05   = fEllP[5]->GetParameter(0);
    eps05e  = fEllP[5]->GetParError(0);
    kn5     = fEllP[5]->GetParameter(2);
    kn5e    = fEllP[5]->GetParError(2);
  }
  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f", eps05,  eps05e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha5, alpha5e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",       kn5,    kn5e));

  cUnfoldDistsBig->cd(2)->SetTickx(0);
  cUnfoldDistsBig->cd(2)->SetTopMargin(0.15);
  double xmin5 = 0;
  double xmax5 = 0.34;
  double ymax5 = 1.9*hFinalUnfold[5]->GetMaximum();

  TGaxis * axEcc5 = new TGaxis(xmin5, ymax5, xmax5, ymax5, xmin5/kn5, xmax5/kn5, 509, "-");
  axEcc5->SetTitle("#epsilon_{2}");
  axEcc5->SetLabelSize(0.055);
  axEcc5->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(2);
  axEcc5->Draw("same");


  cUnfoldDistsBig->cd(3);
  cUnfoldDistsBig->cd(3)->SetLogy();
  hFinalUnfoldSys[7]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[7]->Draw("same");
  fBG[7]->Draw("same");
  fEllP[7]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[7], cent_max[7], "%"));
  double alpha7;
  double alpha7e;
  double eps07;
  double eps07e;
  double kn7;
  double kn7e;
  if(russianFits){
    alpha7  = fEllP[7]->GetParameter(0);
    alpha7e = fEllP[7]->GetParError(0);
    eps07   = fEllP[7]->GetParameter(2);
    eps07e  = fEllP[7]->GetParError(2);
    kn7     = fEllP[7]->GetParameter(1);
    kn7e    = fEllP[7]->GetParError(1);
  }
  else{
    alpha7  = fEllP[7]->GetParameter(1);
    alpha7e = fEllP[7]->GetParError(1);
    eps07   = fEllP[7]->GetParameter(0);
    eps07e  = fEllP[7]->GetParError(0);
    kn7     = fEllP[7]->GetParameter(2);
    kn7e    = fEllP[7]->GetParError(2);
  }
  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps07,  eps07e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha7, alpha7e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn7,    kn7e));

  cUnfoldDistsBig->cd(3)->SetTickx(0);
  cUnfoldDistsBig->cd(3)->SetTopMargin(0.15);
  double xmin7 = 0;
  double xmax7 = 0.34;
  double ymax7 = 1.9*hFinalUnfold[7]->GetMaximum();

  TGaxis * axEcc7 = new TGaxis(xmin7, ymax7, xmax7, ymax7, xmin7/kn7, xmax7/kn7, 509, "-");
  axEcc7->SetTitle("#epsilon_{2}");
  axEcc7->SetLabelSize(0.055);
  axEcc7->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(3);
  axEcc7->Draw("same");

  cUnfoldDistsBig->cd(4);
  cUnfoldDistsBig->cd(4)->SetLogy();
  hFinalUnfoldSys[9]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[9]->Draw("same");
  fBG[9]->Draw("same");
  fEllP[9]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[9], cent_max[9], "%"));
  double alpha9;
  double alpha9e;
  double eps09;
  double eps09e;
  double kn9;
  double kn9e;
  if(russianFits){
    alpha9  = fEllP[9]->GetParameter(0);
    alpha9e = fEllP[9]->GetParError(0);
    eps09   = fEllP[9]->GetParameter(2);
    eps09e  = fEllP[9]->GetParError(2);
    kn9     = fEllP[9]->GetParameter(1);
    kn9e    = fEllP[9]->GetParError(1);
  }
  else{
    alpha9  = fEllP[9]->GetParameter(1);
    alpha9e = fEllP[9]->GetParError(1);
    eps09   = fEllP[9]->GetParameter(0);
    eps09e  = fEllP[9]->GetParError(0);
    kn9     = fEllP[9]->GetParameter(2);
    kn9e    = fEllP[9]->GetParError(2);
  }
  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps09,  eps09e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha9, alpha9e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn9,    kn9e));

  cUnfoldDistsBig->cd(4)->SetTickx(0);
  cUnfoldDistsBig->cd(4)->SetTopMargin(0.15);
  double xmin9 = 0;
  double xmax9 = 0.34;
  double ymax9 = 1.9*hFinalUnfold[9]->GetMaximum();

  TGaxis * axEcc9 = new TGaxis(xmin9, ymax9, xmax9, ymax9, xmin9/kn9, xmax9/kn9, 509, "-");
  axEcc9->SetTitle("#epsilon_{2}");
  axEcc9->SetLabelSize(0.055);
  axEcc9->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(4);
  axEcc9->Draw("same");

  cUnfoldDistsBig->cd(5);
  cUnfoldDistsBig->cd(5)->SetLogy();
  hFinalUnfoldSys[11]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[11]->Draw("same");
  fBG[11]->Draw("same");
  fEllP[11]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[11], cent_max[11], "%"));
  double alpha11;
  double alpha11e;
  double eps011;
  double eps011e;
  double kn11;
  double kn11e;
  if(russianFits){
    alpha11  = fEllP[11]->GetParameter(0);
    alpha11e = fEllP[11]->GetParError(0);
    eps011   = fEllP[11]->GetParameter(2);
    eps011e  = fEllP[11]->GetParError(2);
    kn11     = fEllP[11]->GetParameter(1);
    kn11e    = fEllP[11]->GetParError(1);
  }
  else{
    alpha11  = fEllP[11]->GetParameter(1);
    alpha11e = fEllP[11]->GetParError(1);
    eps011   = fEllP[11]->GetParameter(0);
    eps011e  = fEllP[11]->GetParError(0);
    kn11     = fEllP[11]->GetParameter(2);
    kn11e    = fEllP[11]->GetParError(2);
  }
  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps011,  eps011e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha11, alpha11e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn11,    kn11e));

  cUnfoldDistsBig->cd(5)->SetTickx(0);
  cUnfoldDistsBig->cd(5)->SetTopMargin(0.15);
  double xmin11 = 0;
  double xmax11 = 0.34;
  double ymax11 = 1.9*hFinalUnfold[11]->GetMaximum();

  TGaxis * axEcc11 = new TGaxis(xmin11, ymax11, xmax11, ymax11, xmin11/kn11, xmax11/kn11, 509, "-");
  axEcc11->SetTitle("#epsilon_{2}");
  axEcc11->SetLabelSize(0.055);
  axEcc11->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(5);
  axEcc11->Draw("same");

  cUnfoldDistsBig->cd(6);
  latex2.DrawLatex(0.2, 0.80, "CMS #it{Preliminary}");
  latex2.DrawLatex(0.2, 0.70, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2.DrawLatex(0.2, 0.60, Form("|#eta| < %.1f", tkEta));
  latex2.DrawLatex(0.2, 0.50, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  leg2->Draw("same");

  cUnfoldDistsBig->Update();
  cUnfoldDistsBig->SaveAs(Form("plots/unfolding/finalUnfFit_v%i.pdf",norder_));

}

#include "TPaveStats.h"
#include "TMath.h"
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
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/ATLAS_PV2.h"

#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;
using namespace atlas_pv2;

int BIN = 8;
/*
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

  double p1 = 2.*eccn*alpha;
  double p2 = pow(1-eccn*eccn, alpha-1.);
  double p3 = pow(1-eccn*e0, -1.-2.*alpha);
  double p4 = pow(1-e0*e0, alpha+0.5);
  double p5 = ROOT::Math::hyperg( 0.5, 1+2.*alpha, 1., 2.*eccn*e0/(eccn*e0-1) );

  double ellp = scale * p1 * p2 * p3 * p4 * p5;
  return ellp;

}
*/
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

void atlasPV2Fits(){

  //double knGuess[NCENT] = {16.2, 16.2, 16.2, 0.8308, 0.5066, 0.3792, 0.3383, 0.3235, 0.2939, 0.28, 0.2423, 0.21};
  //double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 511.2, 147.5, 63.9, 40.7, 30.7, 21.2, 20., 10.0, 6.5};
  //double e0Guess[NCENT] = {0.0, 0.0, 0.0, 0.076, 0.1441, 0.2109, 0.2482, 0.2673, 0.2954, 0.3, 0.3303, 0.33};
  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.441, 0.394, 0.392, 0.364, 0.352, 0.344, 0.315, 0.280, 0.260};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 91.1,  56.3,  45.0,  30.2,  23.3,  18.6,  12.7,  8.1,   6.1};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.179, 0.227, 0.251, 0.287, 0.307, 0.314, 0.333, 0.0,   0.0};
  double vnmax[NCENT]   = {0.3,   0.3,   0.3,   0.3,   0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3};
  
  int norder_  = 2;
  double tkEta = 1.0;


  TLatex latex;
  TLatex latex2;

  bool dosys     = 1;
  bool ATLAS     = 1;

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

  TPaveStats * s[NCENT];

  //
  // MAIN
  //

  latex.SetNDC();
  latex2.SetNDC();
  //latex2.SetTextFont(43);
  latex2.SetTextSize(0.08);

  setTDRStyle();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");

  for(int icent = 0; icent < NCENT; icent++){

    fitKn[icent]    = -1;
    fitAlpha[icent] = -1;
    fitE0[icent]    = -1;

    if( icent < 3 || icent > 9 ) continue;
    std::cout<<Form("-------- Centbin %i --------", icent)<<std::endl;

    fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP[icent]->SetLineColor(2);
    fEllP[icent]->SetLineWidth(2);

    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP[icent]->SetParLimits(0, 0., 1.);
    fEllP[icent]->SetParLimits(2, 0., 1.);
    fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], TMath::MaxElement(NATLAS[icent], ATLAS_PV2[icent]->GetY()));
    fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
    ATLAS_PV2[icent]->Fit( Form("fEllP_c%i", icent), "L0", "", 0.0, vnmax[icent]);

    fitKn[icent]    = fEllP[icent]->GetParameter(2);
    fitAlpha[icent] = fEllP[icent]->GetParameter(1);
    fitE0[icent]    = fEllP[icent]->GetParameter(0);

    fitKn_err[icent]    = fEllP[icent]->GetParError(2);
    fitAlpha_err[icent] = fEllP[icent]->GetParError(1);
    fitE0_err[icent]    = fEllP[icent]->GetParError(0);

    ATLAS_PV2[icent]->GetXaxis()->SetNdivisions(508);
    ATLAS_PV2[icent]->GetXaxis()->SetTitle("v_{2}");
    ATLAS_PV2[icent]->GetYaxis()->SetTitle("p(v_{2}), p(#epsilon_{2})");

  } // End cent loop

  //-- Make Graphs
  grFitKn    = new TGraphErrors(NCENT, centBinCenter, fitKn,    CERR, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, centBinCenter, fitAlpha, CERR, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, centBinCenter, fitE0,    CERR, fitE0_err);

  formatGraph(grFitKn,    "Centrality %", 0., 0.75,  "k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "Centrality %", 0., 120,  "#alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "Centrality %", 0., 0.4,  "#epsilon_{0}", 4, 34, "grFitE0");

  //-- ATLAS
  grATLASKn    = new TGraphErrors("k_n_ATLAS_green.txt", "%lg %lg %lg");
  grATLASAlpha = new TGraphErrors("alpha_ATLAS_green.txt", "%lg %lg %lg");
  grATLASEcc0  = new TGraphErrors("ecc0_ATLAS_green.txt", "%lg %lg %lg");

  formatGraph(grATLASKn,    "Centrality %", 0.1, 1.0,  "k_{n}",        1, 24, "grATLASKn");
  formatGraph(grATLASAlpha, "Centrality %", 0.1, 120,  "#alpha",       2, 25, "grATLASAlpha");
  formatGraph(grATLASEcc0,  "Centrality %", 0.1, 0.5,  "#epsilon_{0}", 4, 26, "grATLASE0");


  TLegend * legKn = new TLegend(0.2, 0.2, 0.65, 0.35);
  legKn->SetBorderSize(0);
  legKn->SetFillStyle(0);
  legKn->AddEntry(grFitKn,   "My ATLAS 2.76 TeV",   "lp");
  legKn->AddEntry(grATLASKn, "Theorists's ATLAS 2.76 TeV", "lp");

  TLegend * legAlpha = new TLegend(0.5, 0.75, 0.95, 0.9);
  legAlpha->SetBorderSize(0);
  legAlpha->SetFillStyle(0);
  legAlpha->AddEntry(grFitAlpha,   "My ATLAS 2.76 TeV",   "lp");
  legAlpha->AddEntry(grATLASAlpha, "Theorists's ATLAS 2.76 TeV", "lp");

  TLegend * legEcc0 = new TLegend(0.2, 0.2, 0.65, 0.35);
  legEcc0->SetBorderSize(0);
  legEcc0->SetFillStyle(0);
  legEcc0->AddEntry(grFitE0,     "My ATLAS 2.76 TeV",   "lp");
  legEcc0->AddEntry(grATLASEcc0, "Theorists's ATLAS 2.76 TeV", "lp");

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
    
  cParmSummary->cd(3);
  grFitAlpha->Draw("ap");
  grATLASAlpha->Draw("psame");
  legAlpha->Draw("same");

  cParmSummary->SaveAs(Form("FitParmSummary_v%i.pdf",norder_));

  TCanvas * cUnfoldDistsBig = new TCanvas("cUnfoldDistsBig", "cUnfoldDistsBig", 1500, 1000);
  cUnfoldDistsBig->Divide(3,2);

  cUnfoldDistsBig->cd(1);
  cUnfoldDistsBig->cd(1)->SetLogy();
  ATLAS_PV2[3]->Draw("ap");
  fEllP[3]->Draw("same");
  cUnfoldDistsBig->Update();
  ATLAS_PV2[3]->GetXaxis()->SetLimits(0, vnmax[3]);
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[3], cent_max[3], "%"));
  double alpha3;
  double alpha3e;
  double eps03;
  double eps03e;
  double kn3;
  double kn3e;
  alpha3  = fEllP[3]->GetParameter(1);
  alpha3e = fEllP[3]->GetParError(1);
  eps03   = fEllP[3]->GetParameter(0);
  eps03e  = fEllP[3]->GetParError(0);
  kn3     = fEllP[3]->GetParameter(2);
  kn3e    = fEllP[3]->GetParError(2);

  //latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f", eps03,  eps03e));
  //latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha3, alpha3e));
  //latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",       kn3,    kn3e));

  cUnfoldDistsBig->cd(1)->SetTopMargin(0.15);
  cUnfoldDistsBig->cd(1)->SetTickx(0);
  double xmin3 = 0;
  double xmax3 = vnmax[3];
  double ymax3 = 1.1*TMath::MaxElement(NATLAS[3], ATLAS_PV2[3]->GetY());

  TGaxis * axEcc3 = new TGaxis(xmin3, ymax3, xmax3, ymax3, xmin3/kn3, xmax3/kn3, 509, "-");
  axEcc3->SetLabelSize(0.055);
  axEcc3->SetTitleSize(0.055);
  axEcc3->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(1);
  axEcc3->Draw("same");

  cUnfoldDistsBig->cd(2);
  cUnfoldDistsBig->cd(2)->SetLogy();
  ATLAS_PV2[5]->Draw("ap");
  fEllP[5]->Draw("same");
  ATLAS_PV2[5]->GetXaxis()->SetLimits(0, vnmax[5]);
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[5], cent_max[5], "%"));
  double alpha5;
  double alpha5e;
  double eps05;
  double eps05e;
  double kn5;
  double kn5e;
  alpha5  = fEllP[5]->GetParameter(1);
  alpha5e = fEllP[5]->GetParError(1);
  eps05   = fEllP[5]->GetParameter(0);
  eps05e  = fEllP[5]->GetParError(0);
  kn5     = fEllP[5]->GetParameter(2);
  kn5e    = fEllP[5]->GetParError(2);

  cUnfoldDistsBig->cd(2)->SetTickx(0);
  cUnfoldDistsBig->cd(2)->SetTopMargin(0.15);
  double xmin5 = 0;
  double xmax5 = vnmax[5];
  double ymax5 = 1.1*TMath::MaxElement(NATLAS[5], ATLAS_PV2[5]->GetY());

  TGaxis * axEcc5 = new TGaxis(xmin5, ymax5, xmax5, ymax5, xmin5/kn5, xmax5/kn5, 509, "-");
  axEcc5->SetTitle("#epsilon_{2}");
  axEcc5->SetLabelSize(0.055);
  axEcc5->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(2);
  axEcc5->Draw("same");


  cUnfoldDistsBig->cd(3);
  cUnfoldDistsBig->cd(3)->SetLogy();
  ATLAS_PV2[7]->Draw("ap");
  fEllP[7]->Draw("same");
  ATLAS_PV2[7]->GetXaxis()->SetLimits(0, vnmax[7]);
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[7], cent_max[7], "%"));
  double alpha7;
  double alpha7e;
  double eps07;
  double eps07e;
  double kn7;
  double kn7e;
  alpha7  = fEllP[7]->GetParameter(1);
  alpha7e = fEllP[7]->GetParError(1);
  eps07   = fEllP[7]->GetParameter(0);
  eps07e  = fEllP[7]->GetParError(0);
  kn7     = fEllP[7]->GetParameter(2);
  kn7e    = fEllP[7]->GetParError(2);

  cUnfoldDistsBig->cd(3)->SetTickx(0);
  cUnfoldDistsBig->cd(3)->SetTopMargin(0.15);
  double xmin7 = 0;
  double xmax7 = vnmax[7];
  double ymax7 = 1.1*TMath::MaxElement(NATLAS[7], ATLAS_PV2[7]->GetY());

  TGaxis * axEcc7 = new TGaxis(xmin7, ymax7, xmax7, ymax7, xmin7/kn7, xmax7/kn7, 509, "-");
  axEcc7->SetTitle("#epsilon_{2}");
  axEcc7->SetLabelSize(0.055);
  axEcc7->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(3);
  axEcc7->Draw("same");

  cUnfoldDistsBig->cd(4);
  cUnfoldDistsBig->cd(4)->SetLogy();
  ATLAS_PV2[9]->Draw("ap");
  fEllP[9]->Draw("same");
  ATLAS_PV2[9]->GetXaxis()->SetLimits(0, vnmax[9]);
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[9], cent_max[9], "%"));
  double alpha9;
  double alpha9e;
  double eps09;
  double eps09e;
  double kn9;
  double kn9e;
  alpha9  = fEllP[9]->GetParameter(1);
  alpha9e = fEllP[9]->GetParError(1);
  eps09   = fEllP[9]->GetParameter(0);
  eps09e  = fEllP[9]->GetParError(0);
  kn9     = fEllP[9]->GetParameter(2);
  kn9e    = fEllP[9]->GetParError(2);

  cUnfoldDistsBig->cd(4)->SetTickx(0);
  cUnfoldDistsBig->cd(4)->SetTopMargin(0.15);
  double xmin9 = 0;
  double xmax9 = vnmax[9];
  double ymax9 = 1.1*TMath::MaxElement(NATLAS[9], ATLAS_PV2[9]->GetY());

  TGaxis * axEcc9 = new TGaxis(xmin9, ymax9, xmax9, ymax9, xmin9/kn9, xmax9/kn9, 509, "-");
  axEcc9->SetTitle("#epsilon_{2}");
  axEcc9->SetLabelSize(0.055);
  axEcc9->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(4);
  axEcc9->Draw("same");
  /*
  cUnfoldDistsBig->cd(5);
  cUnfoldDistsBig->cd(5)->SetLogy();
  hFinalUnfoldSys[11]->Draw("e2");
  setex2->Draw();
  ATLAS_PV2[11]->Draw("same");
  fEllP[11]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[11], cent_max[11], "%"));
  double alpha11;
  double alpha11e;
  double eps011;
  double eps011e;
  double kn11;
  double kn11e;
  alpha11  = fEllP[11]->GetParameter(1);
  alpha11e = fEllP[11]->GetParError(1);
  eps011   = fEllP[11]->GetParameter(0);
  eps011e  = fEllP[11]->GetParError(0);
  kn11     = fEllP[11]->GetParameter(2);
  kn11e    = fEllP[11]->GetParError(2);

  cUnfoldDistsBig->cd(5)->SetTickx(0);
  cUnfoldDistsBig->cd(5)->SetTopMargin(0.15);
  double xmin11 = 0;
  double xmax11 = vnmax[11]-binw;
  double ymax11 = 1.9*hFinalUnfold[11]->GetMaximum();

  TGaxis * axEcc11 = new TGaxis(xmin11, ymax11, xmax11, ymax11, xmin11/kn11, xmax11/kn11, 509, "-");
  axEcc11->SetTitle("#epsilon_{2}");
  axEcc11->SetLabelSize(0.055);
  axEcc11->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(5);
  axEcc11->Draw("same");
  */
  cUnfoldDistsBig->Update();
  cUnfoldDistsBig->SaveAs(Form("finalUnfFit_v%i.pdf",norder_));

}

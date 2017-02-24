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

  double knGuess[NCENT] = {16.2, 16.2, 16.2, 0.8308, 0.5066, 0.3792, 0.3383, 0.3235, 0.2939, 0.28, 0.2423, 0.21};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 511.2, 147.5, 63.9, 40.7, 30.7, 21.2, 20., 10.0, 6.5};
  double e0Guess[NCENT] = {0.0, 0.0, 0.0, 0.076, 0.1441, 0.2109, 0.2482, 0.2673, 0.2954, 0.3, 0.3303, 0.33};
  double vnmax[NCENT] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.25, 0.2, 0.2};

  int norder_   = 2;
  double tketa_ = 2.4;


  TLatex latex;
  TLatex latex2;

  //-- Unfold histos
  TFile * fFinalUnf;
  TH1D * hFinalUnfold[NCENT];
  TH1D * hFinalUnfoldStat[NCENT];

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


  //
  // MAIN
  //

  latex.SetNDC();
  latex2.SetNDC();
  //latex2.SetTextFont(43);
  latex2.SetTextSize(0.08);

  setTDRStyle();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");

  fFinalUnf = new TFile( Form("systematicStudies/SysUnfoldDistns_v%i.root", norder_) );

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
    fBG[icent]->SetParameters(1., BGmu, BGDelta);

    fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP[icent]->SetLineColor(2);
    fEllP[icent]->SetLineWidth(2);
    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP[icent]->SetParLimits(0, 0., 1.);
    //fEllP[icent]->SetParLimits(2, 0., 10.);
    fEllP[icent]->SetParLimits(2, 0., 1.);
    double ag = 1./2./BGDelta/BGDelta;
    //fEllP[icent]->SetParameters(1., BGmu, ag, knGuess[icent]);
    fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], 1.0);
    fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");

    hFinalUnfold[icent]->Scale(1./hFinalUnfold[icent]->Integral());
    hFinalUnfold[icent]->SetMarkerColor(1);
    hFinalUnfold[icent]->SetLineColor(1);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->SetStats(0);
    hFinalUnfold[icent]->GetXaxis()->SetRange(1, 85);
    hFinalUnfold[icent]->GetXaxis()->SetTitle("v_{2}");
    hFinalUnfold[icent]->GetYaxis()->SetTitle("p(#epsilon_{2}), p(v_{2})");
    hFinalUnfold[icent]->Fit( Form("fBG_c%i", icent), "0", "", 0.0, 0.3);
    //hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "B0", "", 0.0, 0.3);
    hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "R0", "", 0.0, vnmax[icent]);

    fitKn[icent]    = fEllP[icent]->GetParameter(2);
    fitAlpha[icent] = fEllP[icent]->GetParameter(1);
    fitE0[icent]    = fEllP[icent]->GetParameter(0);

    fitKn_err[icent]    = fEllP[icent]->GetParError(2);
    fitAlpha_err[icent] = fEllP[icent]->GetParError(1);
    fitE0_err[icent]    = fEllP[icent]->GetParError(0);

  } // End cent loop

  //-- Make Graphs
  double c_err[NCENT];
  for(int icent = 0; icent < NCENT; icent++) c_err[icent] = 0;

  grFitKn    = new TGraphErrors(NCENT, Npart, fitKn,    c_err, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, Npart, fitAlpha, c_err, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, Npart, fitE0,    c_err, fitE0_err);

  formatGraph(grFitKn,    "N_{Part}", 0.1, 1.0,  "Fit k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "N_{Part}", 0.1, 700,  "Fit #alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "N_{Part}", 0.1, 0.5,  "Fit #epsilon_{0}", 4, 34, "grFitE0");


  //-- Fit parms vs npart
  TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 1500, 500);
  cParmSummary->Divide(3,1);
  cParmSummary->cd(1);
  grFitKn->Draw("ap");
  cParmSummary->cd(2);
  grFitAlpha->Draw("ap");
  cParmSummary->cd(3);
  grFitE0->Draw("ap");
  latex.DrawLatex(0.45, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.45, 0.82, Form("|#eta| < %.1f", tketa_));
  latex.DrawLatex(0.45, 0.76, Form("%.2f < p_{T} < %.2f GeV/c", pt_min[0], pt_max[NPT-1]));
  cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_v%i.pdf",norder_));


  //-- Fit Summaries
  TLegend * leg2 = new TLegend(0.1730, 0.3718, 0.8615, 0.5666);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(fEllP[3], "Elliptic Power", "l");
  leg2->AddEntry(fBG[3], "Bessel Gaussian", "l");

  TCanvas * cUnfoldDistsBig = new TCanvas("cUnfoldDistsBig", "cUnfoldDistsBig", 1500, 1000);
  cUnfoldDistsBig->Divide(3,2);

  cUnfoldDistsBig->cd(1);
  cUnfoldDistsBig->cd(1)->SetLogy();
  hFinalUnfold[3]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[3]->Draw("same");
  fBG[3]->Draw("same");
  fEllP[3]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[3], cent_max[3], "%"));
  double alpha3  = fEllP[3]->GetParameter(1);
  double alpha3e = fEllP[3]->GetParError(1);
  double eps03   = fEllP[3]->GetParameter(0);
  double eps03e  = fEllP[3]->GetParError(0);
  double kn3     = fEllP[3]->GetParameter(2);
  double kn3e    = fEllP[3]->GetParError(2);
  latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f", eps03,  eps03e));
  latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha3, alpha3e));
  latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",       kn3,    kn3e));

  cUnfoldDistsBig->cd(1)->SetTickx(0);
  cUnfoldDistsBig->cd(1)->SetTopMargin(0.1);
  double xmin3 = 0;
  double xmax3 = 0.34;
  double ymax3 = 1.9*hFinalUnfold[3]->GetMaximum();

  TGaxis * axEcc3 = new TGaxis(xmin3, ymax3, xmax3, ymax3, xmin3/kn3, xmax3/kn3, 509, "-");
  axEcc3->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(1);
  axEcc3->Draw("same");

  cUnfoldDistsBig->cd(2);
  cUnfoldDistsBig->cd(2)->SetLogy();
  hFinalUnfold[5]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[5]->Draw("same");
  fBG[5]->Draw("same");
  fEllP[5]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[5], cent_max[5], "%"));
  double alpha5  = fEllP[5]->GetParameter(1);
  double alpha5e = fEllP[5]->GetParError(1);
  double eps05   = fEllP[5]->GetParameter(0);
  double eps05e  = fEllP[5]->GetParError(0);
  double kn5     = fEllP[5]->GetParameter(2);
  double kn5e    = fEllP[5]->GetParError(2);
  latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f", eps05,  eps05e));
  latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha5, alpha5e));
  latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",       kn5,    kn5e));

  cUnfoldDistsBig->cd(2)->SetTickx(0);
  cUnfoldDistsBig->cd(2)->SetTopMargin(0.1);
  double xmin5 = 0;
  double xmax5 = 0.34;
  double ymax5 = 1.9*hFinalUnfold[5]->GetMaximum();

  TGaxis * axEcc5 = new TGaxis(xmin5, ymax5, xmax5, ymax5, xmin5/kn5, xmax5/kn5, 509, "-");
  axEcc5->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(2);
  axEcc5->Draw("same");


  cUnfoldDistsBig->cd(3);
  cUnfoldDistsBig->cd(3)->SetLogy();
  hFinalUnfold[7]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[7]->Draw("same");
  fBG[7]->Draw("same");
  fEllP[7]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[7], cent_max[7], "%"));
  double alpha7  = fEllP[7]->GetParameter(1);
  double alpha7e = fEllP[7]->GetParError(1);
  double eps07   = fEllP[7]->GetParameter(0);
  double eps07e  = fEllP[7]->GetParError(0);
  double kn7     = fEllP[7]->GetParameter(2);
  double kn7e    = fEllP[7]->GetParError(2);
  latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps07,  eps07e));
  latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha7, alpha7e));
  latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn7,    kn7e));

  cUnfoldDistsBig->cd(3)->SetTickx(0);
  cUnfoldDistsBig->cd(3)->SetTopMargin(0.1);
  double xmin7 = 0;
  double xmax7 = 0.34;
  double ymax7 = 1.9*hFinalUnfold[7]->GetMaximum();

  TGaxis * axEcc7 = new TGaxis(xmin7, ymax7, xmax7, ymax7, xmin7/kn7, xmax7/kn7, 509, "-");
  axEcc7->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(3);
  axEcc7->Draw("same");

  cUnfoldDistsBig->cd(4);
  cUnfoldDistsBig->cd(4)->SetLogy();
  hFinalUnfold[9]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[9]->Draw("same");
  fBG[9]->Draw("same");
  fEllP[9]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[9], cent_max[9], "%"));
  double alpha9  = fEllP[9]->GetParameter(1);
  double alpha9e = fEllP[9]->GetParError(1);
  double eps09   = fEllP[9]->GetParameter(0);
  double eps09e  = fEllP[9]->GetParError(0);
  double kn9     = fEllP[9]->GetParameter(2);
  double kn9e    = fEllP[9]->GetParError(2);
  latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps09,  eps09e));
  latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha9, alpha9e));
  latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn9,    kn9e));

  cUnfoldDistsBig->cd(4)->SetTickx(0);
  cUnfoldDistsBig->cd(4)->SetTopMargin(0.1);
  double xmin9 = 0;
  double xmax9 = 0.34;
  double ymax9 = 1.9*hFinalUnfold[9]->GetMaximum();

  TGaxis * axEcc9 = new TGaxis(xmin9, ymax9, xmax9, ymax9, xmin9/kn9, xmax9/kn9, 509, "-");
  axEcc9->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(4);
  axEcc9->Draw("same");

  cUnfoldDistsBig->cd(5);
  cUnfoldDistsBig->cd(5)->SetLogy();
  hFinalUnfold[11]->Draw("e2");
  setex2->Draw();
  hFinalUnfoldStat[11]->Draw("same");
  fBG[11]->Draw("same");
  fEllP[11]->Draw("same");
  latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[11], cent_max[11], "%"));
  double alpha11  = fEllP[11]->GetParameter(1);
  double alpha11e = fEllP[11]->GetParError(1);
  double eps011   = fEllP[11]->GetParameter(0);
  double eps011e  = fEllP[11]->GetParError(0);
  double kn11     = fEllP[11]->GetParameter(2);
  double kn11e    = fEllP[11]->GetParError(2);
  latex.DrawLatex(0.2, 0.38, Form("#epsilon_{0} = %.4f #pm %.4f",  eps011,  eps011e));
  latex.DrawLatex(0.2, 0.32,  Form("#alpha = %.1f #pm %.1f",       alpha11, alpha11e));
  latex.DrawLatex(0.2, 0.26,  Form("k_{n} = %.3f #pm %.3f",        kn11,    kn11e));

  cUnfoldDistsBig->cd(5)->SetTickx(0);
  cUnfoldDistsBig->cd(5)->SetTopMargin(0.1);
  double xmin11 = 0;
  double xmax11 = 0.34;
  double ymax11 = 1.9*hFinalUnfold[11]->GetMaximum();

  TGaxis * axEcc11 = new TGaxis(xmin11, ymax11, xmax11, ymax11, xmin11/kn11, xmax11/kn11, 509, "-");
  axEcc11->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(5);
  axEcc11->Draw("same");

  cUnfoldDistsBig->cd(6);
  latex2.DrawLatex(0.2, 0.78, "CMS #it{Preliminary}");
  latex2.DrawLatex(0.2, 0.69, Form("|#eta| < %.1f", tketa_));
  latex2.DrawLatex(0.2, 0.60, Form("%.2f < p_{T} < %.2f GeV/c", pt_min[0], pt_max[NPT-1]));
  leg2->Draw("same");

  cUnfoldDistsBig->Update();
  cUnfoldDistsBig->SaveAs(Form("plots/unfolding/finalUnfFit_v%i.pdf",norder_));

}

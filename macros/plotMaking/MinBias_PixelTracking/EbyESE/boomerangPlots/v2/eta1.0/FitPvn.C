#include "TMinuit.h"
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

int BIN = 1;
string fname = "EllPFits.root";
/*
double pEllP(double * x, double * par){

  //-- [0] = e0
  //-- [1] = alpha
  //-- [2] = kn
  //-- [3] = Scale 
  //--  x  = vn

  double e0    = par[0];
  double alpha = par[1];
  double kn    = par[2];
  double scale = par[3];
  double eccn  = x[0] / kn;

  double p1 = 2.*eccn*alpha;
  double p2 = pow(1-eccn*eccn, alpha-1.);
  double p3 = pow(1-eccn*e0, -1.-2.*alpha);
  double p4 = pow(1-e0*e0, alpha+0.5);
  double p5 = ROOT::Math::hyperg( 0.5, 1+2.*alpha, 1., 2.*eccn*e0/(eccn*e0-1.) );

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
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 1; i <= nbin; i++){
    double phi = (2*(double)i-1.)*pi/(2.*(double)nbin);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  double ellp = p1 * integ;

  return ellp;

}

void FitPvn(){

  bool dosys      = 0;
  bool ATLAS      = 0;
  bool fixKn      = 0;
  bool fixAlpha   = 0;
  bool moscowFits = 0;
  bool contours   = 1;
  // kn 0.40
  // al 71.1
  // e0 .17
  //-- Free kn
  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.48,    0.38, 0.37, 0.37, 0.34, 0.32, 0.29, 0.28, 0.28};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 90.0,    63.0, 47.0, 41.0, 30.0, 20.0, 12.0, 11.0, 8.0};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.13,    0.19, 0.22, 0.23, 0.26, 0.29, 0.30, 0.31, 0.30};
  //-- ATLAS KN
  //double knGuess[NCENT] = {16.2,  16.2,  16.2,     0.441, 0.394, 0.392, 0.364, 0.352, 0.344, 0.315, 0.280, 0.260};
  //double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5,    44.0,  35.1,  28.3,  23.4,  19.1,  16.6,  14.6,  13.6,  13.8};
  //double e0Guess[NCENT] = {0.0,   0.0,   0.0,      0.215, 0.248, 0.272, 0.288, 0.294, 0.296, 0.290, 0.281, 0.271};
  //-- kn = 0.47
  //double knGuess[NCENT] = {16.2,  16.2,  16.2,     0.47,  0.47,  0.47,  0.47,  0.47,  0.47,  0.47,  0.47,  0.47};
  //double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5,    119.5, 98.4,  81.7,  69.6,  58.2,  51.6,  46.3,  43.5,  43.2};
  //double e0Guess[NCENT] = {0.0,   0.0,   0.0,      0.138, 0.159, 0.174, 0.185, 0.190, 0.192, 0.190, 0.186, 0.179};
  //-- alpha = (Npart-1)/2
  //double knGuess[NCENT] = {16.2,  16.2,  16.2,     0.472, 0.478, 0.480, 0.475, 0.471, 0.453, 0.432, 0.400, 0.360};
  //double e0Guess[NCENT] = {0.0,   0.0,   0.0,      0.137, 0.156, 0.171, 0.183, 0.190, 0.199, 0.207, 0.217, 0.232};
  //double alGuess[NCENT];
  //for(int icent = 0; icent < NCENT; icent++) alGuess[icent] = (Npart[icent]-1)/2.;

  double vnmax[NCENT]   = {0.25,   0.25,   0.25,   0.25,   0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.23};
  //double vnmax[NCENT]   = {0.25,   0.25,   0.25,   0.25,   0.25,  0.25,  0.25,  0.27,  0.27,  0.27,  0.27,  0.23};

  if( fixKn ){
    for(int icent = 0; icent < NCENT; icent++) alGuess[icent] = (Npart[icent]-1)/2.;
    e0Guess[8]  = 0.2;
    e0Guess[9]  = 0.2;
    e0Guess[10] = 0.2;
    e0Guess[11] = 0.18;
  }
  if( fixAlpha ){
    for(int icent = 0; icent < NCENT; icent++){
      alGuess[icent] = (Npart[icent]-1)/2.;
      knGuess[icent] = 0.4;
      e0Guess[icent] = 0.2;
    }
    e0Guess[11] = 0.18;
  }

  int norder_  = 2;
  double tkEta = 1.0;

  TLatex latex;
  TLatex latex2;
  TLatex latex3;

  TFile * fMoscowFits;

  //-- Unfold histos
  TFile * fFinalUnf;
  TH1D * hFinalUnfold[NCENT];
  TH1D * hFinalUnfoldStat[NCENT];
  TH1D * hFinalUnfoldSys[NCENT];
  TH1D * hNormFactor;

  TF1 * fBG[NCENT];
  TF1 * fEllP[NCENT];

  TFile * fBottomLine;
  TGraph * grChi2EllP;
  TGraph * grChi2BG;

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

  //-- Systematics
  TFile * fSys;
  TH1D * SmoothSysTotKn;
  TH1D * SmoothSysTotAlpha;
  TH1D * SmoothSysTotE0;

  double fitKn_syserr[NCENT];
  double fitAlpha_syserr[NCENT];
  double fitE0_syserr[NCENT];

  TGraphErrors * grFitKnSys;
  TGraphErrors * grFitAlphaSys;
  TGraphErrors * grFitE0Sys;

  //-- ATLAS
  TGraphErrors * grATLASKn;
  TGraphErrors * grATLASAlpha;
  TGraphErrors * grATLASEcc0;

  //-- Npart scaling
  TGraph * grNpartScale;
  double npartScale[NCENT];

  TFile * fOut;

  //-- Contours
  const int Nsig = 3;
  TGraph * grE0Kn[Nsig][NCENT];
  TGraph * grAlphaKn[Nsig][NCENT];
  TGraph * grE0Alpha[Nsig][NCENT];
  int contCol[Nsig] = {1, 2, 4};

  TH1D * hDummy[NCENT];

  //
  // MAIN
  //
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(20);
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(23);
  latex3.SetNDC();
  latex3.SetTextFont(43);
  latex3.SetTextSize(26);

  setTDRStyle();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");

  fFinalUnf    = new TFile( Form("systematicStudies/SysUnfoldDistns_v%i.root", norder_) );
  fMoscowFits = new TFile("MoscowEllpFits.root");

  hNormFactor = (TH1D*) fFinalUnf->Get("hNormFactor");

  //-- Systematics
  fSys = new TFile("systematicStudies/SmoothSysTot.root");
  SmoothSysTotKn    = (TH1D*) fSys->Get("SmoothSysTotKn");
  SmoothSysTotAlpha = (TH1D*) fSys->Get("SmoothSysTotAlpha");
  SmoothSysTotE0    = (TH1D*) fSys->Get("SmoothSysTotE0");

  fOut = new TFile(fname.data(), "recreate");
  fBottomLine = new TFile("SmearSpaceChi2.root");
  grChi2EllP = (TGraph*) fBottomLine->Get("grEllPChi2");
  grChi2BG   = (TGraph*) fBottomLine->Get("grBGChi2");

  for(int icent = 0; icent < NCENT; icent++){

    //if(icent != BIN) continue;

    fitKn[icent]    = -1;
    fitAlpha[icent] = -1;
    fitE0[icent]    = -1;

    hDummy[icent] = new TH1D(Form("hDummy_c%i", icent), Form("hDummy_c%i",icent), 152, 0-binw, vnMax[norder_]);
    hDummy[icent]->GetXaxis()->SetTitle("v_{2}");
    hDummy[icent]->GetYaxis()->SetTitle("p(v_{2})");

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
    hFinalUnfoldSys[icent]->GetXaxis()->SetNdivisions(507);
    //if(icent == 3) hFinalUnfoldSys[icent]->GetXaxis()->SetRange(1, hFinalUnfoldSys[icent]->FindBin(0.2)-2);
    //else if(icent == 5) hFinalUnfoldSys[icent]->GetXaxis()->SetRange(1, hFinalUnfoldSys[icent]->FindBin(0.25)-2);
    hDummy[icent]->GetXaxis()->SetRange(1, hFinalUnfoldSys[icent]->FindBin(0.28));
    hFinalUnfoldSys[icent]->GetXaxis()->SetRange(hFinalUnfoldSys[icent]->FindBin(-0.1), hFinalUnfoldSys[icent]->FindBin(0.28));
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitle("v_{2}");
    hFinalUnfoldSys[icent]->GetYaxis()->SetTitle("p(v_{2})");


    //if( icent != 3 && icent != 5 && icent !=7 ) continue;
    if( icent < 3 ) continue;
    std::cout<<Form("-------- Centbin %i --------", icent)<<std::endl;

    if(moscowFits){
      hFinalUnfold[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfold[icent]->GetMaximum());
      hFinalUnfoldStat[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfoldStat[icent]->GetMaximum());
      hFinalUnfoldSys[icent]->Scale(fEllP[icent]->GetMaximum()/hFinalUnfoldSys[icent]->GetMaximum());
    }
    else{
      //-- Unnormalize
      hFinalUnfold[icent]->Scale(1./hNormFactor->GetBinContent(icent+1));
      hFinalUnfoldStat[icent]->Scale(1./hNormFactor->GetBinContent(icent+1));
      hFinalUnfoldSys[icent]->Scale(1./hNormFactor->GetBinContent(icent+1));
    }

    fBG[icent] = new TF1(Form("fBG_c%i", icent), "[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))", 0.0, 0.4);
    fBG[icent]->SetParLimits(1, 0., 1.);
    fBG[icent]->SetParLimits(2, 0., 1.);
    fBG[icent]->SetLineColor(4);
    fBG[icent]->SetLineWidth(2);

    EbyECumu cumu(hFinalUnfold[icent]);
    double vn2 = cumu.GetCumu_vn2();
    double vn4 = cumu.GetCumu_vn4();
    double BGmu    = vn4;
    double BGDelta = sqrt( (vn2*vn2 - vn4*vn4)/2. ) ;
    if(icent != 6){
      BGmu = hFinalUnfold[icent]->GetMean();
      BGDelta = hFinalUnfold[icent]->GetRMS();
    }
    fBG[icent]->SetParameters(hFinalUnfold[icent]->GetMaximum(), BGmu, BGDelta);

    if(moscowFits) fEllP[icent] = (TF1*) fMoscowFits->Get( Form("elliptic%i", icent) );
    else            fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, 0.26, 4);
    fEllP[icent]->SetLineColor(2);
    fEllP[icent]->SetLineWidth(2);

    if(!moscowFits){
      //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
      fEllP[icent]->SetParLimits(0, 0., 1.);
      fEllP[icent]->SetParLimits(2, 0., 1.);
      fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hFinalUnfold[icent]->GetMaximum());
      if( fixKn ) fEllP[icent]->FixParameter(2, knGuess[icent]);
      if( fixAlpha ) fEllP[icent]->FixParameter(1, alGuess[icent]);
      fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");
    }

    hFinalUnfold[icent]->Fit( Form("fBG_c%i", icent), "BL0", "", 0.0, vnmax[icent]);
    if(!moscowFits){
      hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "L0", "", 0.0, vnmax[icent]);
      //-- contours
      if( contours ){
	for(int isig = 0; isig < Nsig; isig++){

	  gMinuit->SetErrorDef( pow(isig+1, 2) );

	  //-- E0 VS kn
	  grE0Kn[isig][icent] = (TGraph*) gMinuit->Contour(50, 2, 0);
	  grE0Kn[isig][icent]->GetXaxis()->SetTitle( Form("k_{%i}", norder_) );
	  grE0Kn[isig][icent]->GetXaxis()->SetNdivisions(508);
	  grE0Kn[isig][icent]->GetYaxis()->SetTitle( "#epsilon_{0}" );
	  grE0Kn[isig][icent]->GetYaxis()->SetNdivisions(508);
	  grE0Kn[isig][icent]->SetLineColor( contCol[isig] );
	  grE0Kn[isig][icent]->SetMarkerColor( contCol[isig] );
	  grE0Kn[isig][icent]->SetFillColor( contCol[isig] );

	  //-- Alpha VS kn
	  grAlphaKn[isig][icent] = (TGraph*) gMinuit->Contour(50, 2, 1);
	  grAlphaKn[isig][icent]->GetXaxis()->SetTitle( Form("k_{%i}", norder_) );
	  grAlphaKn[isig][icent]->GetXaxis()->SetNdivisions(508);
	  grAlphaKn[isig][icent]->GetYaxis()->SetTitle( "#alpha" );
	  grAlphaKn[isig][icent]->GetYaxis()->SetNdivisions(508);
	  grAlphaKn[isig][icent]->SetLineColor( contCol[isig] );     
	  grAlphaKn[isig][icent]->SetMarkerColor( contCol[isig] );
	  grAlphaKn[isig][icent]->SetFillColor( contCol[isig] );
	  
	  //-- E0 VS Alpha
	  grE0Alpha[isig][icent] = (TGraph*) gMinuit->Contour(50, 1, 0);
	  grE0Alpha[isig][icent]->GetXaxis()->SetTitle( "#alpha" );
	  grE0Alpha[isig][icent]->GetXaxis()->SetNdivisions(508);
	  grE0Alpha[isig][icent]->GetYaxis()->SetTitle( "#epsilon_{0}" );
	  grE0Alpha[isig][icent]->GetYaxis()->SetNdivisions(508);
	  grE0Alpha[isig][icent]->SetLineColor( contCol[isig] );     
	  grE0Alpha[isig][icent]->SetMarkerColor( contCol[isig] );
	  grE0Alpha[isig][icent]->SetFillColor( contCol[isig] );
	  
	  fOut->cd();
	  grE0Kn[isig][icent]->Write( Form( "grE0Kn_%is_c%i", isig+1, icent)  );
	  grAlphaKn[isig][icent]->Write( Form("grAlphaKn_%is_c%i", isig+1, icent) );
	  grE0Alpha[isig][icent]->Write( Form("grE0Alpha_%is_c%i", isig+1, icent) );
  
	}
      }
    }

    //-- Increase hist maximums post-fit.
    if(icent == 3 || icent == 5){
      hFinalUnfold[icent]->SetMaximum( 6e6 );
      hFinalUnfoldStat[icent]->SetMaximum( 6e6 );
      hFinalUnfoldSys[icent]->SetMaximum( 6e6 );
      hDummy[icent]->SetMaximum( 6e6 );

    }
    else{
      hFinalUnfold[icent]->SetMaximum( 5e5 );
      hFinalUnfoldStat[icent]->SetMaximum( 5e5 );
      hFinalUnfoldSys[icent]->SetMaximum( 5e5 );
      hDummy[icent]->SetMaximum( 5e5 );
    }

    hFinalUnfold[icent]->SetMinimum( 12 );
    hFinalUnfoldStat[icent]->SetMinimum( 12 );
    hFinalUnfoldSys[icent]->SetMinimum( 12 );
    hDummy[icent]->SetMinimum( 12 );

    if(icent == 3){
      hFinalUnfoldSys[icent]->GetYaxis()->SetLabelSize(0.07);
      hFinalUnfoldSys[icent]->GetYaxis()->SetTitleSize(0.07);
      hFinalUnfoldSys[icent]->GetYaxis()->SetTitleOffset(1.04);

      hDummy[icent]->GetYaxis()->SetLabelSize(0.07);
      hDummy[icent]->GetYaxis()->SetTitleSize(0.07);
      hDummy[icent]->GetYaxis()->SetTitleOffset(1.04);
    }
    if(icent == 7){
      hFinalUnfoldSys[icent]->GetYaxis()->SetLabelSize(0.06);
      hFinalUnfoldSys[icent]->GetXaxis()->SetTitleOffset(1.02);
      hFinalUnfoldSys[icent]->GetXaxis()->SetTitleSize(0.07);

      hDummy[icent]->GetYaxis()->SetLabelSize(0.06);
      hDummy[icent]->GetXaxis()->SetTitleOffset(1.02);
      hDummy[icent]->GetXaxis()->SetTitleSize(0.07);
    }
    if(icent == 9 || icent == 11){
      hFinalUnfoldSys[icent]->GetXaxis()->SetLabelSize(0.07);
      hFinalUnfoldSys[icent]->GetXaxis()->SetLabelOffset(-0.002);
      hFinalUnfoldSys[icent]->GetXaxis()->SetTitleSize(0.08);

      hDummy[icent]->GetXaxis()->SetLabelSize(0.07);
      hDummy[icent]->GetXaxis()->SetLabelOffset(-0.002);
      hDummy[icent]->GetXaxis()->SetTitleSize(0.08);
    }


    fOut->cd();
    fEllP[icent]->Write();
    fBG[icent]->Write();

    if(moscowFits){
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

      if(dosys){
	fitKn_syserr[icent]    = SmoothSysTotKn->GetBinContent(icent+1);
	fitAlpha_syserr[icent] = SmoothSysTotAlpha->GetBinContent(icent+1);
	fitE0_syserr[icent]    = SmoothSysTotE0->GetBinContent(icent+1);

	fitKn_syserr[icent]    *= fitKn[icent];
	fitAlpha_syserr[icent] *= fitAlpha[icent];
	fitE0_syserr[icent]    *= fitE0[icent];
      }


    }

    npartScale[icent] = (Npart[icent] - 1.)/2.;

  } // End cent loop

  //-- Make Graphs
  double c_err[NCENT];
  for(int icent = 0; icent < NCENT; icent++) c_err[icent] = 0;

  //grFitKn    = new TGraphErrors(NCENT, Npart, fitKn,    c_err, fitKn_err);
  //grFitAlpha = new TGraphErrors(NCENT, Npart, fitAlpha, c_err, fitAlpha_err);
  //grFitE0    = new TGraphErrors(NCENT, Npart, fitE0,    c_err, fitE0_err);

  //formatGraph(grFitKn,    "N_{Part}", 0.1, 1.0,  "Fit k_{n}",        1, 20, "grFitKn");
  //formatGraph(grFitAlpha, "N_{Part}", 0.1, 120,  "Fit #alpha",       2, 21, "grFitAlpha");
  //formatGraph(grFitE0,    "N_{Part}", 0.1, 0.5,  "Fit #epsilon_{0}", 4, 34, "grFitE0");

  grFitKn    = new TGraphErrors(NCENT, centBinCenter, fitKn,    c_err, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, centBinCenter, fitAlpha, c_err, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, centBinCenter, fitE0,    c_err, fitE0_err);

  formatGraph(grFitKn,    "Centrality %", 0., 0.75,  "k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "Centrality %", 0., 130,  "#alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "Centrality %", 0., 0.5,  "#epsilon_{0}", 4, 34, "grFitE0");

  fOut->cd();
  grFitKn->Write();
  grFitAlpha->Write();
  grFitE0->Write();

  grFitKnSys    = 0;
  grFitAlphaSys = 0;
  grFitE0Sys    = 0;

  //-- Npart scale graph
  grNpartScale = new TGraph(NCENT, centBinCenter, npartScale);
  formatGraph(grNpartScale, "Centrality %", 0., 120, "#alpha", 1, 20, "grNpartScale");

  if(dosys){
    grFitKnSys    = new TGraphErrors(NCENT, centBinCenter, fitKn,    centBinErr, fitKn_syserr);
    grFitAlphaSys = new TGraphErrors(NCENT, centBinCenter, fitAlpha, centBinErr, fitAlpha_syserr);
    grFitE0Sys    = new TGraphErrors(NCENT, centBinCenter, fitE0,    centBinErr, fitE0_syserr);

    formatGraph(grFitKnSys,    "Centrality %", 0., 0.75, "k_{n}",        1, 20, "grFitKnSys");
    formatGraph(grFitAlphaSys, "Centrality %", 0., 120, "#alpha",       2, 21, "grFitAlphaSys");
    formatGraph(grFitE0Sys,    "Centrality %", 0., 0.4, "#epsilon_{0}", 4, 34, "grFitE0Sys");

    grFitKnSys->SetFillColor(17);
    grFitAlphaSys->SetFillColor(17);
    grFitE0Sys->SetFillColor(17);
  }

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
  legKn->AddEntry(grFitKn,   "CMS",   "lp");
  legKn->AddEntry(grATLASKn, "ATLAS 2.76 TeV", "lp");

  TLegend * legAlpha = new TLegend(0.5, 0.75, 0.95, 0.9);
  legAlpha->SetBorderSize(0);
  legAlpha->SetFillStyle(0);
  legAlpha->AddEntry(grFitAlpha,   "CMS",   "lp");
  legAlpha->AddEntry(grATLASAlpha, "ATLAS 2.76 TeV", "lp");

  TLegend * legEcc0 = new TLegend(0.2, 0.2, 0.65, 0.35);
  legEcc0->SetBorderSize(0);
  legEcc0->SetFillStyle(0);
  legEcc0->AddEntry(grFitE0,     "CMS",   "lp");
  legEcc0->AddEntry(grATLASEcc0, "ATLAS 2.76 TeV", "lp");

  TLegend * legFit = new TLegend(0.18, 0.20, 0.63, 0.35);
  legFit->SetBorderSize(0);
  legFit->SetFillStyle(0);

  //-- Fit parms vs npart
  if(dosys){
    TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 1000, 1000);
    cParmSummary->Divide(2,2);

    cParmSummary->cd(1);
    grFitKnSys->Draw("apE2");
    grFitKn->Draw("psame");
    latex.DrawLatex(0.2, 0.76, Form("|#eta| < %.1f", tkEta));
    latex.DrawLatex(0.2, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
    if(ATLAS){
      grATLASKn->Draw("psame");
      legKn->Draw("same");
    }

    cParmSummary->cd(2);
    grFitE0Sys->Draw("apE2");
    grFitE0->Draw("psame");
    if(ATLAS){
      grATLASEcc0->Draw("psame");
      legEcc0->Draw("same");
    }

    cParmSummary->cd(3);
    grFitAlphaSys->Draw("apE2");
    grFitAlpha->Draw("psame");
    if(ATLAS){
      grATLASAlpha->Draw("psame");
      legAlpha->Draw("same");
    }

    cParmSummary->cd(4);
    cParmSummary->cd(4)->SetLogy();
    grChi2EllP->Draw("ap");
    grChi2BG->Draw("psame");
    legFit->Draw();
    //lone->Draw("same");

    cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_v%i.pdf",norder_));
  }
  else{
    TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 1500, 500);
    cParmSummary->Divide(3,1);

    cParmSummary->cd(1);
    grFitKn->Draw("ap");
    grATLASKn->Draw("psame");
    legKn->Draw("same");
    latex.DrawLatex(0.2, 0.88, "CMS #it{Preliminary}");
    latex.DrawLatex(0.2, 0.82, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    latex.DrawLatex(0.2, 0.76, Form("|#eta| < %.1f", tkEta));
    latex.DrawLatex(0.2, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

    cParmSummary->cd(2);
    grFitE0->Draw("ap");
    grATLASEcc0->Draw("psame");
    legEcc0->Draw("same");
    
    cParmSummary->cd(3);
    grFitAlpha->Draw("ap");
    grATLASAlpha->Draw("psame");
    legAlpha->Draw("same");

    cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_v%i.pdf",norder_));

    TLegend * legNpart = new TLegend(0.2, 0.2, 0.5, 0.4);
    legInit( legNpart );
    legNpart->AddEntry(grFitAlpha, "Fit #alpha", "p");
    legNpart->AddEntry(grNpartScale, "(Npart - 1)/2", "p");

    TCanvas * cNpartScale = new TCanvas("cNpartScale", "cNpartScale", 500, 500);
    cNpartScale->cd();
    grFitAlpha->Draw("ap");
    grNpartScale->Draw("psame");
    legNpart->Draw("same");
    cNpartScale->SaveAs("plots/unfolding/cNpartScale.pdf");

  }


  //-- Fit Summaries
  grChi2EllP->GetXaxis()->SetLimits(10,62.5);
  TLine * lone = new TLine(grChi2EllP->GetXaxis()->GetXmin(), 1., grChi2EllP->GetXaxis()->GetXmax(), 1.);
  lone->SetLineWidth(2);
  lone->SetLineStyle(2);

  TLegend * leg1 = new TLegend(0.04, 0.73, 0.39, 0.89);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(fBG[3],   "Bessel-Gaussian (v_{n}^{RP}, #delta_{v_{n}})", "l");
  leg1->AddEntry(fEllP[3], "Elliptic Power (#epsilon_{0}, #alpha, k_{n})", "l");

  TLegend * leg2 = new TLegend(0.27, 0.78, 0.61, 0.91);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(grChi2BG,   "Bessel-Gaussian (v_{n}^{RP}, #delta_{v_{n}})", "p");
  leg2->AddEntry(grChi2EllP, "Elliptic Power (#epsilon_{0}, #alpha, k_{n})", "p");

  TCanvas * cUnfoldDistsBig = new TCanvas("cUnfoldDistsBig", "cUnfoldDistsBig", 1500, 1000);
  cUnfoldDistsBig->Divide(3,2,0,0);
  cUnfoldDistsBig->SetFillColor(0);
  cUnfoldDistsBig->SetFrameFillStyle(0);

  cUnfoldDistsBig->cd(1);
  cUnfoldDistsBig->cd(1)->SetLogy();
  cUnfoldDistsBig->cd(1)->SetTopMargin(0.07);
  hDummy[3]->Draw();
  hFinalUnfoldSys[3]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[3]->Draw("same");
  fBG[3]->Draw("same");
  fEllP[3]->Draw("same");
  double alpha3;
  double alpha3e;
  double eps03;
  double eps03e;
  double kn3;
  double kn3e;
  if(moscowFits){
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
  double bgmu3  = fBG[3]->GetParameter(1);
  double bgmu3e = fBG[3]->GetParError(1);
  double bgd3   = fBG[3]->GetParameter(2);
  double bgd3e  = fBG[3]->GetParError(2);

  latex3.DrawLatex(0.19, 0.945, "#bf{CMS} #it{Preliminary}");
  latex3.DrawLatex(0.69, 0.945, "PbPb 5.02 TeV");
  latex3.DrawLatex(0.57, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.81, 0.76, Form("|#eta| < %.1f", tkEta));

  latex3.DrawLatex(0.25, 0.46, Form("#bf{Cent. %i - %i%s}",                cent_min[3], cent_max[3], "%") );
  latex3.DrawLatex(0.25, 0.38, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu3,  100.*bgmu3e, "%") );
  latex3.DrawLatex(0.25, 0.30, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd3,   100.*bgd3e,  "%") );
  latex3.DrawLatex(0.25, 0.22, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps03,  100.*eps03e, "%") );
  latex3.DrawLatex(0.25, 0.14, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn3,    100.*kn3e,   "%") );
  latex3.DrawLatex(0.25, 0.06, Form("#alpha = %.1f #pm %.1f",              alpha3,      alpha3e) );

  //cUnfoldDistsBig->cd(1)->SetTopMargin(0.15);
  //cUnfoldDistsBig->cd(1)->SetTickx(0);
  double xmin3 = 0;
  double xmax3 = 0.2-binw;
  double ymax3 = 1.9*hFinalUnfold[3]->GetMaximum();

  TGaxis * axEcc3 = new TGaxis(xmin3, ymax3, xmax3, ymax3, xmin3/kn3, xmax3/kn3, 508, "-");
  axEcc3->SetLabelSize(0.055);
  axEcc3->SetTitleSize(0.055);
  axEcc3->SetTitle("#epsilon_{2}");
  cUnfoldDistsBig->cd(1);
  //axEcc3->Draw("same");

  cUnfoldDistsBig->cd(2);
  cUnfoldDistsBig->cd(2)->SetLogy();
  cUnfoldDistsBig->cd(2)->SetTopMargin(0.07);
  hDummy[5]->Draw();
  hFinalUnfoldSys[5]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[5]->Draw("same");
  fBG[5]->Draw("same");
  fEllP[5]->Draw("same");
  leg1->Draw("same");
  leg1->SetTextFont(43);
  leg1->SetTextSize(26);
  double alpha5;
  double alpha5e;
  double eps05;
  double eps05e;
  double kn5;
  double kn5e;
  if(moscowFits){
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
  double bgmu5  = fBG[5]->GetParameter(1);
  double bgmu5e = fBG[5]->GetParError(1);
  double bgd5   = fBG[5]->GetParameter(2);
  double bgd5e  = fBG[5]->GetParError(2);

  latex3.DrawLatex(0.09, 0.46, Form("#bf{Cent. %i - %i%s}",                cent_min[5], cent_max[5], "%") );
  latex3.DrawLatex(0.09, 0.38, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu5,  100.*bgmu5e, "%") );
  latex3.DrawLatex(0.09, 0.30, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd5,   100.*bgd5e,  "%") );
  latex3.DrawLatex(0.09, 0.22, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps05,  100.*eps05e, "%") );
  latex3.DrawLatex(0.09, 0.14, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn5,    100.*kn5e,   "%") );
  latex3.DrawLatex(0.09, 0.06, Form("#alpha = %.1f #pm %.1f",              alpha5,      alpha5e) );


  //cUnfoldDistsBig->cd(2)->SetTickx(0);
  //cUnfoldDistsBig->cd(2)->SetTopMargin(0.15);
  double xmin5 = 0;
  double xmax5 = 0.25-binw;
  double ymax5 = 1.9*hFinalUnfold[5]->GetMaximum();

  TGaxis * axEcc5 = new TGaxis(xmin5, ymax5, xmax5, ymax5, xmin5/kn5, xmax5/kn5, 508, "-");
  axEcc5->SetTitle("#epsilon_{2}");
  axEcc5->SetLabelSize(0.055);
  axEcc5->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(2);
  //axEcc5->Draw("same");

  cUnfoldDistsBig->cd(4);
  cUnfoldDistsBig->cd(4)->SetLogy();
  cUnfoldDistsBig->cd(4)->SetBottomMargin(0.22);
  hDummy[7]->Draw();
  hFinalUnfoldSys[7]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[7]->Draw("same");
  fBG[7]->Draw("same");
  fEllP[7]->Draw("same");
  double alpha7;
  double alpha7e;
  double eps07;
  double eps07e;
  double kn7;
  double kn7e;
  if(moscowFits){
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
  double bgmu7  = fBG[7]->GetParameter(1);
  double bgmu7e = fBG[7]->GetParError(1);
  double bgd7   = fBG[7]->GetParameter(2);
  double bgd7e  = fBG[7]->GetParError(2);

  latex3.DrawLatex(0.25, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[7], cent_max[7], "%") );
  latex3.DrawLatex(0.25, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu7,  100.*bgmu7e, "%") );
  latex3.DrawLatex(0.25, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd7,   100.*bgd7e,  "%") );
  latex3.DrawLatex(0.25, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps07,  100.*eps07e, "%") );
  latex3.DrawLatex(0.25, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn7,    100.*kn7e,   "%") );
  latex3.DrawLatex(0.25, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha7,      alpha7e) );

  //cUnfoldDistsBig->cd(4)->SetTickx(0);
  //cUnfoldDistsBig->cd(4)->SetTopMargin(0.15);
  double xmin7 = 0;
  double xmax7 = 0.3-binw;
  double ymax7 = 1.9*hFinalUnfold[7]->GetMaximum();

  TGaxis * axEcc7 = new TGaxis(xmin7, ymax7, xmax7, ymax7, xmin7/kn7, xmax7/kn7, 508, "-");
  axEcc7->SetTitle("#epsilon_{2}");
  axEcc7->SetLabelSize(0.055);
  axEcc7->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(4);
  cUnfoldDistsBig->cd(4)->RedrawAxis();
  //axEcc7->Draw("same");

  cUnfoldDistsBig->cd(5);
  cUnfoldDistsBig->cd(5)->SetLogy();
  cUnfoldDistsBig->cd(5)->SetBottomMargin(0.22);
  hDummy[9]->Draw();
  hFinalUnfoldSys[9]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[9]->Draw("same");
  fBG[9]->Draw("same");
  fEllP[9]->Draw("same");
  double alpha9;
  double alpha9e;
  double eps09;
  double eps09e;
  double kn9;
  double kn9e;
  if(moscowFits){
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
  double bgmu9  = fBG[9]->GetParameter(1);
  double bgmu9e = fBG[9]->GetParError(1);
  double bgd9   = fBG[9]->GetParameter(2);
  double bgd9e  = fBG[9]->GetParError(2);


  latex3.DrawLatex(0.09, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[9], cent_max[9], "%") );
  latex3.DrawLatex(0.09, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu9,  100.*bgmu9e, "%") );
  latex3.DrawLatex(0.09, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd9,   100.*bgd9e,  "%") );
  latex3.DrawLatex(0.09, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps09,  100.*eps09e, "%") );
  latex3.DrawLatex(0.09, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn9,    100.*kn9e,   "%") );
  latex3.DrawLatex(0.09, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha9,      alpha9e) );

  //cUnfoldDistsBig->cd(5)->SetTickx(0);
  //cUnfoldDistsBig->cd(5)->SetTopMargin(0.15);
  double xmin9 = 0;
  double xmax9 = 0.3-binw;
  double ymax9 = 1.9*hFinalUnfold[9]->GetMaximum();

  TGaxis * axEcc9 = new TGaxis(xmin9, ymax9, xmax9, ymax9, xmin9/kn9, xmax9/kn9, 508, "-");
  axEcc9->SetTitle("#epsilon_{2}");
  axEcc9->SetLabelSize(0.055);
  axEcc9->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(5);
  //axEcc9->Draw("same");

  cUnfoldDistsBig->cd(6);
  cUnfoldDistsBig->cd(6)->SetLogy();
  cUnfoldDistsBig->cd(6)->SetBottomMargin(0.22);
  hDummy[11]->Draw();
  hFinalUnfoldSys[11]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[11]->Draw("same");
  fBG[11]->Draw("same");
  fEllP[11]->Draw("same");
  double alpha11;
  double alpha11e;
  double eps011;
  double eps011e;
  double kn11;
  double kn11e;
  if(moscowFits){
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
  double bgmu11  = fBG[11]->GetParameter(1);
  double bgmu11e = fBG[11]->GetParError(1);
  double bgd11   = fBG[11]->GetParameter(2);
  double bgd11e  = fBG[11]->GetParError(2);

  latex3.DrawLatex(0.09, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[11], cent_max[11], "%") );
  latex3.DrawLatex(0.09, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu11,  100.*bgmu11e, "%") );
  latex3.DrawLatex(0.09, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd11,   100.*bgd11e,  "%") );
  latex3.DrawLatex(0.09, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps011,  100.*eps011e, "%") );
  latex3.DrawLatex(0.09, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn11,    100.*kn11e,   "%") );
  latex3.DrawLatex(0.09, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha11,      alpha11e) );

  //cUnfoldDistsBig->cd(6)->SetTickx(0);
  //cUnfoldDistsBig->cd(6)->SetTopMargin(0.15);
  //cUnfoldDistsBig->cd(6)->SetLeftMargin(0.19);
  //cUnfoldDistsBig->cd(6)->SetRightMargin(0.07);
  double xmin11 = 0;
  double xmax11 = 0.3-binw;
  double ymax11 = hFinalUnfold[11]->GetMaximum();

  TGaxis * axEcc11 = new TGaxis(xmin11, ymax11, xmax11, ymax11, xmin11/kn11, xmax11/kn11, 508, "-");
  axEcc11->SetTitle("#epsilon_{2}");
  axEcc11->SetLabelSize(0.055);
  axEcc11->SetTitleSize(0.055);
  cUnfoldDistsBig->cd(6);
  //axEcc11->Draw("same");

  cUnfoldDistsBig->cd(3);
  cUnfoldDistsBig->cd(3)->SetFillStyle(4000);
  cUnfoldDistsBig->cd(3)->SetTopMargin(0.07);
  cUnfoldDistsBig->cd(3)->SetLeftMargin(0.25);
  cUnfoldDistsBig->cd(3)->SetBottomMargin(0.2);
  grChi2EllP->GetYaxis()->SetRangeUser(0.09,15);
  grChi2EllP->GetYaxis()->SetLabelSize(0.06);
  grChi2EllP->GetXaxis()->SetLabelSize(0.06);
  grChi2BG->SetMarkerColor(4);
  grChi2EllP->Draw("ap");
  grChi2BG->Draw("psame");
  lone->Draw("same");
  latex2.DrawLatex(0.25, 0.95, "#bf{CMS} #it{Preliminary}");
  latex2.DrawLatex(0.68, 0.95, "PbPb 5.02 TeV");
  //latex2.DrawLatex(0.29, 0.73, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  //latex2.DrawLatex(0.29, 0.66, Form("|#eta| < %.1f", tkEta));

  leg2->Draw("same");

  leg2->SetTextFont(43);
  leg2->SetTextSize(23);

  cUnfoldDistsBig->Update();
  cUnfoldDistsBig->SaveAs(Form("plots/unfolding/finalUnfFit_v%i.pdf",norder_));

  //-- Standalone plot for Yen-Jie:
  /*
  TLegend * leg3 = new TLegend(0.22, 0.73, 0.57, 0.89);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->AddEntry(fBG[3],   "Bessel-Gaussian (v_{n}^{RP}, #delta_{v_{n}})", "l");
  leg3->AddEntry(fEllP[3], "Elliptic Power (#epsilon_{0}, #alpha, k_{n})", "l");


  TCanvas * c = new TCanvas("c", "c", 500, 500);

  c->cd();
  c->SetTopMargin(0.1);
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.1);
  c->SetLogy();

  hDummy[9]->GetXaxis()->SetLabelSize(0.05);
  hDummy[9]->GetXaxis()->SetLabelOffset(0.001);
  hDummy[9]->SetMaximum( 2e6 );
  hDummy[9]->Draw();
  hFinalUnfoldSys[9]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[9]->Draw("same");
  fBG[9]->Draw("same");
  fEllP[9]->Draw("same");
  leg3->Draw("same");
  leg3->SetTextFont(43);
  leg3->SetTextSize(26);

  latex.DrawLatex(0.20, 0.913, "#bf{CMS} #it{Preliminary}");
  latex.DrawLatex(0.64, 0.913, "PbPb 5.02 TeV");

  latex3.DrawLatex(0.25, 0.55, Form("#bf{Cent. %i - %i%s}",                cent_min[9], cent_max[9], "%") );
  latex3.DrawLatex(0.25, 0.48, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu9,  100.*bgmu9e, "%") );
  latex3.DrawLatex(0.25, 0.41, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd9,   100.*bgd9e,  "%") );
  latex3.DrawLatex(0.25, 0.34, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps09,  100.*eps09e, "%") );
  latex3.DrawLatex(0.25, 0.27, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn9,    100.*kn9e,   "%") );
  latex3.DrawLatex(0.25, 0.20, Form("#alpha = %.1f #pm %.1f",              alpha9,      alpha9e) );

  c->Update();
  c->SaveAs("pv2Fit_ForYenJie.pdf");
  */
  /*
  TCanvas * cUnfoldDistsBig_NoFit = new TCanvas("cUnfoldDistsBig_NoFit", "cUnfoldDistsBig_NoFit", 1500, 1000);
  cUnfoldDistsBig_NoFit->Divide(3,2,0,0);
  cUnfoldDistsBig_NoFit->SetFillColor(0);
  cUnfoldDistsBig_NoFit->SetFrameFillStyle(0);

  cUnfoldDistsBig_NoFit->cd(1);
  cUnfoldDistsBig_NoFit->cd(1)->SetLogy();
  cUnfoldDistsBig_NoFit->cd(1)->SetTopMargin(0.07);
  hDummy[3]->Draw();
  hFinalUnfoldSys[3]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[3]->Draw("same");

  latex3.DrawLatex(0.25, 0.06, Form("#bf{Cent. %i - %i%s}", cent_min[3], cent_max[3], "%") );
  latex3.DrawLatex(0.19, 0.945, "#bf{CMS} #it{Preliminary}");
  latex3.DrawLatex(0.69, 0.945, "PbPb 5.02 TeV");
  latex3.DrawLatex(0.57, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
  latex3.DrawLatex(0.81, 0.76, Form("|#eta| < %.1f", tkEta));

  cUnfoldDistsBig_NoFit->cd(2);
  cUnfoldDistsBig_NoFit->cd(2)->SetLogy();
  cUnfoldDistsBig_NoFit->cd(2)->SetTopMargin(0.07);
  hDummy[5]->Draw();
  hFinalUnfoldSys[5]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[5]->Draw("same");
  latex3.DrawLatex(0.09, 0.06, Form("#bf{Cent. %i - %i%s}", cent_min[5], cent_max[5], "%") );

  cUnfoldDistsBig_NoFit->cd(4);
  cUnfoldDistsBig_NoFit->cd(4)->SetLogy();
  cUnfoldDistsBig_NoFit->cd(4)->SetBottomMargin(0.22);
  hDummy[7]->Draw();
  hFinalUnfoldSys[7]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[7]->Draw("same");
  latex3.DrawLatex(0.25, 0.28, Form("#bf{Cent. %i - %i%s}", cent_min[7], cent_max[7], "%") );

  cUnfoldDistsBig_NoFit->cd(5);
  cUnfoldDistsBig_NoFit->cd(5)->SetLogy();
  cUnfoldDistsBig_NoFit->cd(5)->SetBottomMargin(0.22);
  hDummy[9]->Draw();
  hFinalUnfoldSys[9]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[9]->Draw("same");
  latex3.DrawLatex(0.09, 0.28, Form("#bf{Cent. %i - %i%s}", cent_min[9], cent_max[9], "%") );

  cUnfoldDistsBig_NoFit->cd(6);
  cUnfoldDistsBig_NoFit->cd(6)->SetLogy();
  cUnfoldDistsBig_NoFit->cd(6)->SetBottomMargin(0.22);
  hDummy[11]->Draw();
  hFinalUnfoldSys[11]->Draw("e2same");
  setex2->Draw();
  hFinalUnfoldStat[11]->Draw("same");
  latex3.DrawLatex(0.09, 0.28, Form("#bf{Cent. %i - %i%s}", cent_min[11], cent_max[11], "%") );

  cUnfoldDistsBig_NoFit->cd(3);
  cUnfoldDistsBig_NoFit->cd(3)->SetFillStyle(4000);
  cUnfoldDistsBig_NoFit->cd(3)->SetTopMargin(0.07);
  cUnfoldDistsBig_NoFit->cd(3)->SetLeftMargin(0.25);
  cUnfoldDistsBig_NoFit->cd(3)->SetBottomMargin(0.2);

  cUnfoldDistsBig_NoFit->SaveAs("cUnfoldDistsBig_NoFit.pdf");
  */
}

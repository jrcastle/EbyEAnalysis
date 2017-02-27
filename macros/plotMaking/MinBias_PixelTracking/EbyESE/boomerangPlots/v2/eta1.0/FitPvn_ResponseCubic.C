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

int BIN      = 11;
bool TEST    = 0;
int Scenario = 1; // [1] -> vn = kn*en + knPr*kn*en^3;  [2] -> vn = kn*en + knPr*en^3

string fname = "EllPFits_Cubic1_Test1.root";

double fEccn(double vn, double kn, double knPr){
  double value = 0.;
  if(Scenario == 1) value = (pow(0.6666666666666666,0.3333333333333333)*kn)/pow(-9*pow(kn,2)*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2)),0.3333333333333333) - pow(-9*pow(kn,2)*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2)),0.3333333333333333)/(pow(2,0.3333333333333333)*pow(3,0.6666666666666666)*kn*knPr);
  if(Scenario == 2) value = (pow(0.6666666666666666,0.3333333333333333)*kn)/pow(-9*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2)),0.3333333333333333) - pow(-9*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2)),0.3333333333333333)/(pow(2,0.3333333333333333)*pow(3,0.6666666666666666)*knPr);
  return value;
}

double DEccnDvn(double vn,double kn, double knPr){
  double value = 0.;
  if(Scenario == 1) value = -(pow(0.6666666666666666,0.3333333333333333)*kn*(-9*pow(kn,2)*pow(knPr,2) + (27*sqrt(3)*pow(kn,4)*pow(knPr,4)*vn)/sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2))))/(3.*pow(-9*pow(kn,2)*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2)),1.3333333333333333)) - (-9*pow(kn,2)*pow(knPr,2) + (27*sqrt(3)*pow(kn,4)*pow(knPr,4)*vn)/sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2)))/(3.*pow(2,0.3333333333333333)*pow(3,0.6666666666666666)*kn*knPr*pow(-9*pow(kn,2)*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,6)*pow(knPr,3) + 27*pow(kn,4)*pow(knPr,4)*pow(vn,2)),0.6666666666666666));
  if(Scenario == 2) value = -(pow(0.6666666666666666,0.3333333333333333)*kn*(-9*pow(knPr,2) + (27*sqrt(3)*pow(knPr,4)*vn)/sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2))))/(3.*pow(-9*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2)),1.3333333333333333)) - (-9*pow(knPr,2) + (27*sqrt(3)*pow(knPr,4)*vn)/sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2)))/(3.*pow(2,0.3333333333333333)*pow(3,0.6666666666666666)*knPr*pow(-9*pow(knPr,2)*vn + sqrt(3)*sqrt(4*pow(kn,3)*pow(knPr,3) + 27*pow(knPr,4)*pow(vn,2)),0.6666666666666666));
  return value;
}


double pEllP(double * x, double * par){
  int nbin = 50;

  //-- [0]  = e0
  //-- [1]  = alpha
  //-- [2]  = kn
  //-- [3]  = knPr
  //-- [3]  = Scale 
  //-- x[0] = vn
  double e0    = par[0];
  double alpha = par[1];
  double kn    = par[2];
  double knPr  = par[3];
  double scale = par[4];
  double vn    = x[0];
  double eccn  = fEccn(vn, kn, knPr);

  double pi = TMath::Pi();
  double p1 = (scale * 2. * alpha * eccn / pi / DEccnDvn(vn, kn, knPr) ) * pow( 1 - e0*e0, alpha + 0.5);
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 1; i <= nbin; i++){
    double phi = (2*(double)i-1.)*pi/(2.*(double)nbin);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  double ellp = p1 * integ;
  return ellp;
}

void FitPvn_ResponseCubic(){

  bool dosys       = 0;
  bool ATLAS       = 0;
  bool PRC_93_2016 = 0;
  bool fixKn       = 0;
  bool fixKnPr     = 0;
  bool fixAlpha    = 0;
  bool russianFits = 0;

  double fixedKnPr = 0.0;
  //-- Free kn, knPR [1]     0      1      2      3     4     5     6     7     8     9     10    11
  double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.48, 0.40, 0.30, 0.37, 0.34, 0.35, 0.29, 0.28, 0.24};
  double knPrGuess[NCENT] = {0.10,  0.10,  0.10,  0.10, 0.12, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10};
  double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 60.0, 60.0, 45.0, 41.0, 30.0, 30.0, 12.0, 11.0, 6.8};
  double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.13, 0.18, 0.22, 0.23, 0.26, 0.24, 0.30, 0.31, 0.31};
  double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25, 0.25, 0.25, 0.25, 0.27, 0.26, 0.27, 0.26, 0.26};
  //-- Fix kn = ATLAS, Fix knPr = 0.1 [1]   0      1      2      3      4      5      6      7      8      9      10     11
  //double knGuess[NCENT]                = {16.2,  16.2,  16.2,  0.368, 0.331, 0.328, 0.303, 0.294, 0.269, 0.265, 0.234, 0.221};
  //double knPrGuess[NCENT]              = {0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10};
  //double alGuess[NCENT]                = {2.8e5, 2.8e5, 2.8e5, 60.0,  60.0,  45.0,  41.0,  30.0,  30.0,  12.0,  11.0,  6.8}; 
  //double e0Guess[NCENT]                = {0.0,   0.0,   0.0,   0.13,  0.18,  0.22,  0.23,  0.26,  0.24,  0.30,  0.31,  0.31};
  //double vnmax[NCENT]                  = {0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.27,  0.26,  0.26,  0.25,  0.24};
  //-- Fix knPr = 0.1 [1]      0      1      2      3     4     5     6     7     8     9     10    11 
  //double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.48, 0.40, 0.30, 0.37, 0.34, 0.35, 0.29, 0.28, 0.24};
  //double knPrGuess[NCENT] = {0.10,  0.10,  0.10,  0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10};
  //double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 60.0, 60.0, 45.0, 41.0, 30.0, 30.0, 12.0, 11.0, 6.8};
  //double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.13, 0.20, 0.22, 0.23, 0.26, 0.24, 0.30, 0.31, 0.31};
  //double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25, 0.25, 0.25, 0.25, 0.27, 0.26, 0.26, 0.25, 0.25}; 
  //-- Free kn, knPR [2]     0      1      2      3     4     5     6     7     8     9     10    11
  //double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.48, 0.40, 0.30, 0.37, 0.34, 0.34, 0.29, 0.28, 0.24};
  //double knPrGuess[NCENT] = {0.10,  0.10,  0.10,  0.08, 0.08, 0.10, 0.04, 0.03, 0.03, 0.05, 0.05, 0.05};
  //double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 60.0, 60.0, 45.0, 38.0, 35.0, 30.0, 25.0, 20.0, 12.0};
  //double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.13, 0.18, 0.22, 0.23, 0.25, 0.25, 0.27, 0.27, 0.29};
  //double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25, 0.25, 0.25, 0.25, 0.27, 0.26, 0.27, 0.26, 0.26};
  //-- Fix knPr = PRC [2]  
  //TF1 * f = new TF1("f", "pol1", 0., 60.);
  //f->SetParameters(0.0262387, 0.00190748);
  //double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.245, 0.240, 0.239, 0.233, 0.224, 0.219, 0.210, 0.202, 0.196};
  //double knPrGuess[NCENT];
  //for(int icent = 0; icent < NCENT; icent++) knPrGuess[icent] = f->Eval(centBinCenter[icent]);
  //double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 30.0,  30.0,  45.0,  41.0,  30.0,  30.0,  12.0,  11.0,  6.8};
  //double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.30,  0.30,  0.22,  0.23,  0.26,  0.24,  0.30,  0.31,  0.31};
  //double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.27,  0.26,  0.28,  0.27,  0.26};
  //-- Fix kn = PRC [2] 
  //TF1 * f = new TF1("f", "pol1", 0., 60.);
  //f->SetParameters(0.0262387, 0.00190748);
  //double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.245, 0.240, 0.239, 0.233, 0.224, 0.219, 0.210, 0.202, 0.196};
  //double knPrGuess[NCENT];
  //for(int icent = 0; icent < NCENT; icent++) knPrGuess[icent] = f->Eval(centBinCenter[icent]);
  //double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 30.0,  30.0,  45.0,  41.0,  30.0,  30.0,  12.0,  11.0,  6.8};
  //double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.30,  0.30,  0.22,  0.23,  0.26,  0.24,  0.30,  0.31,  0.31};
  //double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.27,  0.26,  0.28,  0.27,  0.26};
  //-- Fix kn, knPr = PRC [2]
  //TF1 * f = new TF1("f", "pol1", 0., 60.);
  //f->SetParameters(0.0262387, 0.00190748);
  //double knGuess[NCENT]   = {16.2,  16.2,  16.2,  0.245, 0.240, 0.239, 0.233, 0.224, 0.219, 0.210, 0.202, 0.196};
  //double knPrGuess[NCENT];
  //for(int icent = 0; icent < NCENT; icent++) knPrGuess[icent] = f->Eval(centBinCenter[icent]); 
  //double alGuess[NCENT]   = {2.8e5, 2.8e5, 2.8e5, 30.0,  30.0,  45.0,  41.0,  30.0,  30.0,  12.0,  11.0,  6.8};
  //double e0Guess[NCENT]   = {0.0,   0.0,   0.0,   0.30,  0.30,  0.22,  0.23,  0.26,  0.24,  0.30,  0.31,  0.31};
  //double vnmax[NCENT]     = {0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.27,  0.26,  0.28,  0.28,  0.28};

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

  TFile * fRussianFits;

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

  //-- Parms vs cent
  double fitKn[NCENT];
  double fitKnPr[NCENT];
  double fitAlpha[NCENT];
  double fitE0[NCENT];

  double fitKn_err[NCENT];
  double fitKnPr_err[NCENT];
  double fitAlpha_err[NCENT];
  double fitE0_err[NCENT];

  TGraphErrors * grFitKn;
  TGraphErrors * grFitKnPr;
  TGraphErrors * grFitAlpha;
  TGraphErrors * grFitE0;

  //-- ATLAS
  TGraphErrors * grATLASKn;
  TGraphErrors * grATLASAlpha;
  TGraphErrors * grATLASEcc0;

  //-- Phys.Rev. C93 (2016) Kns
  TGraph * grPRC_93_2016Kn;
  TGraph * grPRC_93_2016KnPr;

  //-- Systematics
  TFile * fSys;
  TH1D * SmoothSysTotKn;
  TH1D * SmoothSysTotKnPr;
  TH1D * SmoothSysTotAlpha;
  TH1D * SmoothSysTotE0;

  double fitKn_syserr[NCENT];
  double fitKnPr_syserr[NCENT];
  double fitAlpha_syserr[NCENT];
  double fitE0_syserr[NCENT];

  TGraphErrors * grFitKnSys;
  TGraphErrors * grFitKnPrSys;
  TGraphErrors * grFitAlphaSys;
  TGraphErrors * grFitE0Sys;

  TFile * fOut;

  //-- Contours
  const int Nsig = 3;
  TGraph * grKnE0[Nsig][NCENT];
  TGraph * grKnAlpha[Nsig][NCENT];
  TGraph * grAlphaE0[Nsig][NCENT];
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
  fRussianFits = new TFile("RussianEllpFits.root");

  hNormFactor = (TH1D*) fFinalUnf->Get("hNormFactor");

  //-- Systematics
  fSys = new TFile("systematicStudies/SmoothSysTot.root");
  SmoothSysTotKn    = (TH1D*) fSys->Get("SmoothSysTotKn");
  SmoothSysTotAlpha = (TH1D*) fSys->Get("SmoothSysTotAlpha");
  SmoothSysTotE0    = (TH1D*) fSys->Get("SmoothSysTotE0");

  fOut = new TFile(fname.data(), "recreate");
  fBottomLine = new TFile("SmearSpaceChi2_Cubic.root");
  grChi2EllP = (TGraph*) fBottomLine->Get("grEllPChi2");
  grChi2BG   = (TGraph*) fBottomLine->Get("grBGChi2");

  for(int icent = 0; icent < NCENT; icent++){

    if(TEST && icent != BIN) continue;

    fitKn[icent]    = -1;
    fitKnPr[icent]  = -1;
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

    if(russianFits){
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

    if(russianFits) fEllP[icent] = (TF1*) fRussianFits->Get( Form("elliptic%i", icent) );
    else            fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, vnmax[icent], 5);
    fEllP[icent]->SetLineColor(2);
    fEllP[icent]->SetLineWidth(2);

    if(!russianFits){
      //-- [0] = ecc0; 
      //-- [1] = Alpha; 
      //-- [2] = kn; 
      //-- [3] = knPr;
      //-- [4] = Scale;
      fEllP[icent]->SetParLimits(0, 0., 1.);
      fEllP[icent]->SetParLimits(2, 0., 1.);
      fEllP[icent]->SetParLimits(3, 0., 1.);
      fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], knGuess[icent], hFinalUnfold[icent]->GetMaximum());
      if( fixKn )    fEllP[icent]->FixParameter(2, knGuess[icent]);
      if( fixKnPr )  fEllP[icent]->FixParameter(3, knPrGuess[icent]);
      if( fixAlpha ) fEllP[icent]->FixParameter(1, alGuess[icent]);
      fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "knPr", "Scale");
    }

    hFinalUnfold[icent]->Fit( Form("fBG_c%i", icent), "BL0", "", 0.0, vnmax[icent]);
    if(!russianFits){
      hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "L0", "", 0.0, vnmax[icent]);

      //-- contours
      /*
      for(int isig = 0; isig < Nsig; isig++){

	gMinuit->SetErrorDef( pow(isig+1, 2) );

	//-- E0 VS kn
	grKnE0[isig][icent] = (TGraph*) gMinuit->Contour(50, 2, 0);
	grKnE0[isig][icent]->GetXaxis()->SetTitle( Form("k_{%i}", norder_) );
	grKnE0[isig][icent]->GetXaxis()->SetNdivisions(508);
	grKnE0[isig][icent]->GetYaxis()->SetTitle( "#epsilon_{0}" );
	grKnE0[isig][icent]->GetYaxis()->SetNdivisions(508);
	grKnE0[isig][icent]->SetLineColor( contCol[isig] );
	grKnE0[isig][icent]->SetMarkerColor( contCol[isig] );
	grKnE0[isig][icent]->SetFillColor( contCol[isig] );

	//-- Alpha VS kn
	grKnAlpha[isig][icent] = (TGraph*) gMinuit->Contour(50, 2, 1);
	grKnAlpha[isig][icent]->GetXaxis()->SetTitle( Form("k_{%i}", norder_) );
	grKnAlpha[isig][icent]->GetXaxis()->SetNdivisions(508);
	grKnAlpha[isig][icent]->GetYaxis()->SetTitle( "#alpha" );
	grKnAlpha[isig][icent]->GetYaxis()->SetNdivisions(508);
	grKnAlpha[isig][icent]->SetLineColor( contCol[isig] );     
	grKnAlpha[isig][icent]->SetMarkerColor( contCol[isig] );
        grKnAlpha[isig][icent]->SetFillColor( contCol[isig] );

	//-- E0 VS Alpha
	grAlphaE0[isig][icent] = (TGraph*) gMinuit->Contour(50, 1, 0);
	grAlphaE0[isig][icent]->GetXaxis()->SetTitle( "#alpha" );
	grAlphaE0[isig][icent]->GetXaxis()->SetNdivisions(508);
	grAlphaE0[isig][icent]->GetYaxis()->SetTitle( "#epsilon_{0}" );
	grAlphaE0[isig][icent]->GetYaxis()->SetNdivisions(508);
	grAlphaE0[isig][icent]->SetLineColor( contCol[isig] );     
	grAlphaE0[isig][icent]->SetMarkerColor( contCol[isig] );
        grAlphaE0[isig][icent]->SetFillColor( contCol[isig] );
	
	fOut->cd();
	//grKnE0[isig][icent]->Write( Form( "grKnE0_%is_c%i", isig+1, icent)  );
	//grKnAlpha[isig][icent]->Write( Form("grKnAlpha_%is_c%i", isig+1, icent) );
	grAlphaE0[isig][icent]->Write( Form("grAlphaE0_%is_c%i", isig+1, icent) );
  
      }
      */
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
      fitKnPr[icent]  = fEllP[icent]->GetParameter(3);
      fitAlpha[icent] = fEllP[icent]->GetParameter(1);
      fitE0[icent]    = fEllP[icent]->GetParameter(0);

      fitKn_err[icent]    = fEllP[icent]->GetParError(2);
      fitKnPr_err[icent]  = fEllP[icent]->GetParError(3);
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

  } // End cent loop

  if( TEST ){
    TCanvas * ccc = new TCanvas("ccc", "ccc", 500, 500);
    ccc->cd();
    ccc->SetLogy();
    hDummy[BIN]->Draw();
    hFinalUnfoldSys[BIN]->Draw("e2same");
    setex2->Draw();
    hFinalUnfoldStat[BIN]->Draw("same");
    fBG[BIN]->Draw("same");
    fEllP[BIN]->Draw("same");
    ccc->SaveAs("TESTCubicRespFit.pdf");
  }

  if( !TEST ){

    //-- Make Graphs

    //-- ATLAS 
    grATLASKn    = 0;
    grATLASAlpha = 0;
    grATLASEcc0  = 0;

    if( ATLAS ){

      if( fixKnPr && fixedKnPr == 0.1 ){
	grATLASKn    = new TGraphErrors("k_n_ATLAS_black.txt", "%lg %lg %lg");
	grATLASAlpha = new TGraphErrors("alpha_ATLAS_black.txt", "%lg %lg %lg");
	grATLASEcc0  = new TGraphErrors("ecc0_ATLAS_black.txt", "%lg %lg %lg");
      }
      else if( fixKnPr && fixedKnPr == 0.0 ){
	grATLASKn    = new TGraphErrors("k_n_ATLAS_green.txt", "%lg %lg %lg");
        grATLASAlpha = new TGraphErrors("alpha_ATLAS_green.txt", "%lg %lg %lg");
        grATLASEcc0  = new TGraphErrors("ecc0_ATLAS_green.txt", "%lg %lg %lg");
      }

      formatGraph(grATLASKn,    "Centrality %", 0.1, 1.0,  "k_{n}",        1, 24, "grATLASKn");
      formatGraph(grATLASAlpha, "Centrality %", 0.1, 120,  "#alpha",       2, 25, "grATLASAlpha");
      formatGraph(grATLASEcc0,  "Centrality %", 0.1, 0.5,  "#epsilon_{0}", 4, 26, "grATLASE0");
    }
    //-- PRC_93_2016
    if( PRC_93_2016 ){
      grPRC_93_2016Kn   = new TGraph("k_n_PRC_93_2016.txt", "%lg %lg");
      grPRC_93_2016KnPr = new TGraph("knPr_PRC_93_2016.txt", "%lg %lg");

      formatGraph(grPRC_93_2016Kn,   "Centrality %", 0.1,   1.0,  "k_{n}",    2, 24, "grPRC_93_2016Kn");
      formatGraph(grPRC_93_2016KnPr, "Centrality %", -0.01, 0.15, "#kappa\'", 4, 22, "grPRC_93_2016KnPr");
    }

    double c_err[NCENT];
    for(int icent = 0; icent < NCENT; icent++) c_err[icent] = 0;

    grFitKn    = new TGraphErrors(NCENT, centBinCenter, fitKn,    c_err, fitKn_err);
    grFitKnPr  = new TGraphErrors(NCENT, centBinCenter, fitKnPr,  c_err, fitKnPr_err);
    grFitAlpha = new TGraphErrors(NCENT, centBinCenter, fitAlpha, c_err, fitAlpha_err);
    grFitE0    = new TGraphErrors(NCENT, centBinCenter, fitE0,    c_err, fitE0_err);

    formatGraph(grFitKn,    "Centrality %", 0.,   0.75, "k_{n}",        1,        20, "grFitKn");
    formatGraph(grFitKnPr,  "Centrality %", -0.01, 0.15, "#kappa\'",    kGreen+3, 22, "grFitKnPr");
    formatGraph(grFitAlpha, "Centrality %", 0.,   130,  "#alpha",       2,        21, "grFitAlpha");
    formatGraph(grFitE0,    "Centrality %", 0.,   0.5,  "#epsilon_{0}", 4,        34, "grFitE0");
    grFitE0->SetMarkerSize(1.5);

    fOut->cd();
    grFitKn->Write();
    grFitKnPr->Write();
    grFitAlpha->Write();
    grFitE0->Write();

    TLegend * legKn = new TLegend(0.2, 0.2, 0.65, 0.35);
    legKn->SetBorderSize(0);
    legKn->SetFillStyle(0);
    legKn->AddEntry(grFitKn, "CMS", "lp");
    if(ATLAS) legKn->AddEntry(grATLASKn, "ATLAS", "lp");
    if(PRC_93_2016) legKn->AddEntry(grPRC_93_2016Kn, "Hydro #eta/s = 1/4#pi", "l");

    TLegend * legKnPr = new TLegend(0.2, 0.2, 0.65, 0.35);
    legKnPr->SetBorderSize(0);
    legKnPr->SetFillStyle(0);
    legKnPr->AddEntry(grFitKnPr, "CMS", "lp");
    if(PRC_93_2016) legKnPr->AddEntry(grPRC_93_2016KnPr, "Hydro #eta/s = 1/4#pi", "l");

    TLegend * legAlpha = new TLegend(0.5, 0.75, 0.95, 0.9);
    legAlpha->SetBorderSize(0);
    legAlpha->SetFillStyle(0);
    legAlpha->AddEntry(grFitAlpha,   "CMS",   "lp");
    if(ATLAS) legAlpha->AddEntry(grATLASAlpha, "ATLAS", "lp");

    TLegend * legEcc0 = new TLegend(0.2, 0.2, 0.65, 0.35);
    legEcc0->SetBorderSize(0);
    legEcc0->SetFillStyle(0);
    legEcc0->AddEntry(grFitE0,     "CMS",   "lp");
    if(ATLAS) legEcc0->AddEntry(grATLASEcc0, "ATLAS", "lp");

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

      cParmSummary->cd(2);
      grFitE0Sys->Draw("apE2");
      grFitE0->Draw("psame");

      cParmSummary->cd(3);
      grFitAlphaSys->Draw("apE2");
      grFitAlpha->Draw("psame");

      cParmSummary->cd(4);
      cParmSummary->cd(4)->SetLogy();
      grChi2EllP->Draw("ap");
      grChi2BG->Draw("psame");
      legFit->Draw();
      //lone->Draw("same");

      cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_Cubic_v%i.pdf",norder_));
    }
    else{
      TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 2000, 500);
      cParmSummary->Divide(4,1);

      cParmSummary->cd(1);
      grFitKn->Draw("ap");
      if(ATLAS) grATLASKn->Draw("psame");
      if(PRC_93_2016) grPRC_93_2016Kn->Draw("lsame");
      if(ATLAS || PRC_93_2016) legKn->Draw("same");
      latex.DrawLatex(0.2, 0.88, "CMS #it{Preliminary}");
      latex.DrawLatex(0.2, 0.82, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
      latex.DrawLatex(0.2, 0.76, Form("|#eta| < %.1f", tkEta));
      latex.DrawLatex(0.2, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));

      cParmSummary->cd(2);
      grFitKnPr->Draw("ap");
      if(PRC_93_2016){
	grPRC_93_2016KnPr->Draw("lsame");
	legKnPr->Draw("same");
      }

      cParmSummary->cd(3);
      grFitE0->Draw("ap");
      if(ATLAS){
        grATLASEcc0->Draw("psame");
        legEcc0->Draw("same");
      } 
    
      cParmSummary->cd(4);
      grFitAlpha->Draw("ap");
      if(ATLAS){
	grATLASAlpha->Draw("psame");
        legAlpha->Draw("same");
      } 

      cParmSummary->SaveAs(Form("plots/unfolding/FitParmSummary_Cubic_v%i.pdf",norder_));

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
    double bgmu3  = fBG[3]->GetParameter(1);
    double bgmu3e = fBG[3]->GetParError(1);
    double bgd3   = fBG[3]->GetParameter(2);
    double bgd3e  = fBG[3]->GetParError(2);

    latex3.DrawLatex(0.19, 0.945, "#bf{CMS} #it{Preliminary}");
    latex3.DrawLatex(0.69, 0.945, "PbPb 5.02 TeV");
    latex3.DrawLatex(0.57, 0.83, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]));
    latex3.DrawLatex(0.81, 0.76, Form("|#eta| < %.1f", tkEta));

    latex3.DrawLatex(0.25, 0.06, Form("#bf{Cent. %i - %i%s}",                cent_min[3], cent_max[3], "%") );
    //latex3.DrawLatex(0.25, 0.46, Form("#bf{Cent. %i - %i%s}",                cent_min[3], cent_max[3], "%") );
    //latex3.DrawLatex(0.25, 0.38, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu3,  100.*bgmu3e, "%") );
    //latex3.DrawLatex(0.25, 0.30, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd3,   100.*bgd3e,  "%") );
    //latex3.DrawLatex(0.25, 0.22, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps03,  100.*eps03e, "%") );
    //latex3.DrawLatex(0.25, 0.14, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn3,    100.*kn3e,   "%") );
    //latex3.DrawLatex(0.25, 0.06, Form("#alpha = %.1f #pm %.1f",              alpha3,      alpha3e) );

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
    double bgmu5  = fBG[5]->GetParameter(1);
    double bgmu5e = fBG[5]->GetParError(1);
    double bgd5   = fBG[5]->GetParameter(2);
    double bgd5e  = fBG[5]->GetParError(2);

    latex3.DrawLatex(0.09, 0.06, Form("#bf{Cent. %i - %i%s}",                cent_min[5], cent_max[5], "%") );
    //latex3.DrawLatex(0.09, 0.46, Form("#bf{Cent. %i - %i%s}",                cent_min[5], cent_max[5], "%") );
    //latex3.DrawLatex(0.09, 0.38, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu5,  100.*bgmu5e, "%") );
    //latex3.DrawLatex(0.09, 0.30, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd5,   100.*bgd5e,  "%") );
    //latex3.DrawLatex(0.09, 0.22, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps05,  100.*eps05e, "%") );
    //latex3.DrawLatex(0.09, 0.14, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn5,    100.*kn5e,   "%") );
    //latex3.DrawLatex(0.09, 0.06, Form("#alpha = %.1f #pm %.1f",              alpha5,      alpha5e) );


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
    double bgmu7  = fBG[7]->GetParameter(1);
    double bgmu7e = fBG[7]->GetParError(1);
    double bgd7   = fBG[7]->GetParameter(2);
    double bgd7e  = fBG[7]->GetParError(2);

    latex3.DrawLatex(0.25, 0.28, Form("#bf{Cent. %i - %i%s}",                cent_min[7], cent_max[7], "%") );
    //latex3.DrawLatex(0.25, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[7], cent_max[7], "%") );
    //latex3.DrawLatex(0.25, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu7,  100.*bgmu7e, "%") );
    //latex3.DrawLatex(0.25, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd7,   100.*bgd7e,  "%") );
    //latex3.DrawLatex(0.25, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps07,  100.*eps07e, "%") );
    //latex3.DrawLatex(0.25, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn7,    100.*kn7e,   "%") );
    //latex3.DrawLatex(0.25, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha7,      alpha7e) );

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
    double bgmu9  = fBG[9]->GetParameter(1);
    double bgmu9e = fBG[9]->GetParError(1);
    double bgd9   = fBG[9]->GetParameter(2);
    double bgd9e  = fBG[9]->GetParError(2);


    latex3.DrawLatex(0.09, 0.28, Form("#bf{Cent. %i - %i%s}",                cent_min[9], cent_max[9], "%") );
    //latex3.DrawLatex(0.09, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[9], cent_max[9], "%") );
    //latex3.DrawLatex(0.09, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu9,  100.*bgmu9e, "%") );
    //latex3.DrawLatex(0.09, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd9,   100.*bgd9e,  "%") );
    //latex3.DrawLatex(0.09, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps09,  100.*eps09e, "%") );
    //latex3.DrawLatex(0.09, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn9,    100.*kn9e,   "%") );
    //latex3.DrawLatex(0.09, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha9,      alpha9e) );

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
    double bgmu11  = fBG[11]->GetParameter(1);
    double bgmu11e = fBG[11]->GetParError(1);
    double bgd11   = fBG[11]->GetParameter(2);
    double bgd11e  = fBG[11]->GetParError(2);

    latex3.DrawLatex(0.09, 0.28, Form("#bf{Cent. %i - %i%s}",                cent_min[11], cent_max[11], "%") );
    //latex3.DrawLatex(0.09, 0.63, Form("#bf{Cent. %i - %i%s}",                cent_min[11], cent_max[11], "%") );
    //latex3.DrawLatex(0.09, 0.56, Form("v_{n}^{RP} = %.1f #pm %.1f (%s)",     100.*bgmu11,  100.*bgmu11e, "%") );
    //latex3.DrawLatex(0.09, 0.49, Form("#delta_{v_{n}} = %.1f #pm %.1f (%s)", 100.*bgd11,   100.*bgd11e,  "%") );
    //latex3.DrawLatex(0.09, 0.42, Form("#epsilon_{0} = %.1f #pm %.1f (%s)",   100.*eps011,  100.*eps011e, "%") );
    //latex3.DrawLatex(0.09, 0.35, Form("k_{n} = %.1f #pm %.1f (%s)",          100.*kn11,    100.*kn11e,   "%") );
    //latex3.DrawLatex(0.09, 0.28, Form("#alpha = %.1f #pm %.1f",              alpha11,      alpha11e) );

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
    cUnfoldDistsBig->SaveAs(Form("plots/unfolding/finalUnfFit_Cubic_v%i.pdf",norder_));

    for(int icent = 3; icent < NCENT; icent++) std::cout << "Cent = " << icent << "\t Chi2/NDF = " << Form("%.1f", fEllP[icent]->GetChisquare()/fEllP[icent]->GetNDF()) << std::endl;
    for(int icent = 3; icent < NCENT; icent++) std::cout << "Cent = " << icent << "\tkn+knPr = " << fEllP[icent]->GetParameter(2) + fEllP[icent]->GetParameter(3) << std::endl; 

  } //-- End if( !TEST )

}

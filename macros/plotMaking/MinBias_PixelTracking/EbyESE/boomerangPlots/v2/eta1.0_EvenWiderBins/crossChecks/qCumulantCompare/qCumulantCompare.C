#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void qCumulantCompare(){

  const int norder_ = 2;

  const int NQCENT = 12;
  int quanCent_min[NQCENT] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  int quanCent_max[NQCENT] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  double quanCentBinCenter[NQCENT]; 

  double g1eMin   = -1.0;
  double g1eMax   = 0.5;
  double ratioMin = 0.9;
  double ratioMax = 1.03;

  double cumuMin = 0.0;
  double cumuMax = 0.15;

  double ratioMinVn8Vn6 = 0.97;
  double ratioMaxVn8Vn6 = 1.002;

  double rMinC = 0.4;
  double rMaxC = 1.6;
  double rMinR = 0.95;
  double rMaxR = 1.05;
  double rMinG = -10;
  double rMaxG = 10.;

  //-- Quan's values
  double quanVn2[NQCENT];
  double quanVn4[NQCENT];
  double quanVn6[NQCENT];
  double quanVn8[NQCENT];

  double quanVn2_err[NQCENT];
  double quanVn4_err[NQCENT];
  double quanVn6_err[NQCENT];
  double quanVn8_err[NQCENT];

  double quanVn6Vn4[NQCENT];
  double quanVn8Vn4[NQCENT];
  double quanVn8Vn6[NQCENT];
  double quanG1E[NQCENT];

  double quanVn6Vn4_err[NQCENT];
  double quanVn8Vn4_err[NQCENT];
  double quanVn8Vn6_err[NQCENT];
  double quanG1E_err[NQCENT];

  TFile * fQuan;
  TGraphErrors * grQuanVn2;
  TGraphErrors * grQuanVn4;
  TGraphErrors * grQuanVn6;
  TGraphErrors * grQuanVn8;

  TGraphErrors * grQuanVn6Vn4;
  TGraphErrors * grQuanVn8Vn4;
  TGraphErrors * grQuanVn8Vn6;
  TGraphErrors * grQuanG1E;

  //-- EbyE values
  TFile * fEbyE;
  TGraphErrors * grEbyEVn2;
  TGraphErrors * grEbyEVn4;
  TGraphErrors * grEbyEVn6;
  TGraphErrors * grEbyEVn8;

  TGraphErrors * grEbyEVn2_sys;
  TGraphErrors * grEbyEVn4_sys;
  TGraphErrors * grEbyEVn6_sys;
  TGraphErrors * grEbyEVn8_sys;

  TGraphErrors * grEbyEVn6Vn4;
  TGraphErrors * grEbyEVn8Vn4;
  TGraphErrors * grEbyEVn8Vn6;
  TGraphErrors * grEbyEG1E;

  TGraphErrors * grEbyEVn6Vn4_sys;
  TGraphErrors * grEbyEVn8Vn4_sys;
  TGraphErrors * grEbyEVn8Vn6_sys;
  TGraphErrors * grEbyEG1E_sys;

  TLatex latex;

  //-- Ratio to EbyE
  TGraphErrors * grVn2_RatioQuanToEbyE;
  TGraphErrors * grVn4_RatioQuanToEbyE;
  TGraphErrors * grVn6_RatioQuanToEbyE;
  TGraphErrors * grVn8_RatioQuanToEbyE;

  TGraphErrors * grVn6Vn4_RatioQuanToEbyE;
  TGraphErrors * grVn8Vn4_RatioQuanToEbyE;
  TGraphErrors * grVn8Vn6_RatioQuanToEbyE;
  TGraphErrors * grG1E_RatioQuanToEbyE;

  double Vn2_RatioQuanToEbyE[NQCENT];
  double Vn4_RatioQuanToEbyE[NQCENT];
  double Vn6_RatioQuanToEbyE[NQCENT];
  double Vn8_RatioQuanToEbyE[NQCENT];

  double Vn2_RatioQuanToEbyE_err[NQCENT];
  double Vn4_RatioQuanToEbyE_err[NQCENT];
  double Vn6_RatioQuanToEbyE_err[NQCENT];
  double Vn8_RatioQuanToEbyE_err[NQCENT];

  double Vn6Vn4_RatioQuanToEbyE[NQCENT];
  double Vn8Vn4_RatioQuanToEbyE[NQCENT];
  double Vn8Vn6_RatioQuanToEbyE[NQCENT];
  double G1E_RatioQuanToEbyE[NQCENT];

  double Vn6Vn4_RatioQuanToEbyE_err[NQCENT];
  double Vn8Vn4_RatioQuanToEbyE_err[NQCENT];
  double Vn8Vn6_RatioQuanToEbyE_err[NQCENT];
  double G1E_RatioQuanToEbyE_err[NQCENT];

  TFile * fSave;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Get EbyE physics
  fEbyE = new TFile("../../systematicStudies/PhysicsResults.root");
  grEbyEVn2 = (TGraphErrors*) fEbyE->Get("grVn2Raw");
  grEbyEVn4 = (TGraphErrors*) fEbyE->Get("grVn4Raw");
  grEbyEVn6 = (TGraphErrors*) fEbyE->Get("grVn6Raw");
  grEbyEVn8 = (TGraphErrors*) fEbyE->Get("grVn8Raw");

  grEbyEVn2_sys = (TGraphErrors*) fEbyE->Get("grVn2RawSys");
  grEbyEVn4_sys = (TGraphErrors*) fEbyE->Get("grVn4RawSys");
  grEbyEVn6_sys = (TGraphErrors*) fEbyE->Get("grVn6RawSys");
  grEbyEVn8_sys = (TGraphErrors*) fEbyE->Get("grVn8RawSys");

  grEbyEVn2->GetYaxis()->SetTitle( Form("v_{%i}{2}", norder_) );
  grEbyEVn4->GetYaxis()->SetTitle( Form("v_{%i}{4}", norder_) );
  grEbyEVn6->GetYaxis()->SetTitle( Form("v_{%i}{6}", norder_) );
  grEbyEVn8->GetYaxis()->SetTitle( Form("v_{%i}{8}", norder_) );

  grEbyEVn2_sys->GetYaxis()->SetTitle( Form("v_{%i}{2}", norder_) );
  grEbyEVn4_sys->GetYaxis()->SetTitle( Form("v_{%i}{4}", norder_) );
  grEbyEVn6_sys->GetYaxis()->SetTitle( Form("v_{%i}{6}", norder_) );
  grEbyEVn8_sys->GetYaxis()->SetTitle( Form("v_{%i}{8}", norder_) );

  grEbyEVn6Vn4 = (TGraphErrors*) fEbyE->Get("grvn6vn4Ratio");
  grEbyEVn8Vn4 = (TGraphErrors*) fEbyE->Get("grvn8vn4Ratio");
  grEbyEVn8Vn6 = (TGraphErrors*) fEbyE->Get("grvn8vn6Ratio");
  grEbyEG1E    = (TGraphErrors*) fEbyE->Get("grGamma1Exp");

  grEbyEVn6Vn4_sys = (TGraphErrors*) fEbyE->Get("grvn6vn4RatioSys");
  grEbyEVn8Vn4_sys = (TGraphErrors*) fEbyE->Get("grvn8vn4RatioSys");
  grEbyEVn8Vn6_sys = (TGraphErrors*) fEbyE->Get("grvn8vn6RatioSys");
  grEbyEG1E_sys    = (TGraphErrors*) fEbyE->Get("grGamma1ExpSys");

  //-- Get Quan's values
  fQuan = new TFile("outGraph_PR_pt1_3_eta1p0.root");
  //fQuan = new TFile("outGraph.root");
  grQuanVn2 = (TGraphErrors*) fQuan->Get("gr_vnCentC_2_0");
  grQuanVn2->SetLineColor(1);
  grQuanVn2->SetMarkerColor(1);
  grQuanVn2->SetMarkerStyle(20);
  grQuanVn2->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grQuanVn4 = (TGraphErrors*) fQuan->Get("gr_vnCentC_2_1");
  grQuanVn4->SetLineColor(1);
  grQuanVn4->SetMarkerColor(1);
  grQuanVn4->SetMarkerStyle(20);
  grQuanVn4->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grQuanVn6 = (TGraphErrors*) fQuan->Get("gr_vnCentC_2_2");
  grQuanVn6->SetLineColor(1);
  grQuanVn6->SetMarkerColor(1);
  grQuanVn6->SetMarkerStyle(20);
  grQuanVn6->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  grQuanVn8 = (TGraphErrors*) fQuan->Get("gr_vnCentC_2_3");
  grQuanVn8->SetLineColor(1);
  grQuanVn8->SetMarkerColor(1);
  grQuanVn8->SetMarkerStyle(20);
  grQuanVn8->GetYaxis()->SetRangeUser(cumuMin, cumuMax);

  double c_err[NQCENT];

  //-- fill Quan Arrays
  for(int icent = 0; icent < NQCENT; icent++){

    c_err[icent] = 0.;
    quanCentBinCenter[icent] = ((double)quanCent_max[icent] + (double)quanCent_min[icent])/2.;

    double vn2 = grQuanVn2->GetY()[icent];
    double vn4 = grQuanVn4->GetY()[icent];
    double vn6 = grQuanVn6->GetY()[icent];
    double vn8 = grQuanVn8->GetY()[icent];

    double vn2_err = grQuanVn2->GetErrorY(icent);
    double vn4_err = grQuanVn4->GetErrorY(icent);
    double vn6_err = grQuanVn6->GetErrorY(icent);
    double vn8_err = grQuanVn8->GetErrorY(icent);


    quanVn2[icent] = vn2;
    quanVn4[icent] = vn4;
    quanVn6[icent] = vn6;
    quanVn8[icent] = vn8;

    quanVn2_err[icent] = vn2_err;
    quanVn4_err[icent] = vn4_err;
    quanVn6_err[icent] = vn6_err;
    quanVn8_err[icent] = vn8_err;

    quanVn6Vn4[icent] = vn6 / vn4;
    quanVn8Vn4[icent] = vn8 / vn4;
    quanVn8Vn6[icent] = vn8 / vn6;
    quanG1E[icent]    = -6.*sqrt(2.)*pow(vn4, 2) * (vn4 - vn6) / pow( pow(vn2, 2) - pow(vn4, 2), 3./2. );

    quanVn6Vn4_err[icent] = 0.;
    quanVn8Vn4_err[icent] = 0.;
    quanVn8Vn6_err[icent] = 0.;
    quanG1E_err[icent]    = 0.;

    //-- Ratio quan / EbyE
    double vn2EbyE    = grEbyEVn2->GetY()[icent];
    double vn4EbyE    = grEbyEVn4->GetY()[icent];
    double vn6EbyE    = grEbyEVn6->GetY()[icent];
    double vn8EbyE    = grEbyEVn8->GetY()[icent];
    double vn6vn4EbyE = grEbyEVn6Vn4->GetY()[icent];
    double vn8vn4EbyE = grEbyEVn8Vn4->GetY()[icent];
    double vn8vn6EbyE = grEbyEVn8Vn6->GetY()[icent];
    double g1eEbyE    = grEbyEG1E->GetY()[icent];

    double vn2EbyE_StatErr    = grEbyEVn2->GetErrorY(icent);
    double vn4EbyE_StatErr    = grEbyEVn4->GetErrorY(icent);
    double vn6EbyE_StatErr    = grEbyEVn6->GetErrorY(icent);
    double vn8EbyE_StatErr    = grEbyEVn8->GetErrorY(icent);
    double vn6vn4EbyE_StatErr = grEbyEVn6Vn4->GetErrorY(icent);
    double vn8vn4EbyE_StatErr = grEbyEVn8Vn4->GetErrorY(icent);
    double vn8vn6EbyE_StatErr = grEbyEVn8Vn6->GetErrorY(icent);
    double g1eEbyE_StatErr    = grEbyEG1E->GetErrorY(icent);

    double vn2EbyE_SysErr    = grEbyEVn2_sys->GetErrorY(icent);
    double vn4EbyE_SysErr    = grEbyEVn4_sys->GetErrorY(icent);
    double vn6EbyE_SysErr    = grEbyEVn6_sys->GetErrorY(icent);
    double vn8EbyE_SysErr    = grEbyEVn8_sys->GetErrorY(icent);
    double vn6vn4EbyE_SysErr = grEbyEVn6Vn4_sys->GetErrorY(icent);
    double vn8vn4EbyE_SysErr = grEbyEVn8Vn4_sys->GetErrorY(icent);
    double vn8vn6EbyE_SysErr = grEbyEVn8Vn6_sys->GetErrorY(icent);
    double g1eEbyE_SysErr    = grEbyEG1E_sys->GetErrorY(icent);

    double vn2EbyE_err    = sqrt( pow(vn2EbyE_StatErr,2) + pow(vn2EbyE_SysErr,2) );
    double vn4EbyE_err    = sqrt( pow(vn4EbyE_StatErr,2) + pow(vn4EbyE_SysErr,2) );
    double vn6EbyE_err    = sqrt( pow(vn6EbyE_StatErr,2) + pow(vn6EbyE_SysErr,2) );
    double vn8EbyE_err    = sqrt( pow(vn8EbyE_StatErr,2) + pow(vn8EbyE_SysErr,2) );
    double vn6vn4EbyE_err = sqrt( pow(vn6vn4EbyE_StatErr,2) + pow(vn6vn4EbyE_SysErr,2) );
    double vn8vn4EbyE_err = sqrt( pow(vn8vn4EbyE_StatErr,2) + pow(vn8vn4EbyE_SysErr,2) );
    double vn8vn6EbyE_err = sqrt( pow(vn8vn6EbyE_StatErr,2) + pow(vn8vn6EbyE_SysErr,2) );
    double g1eEbyE_err    = sqrt( pow(g1eEbyE_StatErr,2) + pow(g1eEbyE_SysErr,2) );

    Vn2_RatioQuanToEbyE[icent] = vn2 / vn2EbyE;
    Vn4_RatioQuanToEbyE[icent] = vn4 / vn4EbyE;
    Vn6_RatioQuanToEbyE[icent] = vn6 / vn6EbyE;
    Vn8_RatioQuanToEbyE[icent] = vn8 / vn8EbyE;

    Vn2_RatioQuanToEbyE_err[icent] = sqrt( pow(vn2_err/vn2EbyE,2) + pow(vn2*vn2EbyE_err/vn2EbyE/vn2EbyE,2) );
    Vn4_RatioQuanToEbyE_err[icent] = sqrt( pow(vn4_err/vn4EbyE,2) + pow(vn4*vn4EbyE_err/vn4EbyE/vn4EbyE,2) );
    Vn6_RatioQuanToEbyE_err[icent] = sqrt( pow(vn6_err/vn6EbyE,2) + pow(vn6*vn6EbyE_err/vn6EbyE/vn6EbyE,2) );
    Vn8_RatioQuanToEbyE_err[icent] = sqrt( pow(vn8_err/vn8EbyE,2) + pow(vn8*vn8EbyE_err/vn8EbyE/vn8EbyE,2) );

    if(vn2 == 0 || vn2EbyE == 0){
      Vn2_RatioQuanToEbyE[icent]     = -1;
      Vn2_RatioQuanToEbyE_err[icent] = 0;
    }
    if(vn4 == 0|| vn4EbyE == 0){
      Vn4_RatioQuanToEbyE[icent]     = -1;
      Vn4_RatioQuanToEbyE_err[icent] = 0;
    }
    if(vn6 == 0|| vn6EbyE == 0){
      Vn6_RatioQuanToEbyE[icent]     = -1;
      Vn6_RatioQuanToEbyE_err[icent] = 0;
    }
    if(vn8 == 0|| vn8EbyE == 0){
      Vn8_RatioQuanToEbyE[icent]     = -1;
      Vn8_RatioQuanToEbyE_err[icent] = 0;
    }

    Vn6Vn4_RatioQuanToEbyE[icent] = quanVn6Vn4[icent] / vn6vn4EbyE;
    Vn8Vn4_RatioQuanToEbyE[icent] = quanVn8Vn4[icent] / vn8vn4EbyE;
    Vn8Vn6_RatioQuanToEbyE[icent] = quanVn8Vn6[icent] / vn8vn6EbyE;
    G1E_RatioQuanToEbyE[icent]    = quanG1E[icent] / g1eEbyE;

    Vn6Vn4_RatioQuanToEbyE_err[icent] = fabs( quanVn6Vn4[icent]*vn6vn4EbyE_err/vn6vn4EbyE/vn6vn4EbyE );
    Vn8Vn4_RatioQuanToEbyE_err[icent] = fabs( quanVn8Vn4[icent]*vn8vn4EbyE_err/vn8vn4EbyE/vn8vn4EbyE );
    Vn8Vn6_RatioQuanToEbyE_err[icent] = fabs( quanVn8Vn6[icent]*vn8vn6EbyE_err/vn8vn6EbyE/vn8vn6EbyE );
    G1E_RatioQuanToEbyE_err[icent]    = fabs( quanG1E[icent]*g1eEbyE_err/g1eEbyE/g1eEbyE );

    if( quanVn6Vn4[icent] <= 0 || vn6vn4EbyE <= 0 ){
      Vn6Vn4_RatioQuanToEbyE[icent]     = -1.;
      Vn6Vn4_RatioQuanToEbyE_err[icent] = 0;
    }
    if(quanVn8Vn4[icent] <= 0 || vn8vn4EbyE <= 0){
      Vn8Vn4_RatioQuanToEbyE[icent]     = -1.;
      Vn8Vn4_RatioQuanToEbyE_err[icent] = 0;
    }
    if(quanVn8Vn6[icent] <= 0 || vn8vn6EbyE <= 0){
      Vn8Vn6_RatioQuanToEbyE[icent]     = -1.;
      Vn8Vn6_RatioQuanToEbyE_err[icent] = 0;
    }
    if(quanG1E[icent] < -100 || g1eEbyE < -100){
      G1E_RatioQuanToEbyE[icent]     = -10000;
      G1E_RatioQuanToEbyE_err[icent] = 0;
    }

  }

  //-- Make Quan Cumu Ratios
  grQuanVn6Vn4 = new TGraphErrors(NQCENT, quanCentBinCenter, quanVn6Vn4, c_err, quanVn6Vn4_err);
  grQuanVn8Vn4 = new TGraphErrors(NQCENT, quanCentBinCenter, quanVn8Vn4, c_err, quanVn8Vn4_err);
  grQuanVn8Vn6 = new TGraphErrors(NQCENT, quanCentBinCenter, quanVn8Vn6, c_err, quanVn8Vn6_err);
  grQuanG1E    = new TGraphErrors(NQCENT, quanCentBinCenter, quanG1E,    c_err, quanG1E_err);

  formatGraph(grQuanVn6Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_), 1, 21, "grQuanVn6Vn4");
  formatGraph(grQuanVn8Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_), 1, 34, "grQuanVn8Vn4");
  formatGraph(grQuanVn8Vn6, "Centrality %", ratioMinVn8Vn6, ratioMaxVn8Vn6, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_), 1, 33, "grQuanVn8Vn6");
  formatGraph(grQuanG1E,    "Centrality %", g1eMin,         g1eMax,         "#gamma_{1}^{Exp}",                              1, 20, "grQuanG1E");

  //-- Make ratio quan / EbyE graphs
  grVn2_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn2_RatioQuanToEbyE, c_err, Vn2_RatioQuanToEbyE_err);
  grVn4_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn4_RatioQuanToEbyE, c_err, Vn4_RatioQuanToEbyE_err);
  grVn6_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn6_RatioQuanToEbyE, c_err, Vn6_RatioQuanToEbyE_err);
  grVn8_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn8_RatioQuanToEbyE, c_err, Vn8_RatioQuanToEbyE_err);

  formatGraph(grVn2_RatioQuanToEbyE, "Centrality %", rMinC, rMaxC, "Ratio: QCumu / EbyE", 1, 21, "grVn2_RatioQuanToEbyE");
  formatGraph(grVn4_RatioQuanToEbyE, "Centrality %", rMinC, rMaxC, "Ratio: QCumu / EbyE", 1, 21, "grVn4_RatioQuanToEbyE");
  formatGraph(grVn6_RatioQuanToEbyE, "Centrality %", rMinC, rMaxC, "Ratio: QCumu / EbyE", 1, 21, "grVn6_RatioQuanToEbyE");
  formatGraph(grVn8_RatioQuanToEbyE, "Centrality %", rMinC, rMaxC, "Ratio: QCumu / EbyE", 1, 21, "grVn8_RatioQuanToEbyE");

  grVn6Vn4_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn6Vn4_RatioQuanToEbyE, c_err, Vn6Vn4_RatioQuanToEbyE_err);
  grVn8Vn4_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn8Vn4_RatioQuanToEbyE, c_err, Vn8Vn4_RatioQuanToEbyE_err);
  grVn8Vn6_RatioQuanToEbyE = new TGraphErrors(NQCENT, quanCentBinCenter, Vn8Vn6_RatioQuanToEbyE, c_err, Vn8Vn6_RatioQuanToEbyE_err);
  grG1E_RatioQuanToEbyE    = new TGraphErrors(NQCENT, quanCentBinCenter, G1E_RatioQuanToEbyE,    c_err, G1E_RatioQuanToEbyE_err);

  formatGraph(grVn6Vn4_RatioQuanToEbyE, "Centrality %", rMinR, rMaxR, "Ratio: QCumu / EbyE", 1, 21, "grVn6Vn4_RatioQuanToEbyE");
  formatGraph(grVn8Vn4_RatioQuanToEbyE, "Centrality %", rMinR, rMaxR, "Ratio: QCumu / EbyE", 1, 21, "grVn8Vn4_RatioQuanToEbyE");
  formatGraph(grVn8Vn6_RatioQuanToEbyE, "Centrality %", rMinR, rMaxR, "Ratio: QCumu / EbyE", 1, 21, "grVn8Vn6_RatioQuanToEbyE");
  formatGraph(grG1E_RatioQuanToEbyE,    "Centrality %", rMinG, rMaxG, "Ratio: QCumu / EbyE", 1, 21, "grG1E_RatioQuanToEbyE");


  //-- Draw cumus
  TLegend * legvn2 = new TLegend(0.4720, 0.1909, 0.9909, 0.3834);
  legvn2->SetBorderSize(0);
  legvn2->SetFillStyle(0);
  legvn2->AddEntry(grEbyEVn2, "EbyE", "lp");
  legvn2->AddEntry(grQuanVn2, "Q-Cumulant", "lp");

  TLegend * legvn4 = new TLegend(0.4720, 0.1909, 0.9909, 0.3834);
  legvn4->SetBorderSize(0);
  legvn4->SetFillStyle(0);
  legvn4->AddEntry(grEbyEVn4, "EbyE", "lp");
  legvn4->AddEntry(grQuanVn4, "Q-Cumulant", "lp");

  TLegend * legvn6 = new TLegend(0.4720, 0.1909, 0.9909, 0.3834);
  legvn6->SetBorderSize(0);
  legvn6->SetFillStyle(0);
  legvn6->AddEntry(grEbyEVn6, "EbyE", "lp");
  legvn6->AddEntry(grQuanVn6, "Q-Cumulant", "lp");

  TLegend * legvn8 = new TLegend(0.4720, 0.1909, 0.9909, 0.3834);
  legvn8->SetBorderSize(0);
  legvn8->SetFillStyle(0);
  legvn8->AddEntry(grEbyEVn8, "EbyE", "lp");
  legvn8->AddEntry(grQuanVn8, "Q-Cumulant", "lp");

  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 2000, 500);
  cCumu->Divide(4, 1);

  cCumu->cd(1);
  grEbyEVn2_sys->Draw("apE2");
  grEbyEVn2->Draw("psame");
  grQuanVn2->Draw("psame");
  legvn2->Draw("same");

  cCumu->cd(2);
  grEbyEVn4_sys->Draw("apE2");
  grEbyEVn4->Draw("psame");
  grQuanVn4->Draw("psame");
  legvn4->Draw("same");

  cCumu->cd(3);
  grEbyEVn6_sys->Draw("apE2");
  grEbyEVn6->Draw("psame");
  grQuanVn6->Draw("psame");
  legvn6->Draw("same");

  cCumu->cd(4);
  grEbyEVn8_sys->Draw("apE2");
  grEbyEVn8->Draw("psame");
  grQuanVn8->Draw("psame");
  legvn8->Draw("same");
  cCumu->SaveAs("CumuRaw_QCumuComp.pdf");

  TCanvas * cCumu_RQtoE = new TCanvas("cCumu_RQtoE", "cCumu_RQtoE", 2000, 500);
  cCumu_RQtoE->Divide(4, 1);
  cCumu_RQtoE->cd(1);
  grVn2_RatioQuanToEbyE->Draw("ap");
  cCumu_RQtoE->cd(2);
  grVn4_RatioQuanToEbyE->Draw("ap");
  cCumu_RQtoE->cd(3);
  grVn6_RatioQuanToEbyE->Draw("ap");
  cCumu_RQtoE->cd(4);
  grVn8_RatioQuanToEbyE->Draw("ap");
  cCumu_RQtoE->SaveAs("CumuRaw_QCumuCompRatio.pdf");


  //-- Draw cumu ratio
  TLegend * leg1 = new TLegend(0.5601, 0.1930, 0.9616, 0.3314);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(grEbyEVn6Vn4, "EbyE", "lp");
  leg1->AddEntry(grQuanVn6Vn4, "Q-Cumulant", "lp");

  TLegend * leg2 =new TLegend(0.5601, 0.1930, 0.9616, 0.3314);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(grEbyEVn8Vn4, "EbyE", "lp");
  leg2->AddEntry(grQuanVn8Vn4, "Q-Cumulant", "lp");

  TLegend * leg3 =new TLegend(0.5601, 0.1930, 0.9616, 0.3314);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->AddEntry(grEbyEVn8Vn6, "EbyE", "lp");
  leg3->AddEntry(grQuanVn8Vn6, "Q-Cumulant", "lp");


  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  cCumuRatio->cd(1);
  grQuanVn6Vn4->Draw("ap");
  grEbyEVn6Vn4_sys->Draw("pE2same");
  grEbyEVn6Vn4->Draw("psame");
  grQuanVn6Vn4->Draw("psame");
  leg1->Draw("same");
  cCumuRatio->cd(2);
  grQuanVn8Vn4->Draw("ap");
  grEbyEVn8Vn4_sys->Draw("pE2same");
  grEbyEVn8Vn4->Draw("psame");
  grQuanVn8Vn4->Draw("psame");
  leg2->Draw("same");
  cCumuRatio->cd(3);
  grQuanVn8Vn6->Draw("ap");
  grEbyEVn8Vn6_sys->Draw("pE2same");
  grEbyEVn8Vn6->Draw("psame");
  grQuanVn8Vn6->Draw("psame");
  leg3->Draw("same");
  cCumuRatio->SaveAs("CumuRatio_QCumuComp.pdf");

  TCanvas * cCumuRatio_RQtoE = new TCanvas("cCumuRatio_RQtoE", "cCumuRatio_RQtoE", 1500, 500);
  cCumuRatio_RQtoE->Divide(3, 1);
  cCumuRatio_RQtoE->cd(1);
  grVn6Vn4_RatioQuanToEbyE->Draw("ap");
  cCumuRatio_RQtoE->cd(2);
  grVn8Vn4_RatioQuanToEbyE->Draw("ap");
  cCumuRatio_RQtoE->cd(3);
  grVn8Vn6_RatioQuanToEbyE->Draw("ap");
  cCumuRatio_RQtoE->SaveAs("CumuRatio_QCumuCompRatio.pdf");


  //-- Draw G1E
  TLegend * legG = new TLegend(0.6754, 0.1928, 0.9637, 0.2924);
  legG->SetBorderSize(0);
  legG->SetFillStyle(0);
  legG->AddEntry(grEbyEG1E, "EbyE", "lp");
  legG->AddEntry(grQuanG1E, "Q-Cumulant", "lp");

  TCanvas * cG1E = new TCanvas("cG1E", "cG1E", 500, 500);
  cG1E->cd();
  grQuanG1E->Draw("ap");
  grEbyEG1E_sys->Draw("psameE2");
  grEbyEG1E->Draw("psame");
  grQuanG1E->Draw("psame");
  legG->Draw("same");
  cG1E->SaveAs("G1E_QCumuComp.pdf");

  TCanvas * cG1E_RQtoE = new TCanvas("cG1E_RQtoE", "cG1E_RQtoE", 500, 500);
  cG1E_RQtoE->cd();
  grG1E_RatioQuanToEbyE->Draw("ap");
  cG1E_RQtoE->SaveAs("G1Eo_QCumuCompRatio.pdf");

  //-- Save Quan's files
  fSave = new TFile("fQuan_pt1_3_eta2p4.root", "recreate");
  fSave->cd();
  grQuanVn6Vn4->Write("grQuanVn6Vn4");
  grQuanVn8Vn4->Write("grQuanVn8Vn4");
  grQuanVn8Vn6->Write("grQuanVn8Vn6");


}

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


void sysTails(){

  const int norder_ = 2;

  double vn6vn4Min = 0.95;
  double vn6vn4Max = 1.03;
  double vn8vn4Min = 0.95;
  double vn8vn4Max = 1.03;
  double vn8vn6Min = 0.99;
  double vn8vn6Max = 1.002;
  double g1eMin    = -1.0;
  double g1eMax    = 0.5;

  double rMin = 0.98;
  double rMax = 1.02;
  double g1rMin = 0.5;
  double g1rMax = 1.5;

  TFile * fSysDist;
  TH1D * hUnfoldDefault[NCENT];
  TH1D * hUnfoldUpperSys[NCENT];
  TH1D * hUnfoldLowerSys[NCENT];

  double vn6vn4Default[NCENT];
  double vn8vn4Default[NCENT];
  double vn8vn6Default[NCENT];
  double g1eDefault[NCENT];

  double vn6vn4UpperSys[NCENT];
  double vn8vn4UpperSys[NCENT];
  double vn8vn6UpperSys[NCENT];
  double g1eUpperSys[NCENT];

  double vn6vn4LowerSys[NCENT];
  double vn8vn4LowerSys[NCENT];
  double vn8vn6LowerSys[NCENT];
  double g1eLowerSys[NCENT];

  double vn6vn4UpperSys_RatioToDefault[NCENT];
  double vn8vn4UpperSys_RatioToDefault[NCENT];
  double vn8vn6UpperSys_RatioToDefault[NCENT];
  double g1eUpperSys_RatioToDefault[NCENT];

  double vn6vn4LowerSys_RatioToDefault[NCENT];
  double vn8vn4LowerSys_RatioToDefault[NCENT];
  double vn8vn6LowerSys_RatioToDefault[NCENT];
  double g1eLowerSys_RatioToDefault[NCENT];

  TGraph * grVn6Vn4Default;
  TGraph * grVn6Vn4UpperSys;
  TGraph * grVn6Vn4LowerSys;
  TGraph * grVn6Vn4UpperSys_RatioToDefault;
  TGraph * grVn6Vn4LowerSys_RatioToDefault;

  TGraph * grVn8Vn4Default;
  TGraph * grVn8Vn4UpperSys;
  TGraph * grVn8Vn4LowerSys;
  TGraph * grVn8Vn4UpperSys_RatioToDefault;
  TGraph * grVn8Vn4LowerSys_RatioToDefault;

  TGraph * grVn8Vn6Default;
  TGraph * grVn8Vn6UpperSys;
  TGraph * grVn8Vn6LowerSys;
  TGraph * grVn8Vn6UpperSys_RatioToDefault;
  TGraph * grVn8Vn6LowerSys_RatioToDefault;

  TGraph * grG1EDefault;
  TGraph * grG1EUpperSys;
  TGraph * grG1ELowerSys;
  TGraph * grG1EUpperSys_RatioToDefault;
  TGraph * grG1ELowerSys_RatioToDefault;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Grab file
  fSysDist = new TFile( Form("../SysUnfoldDistns_v%i.root", norder_) );

  //-- Get hists
  for(int icent = 0; icent < NCENT; icent++){

    hUnfoldDefault[icent]  = (TH1D*) fSysDist->Get( Form("hFinalUnfoldStatAndSys_c%i", icent) );
    hUnfoldUpperSys[icent] = (TH1D*) hUnfoldDefault[icent]->Clone( Form("hUnfoldUpperSys_c%i", icent) );
    hUnfoldLowerSys[icent] = (TH1D*) hUnfoldDefault[icent]->Clone( Form("hUnfoldLowerSys_c%i",icent) );

    hUnfoldDefault[icent]  -> SetLineColor(1);
    hUnfoldUpperSys[icent] -> SetLineColor(2);
    hUnfoldLowerSys[icent] -> SetLineColor(4);

    hUnfoldDefault[icent]  -> SetMarkerColor(1);
    hUnfoldUpperSys[icent] -> SetMarkerColor(2);
    hUnfoldLowerSys[icent] -> SetMarkerColor(4);

    hUnfoldDefault[icent]  -> GetXaxis() -> SetRange(1, hUnfoldDefault[icent]->FindBin(0.2));
    hUnfoldUpperSys[icent] -> GetXaxis() -> SetRange(1, hUnfoldDefault[icent]->FindBin(0.2));
    hUnfoldLowerSys[icent] -> GetXaxis() -> SetRange(1, hUnfoldDefault[icent]->FindBin(0.2));

    //-- Alter the upper and lower histos bin contents
    hUnfoldUpperSys[icent]->Reset();
    hUnfoldLowerSys[icent]->Reset();

    for(int i = 1; i <= NBins; i++){
      double binc = hUnfoldDefault[icent]->GetBinContent(i);
      double bine = hUnfoldDefault[icent]->GetBinError(i);
      if( binc == 0 ) continue;

      hUnfoldUpperSys[icent]->SetBinContent(i, binc+bine);
      hUnfoldLowerSys[icent]->SetBinContent(i, binc-bine);
      hUnfoldUpperSys[icent]->SetBinError(i, 0);
      hUnfoldLowerSys[icent]->SetBinError(i, 0);
    } //-- End bin loop

    //-- Raw Values
    //-- Default
    EbyECumu cumuDefault(hUnfoldDefault[icent]);
    double vn4Default  = cumuDefault.GetCumu_vn4();
    double vn6Default  = cumuDefault.GetCumu_vn6();
    double vn8Default  = cumuDefault.GetCumu_vn8();
    double g1exDefault = cumuDefault.GetGamma1Exp();

    g1eDefault[icent] = g1exDefault;
    if(vn4Default == 0 || vn6Default == 0) vn6vn4Default[icent] = -1000.;
    else                                   vn6vn4Default[icent] = vn6Default / vn4Default;
    if(vn4Default == 0 || vn8Default == 0) vn8vn4Default[icent] = -1000.;
    else                                   vn8vn4Default[icent] = vn8Default / vn4Default;
    if(vn6Default == 0 || vn8Default == 0) vn8vn6Default[icent] = -1000.;
    else                                   vn8vn6Default[icent] = vn8Default / vn6Default;

    //-- UpperSys
    EbyECumu cumuUpperSys(hUnfoldUpperSys[icent]);
    double vn4UpperSys  = cumuUpperSys.GetCumu_vn4();
    double vn6UpperSys  = cumuUpperSys.GetCumu_vn6();
    double vn8UpperSys  = cumuUpperSys.GetCumu_vn8();
    double g1exUpperSys = cumuUpperSys.GetGamma1Exp();

    g1eUpperSys[icent] = g1exUpperSys;
    if(vn4UpperSys == 0 || vn6UpperSys == 0) vn6vn4UpperSys[icent] = -1000.;
    else                                     vn6vn4UpperSys[icent] = vn6UpperSys / vn4UpperSys;
    if(vn4UpperSys == 0 || vn8UpperSys == 0) vn8vn4UpperSys[icent] = -1000.;
    else                                     vn8vn4UpperSys[icent] = vn8UpperSys / vn4UpperSys;
    if(vn6UpperSys == 0 || vn8UpperSys == 0) vn8vn6UpperSys[icent] = -1000.;
    else                                     vn8vn6UpperSys[icent] = vn8UpperSys / vn6UpperSys;

    //-- LowerSys
    EbyECumu cumuLowerSys(hUnfoldLowerSys[icent]);
    double vn4LowerSys  = cumuLowerSys.GetCumu_vn4();
    double vn6LowerSys  = cumuLowerSys.GetCumu_vn6();
    double vn8LowerSys  = cumuLowerSys.GetCumu_vn8();
    double g1exLowerSys = cumuLowerSys.GetGamma1Exp();

    g1eLowerSys[icent] = g1exLowerSys;
    if(vn4LowerSys == 0 || vn6LowerSys == 0) vn6vn4LowerSys[icent] = -1000.;
    else                                     vn6vn4LowerSys[icent] = vn6LowerSys / vn4LowerSys;
    if(vn4LowerSys == 0 || vn8LowerSys == 0) vn8vn4LowerSys[icent] = -1000.;
    else                                     vn8vn4LowerSys[icent] = vn8LowerSys / vn4LowerSys;
    if(vn6LowerSys == 0 || vn8LowerSys == 0) vn8vn6LowerSys[icent] = -1000.;
    else                                     vn8vn6LowerSys[icent] = vn8LowerSys / vn6LowerSys;

    //-- Ratio to Default
    //-- UpperSys
    if(g1exUpperSys == -10000. || g1exDefault == -10000.) g1eUpperSys_RatioToDefault[icent] = -10000.;
    else                                                  g1eUpperSys_RatioToDefault[icent] = g1eUpperSys[icent] / g1eDefault[icent];
    if(vn4Default == 0 || vn6Default == 0 || vn4UpperSys == 0 || vn6UpperSys == 0) vn6vn4UpperSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn6vn4UpperSys_RatioToDefault[icent] = vn6vn4UpperSys[icent] / vn6vn4Default[icent];
    if(vn4Default == 0 || vn8Default == 0 || vn4UpperSys == 0 || vn8UpperSys == 0) vn8vn4UpperSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn8vn4UpperSys_RatioToDefault[icent] = vn8vn4UpperSys[icent] / vn8vn4Default[icent];
    if(vn6Default == 0 || vn8Default == 0 || vn6UpperSys == 0 || vn8UpperSys == 0) vn8vn6UpperSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn8vn6UpperSys_RatioToDefault[icent] = vn8vn6UpperSys[icent] / vn8vn6Default[icent];

    //-- LowerSys
    if(g1exLowerSys == -10000. || g1exDefault == -10000.) g1eLowerSys_RatioToDefault[icent] = -10000.;
    else                                                  g1eLowerSys_RatioToDefault[icent] = g1eLowerSys[icent] / g1eDefault[icent];
    if(vn4Default == 0 || vn6Default == 0 || vn4LowerSys == 0 || vn6LowerSys == 0) vn6vn4LowerSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn6vn4LowerSys_RatioToDefault[icent] = vn6vn4LowerSys[icent] / vn6vn4Default[icent];
    if(vn4Default == 0 || vn8Default == 0 || vn4LowerSys == 0 || vn8LowerSys == 0) vn8vn4LowerSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn8vn4LowerSys_RatioToDefault[icent] = vn8vn4LowerSys[icent] / vn8vn4Default[icent];
    if(vn6Default == 0 || vn8Default == 0 || vn6LowerSys == 0 || vn8LowerSys == 0) vn8vn6LowerSys_RatioToDefault[icent] = -1000.;
    else                                                                           vn8vn6LowerSys_RatioToDefault[icent] = vn8vn6LowerSys[icent] / vn8vn6Default[icent];

  } //-- End cent loop

  //-- Make TGraphs
  grVn6Vn4Default                 = new TGraph(NCENT, centBinCenter, vn6vn4Default);
  grVn6Vn4UpperSys                = new TGraph(NCENT, centBinCenter, vn6vn4UpperSys);
  grVn6Vn4LowerSys                = new TGraph(NCENT, centBinCenter, vn6vn4LowerSys);
  grVn6Vn4UpperSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn6vn4UpperSys_RatioToDefault);
  grVn6Vn4LowerSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn6vn4LowerSys_RatioToDefault);

  formatGraph(grVn6Vn4Default,                 "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_),       4, 21, "grVn6Vn4Default");
  formatGraph(grVn6Vn4UpperSys,                "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_),       4, 24, "grVn6Vn4UpperSys");
  formatGraph(grVn6Vn4LowerSys,                "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_),       4, 26, "grVn6Vn4LowerSys");
  formatGraph(grVn6Vn4UpperSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 4, 24, "grVn6Vn4UpperSys_RatioToDefault");
  formatGraph(grVn6Vn4LowerSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{6} / v_{%i}{4} Ratio", norder_, norder_), 4, 26, "grVn6Vn4LowerSys_RatioToDefault");

  grVn8Vn4Default                 = new TGraph(NCENT, centBinCenter, vn8vn4Default);
  grVn8Vn4UpperSys                = new TGraph(NCENT, centBinCenter, vn8vn4UpperSys);
  grVn8Vn4LowerSys                = new TGraph(NCENT, centBinCenter, vn8vn4LowerSys);
  grVn8Vn4UpperSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn8vn4UpperSys_RatioToDefault);
  grVn8Vn4LowerSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn8vn4LowerSys_RatioToDefault);

  formatGraph(grVn8Vn4Default,                 "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_),       kGreen+2, 34, "grVn8Vn4Default");
  formatGraph(grVn8Vn4UpperSys,                "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_),       kGreen+2, 25, "grVn8Vn4UpperSys");
  formatGraph(grVn8Vn4LowerSys,                "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_),       kGreen+2, 27, "grVn8Vn4LowerSys");
  formatGraph(grVn8Vn4UpperSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), kGreen+2, 25, "grVn8Vn4UpperSys_RatioToDefault");
  formatGraph(grVn8Vn4LowerSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{8} / v_{%i}{4} Ratio", norder_, norder_), kGreen+2, 27, "grVn8Vn4LowerSys_RatioToDefault");

  grVn8Vn6Default                 = new TGraph(NCENT, centBinCenter, vn8vn6Default);
  grVn8Vn6UpperSys                = new TGraph(NCENT, centBinCenter, vn8vn6UpperSys);
  grVn8Vn6LowerSys                = new TGraph(NCENT, centBinCenter, vn8vn6LowerSys);
  grVn8Vn6UpperSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn8vn6UpperSys_RatioToDefault);
  grVn8Vn6LowerSys_RatioToDefault = new TGraph(NCENT, centBinCenter, vn8vn6LowerSys_RatioToDefault);

  formatGraph(grVn8Vn6Default,                 "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_),       kViolet-1, 33, "grVn8Vn6Default");
  formatGraph(grVn8Vn6UpperSys,                "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_),       kViolet-1, 28, "grVn8Vn6UpperSys");
  formatGraph(grVn8Vn6LowerSys,                "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_),       kViolet-1, 32, "grVn8Vn6LowerSys");
  formatGraph(grVn8Vn6UpperSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kViolet-1, 28, "grVn8Vn6UpperSys_RatioToDefault");
  formatGraph(grVn8Vn6LowerSys_RatioToDefault, "Centrality %", rMin,      rMax,      Form("v_{%i}{8} / v_{%i}{6} Ratio", norder_, norder_), kViolet-1, 32, "grVn8Vn6LowerSys_RatioToDefault");

  grG1EDefault                 = new TGraph(NCENT, centBinCenter, g1eDefault);
  grG1EUpperSys                = new TGraph(NCENT, centBinCenter, g1eUpperSys);
  grG1ELowerSys                = new TGraph(NCENT, centBinCenter, g1eLowerSys);
  grG1EUpperSys_RatioToDefault = new TGraph(NCENT, centBinCenter, g1eUpperSys_RatioToDefault);
  grG1ELowerSys_RatioToDefault = new TGraph(NCENT, centBinCenter, g1eLowerSys_RatioToDefault);

  formatGraph(grG1EDefault,                 "Centrality %", g1eMin, g1eMax, "#gamma_{1}^{exp}",       2, 20, "grG1EDefault");
  formatGraph(grG1EUpperSys,                "Centrality %", g1eMin, g1eMax, "#gamma_{1}^{exp}",       2, 25, "grG1EUpperSys");
  formatGraph(grG1ELowerSys,                "Centrality %", g1eMin, g1eMax, "#gamma_{1}^{exp}",       2, 26, "grG1ELowerSys");
  formatGraph(grG1EUpperSys_RatioToDefault, "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio", 2, 25, "grG1EUpperSys_RatioToDefault");
  formatGraph(grG1ELowerSys_RatioToDefault, "Centrality %", g1rMin, g1rMax, "#gamma_{1}^{exp} Ratio", 2, 26, "grG1ELowerSys_RatioToDefault");

  //-- Draw
  //--Distns
  TCanvas * cDistn = new TCanvas("cDistn", "cDistn", 500, 500);
  cDistn->cd();
  cDistn->SetLogy();
  hUnfoldDefault[1]->Draw();
  hUnfoldUpperSys[1]->Draw("same");
  hUnfoldLowerSys[1]->Draw("same");
  latex.DrawLatex(0.6, 0.8, "Cent 5-10%");
  cDistn->SaveAs("cDistn.pdf");

  //-- cumuRatio
  TLegend * legvn6vn4 = new TLegend(0.6538, 0.7421, 0.9529, 0.9172);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grVn6Vn4Default,  "Default",  "lp");
  legvn6vn4->AddEntry(grVn6Vn4UpperSys, "UpperSys", "lp");
  legvn6vn4->AddEntry(grVn6Vn4LowerSys, "LowerSys", "lp");

  TLegend * legvn8vn4 = new TLegend(0.6538, 0.7421, 0.9529, 0.9172);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grVn8Vn4Default,  "Default",  "lp");
  legvn8vn4->AddEntry(grVn8Vn4UpperSys, "UpperSys", "lp");
  legvn8vn4->AddEntry(grVn8Vn4LowerSys, "LowerSys", "lp");

  TLegend * legvn8vn6 = new TLegend(0.6520, 0.1973, 0.9520, 0.3703);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grVn8Vn6Default,  "Default",  "lp");
  legvn8vn6->AddEntry(grVn8Vn6UpperSys, "UpperSys", "lp");
  legvn8vn6->AddEntry(grVn8Vn6LowerSys, "LowerSys", "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  cCumuRatio->cd(1);
  grVn6Vn4Default->Draw("ap");
  grVn6Vn4UpperSys->Draw("psame");
  grVn6Vn4LowerSys->Draw("psame");
  legvn6vn4->Draw("same");
  cCumuRatio->cd(2);
  grVn8Vn4Default->Draw("ap");
  grVn8Vn4UpperSys->Draw("psame");
  grVn8Vn4LowerSys->Draw("psame");
  legvn8vn4->Draw("same");
  cCumuRatio->cd(3);
  grVn8Vn6Default->Draw("ap");
  grVn8Vn6UpperSys->Draw("psame");
  grVn8Vn6LowerSys->Draw("psame");
  legvn8vn6->Draw("same");
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- cumuRatio To Default
  TLegend * legvn6vn4R = new TLegend(0.4837, 0.8070, 0.9529, 0.9172);
  legvn6vn4R->SetBorderSize(0);
  legvn6vn4R->SetFillStyle(0);
  legvn6vn4R->AddEntry(grVn6Vn4UpperSys_RatioToDefault, "UpperSys / Default", "lp");
  legvn6vn4R->AddEntry(grVn6Vn4LowerSys_RatioToDefault, "LowerSys / Default", "lp");

  TLegend * legvn8vn4R = new TLegend(0.4837, 0.8070, 0.9529, 0.9172);
  legvn8vn4R->SetBorderSize(0);
  legvn8vn4R->SetFillStyle(0);
  legvn8vn4R->AddEntry(grVn8Vn4UpperSys_RatioToDefault, "UpperSys / Default", "lp");
  legvn8vn4R->AddEntry(grVn8Vn4LowerSys_RatioToDefault, "LowerSys / Default", "lp");

  TLegend * legvn8vn6R = new TLegend(0.4837, 0.8070, 0.9529, 0.9172);
  legvn8vn6R->SetBorderSize(0);
  legvn8vn6R->SetFillStyle(0);
  legvn8vn6R->AddEntry(grVn8Vn6UpperSys_RatioToDefault, "UpperSys / Default", "lp");
  legvn8vn6R->AddEntry(grVn8Vn6LowerSys_RatioToDefault, "LowerSys / Default", "lp");

  TCanvas * cCumuRatio_RatToDef = new TCanvas("cCumuRatio_RatToDef", "cCumuRatio_RatToDef", 1500, 500);
  cCumuRatio_RatToDef->Divide(3,1);
  cCumuRatio_RatToDef->cd(1);
  grVn6Vn4UpperSys_RatioToDefault->Draw("ap");
  grVn6Vn4LowerSys_RatioToDefault->Draw("psame");
  legvn6vn4R->Draw("same");
  cCumuRatio_RatToDef->cd(2);
  grVn8Vn4UpperSys_RatioToDefault->Draw("ap");
  grVn8Vn4LowerSys_RatioToDefault->Draw("psame");
  legvn8vn4R->Draw("same");
  cCumuRatio_RatToDef->cd(3);
  grVn8Vn6UpperSys_RatioToDefault->Draw("ap");
  grVn8Vn6LowerSys_RatioToDefault->Draw("psame");
  legvn8vn6R->Draw("same");
  cCumuRatio_RatToDef->SaveAs("cCumuRatio_RatToDef.pdf");

  //-- g1e
  TLegend * legG1E = new TLegend(0.6538, 0.7421, 0.9529, 0.9172);
  legG1E->SetBorderSize(0);
  legG1E->SetFillStyle(0);
  legG1E->AddEntry(grG1EDefault,  "Default", "lp");
  legG1E->AddEntry(grG1EUpperSys, "UpperSys", "lp");
  legG1E->AddEntry(grG1ELowerSys, "LowerSys", "lp");

  TLegend * legG1ER = new TLegend(0.6520, 0.1973, 0.9520, 0.3703);
  legG1ER->SetBorderSize(0);
  legG1ER->SetFillStyle(0);
  legG1ER->AddEntry(grG1EUpperSys_RatioToDefault, "UpperSys / Default", "lp");
  legG1ER->AddEntry(grG1ELowerSys_RatioToDefault, "LowerSys / Default", "lp");

  TCanvas* cG1e = new TCanvas("cG1e", "cG1e", 1000, 500);
  cG1e->Divide(2,1);
  cG1e->cd(1);
  grG1EDefault->Draw("ap");
  grG1EUpperSys->Draw("psame");
  grG1ELowerSys->Draw("psame");
  legG1E->Draw("same");
  cG1e->cd(2);
  grG1EUpperSys_RatioToDefault->Draw("ap");
  grG1ELowerSys_RatioToDefault->Draw("psame");
  legG1ER->Draw("same");
  cG1e->SaveAs("cG1e.pdf");




}

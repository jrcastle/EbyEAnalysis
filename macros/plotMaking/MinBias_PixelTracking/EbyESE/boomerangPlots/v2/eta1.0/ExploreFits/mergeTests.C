#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

/*

    Assume a linear flow response: vn = kn * eccn
      1. (kn, alpha, ecc0) free
      2. (alpha, ecc0) free; kn fixed at ATLAS values

    Assume 2015 cubic flow response (Cubic1): vn = kn*eccn + kn*knPr*pow(eccn, 3)
      1. (kn, knPr, alpha, ecc0) free 
      2. (kn, alpha, ecc0) free; knPr fixed at 0.1
      3. (alpha, ecc0) free; kn fixed at ATLAS values, knPr fixed at 0.1

    Assume 2016 cubic flow response (Cubic1): vn = kn*eccn + knPr*pow(eccn, 3)  
      1. (kn, knPr, alpha, ecc0) free
      2. (kn, alpha, ecc0) free; knPr fixed at PRC hydro
      3. (knPr, alpha, ecc0) free; kn fixed at PRC hydro  
      4. (alpha, ecc0) free; kn and knPr fixed at PRC hydro

*/

void mergeTests(){

  double knMin    = 0.;
  double knMax    = 0.75;
  double knPrMin  = -0.01;
  double knPrMax  = 0.15;
  double alphaMin = 0.;
  double alphaMax = 130.;
  double e0Min    = 0.;
  double e0Max    = 0.5;



  //-- Linear, test 1
  TFile * fLT1;
  TGraphErrors * grFitKn_LT1;
  TGraphErrors * grFitAlpha_LT1;
  TGraphErrors * grFitE0_LT1;

  //-- Linear, test 2
  TFile * fLT2;
  TGraphErrors * grFitKn_LT2;
  TGraphErrors * grFitAlpha_LT2;
  TGraphErrors * grFitE0_LT2;

  //-- Cubic1, test 1 
  TFile * fC1T1;
  TGraphErrors * grFitKn_C1T1;
  TGraphErrors * grFitKnPr_C1T1;
  TGraphErrors * grFitAlpha_C1T1;
  TGraphErrors * grFitE0_C1T1;

  //-- Cubic1, test 2
  TFile * fC1T2;
  TGraphErrors * grFitKn_C1T2;
  TGraphErrors * grFitKnPr_C1T2;
  TGraphErrors * grFitAlpha_C1T2;
  TGraphErrors * grFitE0_C1T2;

  //-- Cubic1, test 3
  TFile * fC1T3;
  TGraphErrors * grFitKn_C1T3;
  TGraphErrors * grFitKnPr_C1T3;
  TGraphErrors * grFitAlpha_C1T3;
  TGraphErrors * grFitE0_C1T3;

  //-- Cubic2, test 1
  TFile * fC2T1;
  TGraphErrors * grFitKn_C2T1;
  TGraphErrors * grFitKnPr_C2T1;
  TGraphErrors * grFitAlpha_C2T1;
  TGraphErrors * grFitE0_C2T1;

  //-- Cubic2, test 2
  TFile * fC2T2;
  TGraphErrors * grFitKn_C2T2;
  TGraphErrors * grFitKnPr_C2T2;
  TGraphErrors * grFitAlpha_C2T2;
  TGraphErrors * grFitE0_C2T2;

  //-- Cubic2, test 3
  TFile * fC2T3;
  TGraphErrors * grFitKn_C2T3;
  TGraphErrors * grFitKnPr_C2T3;
  TGraphErrors * grFitAlpha_C2T3;
  TGraphErrors * grFitE0_C2T3;

  //-- Cubic2, test 4
  TFile * fC2T4;
  TGraphErrors * grFitKn_C2T4;
  TGraphErrors * grFitKnPr_C2T4;
  TGraphErrors * grFitAlpha_C2T4;
  TGraphErrors * grFitE0_C2T4;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Grab data and format
  //-- Linear, test 1
  fLT1 = new TFile("EllPFits_Linear_Test1.root");
  grFitKn_LT1    = (TGraphErrors*) fLT1->Get("grFitKn");
  grFitAlpha_LT1 = (TGraphErrors*) fLT1->Get("grFitAlpha");
  grFitE0_LT1    = (TGraphErrors*) fLT1->Get("grFitE0");

  formatGraph(grFitKn_LT1,    "Centrality %", knMin,    knMax,    "k_{n}",        1,        20, "grFitKn_LT1");
  formatGraph(grFitAlpha_LT1, "Centrality %", alphaMin, alphaMax, "#alpha",       1,        20, "grFitAlpha_LT1");
  formatGraph(grFitE0_LT1,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", 1,        20, "grFitE0_LT1");

  //-- Linear, test 2
  fLT2 = new TFile("EllPFits_Linear_Test2.root");
  grFitKn_LT2    = (TGraphErrors*) fLT2->Get("grFitKn");
  grFitAlpha_LT2 = (TGraphErrors*) fLT2->Get("grFitAlpha");
  grFitE0_LT2    = (TGraphErrors*) fLT2->Get("grFitE0");

  formatGraph(grFitKn_LT2,    "Centrality %", knMin,    knMax,    "k_{n}",        34,        21, "grFitKn_LT2");
  formatGraph(grFitAlpha_LT2, "Centrality %", alphaMin, alphaMax, "#alpha",       34,        21, "grFitAlpha_LT2");
  formatGraph(grFitE0_LT2,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", 34,        21, "grFitE0_LT2");

  //-- Cubic1, test 1
  fC1T1 = new TFile("EllPFits_Cubic1_Test1.root");
  grFitKn_C1T1    = (TGraphErrors*) fC1T1->Get("grFitKn");
  grFitKnPr_C1T1  = (TGraphErrors*) fC1T1->Get("grFitKnPr");
  grFitAlpha_C1T1 = (TGraphErrors*) fC1T1->Get("grFitAlpha");
  grFitE0_C1T1    = (TGraphErrors*) fC1T1->Get("grFitE0");
 
  formatGraph(grFitKn_C1T1,    "Centrality %", knMin,    knMax,    "k_{n}",        2,        20, "grFitKn_C1T1");
  formatGraph(grFitKnPr_C1T1,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     2,        20, "grFitKnPr_C1T1");
  formatGraph(grFitAlpha_C1T1, "Centrality %", alphaMin, alphaMax, "#alpha",       2,        20, "grFitAlpha_C1T1");
  formatGraph(grFitE0_C1T1,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", 2,        20, "grFitE0_C1T1");

  //-- Cubic1, test 2
  fC1T2 = new TFile("EllPFits_Cubic1_Test2.root");
  grFitKn_C1T2    = (TGraphErrors*) fC1T2->Get("grFitKn");
  grFitKnPr_C1T2  = (TGraphErrors*) fC1T2->Get("grFitKnPr");
  grFitAlpha_C1T2 = (TGraphErrors*) fC1T2->Get("grFitAlpha");
  grFitE0_C1T2    = (TGraphErrors*) fC1T2->Get("grFitE0");

  formatGraph(grFitKn_C1T2,    "Centrality %", knMin,    knMax,    "k_{n}",        kRed+1,        22, "grFitKn_C1T2");
  formatGraph(grFitKnPr_C1T2,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     kRed+1,        22, "grFitKnPr_C1T2");
  formatGraph(grFitAlpha_C1T2, "Centrality %", alphaMin, alphaMax, "#alpha",       kRed+1,        22, "grFitAlpha_C1T2");
  formatGraph(grFitE0_C1T2,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", kRed+1,        22, "grFitE0_C1T2");

  //-- Cubic1, test 3
  fC1T3 = new TFile("EllPFits_Cubic1_Test3.root");
  grFitKn_C1T3    = (TGraphErrors*) fC1T3->Get("grFitKn");
  grFitKnPr_C1T3  = (TGraphErrors*) fC1T3->Get("grFitKnPr");
  grFitAlpha_C1T3 = (TGraphErrors*) fC1T3->Get("grFitAlpha");
  grFitE0_C1T3    = (TGraphErrors*) fC1T3->Get("grFitE0");

  formatGraph(grFitKn_C1T3,    "Centrality %", knMin,    knMax,    "k_{n}",        kRed+2,        21, "grFitKn_C1T3");
  formatGraph(grFitKnPr_C1T3,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     kRed+2,        21, "grFitKnPr_C1T3");
  formatGraph(grFitAlpha_C1T3, "Centrality %", alphaMin, alphaMax, "#alpha",       kRed+2,        21, "grFitAlpha_C1T3");
  formatGraph(grFitE0_C1T3,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", kRed+2,        21, "grFitE0_C1T3");

  //-- Cubic2, test 1
  fC2T1 = new TFile("EllPFits_Cubic2_Test1.root");
  grFitKn_C2T1    = (TGraphErrors*) fC2T1->Get("grFitKn");
  grFitKnPr_C2T1  = (TGraphErrors*) fC2T1->Get("grFitKnPr");
  grFitAlpha_C2T1 = (TGraphErrors*) fC2T1->Get("grFitAlpha");
  grFitE0_C2T1    = (TGraphErrors*) fC2T1->Get("grFitE0");

  formatGraph(grFitKn_C2T1,    "Centrality %", knMin,    knMax,    "k_{n}",        4,        20, "grFitKn_C2T1");
  formatGraph(grFitKnPr_C2T1,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     4,        20, "grFitKnPr_C2T1");
  formatGraph(grFitAlpha_C2T1, "Centrality %", alphaMin, alphaMax, "#alpha",       4,        20, "grFitAlpha_C2T1");
  formatGraph(grFitE0_C2T1,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", 4,        20, "grFitE0_C2T1");

  //-- Cubic2, test 2
  fC2T2 = new TFile("EllPFits_Cubic2_Test2.root");
  grFitKn_C2T2    = (TGraphErrors*) fC2T2->Get("grFitKn");
  grFitKnPr_C2T2  = (TGraphErrors*) fC2T2->Get("grFitKnPr");
  grFitAlpha_C2T2 = (TGraphErrors*) fC2T2->Get("grFitAlpha");
  grFitE0_C2T2    = (TGraphErrors*) fC2T2->Get("grFitE0");

  formatGraph(grFitKn_C2T2,    "Centrality %", knMin,    knMax,    "k_{n}",        38,        23, "grFitKn_C2T2");
  formatGraph(grFitKnPr_C2T2,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     38,        23, "grFitKnPr_C2T2");
  formatGraph(grFitAlpha_C2T2, "Centrality %", alphaMin, alphaMax, "#alpha",       38,        23, "grFitAlpha_C2T2");
  formatGraph(grFitE0_C2T2,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", 38,        23, "grFitE0_C2T2");

  //-- Cubic2, test 3
  fC2T3 = new TFile("EllPFits_Cubic2_Test3.root");
  grFitKn_C2T3    = (TGraphErrors*) fC2T3->Get("grFitKn");
  grFitKnPr_C2T3  = (TGraphErrors*) fC2T3->Get("grFitKnPr");
  grFitAlpha_C2T3 = (TGraphErrors*) fC2T3->Get("grFitAlpha");
  grFitE0_C2T3    = (TGraphErrors*) fC2T3->Get("grFitE0");

  formatGraph(grFitKn_C2T3,    "Centrality %", knMin,    knMax,    "k_{n}",        kBlue+2,        22, "grFitKn_C2T3");
  formatGraph(grFitKnPr_C2T3,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     kBlue+2,        22, "grFitKnPr_C2T3");
  formatGraph(grFitAlpha_C2T3, "Centrality %", alphaMin, alphaMax, "#alpha",       kBlue+2,        22, "grFitAlpha_C2T3");
  formatGraph(grFitE0_C2T3,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", kBlue+2,        22, "grFitE0_C2T3");

  //-- Cubic2, test 4
  fC2T4 = new TFile("EllPFits_Cubic2_Test4.root");
  grFitKn_C2T4    = (TGraphErrors*) fC2T4->Get("grFitKn");
  grFitKnPr_C2T4  = (TGraphErrors*) fC2T4->Get("grFitKnPr");
  grFitAlpha_C2T4 = (TGraphErrors*) fC2T4->Get("grFitAlpha");
  grFitE0_C2T4    = (TGraphErrors*) fC2T4->Get("grFitE0");

  formatGraph(grFitKn_C2T4,    "Centrality %", knMin,    knMax,    "k_{n}",        kCyan+2,        21, "grFitKn_C2T4");
  formatGraph(grFitKnPr_C2T4,  "Centrality %", knPrMin,  knPrMax,  "#kappa\'",     kCyan+2,        21, "grFitKnPr_C2T4");
  formatGraph(grFitAlpha_C2T4, "Centrality %", alphaMin, alphaMax, "#alpha",       kCyan+2,        21, "grFitAlpha_C2T4");
  formatGraph(grFitE0_C2T4,    "Centrality %", e0Min,    e0Max,    "#epsilon_{0}", kCyan+2,        21, "grFitE0_C2T4");

  //-- Reset Marker sizes
  grFitE0_C1T1->SetMarkerSize(1.);
  grFitE0_C1T2->SetMarkerSize(1.);
  grFitE0_C1T3->SetMarkerSize(1.);
  grFitE0_C2T1->SetMarkerSize(1.);
  grFitE0_C2T2->SetMarkerSize(1.);
  grFitE0_C2T3->SetMarkerSize(1.);
  grFitE0_C2T4->SetMarkerSize(1.);

  //-- Plot Number 1, show free parm fits from different response scenarios
  TLegend * leg1 = new TLegend(0.20, 0.72, 0.60, 0.92);
  legInit( leg1 );
  leg1->AddEntry(grFitKn_LT1,  "Linear",     "lp");
  leg1->AddEntry(grFitKn_C1T1, "2015 Cubic", "lp");
  leg1->AddEntry(grFitKn_C2T1, "2016 Cubic", "lp");

  TCanvas * c1 = new TCanvas("c1", "c1", 2000, 500);
  c1->Divide(4,1);

  c1->cd(1);
  grFitKn_LT1->Draw("ap");
  grFitKn_C1T1->Draw("psame");
  grFitKn_C2T1->Draw("psame");
  leg1->Draw("same");

  c1->cd(2);
  grFitKnPr_C1T1->Draw("ap");
  grFitKnPr_C2T1->Draw("psame");

  c1->cd(3);
  grFitE0_LT1->Draw("ap");
  grFitE0_C1T1->Draw("psame");
  grFitE0_C2T1->Draw("psame");

  c1->cd(4);
  grFitAlpha_LT1->Draw("ap");
  grFitAlpha_C1T1->Draw("psame");
  grFitAlpha_C2T1->Draw("psame");

  leg1->SetTextFont(43);
  leg1->SetTextSize(23);

  c1->Update();
  c1->SaveAs( "FreeParms.pdf");

  //-- Plot Number 2, show linear test results
  TLegend * leg2 = new TLegend(0.20, 0.72, 0.60, 0.92);
  legInit( leg2 );
  leg2->AddEntry(grFitKn_LT1, "All Free",      "lp");
  leg2->AddEntry(grFitKn_LT2, "k_{n} = ATLAS", "lp");

  TGraphErrors * g = (TGraphErrors*) grFitKnPr_C1T1->Clone("g");
  for(int i = 0; i < NCENT; i++) g->SetPoint(i, centBinCenter[i], -1);
  formatGraph(g, "Centrality %", knPrMin,  knPrMax, "#kappa\'", 4, 21, "g");

  TCanvas * c2 = new TCanvas("c2", "c2", 2000, 500);
  c2->Divide(4,1);

  c2->cd(1);
  grFitKn_LT1->Draw("ap");
  grFitKn_LT2->Draw("psame");
  leg2->Draw("same");

  c2->cd(2);
  g->Draw("ap");

  c2->cd(3);
  grFitE0_LT1->Draw("ap");
  grFitE0_LT2->Draw("psame");

  c2->cd(4);
  grFitAlpha_LT1->Draw("ap");
  grFitAlpha_LT2->Draw("psame");

  leg2->SetTextFont(43);
  leg2->SetTextSize(23);

  c2->Update();
  c2->SaveAs( "LinearTests.pdf");

  //-- Plot Number 3, show 2015 cubic test results
  TLegend * leg3 = new TLegend(0.20, 0.72, 0.60, 0.92);
  legInit( leg3 );
  leg3->AddEntry(grFitKn_C1T1, "All Free",                      "lp");
  leg3->AddEntry(grFitKn_C1T2, "#kappa\' = 0.1",                "lp");
  leg3->AddEntry(grFitKn_C1T3, "k_{n} = ATLAS; #kappa\' = 0.1", "lp");

  TCanvas * c3 = new TCanvas("c3", "c3", 2000, 500);
  c3->Divide(4,1);

  c3->cd(1);
  grFitKn_C1T1->Draw("ap");
  grFitKn_C1T2->Draw("psame");
  grFitKn_C1T3->Draw("psame");
  leg3->Draw("same");

  c3->cd(2);
  grFitKnPr_C1T1->Draw("ap");
  grFitKnPr_C1T2->Draw("ap");
  grFitKnPr_C1T3->Draw("ap");

  c3->cd(3);
  grFitE0_C1T1->Draw("ap");
  grFitE0_C1T2->Draw("psame");
  grFitE0_C1T3->Draw("psame");

  c3->cd(4);
  grFitAlpha_C1T1->Draw("ap");
  grFitAlpha_C1T2->Draw("psame");
  grFitAlpha_C1T3->Draw("psame");

  leg3->SetTextFont(43);
  leg3->SetTextSize(23);

  c3->Update();
  c3->SaveAs( "Cubic1Tests.pdf");

  //-- Plot Number 4, show 2016 cubic test results
  TLegend * leg4 = new TLegend(0.20, 0.72, 0.60, 0.92);
  legInit( leg4 );
  leg4->AddEntry(grFitKn_C2T1, "All Free",                        "lp");
  leg4->AddEntry(grFitKn_C2T2, "#kappa\' = Hydro",                "lp");
  leg4->AddEntry(grFitKn_C2T3, "k_{n} = Hydro",                   "lp");
  leg4->AddEntry(grFitKn_C2T4, "k_{n} = Hydro; #kappa\' = Hydro", "lp");

  TCanvas * c4 = new TCanvas("c4", "c4", 2000, 500);
  c4->Divide(4,1);

  c4->cd(1);
  grFitKn_C2T1->Draw("ap");
  grFitKn_C2T2->Draw("psame");
  grFitKn_C2T3->Draw("psame");
  grFitKn_C2T4->Draw("psame");
  leg4->Draw("same");

  c4->cd(2);
  grFitKnPr_C2T1->Draw("ap");
  grFitKnPr_C2T2->Draw("psame");
  grFitKnPr_C2T3->Draw("psame");
  grFitKnPr_C2T4->Draw("psame");

  c4->cd(3);
  grFitE0_C2T1->Draw("ap");
  grFitE0_C2T2->Draw("psame");
  grFitE0_C2T3->Draw("psame");
  grFitE0_C2T4->Draw("psame");

  c4->cd(4);
  grFitAlpha_C2T1->Draw("ap");
  grFitAlpha_C2T2->Draw("psame");
  grFitAlpha_C2T3->Draw("psame");
  grFitAlpha_C2T4->Draw("psame");

  leg4->SetTextFont(43);
  leg4->SetTextSize(23);

  c4->Update();
  c4->SaveAs( "Cubic2Tests.pdf");



}

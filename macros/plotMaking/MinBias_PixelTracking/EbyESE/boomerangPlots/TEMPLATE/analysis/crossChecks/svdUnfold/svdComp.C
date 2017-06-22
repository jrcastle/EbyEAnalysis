#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void svdComp(int n = 2){

  const int norder_ = n;
  double ratioMin = 0.95;
  double ratioMax = 1.05;
  double ratioMinVn8Vn6 = 0.98;
  double ratioMaxVn8Vn6 = 1.003;
  double ratioMinG1E = 0.5;
  double ratioMaxG1E = 3.2;

  double histMin = -0.2;
  double histMax = 1.2;

  //-- D'Agostini
  TFile * fDAgUnf;
  TH1D * hUnfoldStatDAg[NCENT];
  TH1D * hUnfoldSysDAg[NCENT];
  TH1D * hUnfoldStatAndSysDAg[NCENT];

  TFile * fDAgPhysics;
  TGraphErrors * grVn2RawDAg;
  TGraphErrors * grVn2RawDAgSys;
  TGraphErrors * grVn4RawDAg;
  TGraphErrors * grVn4RawDAgSys;
  TGraphErrors * grVn6RawDAg;
  TGraphErrors * grVn6RawDAgSys;
  TGraphErrors * grVn8RawDAg;
  TGraphErrors * grVn8RawDAgSys;
  TGraphErrors * grvn6vn4RatioDAg;
  TGraphErrors * grvn6vn4RatioDAgSys;
  TGraphErrors * grvn8vn4RatioDAg;
  TGraphErrors * grvn8vn4RatioDAgSys;
  TGraphErrors * grvn8vn6RatioDAg;
  TGraphErrors * grvn8vn6RatioDAgSys;
  TGraphErrors * grGamma1ExpDAg;
  TGraphErrors * grGamma1ExpDAgSys;

  //-- SVD
  TFile * fSVDPhysics;
  TH1D * hUnfoldSVD[NCENT];
  TGraphErrors * grVn2RawSVD;
  TGraphErrors * grVn4RawSVD;
  TGraphErrors * grVn6RawSVD;
  TGraphErrors * grVn8RawSVD;
  TGraphErrors * grvn6vn4RatioSVD;
  TGraphErrors * grvn8vn4RatioSVD;
  TGraphErrors * grvn8vn6RatioSVD;
  TGraphErrors * grGamma1ExpSVD;

  //-- Unfold ratios
  TCanvas * cUnfold;
  TCanvas * cUnfoldRatio_SVD_DAg;
  TH1D * hUnfoldRatio_SVD_DAg[NCENT];

  //-- Physics ratios
  double vn2SVD_DAg[NCENT];
  double vn4SVD_DAg[NCENT];
  double vn6SVD_DAg[NCENT];
  double vn8SVD_DAg[NCENT];
  double vn6vn4SVD_DAg[NCENT];
  double vn8vn4SVD_DAg[NCENT];
  double vn8vn6SVD_DAg[NCENT];
  double g1eSVD_DAg[NCENT];

  double vn2SVD_DAg_err[NCENT];
  double vn4SVD_DAg_err[NCENT];
  double vn6SVD_DAg_err[NCENT];
  double vn8SVD_DAg_err[NCENT];
  double vn6vn4SVD_DAg_err[NCENT];
  double vn8vn4SVD_DAg_err[NCENT];
  double vn8vn6SVD_DAg_err[NCENT];
  double g1eSVD_DAg_err[NCENT];

  TGraphErrors * grVn2RawRatio_SVD_DAg;
  TGraphErrors * grVn4RawRatio_SVD_DAg;
  TGraphErrors * grVn6RawRatio_SVD_DAg;
  TGraphErrors * grVn8RawRatio_SVD_DAg;
  TGraphErrors * grvn6vn4Ratio_SVD_DAg;
  TGraphErrors * grvn8vn4Ratio_SVD_DAg;
  TGraphErrors * grvn8vn6Ratio_SVD_DAg;
  TGraphErrors * grGamma1ExpRatio_SVD_DAg;

  TLatex latex;
  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  //-- Get D'Agostini objects
  fDAgUnf = new TFile( Form("../../systematicStudies/SysUnfoldDistns_v%i.root", norder_) );
  for(int icent = 0; icent < NCENT; icent++){
    hUnfoldStatDAg[icent]       = (TH1D*) fDAgUnf->Get( Form("hFinalUnfoldStat_c%i", icent) );
    hUnfoldStatDAg[icent]->SetLineColor(1);
    hUnfoldStatDAg[icent]->SetMarkerColor(1);
    hUnfoldSysDAg[icent]        = (TH1D*) fDAgUnf->Get( Form("hFinalUnfoldSys_c%i", icent) );
    hUnfoldStatAndSysDAg[icent] = (TH1D*) fDAgUnf->Get( Form("hFinalUnfoldStatAndSys_c%i", icent) );
  }

  fDAgPhysics = new TFile( "../../systematicStudies/PhysicsResults.root" );
  grVn2RawDAg         = (TGraphErrors*) fDAgPhysics->Get( "grVn2Raw" );
  grVn2RawDAgSys      = (TGraphErrors*) fDAgPhysics->Get( "grVn2RawSys" );
  grVn4RawDAg         = (TGraphErrors*) fDAgPhysics->Get( "grVn4Raw" );
  grVn4RawDAgSys      = (TGraphErrors*) fDAgPhysics->Get( "grVn4RawSys" );
  grVn6RawDAg         = (TGraphErrors*) fDAgPhysics->Get( "grVn6Raw" );
  grVn6RawDAgSys      = (TGraphErrors*) fDAgPhysics->Get( "grVn6RawSys" );
  grVn8RawDAg         = (TGraphErrors*) fDAgPhysics->Get( "grVn8Raw" );
  grVn8RawDAgSys      = (TGraphErrors*) fDAgPhysics->Get( "grVn8RawSys" );
  grvn6vn4RatioDAg    = (TGraphErrors*) fDAgPhysics->Get( "grvn6vn4Ratio" );
  grvn6vn4RatioDAgSys = (TGraphErrors*) fDAgPhysics->Get( "grvn6vn4RatioSys" );
  grvn8vn4RatioDAg    = (TGraphErrors*) fDAgPhysics->Get( "grvn8vn4Ratio" );
  grvn8vn4RatioDAgSys = (TGraphErrors*) fDAgPhysics->Get( "grvn8vn4RatioSys" );
  grvn8vn6RatioDAg    = (TGraphErrors*) fDAgPhysics->Get( "grvn8vn6Ratio" );
  grvn8vn6RatioDAgSys = (TGraphErrors*) fDAgPhysics->Get( "grvn8vn6RatioSys" );
  grGamma1ExpDAg      = (TGraphErrors*) fDAgPhysics->Get( "grGamma1Exp" );
  grGamma1ExpDAgSys   = (TGraphErrors*) fDAgPhysics->Get( "grGamma1ExpSys" );

  grVn2RawDAgSys->GetYaxis()->SetTitle( Form("v_{%i}{2}", norder_) );
  grVn4RawDAgSys->GetYaxis()->SetTitle( Form("v_{%i}{4}", norder_) );
  grVn6RawDAgSys->GetYaxis()->SetTitle( Form("v_{%i}{6}", norder_) );
  grVn8RawDAgSys->GetYaxis()->SetTitle( Form("v_{%i}{8}", norder_) );

  //-- Get SVD objects
  fSVDPhysics = new TFile( "../../systematicStudies/SVDPhysics.root" );

  for(int icent = 0; icent < NCENT; icent++){
    hUnfoldSVD[icent] = (TH1D*) fSVDPhysics->Get( Form("hrecokreg4_c%i", icent) );
    hUnfoldSVD[icent]->SetLineColor(2);
    hUnfoldSVD[icent]->SetMarkerColor(2);
    hUnfoldSVD[icent]->Scale( 1./hUnfoldSVD[icent]->Integral() );

    //-- Ratio to D'Agostini
    hUnfoldRatio_SVD_DAg[icent] = (TH1D*) hUnfoldSVD[icent]->Clone( Form("hUnfoldRatio_SVD_DAg_c%i", icent) );
    hUnfoldRatio_SVD_DAg[icent]->Divide( hUnfoldStatAndSysDAg[icent] );
    hUnfoldRatio_SVD_DAg[icent]->SetMinimum(histMin);
    hUnfoldRatio_SVD_DAg[icent]->SetMaximum(histMax);

  }

  grVn2RawSVD      = (TGraphErrors*) fSVDPhysics->Get( "grVn2" );
  grVn4RawSVD      = (TGraphErrors*) fSVDPhysics->Get( "grVn4" );
  grVn6RawSVD      = (TGraphErrors*) fSVDPhysics->Get( "grVn6" );
  grVn8RawSVD      = (TGraphErrors*) fSVDPhysics->Get( "grVn8" );
  grvn6vn4RatioSVD = (TGraphErrors*) fSVDPhysics->Get( "grVn6Vn4" );
  grvn8vn4RatioSVD = (TGraphErrors*) fSVDPhysics->Get( "grVn8Vn4" );
  grvn8vn6RatioSVD = (TGraphErrors*) fSVDPhysics->Get( "grVn8Vn6" );
  grGamma1ExpSVD   = (TGraphErrors*) fSVDPhysics->Get( "grG1E" );

  //-- Fill Ratio Arrays
  for(int icent = 0; icent < NCENT; icent++){

    //-- D'Agostini
    double vn2D    = grVn2RawDAg->GetY()[icent];
    double vn4D    = grVn4RawDAg->GetY()[icent];
    double vn6D    = grVn6RawDAg->GetY()[icent];
    double vn8D    = grVn8RawDAg->GetY()[icent];
    double vn6vn4D = grvn6vn4RatioDAg->GetY()[icent];
    double vn8vn4D = grvn8vn4RatioDAg->GetY()[icent];
    double vn8vn6D = grvn8vn6RatioDAg->GetY()[icent];
    double g1eD    = grGamma1ExpDAg->GetY()[icent];

    double vn2D_Ste    = grVn2RawDAg->GetErrorY(icent);
    double vn4D_Ste    = grVn4RawDAg->GetErrorY(icent);
    double vn6D_Ste    = grVn6RawDAg->GetErrorY(icent);
    double vn8D_Ste    = grVn8RawDAg->GetErrorY(icent);
    double vn6vn4D_Ste = grvn6vn4RatioDAg->GetErrorY(icent);
    double vn8vn4D_Ste = grvn8vn4RatioDAg->GetErrorY(icent);
    double vn8vn6D_Ste = grvn8vn6RatioDAg->GetErrorY(icent);
    double g1eD_Ste    = grGamma1ExpDAg->GetErrorY(icent);

    double vn2D_Sye    = grVn2RawDAgSys->GetErrorY(icent);
    double vn4D_Sye    = grVn4RawDAgSys->GetErrorY(icent);
    double vn6D_Sye    = grVn6RawDAgSys->GetErrorY(icent);
    double vn8D_Sye    = grVn8RawDAgSys->GetErrorY(icent);
    double vn6vn4D_Sye = grvn6vn4RatioDAgSys->GetErrorY(icent);
    double vn8vn4D_Sye = grvn8vn4RatioDAgSys->GetErrorY(icent);
    double vn8vn6D_Sye = grvn8vn6RatioDAgSys->GetErrorY(icent);
    double g1eD_Sye    = grGamma1ExpDAgSys->GetErrorY(icent);

    double vn2D_Tote    = sqrt( pow(vn2D_Ste, 2) + pow(vn2D_Sye, 2) );
    double vn4D_Tote    = sqrt( pow(vn4D_Ste, 2) + pow(vn4D_Sye, 2) );
    double vn6D_Tote    = sqrt( pow(vn6D_Ste, 2) + pow(vn6D_Sye, 2) );
    double vn8D_Tote    = sqrt( pow(vn8D_Ste, 2) + pow(vn8D_Sye, 2) );
    double vn6vn4D_Tote = sqrt( pow(vn6vn4D_Ste, 2) + pow(vn6vn4D_Sye, 2) );
    double vn8vn4D_Tote = sqrt( pow(vn8vn4D_Ste, 2) + pow(vn8vn4D_Sye, 2) );
    double vn8vn6D_Tote = sqrt( pow(vn8vn6D_Ste, 2) + pow(vn8vn6D_Sye, 2) );
    double g1eD_Tote    = sqrt( pow(g1eD_Ste, 2) + pow(g1eD_Sye, 2) );

    //-- SVD
    double vn2S    = grVn2RawSVD->GetY()[icent];
    double vn4S    = grVn4RawSVD->GetY()[icent];
    double vn6S    = grVn6RawSVD->GetY()[icent];
    double vn8S    = grVn8RawSVD->GetY()[icent];
    double vn6vn4S = grvn6vn4RatioSVD->GetY()[icent];
    double vn8vn4S = grvn8vn4RatioSVD->GetY()[icent];
    double vn8vn6S = grvn8vn6RatioSVD->GetY()[icent];
    double g1eS    = grGamma1ExpSVD->GetY()[icent];

    double vn2S_Ste    = grVn2RawSVD->GetErrorY(icent);
    double vn4S_Ste    = grVn4RawSVD->GetErrorY(icent);
    double vn6S_Ste    = grVn6RawSVD->GetErrorY(icent);
    double vn8S_Ste    = grVn8RawSVD->GetErrorY(icent);
    double vn6vn4S_Ste = grvn6vn4RatioSVD->GetErrorY(icent);
    double vn8vn4S_Ste = grvn8vn4RatioSVD->GetErrorY(icent);
    double vn8vn6S_Ste = grvn8vn6RatioSVD->GetErrorY(icent);
    double g1eS_Ste    = grGamma1ExpSVD->GetErrorY(icent);

    //-- Ratios
    double vn2SD;
    double vn4SD;
    double vn6SD;
    double vn8SD;
    double vn6vn4SD;
    double vn8vn4SD;
    double vn8vn6SD;
    double g1eSD;

    if(vn2D <= 0 || vn2S <= 0) vn2SD = -1;
    else                       vn2SD = vn2S / vn2D;
    if(vn4D <= 0 || vn4S <= 0) vn4SD = -1;
    else                       vn4SD = vn4S / vn4D;
    if(vn6D <= 0 || vn6S <= 0) vn6SD = -1;
    else                       vn6SD = vn6S / vn6D;
    if(vn8D <= 0 || vn8S <= 0) vn8SD = -1;
    else                       vn8SD = vn8S / vn8D;

    if(vn6vn4D <= 0 || vn6vn4S <= 0) vn6vn4SD = -1;
    else                             vn6vn4SD = vn6vn4S / vn6vn4D;
    if(vn8vn4D <= 0 || vn8vn4S <= 0) vn8vn4SD =-1;
    else                             vn8vn4SD =vn8vn4S/ vn8vn4D;
    if(vn8vn6D <= 0 || vn8vn6S <= 0) vn8vn6SD =-1;
    else                             vn8vn6SD =vn8vn6S/ vn8vn6D;
    if(g1eS < -100. || g1eD < -100.) g1eSD = -10000.;
    else                             g1eSD = g1eS / g1eD;

    double vn2SDe    = sqrt( pow(vn2S_Ste/vn2D, 2) + pow(vn2S*vn2D_Tote/vn2D/vn2D, 2) );
    double vn4SDe    = sqrt( pow(vn4S_Ste/vn4D, 2) + pow(vn4S*vn4D_Tote/vn4D/vn4D, 2) );
    double vn6SDe    = sqrt( pow(vn6S_Ste/vn6D, 2) + pow(vn6S*vn6D_Tote/vn6D/vn6D, 2) );
    double vn8SDe    = sqrt( pow(vn8S_Ste/vn8D, 2) + pow(vn8S*vn8D_Tote/vn8D/vn8D, 2) );
    double vn6vn4SDe = sqrt( pow(vn6vn4S_Ste/vn6vn4D, 2) + pow(vn6vn4S*vn6vn4D_Tote/vn6vn4D/vn6vn4D, 2) );
    double vn8vn4SDe = sqrt( pow(vn8vn4S_Ste/vn8vn4D, 2) + pow(vn8vn4S*vn8vn4D_Tote/vn8vn4D/vn8vn4D, 2) );
    double vn8vn6SDe = sqrt( pow(vn8vn6S_Ste/vn8vn6D, 2) + pow(vn8vn6S*vn8vn6D_Tote/vn8vn6D/vn8vn6D, 2) );
    double g1eSDe    = sqrt( pow(g1eS_Ste/g1eD, 2) + pow(g1eS*g1eD_Tote/g1eD/g1eD, 2) );

    vn2SVD_DAg[icent]    = vn2SD;
    vn4SVD_DAg[icent]    = vn4SD;
    vn6SVD_DAg[icent]    = vn6SD;
    vn8SVD_DAg[icent]    = vn8SD;
    vn6vn4SVD_DAg[icent] = vn6vn4SD;
    vn8vn4SVD_DAg[icent] = vn8vn4SD;
    vn8vn6SVD_DAg[icent] = vn8vn6SD;
    g1eSVD_DAg[icent]    = g1eSD;

    vn2SVD_DAg_err[icent]    = vn2SDe;
    vn4SVD_DAg_err[icent]    = vn4SDe;
    vn6SVD_DAg_err[icent]    = vn6SDe;
    vn8SVD_DAg_err[icent]    = vn8SDe;
    vn6vn4SVD_DAg_err[icent] = vn6vn4SDe;
    vn8vn4SVD_DAg_err[icent] = vn8vn4SDe;
    vn8vn6SVD_DAg_err[icent] = vn8vn6SDe;
    g1eSVD_DAg_err[icent]    = g1eSDe;

  }

  //-- TGraph time
  grVn2RawRatio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn2SVD_DAg,    CERR, vn2SVD_DAg_err);
  grVn4RawRatio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn4SVD_DAg,    CERR, vn4SVD_DAg_err);
  grVn6RawRatio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn6SVD_DAg,    CERR, vn6SVD_DAg_err);
  grVn8RawRatio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn8SVD_DAg,    CERR, vn8SVD_DAg_err);
  grvn6vn4Ratio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn6vn4SVD_DAg, CERR, vn6vn4SVD_DAg_err);
  grvn8vn4Ratio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn8vn4SVD_DAg, CERR, vn8vn4SVD_DAg_err);
  grvn8vn6Ratio_SVD_DAg    = new TGraphErrors(NCENT, centBinCenter, vn8vn6SVD_DAg, CERR, vn8vn6SVD_DAg_err);
  grGamma1ExpRatio_SVD_DAg = new TGraphErrors(NCENT, centBinCenter, g1eSVD_DAg,    CERR, g1eSVD_DAg_err);

  formatGraph(grVn2RawRatio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{2} Ratio: SVD/DAg", norder_),                    9,         20, "grVn2RawRatio_SVD_DAg");
  formatGraph(grVn4RawRatio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{4} Ratio: SVD/DAg", norder_),                    kSpring+4, 20, "grVn4RawRatio_SVD_DAg");
  formatGraph(grVn6RawRatio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{6} Ratio: SVD/DAg", norder_),                    6,         20, "grVn6RawRatio_SVD_DAg");
  formatGraph(grVn8RawRatio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{8} Ratio: SVD/DAg", norder_),                    kOrange+7, 20, "grVn8RawRatio_SVD_DAg");
  formatGraph(grvn6vn4Ratio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{6}/v_{%i}{4} Ratio: SVD/DAg", norder_, norder_), 4,         20, "grvn6vn4Ratio_SVD_DAg");
  formatGraph(grvn8vn4Ratio_SVD_DAg,    "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{8}/v_{%i}{4} Ratio: SVD/DAg", norder_, norder_), kGreen+2,  20, "grvn8vn4Ratio_SVD_DAg");
  formatGraph(grvn8vn6Ratio_SVD_DAg,    "Centrality %", ratioMinVn8Vn6, ratioMaxVn8Vn6, Form("v_{%i}{8}/v_{%i}{6} Ratio: SVD/DAg", norder_, norder_), kViolet-1, 20, "grvn8vn6Ratio_SVD_DAg");
  formatGraph(grGamma1ExpRatio_SVD_DAg, "Centrality %", ratioMinG1E,    ratioMaxG1E,    "#gamma_{1}^{exp} Ratio: SVD/DAg",                            2,         20, "grGamma1ExpRatio_SVD_DAg");


  //-- DRAW

  TLine * lOne = new TLine(0, 1.0, grVn2RawRatio_SVD_DAg->GetXaxis()->GetXmax(), 1.0);
  lOne->SetLineStyle(2);
  lOne->SetLineWidth(2);

  //-- cumus
  TLegend * legvn2 = new TLegend(0.4939, 0.2013, 0.8977, 0.3675);
  legvn2->SetBorderSize(0);
  legvn2->SetFillStyle(0);
  legvn2->AddEntry(grVn2RawDAg, "D'Agostini", "lp");
  legvn2->AddEntry(grVn2RawSVD, "SVD",        "lp");

  TLegend * legvn4 = new TLegend(0.4939, 0.2013, 0.8977, 0.3675);
  legvn4->SetBorderSize(0);
  legvn4->SetFillStyle(0);
  legvn4->AddEntry(grVn4RawDAg, "D'Agostini", "lp");
  legvn4->AddEntry(grVn4RawSVD, "SVD",        "lp");

  TLegend * legvn6 = new TLegend(0.4939, 0.2013, 0.8977, 0.3675);
  legvn6->SetBorderSize(0);
  legvn6->SetFillStyle(0);
  legvn6->AddEntry(grVn6RawDAg, "D'Agostini", "lp");
  legvn6->AddEntry(grVn6RawSVD, "SVD",        "lp");

  TLegend * legvn8 = new TLegend(0.4939, 0.2013, 0.8977, 0.3675);
  legvn8->SetBorderSize(0);
  legvn8->SetFillStyle(0);
  legvn8->AddEntry(grVn8RawDAg, "D'Agostini", "lp");
  legvn8->AddEntry(grVn8RawSVD, "SVD",        "lp");


  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 2000, 1000);
  cCumu->Divide(4, 2);

  cCumu->cd(1);
  grVn2RawDAgSys->Draw("apE2");
  grVn2RawDAg->Draw("psame");
  grVn2RawSVD->Draw("psame");
  legvn2->Draw("same");

  cCumu->cd(2);
  grVn4RawDAgSys->Draw("apE2");
  grVn4RawDAg->Draw("psame");
  grVn4RawSVD->Draw("psame");
  legvn4->Draw("same");

  cCumu->cd(3);
  grVn6RawDAgSys->Draw("apE2");
  grVn6RawDAg->Draw("psame");
  grVn6RawSVD->Draw("psame");
  legvn6->Draw("same");

  cCumu->cd(4);
  grVn8RawDAgSys->Draw("apE2");
  grVn8RawDAg->Draw("psame");
  grVn8RawSVD->Draw("psame");
  legvn8->Draw("same");

  cCumu->cd(5);
  grVn2RawRatio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumu->cd(6);
  grVn4RawRatio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumu->cd(7);
  grVn6RawRatio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumu->cd(8);
  grVn8RawRatio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumu->SaveAs("cumu.pdf");

  //-- cumu Ratios
  TLegend * legvn6vn4 = new TLegend(0.1846, 0.1898, 0.5894, 0.3560);
  legvn6vn4->SetBorderSize(0);
  legvn6vn4->SetFillStyle(0);
  legvn6vn4->AddEntry(grvn6vn4RatioDAg, "D'Agostini", "lp");
  legvn6vn4->AddEntry(grvn6vn4RatioSVD, "SVD",        "lp");
 
  TLegend * legvn8vn4 = new TLegend(0.1846, 0.1898, 0.5894, 0.3560);
  legvn8vn4->SetBorderSize(0);
  legvn8vn4->SetFillStyle(0);
  legvn8vn4->AddEntry(grvn8vn4RatioDAg, "D'Agostini", "lp");
  legvn8vn4->AddEntry(grvn8vn4RatioSVD, "SVD",        "lp");

  TLegend * legvn8vn6 = new TLegend(0.1846, 0.1898, 0.5894, 0.3560);
  legvn8vn6->SetBorderSize(0);
  legvn8vn6->SetFillStyle(0);
  legvn8vn6->AddEntry(grvn8vn6RatioDAg, "D'Agostini", "lp");
  legvn8vn6->AddEntry(grvn8vn6RatioSVD, "SVD",        "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 1000);
  cCumuRatio->Divide(3, 2);

  cCumuRatio->cd(1);
  grvn6vn4RatioDAgSys->Draw("apE2");
  grvn6vn4RatioDAg->Draw("psame");
  grvn6vn4RatioSVD->Draw("psame");
  legvn6vn4->Draw("same");

  cCumuRatio->cd(2);
  grvn8vn4RatioDAgSys->Draw("apE2");
  grvn8vn4RatioDAg->Draw("psame");
  grvn8vn4RatioSVD->Draw("psame");
  legvn8vn4->Draw("same");

  cCumuRatio->cd(3);
  grvn8vn6RatioDAgSys->Draw("apE2");
  grvn8vn6RatioDAg->Draw("psame");
  grvn8vn6RatioSVD->Draw("psame");
  legvn8vn6->Draw("same");

  cCumuRatio->cd(4);
  grvn6vn4Ratio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumuRatio->cd(5);
  grvn8vn4Ratio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumuRatio->cd(6);
  grvn8vn6Ratio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cCumuRatio->SaveAs("cumuRatio.pdf");

  //-- Gamma1Exp
  TLegend * legg1e = new TLegend(0.1846, 0.1898, 0.5894, 0.3560);
  legg1e->SetBorderSize(0);
  legg1e->SetFillStyle(0);
  legg1e->AddEntry(grGamma1ExpDAg, "D'Agostini", "lp");
  legg1e->AddEntry(grGamma1ExpSVD, "SVD",        "lp");

  TCanvas * cG1E = new TCanvas("cG1E", "cG1E", 1000, 500);
  cG1E->Divide(2, 1);

  cG1E->cd(1);
  grGamma1ExpDAgSys->Draw("apE2");
  grGamma1ExpDAg->Draw("psame");
  grGamma1ExpSVD->Draw("psame");
  legg1e->Draw("same");

  cG1E->cd(2);
  grGamma1ExpRatio_SVD_DAg->Draw("ap");
  lOne->Draw("same");

  cG1E->SaveAs("G1E.pdf");

  //-- Unfold Distns

  TLegend * legUnf = new TLegend(0.6, 0.6, 0.9, 0.8);
  legInit(legUnf);
  legUnf->AddEntry(hUnfoldStatDAg[0], "D'Agostini", "lp");
  legUnf->AddEntry(hUnfoldSVD[0],     "SVD",        "lp");

  cUnfold = new TCanvas("cUnfold", "cUnfold", 2000, 1500);
  cUnfold->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    cUnfold->cd(icent+1);
    cUnfold->cd(icent+1)->SetLogy();
    hUnfoldSysDAg[icent]->Draw("e2");
    setex2->Draw();
    if(icent == 0) legUnf->Draw("same");
  }
  for(int icent= 0; icent < NCENT; icent++){
    cUnfold->cd(icent+1);
    hUnfoldStatDAg[icent]->Draw("same");
    hUnfoldSVD[icent]->Draw("same");
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
  }
  cUnfold->SaveAs("UnfoldCompare.pdf");

  //-- Unfold Ratio SVD/DAg
  cUnfoldRatio_SVD_DAg = new TCanvas("cUnfoldRatio_SVD_DAg", "cUnfoldRatio_SVD_DAg", 2000, 1500);
  cUnfoldRatio_SVD_DAg->Divide(4,3);
  for(int icent= 0; icent < NCENT; icent++){
    cUnfoldRatio_SVD_DAg->cd(icent+1);
    hUnfoldRatio_SVD_DAg[icent]->Draw();
    latex.DrawLatex(0.67, 0.88, Form( "Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );
  }
  cUnfoldRatio_SVD_DAg->SaveAs("UnfoldRatioCompare.pdf");

}

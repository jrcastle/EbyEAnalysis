#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void svdUnfoldResults(){

  int norder_ = 2;
  double tkEta = 2.4;

  double cumuMin = 0.0;
  double cumuMax = 0.15;

  double g1eMin   = -1.0;
  double g1eMax   = 0.5;

  double ratioMin       = 0.9;
  double ratioMax       = 1.03;
  double ratioMinVn8Vn6 = 0.97;
  double ratioMaxVn8Vn6 = 1.002;

  TFile * fSVD;
  TH1D * finalUnf[NCENT];

  TFile * fStatSVD;
  TH1D * hVarianceOfMean_Vn2;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
  TH1D * hVarianceOfMean_Gamma1Exp;
  TH1D * hVarianceOfMean_Vn6Vn4;
  TH1D * hVarianceOfMean_Vn8Vn4;
  TH1D * hVarianceOfMean_Vn8Vn6;

  double Vn2[NCENT];
  double Vn4[NCENT];
  double Vn6[NCENT];
  double Vn8[NCENT];

  double Vn6Vn4[NCENT];
  double Vn8Vn4[NCENT];
  double Vn8Vn6[NCENT];
  double G1E[NCENT];

  double Vn2_err[NCENT];
  double Vn4_err[NCENT];
  double Vn6_err[NCENT];
  double Vn8_err[NCENT];

  double Vn6Vn4_err[NCENT];
  double Vn8Vn4_err[NCENT];
  double Vn8Vn6_err[NCENT];
  double G1E_err[NCENT];

  TGraphErrors * grVn2;
  TGraphErrors * grVn4;
  TGraphErrors * grVn6;
  TGraphErrors * grVn8;
  TGraphErrors * grVn6Vn4;
  TGraphErrors * grVn8Vn4;
  TGraphErrors * grVn8Vn6;
  TGraphErrors * grG1E;

  TFile * fOut;

  //
  // MAIN
  //
  setTDRStyle();

  //-- Get Stat Errors
  fStatSVD = new TFile( Form("../../../statErrorHandle/v%i/eta%.1f/StatisticalUncertaintiesSVD_v%i.root", norder_, tkEta, norder_)  );
  hVarianceOfMean_Vn2       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn2");
  hVarianceOfMean_Vn4       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn4");
  hVarianceOfMean_Vn6       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn6");
  hVarianceOfMean_Vn8       = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8");
  hVarianceOfMean_Gamma1Exp = (TH1D*) fStatSVD->Get("hVarianceOfMean_Gamma1Exp");
  hVarianceOfMean_Vn6Vn4    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn6Vn4");
  hVarianceOfMean_Vn8Vn4    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8Vn4");
  hVarianceOfMean_Vn8Vn6    = (TH1D*) fStatSVD->Get("hVarianceOfMean_Vn8Vn6");

  //-- Get unfold file
  fSVD = new TFile("../UnfoldResults/dataResp/data2_svd.root");

  for(int icent = 0; icent < NCENT; icent++){

    finalUnf[icent] = (TH1D*) fSVD->Get( Form("hrecokreg4_c%i", icent) ); 
    finalUnf[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    finalUnf[icent]->SetLineColor(4);
    finalUnf[icent]->SetMarkerColor(4);

    EbyECumu cumu( finalUnf[icent] );
    double vn2 = cumu.GetCumu_vn2();
    double vn4 = cumu.GetCumu_vn4();
    double vn6 = cumu.GetCumu_vn6();
    double vn8 = cumu.GetCumu_vn8();
    double g1e = cumu.GetGamma1Exp();

    double vn6vn4;
    double vn8vn4;
    double vn8vn6;

    if(vn4 == 0 || vn6 == 0) vn6vn4 = -1;
    else                     vn6vn4 = vn6 / vn4;
    if(vn4 == 0 || vn8 == 0) vn8vn4 = -1;
    else                     vn8vn4 = vn8 / vn4;
    if(vn6 == 0 || vn8 == 0) vn8vn6 = -1;
    else                     vn8vn6 = vn8 / vn6;

    double vn2e    = sqrt( hVarianceOfMean_Vn2->GetBinContent(icent+1) );
    double vn4e    = sqrt( hVarianceOfMean_Vn4->GetBinContent(icent+1) );
    double vn6e    = sqrt( hVarianceOfMean_Vn6->GetBinContent(icent+1) );
    double vn8e    = sqrt( hVarianceOfMean_Vn8->GetBinContent(icent+1) );
    double vn6vn4e = sqrt( hVarianceOfMean_Vn6Vn4->GetBinContent(icent+1) );
    double vn8vn4e = sqrt( hVarianceOfMean_Vn8Vn4->GetBinContent(icent+1) );
    double vn8vn6e = sqrt( hVarianceOfMean_Vn8Vn6->GetBinContent(icent+1) );
    double g1ee    = sqrt( hVarianceOfMean_Gamma1Exp->GetBinContent(icent+1) );

    Vn2[icent] = vn2;
    Vn4[icent] = vn4;
    Vn6[icent] = vn6;
    Vn8[icent] = vn8;

    Vn6Vn4[icent] = vn6vn4;
    Vn8Vn4[icent] = vn8vn4;
    Vn8Vn6[icent] = vn8vn6;
    G1E[icent]    = g1e;

    Vn2_err[icent] = vn2e;
    Vn4_err[icent] = vn4e;
    Vn6_err[icent] = vn6e;
    Vn8_err[icent] = vn8e;

    Vn6Vn4_err[icent] = vn6vn4e;
    Vn8Vn4_err[icent] = vn8vn4e;
    Vn8Vn6_err[icent] = vn8vn6e;
    G1E_err[icent]    = g1ee;

  }

  //- Tgraph time
  grVn2    = new TGraphErrors(NCENT, centBinCenter, Vn2,    CERR, Vn2_err);
  grVn4    = new TGraphErrors(NCENT, centBinCenter, Vn4,    CERR, Vn4_err);
  grVn6    = new TGraphErrors(NCENT, centBinCenter, Vn6,    CERR, Vn6_err);
  grVn8    = new TGraphErrors(NCENT, centBinCenter, Vn8,    CERR, Vn8_err);
  grVn6Vn4 = new TGraphErrors(NCENT, centBinCenter, Vn6Vn4, CERR, Vn6Vn4_err);
  grVn8Vn4 = new TGraphErrors(NCENT, centBinCenter, Vn8Vn4, CERR, Vn8Vn4_err);
  grVn8Vn6 = new TGraphErrors(NCENT, centBinCenter, Vn8Vn6, CERR, Vn8Vn6_err);
  grG1E    = new TGraphErrors(NCENT, centBinCenter, G1E,    CERR, G1E_err);

  formatGraph(grVn2,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{2}", norder_),                      9,         24, "grVn2");
  formatGraph(grVn4,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{4}", norder_),                      kSpring+4, 24, "grVn4");
  formatGraph(grVn6,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{6}", norder_),                      6,         24, "grVn6");
  formatGraph(grVn8,    "Centrality %", cumuMin,        cumuMax,        Form("v_{%i}{8}", norder_),                      kOrange+7, 24, "grVn8");
  formatGraph(grVn6Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_), 4,         24, "grVn6Vn4");
  formatGraph(grVn8Vn4, "Centrality %", ratioMin,       ratioMax,       Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_), kGreen+2,  24, "grVn8Vn4");
  formatGraph(grVn8Vn6, "Centrality %", ratioMinVn8Vn6, ratioMaxVn8Vn6, Form("v_{%i}{8} / v_{%i}{6}", norder_, norder_), kViolet-1, 24, "grVn8Vn6");
  formatGraph(grG1E,    "Centrality %", g1eMin,         g1eMax,         "#gamma_{1}^{Exp}",                              2,         24, "grG1E");

  //Save to a file
  fOut = new TFile("SVDPhysics.root", "recreate");
  grVn2->Write();
  grVn4->Write();
  grVn6->Write();
  grVn8->Write();
  grVn6Vn4->Write();
  grVn8Vn4->Write();
  grVn8Vn6->Write();
  grG1E->Write();
  for(int icent = 0; icent < NCENT; icent++) finalUnf[icent]->Write();

}

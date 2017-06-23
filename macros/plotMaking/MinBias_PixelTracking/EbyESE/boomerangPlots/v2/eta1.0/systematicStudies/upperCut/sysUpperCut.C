#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void FIXUNFOLD(TH1D * h, double sig){
  double upBound = h->GetMean() + sig*h->GetRMS();
  int nb = h->GetNbinsX();
  for(int ibin = 1; ibin < nb; ibin++){
    double bc = h->GetBinCenter(ibin);
    if( bc > upBound ){
      h->SetBinContent(ibin, 0);
      h->SetBinError(ibin, 0);
    }
  }
}

void sysUpperCut(){

  const int norder_ = 2;
  bool compATLAS    = 1;

  double vnCumuMin = 0.0;
  double vnCumuMax = 0.15;
  double g1eMin    = -1.0;
  double g1eMax    = 0.5;
  double vn6vn4Min = 0.95;
  double vn6vn4Max = 1.05;
  double vn8vn4Min = 0.95;
  double vn8vn4Max = 1.05;
  double vn8vn6Min = 0.99;
  double vn8vn6Max = 1.01;

  double rMinCumu = 0.98;
  double rMaxCumu = 1.02;
  double rMin6484 = 0.97;
  double rMax6484 = 1.02;
  double rMin86   = 0.997;
  double rMax86   = 1.003;
  double rMing1e  = 0.0;
  double rMaxg1e  = 3.0;

  TFile * fAna;
  TH1D * hObs[NCENT];

  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  TH1D * hFinalUnfold_Default[NCENT];
  TH1D * hFinalUnfold_Cut3p5Sigma[NCENT];
  TH1D * hFinalUnfold_Cut3p8Sigma[NCENT];
  TH1D * hFinalUnfold_Cut4Sigma[NCENT];
  TH1D * hFinalUnfold_Cut4p2Sigma[NCENT];
  TH1D * hFinalUnfold_Cut4p5Sigma[NCENT];

  TLatex latex;

  //-- Vectors galore!
  //-- 3.5 sigma
  std::vector<double> vn2_3p5Sigma;
  std::vector<double> vn2Cent_3p5Sigma;
  std::vector<double> vn4_3p5Sigma;
  std::vector<double> vn4Cent_3p5Sigma;
  std::vector<double> vn6_3p5Sigma;
  std::vector<double> vn6Cent_3p5Sigma;
  std::vector<double> vn8_3p5Sigma;
  std::vector<double> vn8Cent_3p5Sigma;

  std::vector<double> g1e_3p5Sigma;
  std::vector<double> g1eCent_3p5Sigma;
  std::vector<double> vn6vn4_3p5Sigma;
  std::vector<double> vn6vn4Cent_3p5Sigma;
  std::vector<double> vn6vn4Npart_3p5Sigma;
  std::vector<double> vn8vn4_3p5Sigma;
  std::vector<double> vn8vn4Cent_3p5Sigma;
  std::vector<double> vn8vn4Npart_3p5Sigma;
  std::vector<double> vn8vn6_3p5Sigma;
  std::vector<double> vn8vn6Cent_3p5Sigma;

  //-- 3.8 Sigma
  std::vector<double> vn2_3p8Sigma;
  std::vector<double> vn2Cent_3p8Sigma;
  std::vector<double> vn4_3p8Sigma;
  std::vector<double> vn4Cent_3p8Sigma;
  std::vector<double> vn6_3p8Sigma;
  std::vector<double> vn6Cent_3p8Sigma;
  std::vector<double> vn8_3p8Sigma;
  std::vector<double> vn8Cent_3p8Sigma;

  std::vector<double> g1e_3p8Sigma;
  std::vector<double> g1eCent_3p8Sigma;
  std::vector<double> vn6vn4_3p8Sigma;
  std::vector<double> vn6vn4Cent_3p8Sigma;
  std::vector<double> vn6vn4Npart_3p8Sigma;
  std::vector<double> vn8vn4_3p8Sigma;
  std::vector<double> vn8vn4Cent_3p8Sigma;
  std::vector<double> vn8vn4Npart_3p8Sigma;
  std::vector<double> vn8vn6_3p8Sigma;
  std::vector<double> vn8vn6Cent_3p8Sigma;

  //-- 4 Sigma
  std::vector<double> vn2_4Sigma;
  std::vector<double> vn2Cent_4Sigma;
  std::vector<double> vn4_4Sigma;
  std::vector<double> vn4Cent_4Sigma;
  std::vector<double> vn6_4Sigma;
  std::vector<double> vn6Cent_4Sigma;
  std::vector<double> vn8_4Sigma;
  std::vector<double> vn8Cent_4Sigma;

  std::vector<double> g1e_4Sigma;
  std::vector<double> g1eCent_4Sigma;
  std::vector<double> vn6vn4_4Sigma;
  std::vector<double> vn6vn4Cent_4Sigma;
  std::vector<double> vn6vn4Npart_4Sigma;
  std::vector<double> vn8vn4_4Sigma;
  std::vector<double> vn8vn4Cent_4Sigma;
  std::vector<double> vn8vn4Npart_4Sigma;
  std::vector<double> vn8vn6_4Sigma;
  std::vector<double> vn8vn6Cent_4Sigma;

  //-- 4.2 Sigma
  std::vector<double> vn2_4p2Sigma;
  std::vector<double> vn2Cent_4p2Sigma;
  std::vector<double> vn4_4p2Sigma;
  std::vector<double> vn4Cent_4p2Sigma;
  std::vector<double> vn6_4p2Sigma;
  std::vector<double> vn6Cent_4p2Sigma;
  std::vector<double> vn8_4p2Sigma;
  std::vector<double> vn8Cent_4p2Sigma;

  std::vector<double> g1e_4p2Sigma;
  std::vector<double> g1eCent_4p2Sigma;
  std::vector<double> vn6vn4_4p2Sigma;
  std::vector<double> vn6vn4Cent_4p2Sigma;
  std::vector<double> vn6vn4Npart_4p2Sigma;
  std::vector<double> vn8vn4_4p2Sigma;
  std::vector<double> vn8vn4Cent_4p2Sigma;
  std::vector<double> vn8vn4Npart_4p2Sigma;
  std::vector<double> vn8vn6_4p2Sigma;
  std::vector<double> vn8vn6Cent_4p2Sigma;

  //-- 4.5 sigma 
  std::vector<double> vn2_4p5Sigma;
  std::vector<double> vn2Cent_4p5Sigma;
  std::vector<double> vn4_4p5Sigma;
  std::vector<double> vn4Cent_4p5Sigma;
  std::vector<double> vn6_4p5Sigma;
  std::vector<double> vn6Cent_4p5Sigma;
  std::vector<double> vn8_4p5Sigma;
  std::vector<double> vn8Cent_4p5Sigma;

  std::vector<double> g1e_4p5Sigma;
  std::vector<double> g1eCent_4p5Sigma;
  std::vector<double> vn6vn4_4p5Sigma;
  std::vector<double> vn6vn4Cent_4p5Sigma;
  std::vector<double> vn6vn4Npart_4p5Sigma;
  std::vector<double> vn8vn4_4p5Sigma;
  std::vector<double> vn8vn4Cent_4p5Sigma;
  std::vector<double> vn8vn4Npart_4p5Sigma;
  std::vector<double> vn8vn6_4p5Sigma;
  std::vector<double> vn8vn6Cent_4p5Sigma;

  TGraph * grVn2_3p5Sigma;
  TGraph * grVn4_3p5Sigma;
  TGraph * grVn6_3p5Sigma;
  TGraph * grVn8_3p5Sigma;
  TGraph * grG1e_3p5Sigma;
  TGraph * grVn6Vn4_3p5Sigma;
  TGraph * grVn8Vn4_3p5Sigma;
  TGraph * grVn8Vn6_3p5Sigma;

  TGraph * grVn2_3p8Sigma;
  TGraph * grVn4_3p8Sigma;
  TGraph * grVn6_3p8Sigma;
  TGraph * grVn8_3p8Sigma;
  TGraph * grG1e_3p8Sigma;
  TGraph * grVn6Vn4_3p8Sigma;
  TGraph * grVn8Vn4_3p8Sigma;
  TGraph * grVn8Vn6_3p8Sigma;

  TGraphErrors * grVn2_4Sigma;
  TGraphErrors * grVn4_4Sigma;
  TGraphErrors * grVn6_4Sigma;
  TGraphErrors * grVn8_4Sigma;
  TGraphErrors * grG1e_4Sigma;
  TGraphErrors * grVn6Vn4_4Sigma;
  TGraphErrors * grVn8Vn4_4Sigma;
  TGraphErrors * grVn8Vn6_4Sigma;

  TGraphErrors * grVn2Sys_4Sigma;
  TGraphErrors * grVn4Sys_4Sigma;
  TGraphErrors * grVn6Sys_4Sigma;
  TGraphErrors * grVn8Sys_4Sigma;
  TGraphErrors * grG1eSys_4Sigma;
  TGraphErrors * grVn6Vn4Sys_4Sigma;
  TGraphErrors * grVn8Vn4Sys_4Sigma;
  TGraphErrors * grVn8Vn6Sys_4Sigma;

  TGraph * grVn2_4p2Sigma;
  TGraph * grVn4_4p2Sigma;
  TGraph * grVn6_4p2Sigma;
  TGraph * grVn8_4p2Sigma;
  TGraph * grG1e_4p2Sigma;
  TGraph * grVn6Vn4_4p2Sigma;
  TGraph * grVn8Vn4_4p2Sigma;
  TGraph * grVn8Vn6_4p2Sigma;

  TGraph * grVn2_4p5Sigma;
  TGraph * grVn4_4p5Sigma;
  TGraph * grVn6_4p5Sigma;
  TGraph * grVn8_4p5Sigma;
  TGraph * grG1e_4p5Sigma;
  TGraph * grVn6Vn4_4p5Sigma;
  TGraph * grVn8Vn4_4p5Sigma;
  TGraph * grVn8Vn6_4p5Sigma;

  //-- 3p5 sigma / Default
  std::vector<double> vn2_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn4_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8_3p5Sigma_RatioTo4Sigma;
  std::vector<double> g1e_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6_3p5Sigma_RatioTo4Sigma;

  std::vector<double> vn2_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn4_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> g1e_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6vn4_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn4_3p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn6_3p5Sigma_RatioTo4Sigma_err;

  std::vector<double> vn2Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn4Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> g1eCent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4Cent_3p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6Cent_3p5Sigma_RatioTo4Sigma;


  //-- 3p8Sigma / Default
  std::vector<double> vn2_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn4_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn6_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8_3p8Sigma_RatioTo4Sigma;
  std::vector<double> g1e_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6_3p8Sigma_RatioTo4Sigma;

  std::vector<double> vn2_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn4_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> g1e_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6vn4_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn4_3p8Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn6_3p8Sigma_RatioTo4Sigma_err;

  std::vector<double> vn2Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn4Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn6Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> g1eCent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4Cent_3p8Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6Cent_3p8Sigma_RatioTo4Sigma;

  //-- 4p2Sigma / Default
  std::vector<double> vn2_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn4_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn6_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8_4p2Sigma_RatioTo4Sigma;
  std::vector<double> g1e_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6_4p2Sigma_RatioTo4Sigma;

  std::vector<double> vn2_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn4_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> g1e_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6vn4_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn4_4p2Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn6_4p2Sigma_RatioTo4Sigma_err;

  std::vector<double> vn2Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn4Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn6Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> g1eCent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4Cent_4p2Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6Cent_4p2Sigma_RatioTo4Sigma;

  //-- 4.5 sigma / Default
  std::vector<double> vn2_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn4_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8_4p5Sigma_RatioTo4Sigma;
  std::vector<double> g1e_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6_4p5Sigma_RatioTo4Sigma;

  std::vector<double> vn2_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn4_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> g1e_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn6vn4_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn4_4p5Sigma_RatioTo4Sigma_err;
  std::vector<double> vn8vn6_4p5Sigma_RatioTo4Sigma_err;

  std::vector<double> vn2Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn4Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> g1eCent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn6vn4Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn4Cent_4p5Sigma_RatioTo4Sigma;
  std::vector<double> vn8vn6Cent_4p5Sigma_RatioTo4Sigma;

  TGraphErrors * grVn2_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn4_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grG1e_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6Vn4_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn4_3p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn6_3p5Sigma_RatioTo4Sigma;

  TGraphErrors * grVn2_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn4_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grG1e_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6Vn4_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn4_3p8Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn6_3p8Sigma_RatioTo4Sigma;

  TGraphErrors * grVn2_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn4_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grG1e_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6Vn4_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn4_4p2Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn6_4p2Sigma_RatioTo4Sigma;

  TGraphErrors * grVn2_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn4_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grG1e_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn6Vn4_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn4_4p5Sigma_RatioTo4Sigma;
  TGraphErrors * grVn8Vn6_4p5Sigma_RatioTo4Sigma;

  TFile * fPhy;
  TFile * fPhyUnf;

  TLine * l3p5s[NCENT];
  TLine * l3p8s[NCENT];
  TLine * l4s[NCENT];
  TLine * l4p2s[NCENT];
  TLine * l4p5s[NCENT];

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //-- Get the output files for the reported distributions
  fAna    = new TFile( "../../AnalyzerResults/CastleEbyE.root" );
  fUnfold = new TFile( Form("../../UnfoldResults/dataResp/data%i.root", norder_) );
  fPhy    = new TFile("../PhysicsResults.root");
  fPhyUnf = new TFile( Form("../SysUnfoldDistns_v%i.root", norder_) );

  //-- Get final physics results
  grVn2_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn2Raw" );
  grVn4_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn4Raw" );
  grVn6_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn6Raw" );
  grVn8_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn8Raw" );
  grG1e_4Sigma    = (TGraphErrors*) fPhy->Get( "grGamma1Exp" );
  grVn6Vn4_4Sigma = (TGraphErrors*) fPhy->Get( "grvn6vn4Ratio" );
  grVn8Vn4_4Sigma = (TGraphErrors*) fPhy->Get( "grvn8vn4Ratio" );
  grVn8Vn6_4Sigma = (TGraphErrors*) fPhy->Get( "grvn8vn6Ratio" );

  grVn2Sys_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn2RawSys" );
  grVn4Sys_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn4RawSys" );
  grVn6Sys_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn6RawSys" );
  grVn8Sys_4Sigma    = (TGraphErrors*) fPhy->Get( "grVn8RawSys" );
  grG1eSys_4Sigma    = (TGraphErrors*) fPhy->Get( "grGamma1ExpSys" );
  grVn6Vn4Sys_4Sigma = (TGraphErrors*) fPhy->Get( "grvn6vn4RatioSys" );
  grVn8Vn4Sys_4Sigma = (TGraphErrors*) fPhy->Get( "grvn8vn4RatioSys" );
  grVn8Vn6Sys_4Sigma = (TGraphErrors*) fPhy->Get( "grvn8vn6RatioSys" );

  //-- Start fetching histograms....
  for(int icent = 0; icent < NCENT; icent++){

    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    for(int i = 0; i < NITER; i++){
      hUnfold[icent][i]     = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i]     = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");

      if( chi2 < 1.2 ){
	hFinalUnfold_Default[icent]     = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Default_c%i", icent) );
	hFinalUnfold_Cut3p5Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut3p5Sigma_c%i", icent) );
	hFinalUnfold_Cut3p8Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut3p8Sigma_c%i", icent) );
	hFinalUnfold_Cut4Sigma[icent]   = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4Sigma_c%i", icent) );
	hFinalUnfold_Cut4p2Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4p2Sigma_c%i", icent) );
	hFinalUnfold_Cut4p5Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4p5Sigma_c%i", icent) );
	break;
      }
      if(i == NITER-1){
	hFinalUnfold_Default[icent]     = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Default_c%i", icent) );
        hFinalUnfold_Cut3p5Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut3p5Sigma_c%i", icent) );
        hFinalUnfold_Cut3p8Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut3p8Sigma_c%i", icent) );
        hFinalUnfold_Cut4Sigma[icent]   = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4Sigma_c%i", icent) );
        hFinalUnfold_Cut4p2Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4p2Sigma_c%i", icent) );
	hFinalUnfold_Cut4p5Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4p5Sigma_c%i", icent) );
        break;
      }

    } //-- End iter loop

  } //-- End cent loop

  //-- Now let's calculate some values!
  for(int icent = 0; icent < NCENT; icent++){

    //-- 3p5 sigma
    FIXUNFOLD( hFinalUnfold_Cut3p5Sigma[icent], 3.5 );
    EbyECumu cumu3p5Sig( hFinalUnfold_Cut3p5Sigma[icent] );
    double vn2_3p5Sig    = cumu3p5Sig.GetCumu_vn2();
    double vn4_3p5Sig    = cumu3p5Sig.GetCumu_vn4();
    double vn6_3p5Sig    = cumu3p5Sig.GetCumu_vn6();
    double vn8_3p5Sig    = cumu3p5Sig.GetCumu_vn8();
    double g1e_3p5Sig    = cumu3p5Sig.GetGamma1Exp();
    double vn6vn4_3p5Sig;
    double vn8vn4_3p5Sig;
    double vn8vn6_3p5Sig;

    if(vn2_3p5Sig != 0){
      vn2_3p5Sigma.push_back(vn2_3p5Sig);
      vn2Cent_3p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p5Sig != 0){
      vn4_3p5Sigma.push_back(vn4_3p5Sig);
      vn4Cent_3p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p5Sig != 0){
      vn6_3p5Sigma.push_back(vn6_3p5Sig);
      vn6Cent_3p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_3p5Sig != 0){
      vn8_3p5Sigma.push_back(vn8_3p5Sig);
      vn8Cent_3p5Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_3p5Sig != -10000){
      g1e_3p5Sigma.push_back(g1e_3p5Sig);
      g1eCent_3p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p5Sig != 0 && vn6_3p5Sig != 0){
      vn6vn4_3p5Sigma.push_back( vn6_3p5Sig/vn4_3p5Sig );
      vn6vn4Cent_3p5Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_3p5Sigma.push_back(Npart[icent]);
    }
    if(vn4_3p5Sig != 0 && vn8_3p5Sig!= 0){
      vn8vn4_3p5Sigma.push_back( vn8_3p5Sig/vn4_3p5Sig );
      vn8vn4Cent_3p5Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_3p5Sigma.push_back(Npart[icent]);
    }
    if(vn6_3p5Sig != 0 && vn8_3p5Sig!= 0){
      vn8vn6_3p5Sigma.push_back( vn8_3p5Sig/vn6_3p5Sig );
      vn8vn6Cent_3p5Sigma.push_back(centBinCenter[icent]);
    }

    //-- 3p8Sigma
    FIXUNFOLD( hFinalUnfold_Cut3p8Sigma[icent], 3.8 );
    EbyECumu cumu3p8Sig( hFinalUnfold_Cut3p8Sigma[icent] );
    double vn2_3p8Sig    = cumu3p8Sig.GetCumu_vn2();
    double vn4_3p8Sig    = cumu3p8Sig.GetCumu_vn4();
    double vn6_3p8Sig    = cumu3p8Sig.GetCumu_vn6();
    double vn8_3p8Sig    = cumu3p8Sig.GetCumu_vn8();
    double g1e_3p8Sig    = cumu3p8Sig.GetGamma1Exp();
    double vn6vn4_3p8Sig;
    double vn8vn4_3p8Sig;
    double vn8vn6_3p8Sig;

    if(vn2_3p8Sig != 0){
      vn2_3p8Sigma.push_back(vn2_3p8Sig);
      vn2Cent_3p8Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p8Sig != 0){
      vn4_3p8Sigma.push_back(vn4_3p8Sig);
      vn4Cent_3p8Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p8Sig != 0){
      vn6_3p8Sigma.push_back(vn6_3p8Sig);
      vn6Cent_3p8Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_3p8Sig != 0){
      vn8_3p8Sigma.push_back(vn8_3p8Sig);
      vn8Cent_3p8Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_3p8Sig != -10000){
      g1e_3p8Sigma.push_back(g1e_3p8Sig);
      g1eCent_3p8Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p8Sig != 0 && vn6_3p8Sig != 0){
      vn6vn4_3p8Sigma.push_back( vn6_3p8Sig/vn4_3p8Sig );
      vn6vn4Cent_3p8Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_3p8Sigma.push_back(Npart[icent]);
    }
    if(vn4_3p8Sig != 0 && vn8_3p8Sig!= 0){
      vn8vn4_3p8Sigma.push_back( vn8_3p8Sig/vn4_3p8Sig );
      vn8vn4Cent_3p8Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_3p8Sigma.push_back(Npart[icent]);
    }
    if(vn6_3p8Sig != 0 && vn8_3p8Sig!= 0){
      vn8vn6_3p8Sigma.push_back( vn8_3p8Sig/vn6_3p8Sig );
      vn8vn6Cent_3p8Sigma.push_back(centBinCenter[icent]);
    }

    //-- 4Sigma
    double vn2_4Sig    = grVn2_4Sigma->GetY()[icent];
    double vn4_4Sig    = grVn4_4Sigma->GetY()[icent];
    double vn6_4Sig    = grVn6_4Sigma->GetY()[icent];
    double vn8_4Sig    = grVn8_4Sigma->GetY()[icent];
    double g1e_4Sig    = grG1e_4Sigma->GetY()[icent];
    double vn6vn4_4Sig = grVn6Vn4_4Sigma->GetY()[icent];
    double vn8vn4_4Sig = grVn8Vn4_4Sigma->GetY()[icent];
    double vn8vn6_4Sig = grVn8Vn6_4Sigma->GetY()[icent];

    double vn2_4Sige    = sqrt( pow(grVn2_4Sigma->GetErrorY(icent), 2) + pow(grVn2Sys_4Sigma->GetErrorY(icent), 2) );
    double vn4_4Sige    = sqrt( pow(grVn4_4Sigma->GetErrorY(icent), 2) + pow(grVn4Sys_4Sigma->GetErrorY(icent), 2) );
    double vn6_4Sige    = sqrt( pow(grVn6_4Sigma->GetErrorY(icent), 2) + pow(grVn6Sys_4Sigma->GetErrorY(icent), 2) );
    double vn8_4Sige    = sqrt( pow(grVn8_4Sigma->GetErrorY(icent), 2) + pow(grVn8Sys_4Sigma->GetErrorY(icent), 2) );
    double g1e_4Sige    = sqrt( pow(grG1e_4Sigma->GetErrorY(icent), 2) + pow(grG1eSys_4Sigma->GetErrorY(icent), 2) );
    double vn6vn4_4Sige = sqrt( pow(grVn6Vn4_4Sigma->GetErrorY(icent), 2) + pow(grVn6Vn4Sys_4Sigma->GetErrorY(icent), 2) );
    double vn8vn4_4Sige = sqrt( pow(grVn8Vn4_4Sigma->GetErrorY(icent), 2) + pow(grVn8Vn4Sys_4Sigma->GetErrorY(icent), 2) );
    double vn8vn6_4Sige = sqrt( pow(grVn8Vn6_4Sigma->GetErrorY(icent), 2) + pow(grVn8Vn6Sys_4Sigma->GetErrorY(icent), 2) );

    if(vn2_4Sig > 0){
      vn2_4Sigma.push_back(vn2_4Sig);
      vn2Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig > 0){
      vn4_4Sigma.push_back(vn4_4Sig);
      vn4Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4Sig > 0){
      vn6_4Sigma.push_back(vn6_4Sig);
      vn6Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4Sig > 0){
      vn8_4Sigma.push_back(vn8_4Sig);
      vn8Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4Sig != -10000){
      g1e_4Sigma.push_back(g1e_4Sig);
      g1eCent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6vn4_4Sig > 0){
      vn6vn4_4Sigma.push_back( vn6vn4_4Sig );
      vn6vn4Cent_4Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_4Sigma.push_back(Npart[icent]);
    }
    if(vn8vn4_4Sig > 0){
      vn8vn4_4Sigma.push_back( vn8vn4_4Sig );
      vn8vn4Cent_4Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_4Sigma.push_back(Npart[icent]);
    }
    if(vn8vn6_4Sig > 0){
      vn8vn6_4Sigma.push_back( vn8vn6_4Sig );
      vn8vn6Cent_4Sigma.push_back(centBinCenter[icent]);
    }

    //-- 4p2Sigma
    FIXUNFOLD( hFinalUnfold_Cut4p2Sigma[icent], 4.2 );
    EbyECumu cumu4p2Sig( hFinalUnfold_Cut4p2Sigma[icent] );
    double vn2_4p2Sig    = cumu4p2Sig.GetCumu_vn2();
    double vn4_4p2Sig    = cumu4p2Sig.GetCumu_vn4();
    double vn6_4p2Sig    = cumu4p2Sig.GetCumu_vn6();
    double vn8_4p2Sig    = cumu4p2Sig.GetCumu_vn8();
    double g1e_4p2Sig    = cumu4p2Sig.GetGamma1Exp();
    double vn6vn4_4p2Sig;
    double vn8vn4_4p2Sig;
    double vn8vn6_4p2Sig;

    if(vn2_4p2Sig != 0){
      vn2_4p2Sigma.push_back(vn2_4p2Sig);
      vn2Cent_4p2Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p2Sig != 0){
      vn4_4p2Sigma.push_back(vn4_4p2Sig);
      vn4Cent_4p2Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p2Sig != 0){
      vn6_4p2Sigma.push_back(vn6_4p2Sig);
      vn6Cent_4p2Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4p2Sig != 0){
      vn8_4p2Sigma.push_back(vn8_4p2Sig);
      vn8Cent_4p2Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4p2Sig != -10000){
      g1e_4p2Sigma.push_back(g1e_4p2Sig);
      g1eCent_4p2Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p2Sig != 0 && vn6_4p2Sig != 0){
      vn6vn4_4p2Sigma.push_back( vn6_4p2Sig/vn4_4p2Sig );
      vn6vn4Cent_4p2Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_4p2Sigma.push_back(Npart[icent]);
    }
    if(vn4_4p2Sig != 0 && vn8_4p2Sig!= 0){
      vn8vn4_4p2Sigma.push_back( vn8_4p2Sig/vn4_4p2Sig );
      vn8vn4Cent_4p2Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_4p2Sigma.push_back(Npart[icent]);
    }
    if(vn6_4p2Sig != 0 && vn8_4p2Sig!= 0){
      vn8vn6_4p2Sigma.push_back( vn8_4p2Sig/vn6_4p2Sig );
      vn8vn6Cent_4p2Sigma.push_back(centBinCenter[icent]);
    }

    //-- 4p5 sigma
    FIXUNFOLD( hFinalUnfold_Cut4p5Sigma[icent], 4.5 );
    EbyECumu cumu4p5Sig( hFinalUnfold_Cut4p5Sigma[icent] );
    double vn2_4p5Sig    = cumu4p5Sig.GetCumu_vn2();
    double vn4_4p5Sig    = cumu4p5Sig.GetCumu_vn4();
    double vn6_4p5Sig    = cumu4p5Sig.GetCumu_vn6();
    double vn8_4p5Sig    = cumu4p5Sig.GetCumu_vn8();
    double g1e_4p5Sig    = cumu4p5Sig.GetGamma1Exp();
    double vn6vn4_4p5Sig;
    double vn8vn4_4p5Sig;
    double vn8vn6_4p5Sig;

    if(vn2_4p5Sig != 0){
      vn2_4p5Sigma.push_back(vn2_4p5Sig);
      vn2Cent_4p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p5Sig != 0){
      vn4_4p5Sigma.push_back(vn4_4p5Sig);
      vn4Cent_4p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p5Sig != 0){
      vn6_4p5Sigma.push_back(vn6_4p5Sig);
      vn6Cent_4p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4p5Sig != 0){
      vn8_4p5Sigma.push_back(vn8_4p5Sig);
      vn8Cent_4p5Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4p5Sig != -10000){
      g1e_4p5Sigma.push_back(g1e_4p5Sig);
      g1eCent_4p5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p5Sig != 0 && vn6_4p5Sig != 0){
      vn6vn4_4p5Sigma.push_back( vn6_4p5Sig/vn4_4p5Sig );
      vn6vn4Cent_4p5Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_4p5Sigma.push_back(Npart[icent]);
    }
    if(vn4_4p5Sig != 0 && vn8_4p5Sig!= 0){
      vn8vn4_4p5Sigma.push_back( vn8_4p5Sig/vn4_4p5Sig );
      vn8vn4Cent_4p5Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_4p5Sigma.push_back(Npart[icent]);
    }
    if(vn6_4p5Sig != 0 && vn8_4p5Sig!= 0){
      vn8vn6_4p5Sigma.push_back( vn8_4p5Sig/vn6_4p5Sig );
      vn8vn6Cent_4p5Sigma.push_back(centBinCenter[icent]);
    }

    //-- 3p5sigma / default
    if(vn2_3p5Sig > 0 && vn2_4Sig > 0){
      vn2_3p5Sigma_RatioTo4Sigma.push_back(vn2_3p5Sig / vn2_4Sig);
      vn2_3p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn2_3p5Sig*vn2_4Sige/vn2_4Sig/vn2_4Sig) );
      vn2Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p5Sig > 0 && vn4_4Sig > 0){
      vn4_3p5Sigma_RatioTo4Sigma.push_back(vn4_3p5Sig / vn4_4Sig);
      vn4_3p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn4_3p5Sig*vn4_4Sige/vn4_4Sig/vn4_4Sig) );
      vn4Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p5Sig > 0 && vn6_4Sig > 0){
      vn6_3p5Sigma_RatioTo4Sigma.push_back(vn6_3p5Sig / vn6_4Sig);
      vn6_3p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn6_3p5Sig*vn6_4Sige/vn6_4Sig/vn6_4Sig) );
      vn6Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_3p5Sig > 0 && vn8_4Sig > 0){
      vn8_3p5Sigma_RatioTo4Sigma.push_back(vn8_3p5Sig / vn8_4Sig);
      vn8_3p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn8_3p5Sig*vn8_4Sige/vn8_4Sig/vn8_4Sig) );
      vn8Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_3p5Sig != -10000 && g1e_4Sig != -10000){
      g1e_3p5Sigma_RatioTo4Sigma.push_back(g1e_3p5Sig / g1e_4Sig);
      g1e_3p5Sigma_RatioTo4Sigma_err.push_back( fabs(g1e_3p5Sig*g1e_4Sige/g1e_4Sig/g1e_4Sig) );
      g1eCent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p5Sig > 0 && vn6_3p5Sig > 0 && vn6vn4_4Sig > 0){
      vn6vn4_3p5Sigma_RatioTo4Sigma.push_back( (vn6_3p5Sig/vn4_3p5Sig) / (vn6_4Sig/vn4_4Sig) );
      vn6vn4_3p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn6_3p5Sig/vn4_3p5Sig)*vn6vn4_4Sige/vn6vn4_4Sig/vn6vn4_4Sig) );
      vn6vn4Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p5Sig > 0 && vn8_3p5Sig > 0 && vn8vn4_4Sig > 0){
      vn8vn4_3p5Sigma_RatioTo4Sigma.push_back( (vn8_3p5Sig/vn4_3p5Sig) / (vn8_4Sig/vn4_4Sig) );
      vn8vn4_3p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_3p5Sig/vn4_3p5Sig)*vn8vn4_4Sige/vn8vn4_4Sig/vn8vn4_4Sig) );
      vn8vn4Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p5Sig > 0 && vn8_3p5Sig > 0 && vn8vn6_4Sig > 0){
      vn8vn6_3p5Sigma_RatioTo4Sigma.push_back( (vn8_3p5Sig/vn6_3p5Sig) / (vn8_4Sig/vn6_4Sig) );
      vn8vn6_3p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_3p5Sig/vn6_3p5Sig)*vn8vn6_4Sige/vn8vn6_4Sig/vn8vn6_4Sig) );
      vn8vn6Cent_3p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }

    //-- 3p8Sigma/Default
    if(vn2_3p8Sig > 0 && vn2_4Sig > 0){
      vn2_3p8Sigma_RatioTo4Sigma.push_back(vn2_3p8Sig / vn2_4Sig);
      vn2_3p8Sigma_RatioTo4Sigma_err.push_back( fabs(vn2_3p8Sig*vn2_4Sige/vn2_4Sig/vn2_4Sig) );
      vn2Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p8Sig > 0 && vn4_4Sig > 0){
      vn4_3p8Sigma_RatioTo4Sigma.push_back(vn4_3p8Sig / vn4_4Sig);
      vn4_3p8Sigma_RatioTo4Sigma_err.push_back( fabs(vn4_3p8Sig*vn4_4Sige/vn4_4Sig/vn4_4Sig) );
      vn4Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p8Sig > 0 && vn6_4Sig > 0){
      vn6_3p8Sigma_RatioTo4Sigma.push_back(vn6_3p8Sig / vn6_4Sig);
      vn6_3p8Sigma_RatioTo4Sigma_err.push_back( fabs(vn6_3p8Sig*vn6_4Sige/vn6_4Sig/vn6_4Sig) );
      vn6Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_3p8Sig > 0 && vn8_4Sig > 0){
      vn8_3p8Sigma_RatioTo4Sigma.push_back(vn8_3p8Sig / vn8_4Sig);
      vn8_3p8Sigma_RatioTo4Sigma_err.push_back( fabs(vn8_3p8Sig*vn8_4Sige/vn8_4Sig/vn8_4Sig) );
      vn8Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_3p8Sig != -10000 && g1e_4Sig != -10000){
      g1e_3p8Sigma_RatioTo4Sigma.push_back(g1e_3p8Sig / g1e_4Sig);
      g1e_3p8Sigma_RatioTo4Sigma_err.push_back( fabs(g1e_3p8Sig*g1e_4Sige/g1e_4Sig/g1e_4Sig) );
      g1eCent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p8Sig > 0 && vn6_3p8Sig > 0 && vn6vn4_4Sig > 0){
      vn6vn4_3p8Sigma_RatioTo4Sigma.push_back( (vn6_3p8Sig/vn4_3p8Sig) / (vn6_4Sig/vn4_4Sig) );
      vn6vn4_3p8Sigma_RatioTo4Sigma_err.push_back( fabs((vn6_3p8Sig/vn4_3p8Sig)*vn6vn4_4Sige/vn6vn4_4Sig/vn6vn4_4Sig) );
      vn6vn4Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3p8Sig > 0 && vn8_3p8Sig > 0 && vn8vn4_4Sig > 0){
      vn8vn4_3p8Sigma_RatioTo4Sigma.push_back( (vn8_3p8Sig/vn4_3p8Sig) / (vn8_4Sig/vn4_4Sig) );
      vn8vn4_3p8Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_3p8Sig/vn4_3p8Sig)*vn8vn4_4Sige/vn8vn4_4Sig/vn8vn4_4Sig) );
      vn8vn4Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3p8Sig > 0 && vn8_3p8Sig > 0 && vn8vn6_4Sig > 0){
      vn8vn6_3p8Sigma_RatioTo4Sigma.push_back( (vn8_3p8Sig/vn6_3p8Sig) / (vn8_4Sig/vn6_4Sig) );
      vn8vn6_3p8Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_3p8Sig/vn6_3p8Sig)*vn8vn6_4Sige/vn8vn6_4Sig/vn8vn6_4Sig) );
      vn8vn6Cent_3p8Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }


    //-- 4p2Sigma/Default
    if(vn2_4p2Sig > 0 && vn2_4Sig > 0){
      vn2_4p2Sigma_RatioTo4Sigma.push_back(vn2_4p2Sig / vn2_4Sig);
      vn2_4p2Sigma_RatioTo4Sigma_err.push_back( fabs(vn2_4p2Sig*vn2_4Sige/vn2_4Sig/vn2_4Sig) );
      vn2Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p2Sig > 0 && vn4_4Sig > 0){
      vn4_4p2Sigma_RatioTo4Sigma.push_back(vn4_4p2Sig / vn4_4Sig);
      vn4_4p2Sigma_RatioTo4Sigma_err.push_back( fabs(vn4_4p2Sig*vn4_4Sige/vn4_4Sig/vn4_4Sig) );
      vn4Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p2Sig > 0 && vn6_4Sig > 0){
      vn6_4p2Sigma_RatioTo4Sigma.push_back(vn6_4p2Sig / vn6_4Sig);
      vn6_4p2Sigma_RatioTo4Sigma_err.push_back( fabs(vn6_4p2Sig*vn6_4Sige/vn6_4Sig/vn6_4Sig) );
      vn6Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4p2Sig > 0 && vn8_4Sig > 0){
      vn8_4p2Sigma_RatioTo4Sigma.push_back(vn8_4p2Sig / vn8_4Sig);
      vn8_4p2Sigma_RatioTo4Sigma_err.push_back( fabs(vn8_4p2Sig*vn8_4Sige/vn8_4Sig/vn8_4Sig) );
      vn8Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4p2Sig != -10000 && g1e_4Sig != -10000){
      g1e_4p2Sigma_RatioTo4Sigma.push_back(g1e_4p2Sig / g1e_4Sig);
      g1e_4p2Sigma_RatioTo4Sigma_err.push_back( fabs(g1e_4p2Sig*g1e_4Sige/g1e_4Sig/g1e_4Sig) );
      g1eCent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p2Sig > 0 && vn6_4p2Sig > 0 && vn6vn4_4Sig > 0){
      vn6vn4_4p2Sigma_RatioTo4Sigma.push_back( (vn6_4p2Sig/vn4_4p2Sig) / (vn6_4Sig/vn4_4Sig) );
      vn6vn4_4p2Sigma_RatioTo4Sigma_err.push_back( fabs((vn6_4p2Sig/vn4_4p2Sig)*vn6vn4_4Sige/vn6vn4_4Sig/vn6vn4_4Sig) );
      vn6vn4Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p2Sig > 0 && vn8_4p2Sig > 0 && vn8vn4_4Sig > 0){
      vn8vn4_4p2Sigma_RatioTo4Sigma.push_back( (vn8_4p2Sig/vn4_4p2Sig) / (vn8_4Sig/vn4_4Sig) );
      vn8vn4_4p2Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_4p2Sig/vn4_4p2Sig)*vn8vn4_4Sige/vn8vn4_4Sig/vn8vn4_4Sig) );
      vn8vn4Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p2Sig > 0 && vn8_4p2Sig > 0 && vn8vn6_4Sig > 0){
      vn8vn6_4p2Sigma_RatioTo4Sigma.push_back( (vn8_4p2Sig/vn6_4p2Sig) / (vn8_4Sig/vn6_4Sig) );
      vn8vn6_4p2Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_4p2Sig/vn6_4p2Sig)*vn8vn6_4Sige/vn8vn6_4Sig/vn8vn6_4Sig) );
      vn8vn6Cent_4p2Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }

    //-- 4p5 sigma / Default
    if(vn2_4p5Sig > 0 && vn2_4Sig > 0){
      vn2_4p5Sigma_RatioTo4Sigma.push_back(vn2_4p5Sig / vn2_4Sig);
      vn2_4p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn2_4p5Sig*vn2_4Sige/vn2_4Sig/vn2_4Sig) );
      vn2Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p5Sig > 0 && vn4_4Sig > 0){
      vn4_4p5Sigma_RatioTo4Sigma.push_back(vn4_4p5Sig / vn4_4Sig);
      vn4_4p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn4_4p5Sig*vn4_4Sige/vn4_4Sig/vn4_4Sig) );
      vn4Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p5Sig > 0 && vn6_4Sig > 0){
      vn6_4p5Sigma_RatioTo4Sigma.push_back(vn6_4p5Sig / vn6_4Sig);
      vn6_4p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn6_4p5Sig*vn6_4Sige/vn6_4Sig/vn6_4Sig) );
      vn6Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4p5Sig > 0 && vn8_4Sig > 0){
      vn8_4p5Sigma_RatioTo4Sigma.push_back(vn8_4p5Sig / vn8_4Sig);
      vn8_4p5Sigma_RatioTo4Sigma_err.push_back( fabs(vn8_4p5Sig*vn8_4Sige/vn8_4Sig/vn8_4Sig) );
      vn8Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4p5Sig != -10000 && g1e_4Sig != -10000){
      g1e_4p5Sigma_RatioTo4Sigma.push_back(g1e_4p5Sig / g1e_4Sig);
      g1e_4p5Sigma_RatioTo4Sigma_err.push_back( fabs(g1e_4p5Sig*g1e_4Sige/g1e_4Sig/g1e_4Sig) );
      g1eCent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p5Sig > 0 && vn6_4p5Sig > 0 && vn6vn4_4Sig > 0){
      vn6vn4_4p5Sigma_RatioTo4Sigma.push_back( (vn6_4p5Sig/vn4_4p5Sig) / (vn6_4Sig/vn4_4Sig) );
      vn6vn4_4p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn6_4p5Sig/vn4_4p5Sig)*vn6vn4_4Sige/vn6vn4_4Sig/vn6vn4_4Sig) );
      vn6vn4Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4p5Sig > 0 && vn8_4p5Sig > 0 && vn8vn4_4Sig > 0){
      vn8vn4_4p5Sigma_RatioTo4Sigma.push_back( (vn8_4p5Sig/vn4_4p5Sig) / (vn8_4Sig/vn4_4Sig) );
      vn8vn4_4p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_4p5Sig/vn4_4p5Sig)*vn8vn4_4Sige/vn8vn4_4Sig/vn8vn4_4Sig) );
      vn8vn4Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4p5Sig > 0 && vn8_4p5Sig > 0 && vn8vn6_4Sig > 0){
      vn8vn6_4p5Sigma_RatioTo4Sigma.push_back( (vn8_4p5Sig/vn6_4p5Sig) / (vn8_4Sig/vn6_4Sig) );
      vn8vn6_4p5Sigma_RatioTo4Sigma_err.push_back( fabs((vn8_4p5Sig/vn6_4p5Sig)*vn8vn6_4Sige/vn8vn6_4Sig/vn8vn6_4Sig) );
      vn8vn6Cent_4p5Sigma_RatioTo4Sigma.push_back(centBinCenter[icent]);
    }

    //-- Line Markers
    hFinalUnfold_Default[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hFinalUnfold_Default[icent]->GetYaxis()->SetTitle( "Events" );
    hFinalUnfold_Default[icent]->SetMinimum(1e-10*hFinalUnfold_Default[icent]->GetMaximum());
    hFinalUnfold_Default[icent]->SetMaximum(10.*hFinalUnfold_Default[icent]->GetMaximum());
    hFinalUnfold_Default[icent]->GetXaxis()->SetRange(1, hFinalUnfold_Default[icent]->FindBin(0.5));
    double min  = hFinalUnfold_Default[icent]->GetMinimum();
    double max  = hFinalUnfold_Default[icent]->GetMaximum();
    double s3p5 = hFinalUnfold_Default[icent]->GetMean() + 3.5*hFinalUnfold_Default[icent]->GetRMS();
    double s3p8 = hFinalUnfold_Default[icent]->GetMean() + 3.8*hFinalUnfold_Default[icent]->GetRMS();
    double s4   = hFinalUnfold_Default[icent]->GetMean() + 4.*hFinalUnfold_Default[icent]->GetRMS();
    double s4p2 = hFinalUnfold_Default[icent]->GetMean() + 4.2*hFinalUnfold_Default[icent]->GetRMS();
    double s4p5 = hFinalUnfold_Default[icent]->GetMean() + 4.5*hFinalUnfold_Default[icent]->GetRMS();

    l3p5s[icent] = new TLine(s3p5, min, s3p5, max);
    l3p8s[icent] = new TLine(s3p8, min, s3p8, max);
    l4s[icent]   = new TLine(s4,   min, s4,   max);
    l4p2s[icent] = new TLine(s4p2, min, s4p2, max);
    l4p5s[icent] = new TLine(s4p5, min, s4p5, max);
    /*
    l3p5s[icent]->SetLineWidth(2);
    l3p8s[icent]->SetLineWidth(2);
    l4s[icent]->SetLineWidth(2);
    l4p2s[icent]->SetLineWidth(2);
    l4p5s[icent]->SetLineWidth(2);
    */
    l3p5s[icent]->SetLineColor(2);
    l3p8s[icent]->SetLineColor(4);
    l4s[icent]->SetLineColor(6);
    l4p2s[icent]->SetLineColor(8);
    l4p5s[icent]->SetLineColor(9);

  } //-- End Cent loop

  //-- Make sweet, sweet TGraphs

  //-- 3p5Sigma
  grVn2_3p5Sigma    = new TGraph(vn2Cent_3p5Sigma.size(), &(vn2Cent_3p5Sigma[0]), &(vn2_3p5Sigma[0]));
  grVn4_3p5Sigma    = new TGraph(vn4Cent_3p5Sigma.size(), &(vn4Cent_3p5Sigma[0]), &(vn4_3p5Sigma[0]));
  grVn6_3p5Sigma    = new TGraph(vn6Cent_3p5Sigma.size(), &(vn6Cent_3p5Sigma[0]), &(vn6_3p5Sigma[0]));
  grVn8_3p5Sigma    = new TGraph(vn8Cent_3p5Sigma.size(), &(vn8Cent_3p5Sigma[0]), &(vn8_3p5Sigma[0]));
  grG1e_3p5Sigma    = new TGraph(g1eCent_3p5Sigma.size(), &(g1eCent_3p5Sigma[0]), &(g1e_3p5Sigma[0]));
  grVn6Vn4_3p5Sigma = new TGraph(vn6vn4Cent_3p5Sigma.size(), &(vn6vn4Cent_3p5Sigma[0]), &(vn6vn4_3p5Sigma[0]));
  grVn8Vn4_3p5Sigma = new TGraph(vn8vn4Cent_3p5Sigma.size(), &(vn8vn4Cent_3p5Sigma[0]), &(vn8vn4_3p5Sigma[0]));
  grVn8Vn6_3p5Sigma = new TGraph(vn8vn6Cent_3p5Sigma.size(), &(vn8vn6Cent_3p5Sigma[0]), &(vn8vn6_3p5Sigma[0]));

  formatGraph(grVn2_3p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         22, "grVn2_3p5Sigma");
  formatGraph(grVn4_3p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 22, "grVn4_3p5Sigma");
  formatGraph(grVn6_3p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         22, "grVn6_3p5Sigma");
  formatGraph(grVn8_3p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 22, "grVn8_3p5Sigma");
  formatGraph(grG1e_3p5Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         22, "grG1e_3p5Sigma");
  formatGraph(grVn6Vn4_3p5Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         22, "grVn6Vn4_3p5Sigma");
  formatGraph(grVn8Vn4_3p5Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  22, "grVn8Vn4_3p5Sigma");
  formatGraph(grVn8Vn6_3p5Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 22, "grVn8Vn6_3p5Sigma");

  //-- 3p5Sigma/Default
  grVn2_3p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn2Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn2Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn2_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn2_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn4_3p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn4Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn4Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn4_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn4_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn6_3p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn6Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn6Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn6_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn6_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn8_3p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn8Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn8Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn8_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8_3p5Sigma_RatioTo4Sigma_err[0]));
  grG1e_3p5Sigma_RatioTo4Sigma    = new TGraphErrors(g1eCent_3p5Sigma_RatioTo4Sigma.size(), &(g1eCent_3p5Sigma_RatioTo4Sigma[0]), &(g1e_3p5Sigma_RatioTo4Sigma[0]), CERR, &(g1e_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn6Vn4_3p5Sigma_RatioTo4Sigma = new TGraphErrors(vn6vn4Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn6vn4Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn6vn4_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn6vn4_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn4_3p5Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn4Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn8vn4Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn8vn4_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn4_3p5Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn6_3p5Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn6Cent_3p5Sigma_RatioTo4Sigma.size(), &(vn8vn6Cent_3p5Sigma_RatioTo4Sigma[0]), &(vn8vn6_3p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn6_3p5Sigma_RatioTo4Sigma_err[0]));

  formatGraph(grVn2_3p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         22, "grVn2_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn4_3p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 22, "grVn4_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn6_3p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         22, "grVn6_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8_3p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 22, "grVn8_3p5Sigma_RatioTo4Sigma");
  formatGraph(grG1e_3p5Sigma_RatioTo4Sigma,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         22, "grG1e_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn6Vn4_3p5Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         22, "grVn6Vn4_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn4_3p5Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  22, "grVn8Vn4_3p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn6_3p5Sigma_RatioTo4Sigma, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 22, "grVn8Vn6_3p5Sigma_RatioTo4Sigma");

  //-- 3p8Sigma
  grVn2_3p8Sigma    = new TGraph(vn2Cent_3p8Sigma.size(), &(vn2Cent_3p8Sigma[0]), &(vn2_3p8Sigma[0]));
  grVn4_3p8Sigma    = new TGraph(vn4Cent_3p8Sigma.size(), &(vn4Cent_3p8Sigma[0]), &(vn4_3p8Sigma[0]));
  grVn6_3p8Sigma    = new TGraph(vn6Cent_3p8Sigma.size(), &(vn6Cent_3p8Sigma[0]), &(vn6_3p8Sigma[0]));
  grVn8_3p8Sigma    = new TGraph(vn8Cent_3p8Sigma.size(), &(vn8Cent_3p8Sigma[0]), &(vn8_3p8Sigma[0]));
  grG1e_3p8Sigma    = new TGraph(g1eCent_3p8Sigma.size(), &(g1eCent_3p8Sigma[0]), &(g1e_3p8Sigma[0]));
  grVn6Vn4_3p8Sigma = new TGraph(vn6vn4Cent_3p8Sigma.size(), &(vn6vn4Cent_3p8Sigma[0]), &(vn6vn4_3p8Sigma[0]));
  grVn8Vn4_3p8Sigma = new TGraph(vn8vn4Cent_3p8Sigma.size(), &(vn8vn4Cent_3p8Sigma[0]), &(vn8vn4_3p8Sigma[0]));
  grVn8Vn6_3p8Sigma = new TGraph(vn8vn6Cent_3p8Sigma.size(), &(vn8vn6Cent_3p8Sigma[0]), &(vn8vn6_3p8Sigma[0]));

  formatGraph(grVn2_3p8Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         25, "grVn2_3p8Sigma");
  formatGraph(grVn4_3p8Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 25, "grVn4_3p8Sigma");
  formatGraph(grVn6_3p8Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         25, "grVn6_3p8Sigma");
  formatGraph(grVn8_3p8Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 25, "grVn8_3p8Sigma");
  formatGraph(grG1e_3p8Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         25, "grG1e_3p8Sigma");
  formatGraph(grVn6Vn4_3p8Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         25, "grVn6Vn4_3p8Sigma");
  formatGraph(grVn8Vn4_3p8Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  25, "grVn8Vn4_3p8Sigma");
  formatGraph(grVn8Vn6_3p8Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 25, "grVn8Vn6_3p8Sigma");

  //-- 3p8Sigma/Default
  grVn2_3p8Sigma_RatioTo4Sigma    = new TGraphErrors(vn2Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn2Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn2_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn2_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn4_3p8Sigma_RatioTo4Sigma    = new TGraphErrors(vn4Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn4Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn4_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn4_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn6_3p8Sigma_RatioTo4Sigma    = new TGraphErrors(vn6Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn6Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn6_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn6_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn8_3p8Sigma_RatioTo4Sigma    = new TGraphErrors(vn8Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn8Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn8_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn8_3p8Sigma_RatioTo4Sigma_err[0]));
  grG1e_3p8Sigma_RatioTo4Sigma    = new TGraphErrors(g1eCent_3p8Sigma_RatioTo4Sigma.size(), &(g1eCent_3p8Sigma_RatioTo4Sigma[0]), &(g1e_3p8Sigma_RatioTo4Sigma[0]), CERR, &(g1e_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn6Vn4_3p8Sigma_RatioTo4Sigma = new TGraphErrors(vn6vn4Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn6vn4Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn6vn4_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn6vn4_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn4_3p8Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn4Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn8vn4Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn8vn4_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn4_3p8Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn6_3p8Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn6Cent_3p8Sigma_RatioTo4Sigma.size(), &(vn8vn6Cent_3p8Sigma_RatioTo4Sigma[0]), &(vn8vn6_3p8Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn6_3p8Sigma_RatioTo4Sigma_err[0]));

  formatGraph(grVn2_3p8Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         25, "grVn2_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn4_3p8Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 25, "grVn4_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn6_3p8Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         25, "grVn6_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn8_3p8Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 25, "grVn8_3p8Sigma_RatioTo4Sigma");
  formatGraph(grG1e_3p8Sigma_RatioTo4Sigma,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         25, "grG1e_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn6Vn4_3p8Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         25, "grVn6Vn4_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn4_3p8Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  25, "grVn8Vn4_3p8Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn6_3p8Sigma_RatioTo4Sigma, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 25, "grVn8Vn6_3p8Sigma_RatioTo4Sigma");

  //-- 4Sigma
  /*
  grVn2_4Sigma    = new TGraph(vn2Cent_4Sigma.size(), &(vn2Cent_4Sigma[0]), &(vn2_4Sigma[0]));
  grVn4_4Sigma    = new TGraph(vn4Cent_4Sigma.size(), &(vn4Cent_4Sigma[0]), &(vn4_4Sigma[0]));
  grVn6_4Sigma    = new TGraph(vn6Cent_4Sigma.size(), &(vn6Cent_4Sigma[0]), &(vn6_4Sigma[0]));
  grVn8_4Sigma    = new TGraph(vn8Cent_4Sigma.size(), &(vn8Cent_4Sigma[0]), &(vn8_4Sigma[0]));
  grG1e_4Sigma    = new TGraph(g1eCent_4Sigma.size(), &(g1eCent_4Sigma[0]), &(g1e_4Sigma[0]));
  grVn6Vn4_4Sigma = new TGraph(vn6vn4Cent_4Sigma.size(), &(vn6vn4Cent_4Sigma[0]), &(vn6vn4_4Sigma[0]));
  grVn8Vn4_4Sigma = new TGraph(vn8vn4Cent_4Sigma.size(), &(vn8vn4Cent_4Sigma[0]), &(vn8vn4_4Sigma[0]));
  grVn8Vn6_4Sigma = new TGraph(vn8vn6Cent_4Sigma.size(), &(vn8vn6Cent_4Sigma[0]), &(vn8vn6_4Sigma[0]));
  */
  formatGraph(grVn2_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         20, "grVn2_4Sigma");
  formatGraph(grVn4_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 20, "grVn4_4Sigma");
  formatGraph(grVn6_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         20, "grVn6_4Sigma");
  formatGraph(grVn8_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 20, "grVn8_4Sigma");
  formatGraph(grG1e_4Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         20, "grG1e_4Sigma");
  formatGraph(grVn6Vn4_4Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         20, "grVn6Vn4_4Sigma");
  formatGraph(grVn8Vn4_4Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  20, "grVn8Vn4_4Sigma");
  formatGraph(grVn8Vn6_4Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 20, "grVn8Vn6_4Sigma");

  //-- 4p2Sigma
  grVn2_4p2Sigma    = new TGraph(vn2Cent_4p2Sigma.size(), &(vn2Cent_4p2Sigma[0]), &(vn2_4p2Sigma[0]));
  grVn4_4p2Sigma    = new TGraph(vn4Cent_4p2Sigma.size(), &(vn4Cent_4p2Sigma[0]), &(vn4_4p2Sigma[0]));
  grVn6_4p2Sigma    = new TGraph(vn6Cent_4p2Sigma.size(), &(vn6Cent_4p2Sigma[0]), &(vn6_4p2Sigma[0]));
  grVn8_4p2Sigma    = new TGraph(vn8Cent_4p2Sigma.size(), &(vn8Cent_4p2Sigma[0]), &(vn8_4p2Sigma[0]));
  grG1e_4p2Sigma    = new TGraph(g1eCent_4p2Sigma.size(), &(g1eCent_4p2Sigma[0]), &(g1e_4p2Sigma[0]));
  grVn6Vn4_4p2Sigma = new TGraph(vn6vn4Cent_4p2Sigma.size(), &(vn6vn4Cent_4p2Sigma[0]), &(vn6vn4_4p2Sigma[0]));
  grVn8Vn4_4p2Sigma = new TGraph(vn8vn4Cent_4p2Sigma.size(), &(vn8vn4Cent_4p2Sigma[0]), &(vn8vn4_4p2Sigma[0]));
  grVn8Vn6_4p2Sigma = new TGraph(vn8vn6Cent_4p2Sigma.size(), &(vn8vn6Cent_4p2Sigma[0]), &(vn8vn6_4p2Sigma[0]));

  formatGraph(grVn2_4p2Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         28, "grVn2_4p2Sigma");
  formatGraph(grVn4_4p2Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 28, "grVn4_4p2Sigma");
  formatGraph(grVn6_4p2Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         28, "grVn6_4p2Sigma");
  formatGraph(grVn8_4p2Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 28, "grVn8_4p2Sigma");
  formatGraph(grG1e_4p2Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         28, "grG1e_4p2Sigma");
  formatGraph(grVn6Vn4_4p2Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         28, "grVn6Vn4_4p2Sigma");
  formatGraph(grVn8Vn4_4p2Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  28, "grVn8Vn4_4p2Sigma");
  formatGraph(grVn8Vn6_4p2Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 28, "grVn8Vn6_4p2Sigma");

  //-- 4p2Sigma/Default
  grVn2_4p2Sigma_RatioTo4Sigma    = new TGraphErrors(vn2Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn2Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn2_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn2_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn4_4p2Sigma_RatioTo4Sigma    = new TGraphErrors(vn4Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn4Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn4_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn4_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn6_4p2Sigma_RatioTo4Sigma    = new TGraphErrors(vn6Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn6Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn6_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn6_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn8_4p2Sigma_RatioTo4Sigma    = new TGraphErrors(vn8Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn8Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn8_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn8_4p2Sigma_RatioTo4Sigma_err[0]));
  grG1e_4p2Sigma_RatioTo4Sigma    = new TGraphErrors(g1eCent_4p2Sigma_RatioTo4Sigma.size(), &(g1eCent_4p2Sigma_RatioTo4Sigma[0]), &(g1e_4p2Sigma_RatioTo4Sigma[0]), CERR, &(g1e_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn6Vn4_4p2Sigma_RatioTo4Sigma = new TGraphErrors(vn6vn4Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn6vn4Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn6vn4_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn6vn4_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn4_4p2Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn4Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn8vn4Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn8vn4_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn4_4p2Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn6_4p2Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn6Cent_4p2Sigma_RatioTo4Sigma.size(), &(vn8vn6Cent_4p2Sigma_RatioTo4Sigma[0]), &(vn8vn6_4p2Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn6_4p2Sigma_RatioTo4Sigma_err[0]));

  formatGraph(grVn2_4p2Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         28, "grVn2_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn4_4p2Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 28, "grVn4_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn6_4p2Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         28, "grVn6_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn8_4p2Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 28, "grVn8_4p2Sigma_RatioTo4Sigma");
  formatGraph(grG1e_4p2Sigma_RatioTo4Sigma,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         28, "grG1e_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn6Vn4_4p2Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         28, "grVn6Vn4_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn4_4p2Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  28, "grVn8Vn4_4p2Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn6_4p2Sigma_RatioTo4Sigma, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 28, "grVn8Vn6_4p2Sigma_RatioTo4Sigma");

  //-- 4p5Sigma
  grVn2_4p5Sigma    = new TGraph(vn2Cent_4p5Sigma.size(), &(vn2Cent_4p5Sigma[0]), &(vn2_4p5Sigma[0]));
  grVn4_4p5Sigma    = new TGraph(vn4Cent_4p5Sigma.size(), &(vn4Cent_4p5Sigma[0]), &(vn4_4p5Sigma[0]));
  grVn6_4p5Sigma    = new TGraph(vn6Cent_4p5Sigma.size(), &(vn6Cent_4p5Sigma[0]), &(vn6_4p5Sigma[0]));
  grVn8_4p5Sigma    = new TGraph(vn8Cent_4p5Sigma.size(), &(vn8Cent_4p5Sigma[0]), &(vn8_4p5Sigma[0]));
  grG1e_4p5Sigma    = new TGraph(g1eCent_4p5Sigma.size(), &(g1eCent_4p5Sigma[0]), &(g1e_4p5Sigma[0]));
  grVn6Vn4_4p5Sigma = new TGraph(vn6vn4Cent_4p5Sigma.size(), &(vn6vn4Cent_4p5Sigma[0]), &(vn6vn4_4p5Sigma[0]));
  grVn8Vn4_4p5Sigma = new TGraph(vn8vn4Cent_4p5Sigma.size(), &(vn8vn4Cent_4p5Sigma[0]), &(vn8vn4_4p5Sigma[0]));
  grVn8Vn6_4p5Sigma = new TGraph(vn8vn6Cent_4p5Sigma.size(), &(vn8vn6Cent_4p5Sigma[0]), &(vn8vn6_4p5Sigma[0]));

  formatGraph(grVn2_4p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         23, "grVn2_4p5Sigma");
  formatGraph(grVn4_4p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 23, "grVn4_4p5Sigma");
  formatGraph(grVn6_4p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         23, "grVn6_4p5Sigma");
  formatGraph(grVn8_4p5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 23, "grVn8_4p5Sigma");
  formatGraph(grG1e_4p5Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         23, "grG1e_4p5Sigma");
  formatGraph(grVn6Vn4_4p5Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         23, "grVn6Vn4_4p5Sigma");
  formatGraph(grVn8Vn4_4p5Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  23, "grVn8Vn4_4p5Sigma");
  formatGraph(grVn8Vn6_4p5Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 23, "grVn8Vn6_4p5Sigma");

  //-- 4p5Sigma/Default
  grVn2_4p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn2Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn2Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn2_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn2_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn4_4p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn4Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn4Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn4_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn4_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn6_4p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn6Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn6Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn6_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn6_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn8_4p5Sigma_RatioTo4Sigma    = new TGraphErrors(vn8Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn8Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn8_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8_4p5Sigma_RatioTo4Sigma_err[0]));
  grG1e_4p5Sigma_RatioTo4Sigma    = new TGraphErrors(g1eCent_4p5Sigma_RatioTo4Sigma.size(), &(g1eCent_4p5Sigma_RatioTo4Sigma[0]), &(g1e_4p5Sigma_RatioTo4Sigma[0]), CERR, &(g1e_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn6Vn4_4p5Sigma_RatioTo4Sigma = new TGraphErrors(vn6vn4Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn6vn4Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn6vn4_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn6vn4_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn4_4p5Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn4Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn8vn4Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn8vn4_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn4_4p5Sigma_RatioTo4Sigma_err[0]));
  grVn8Vn6_4p5Sigma_RatioTo4Sigma = new TGraphErrors(vn8vn6Cent_4p5Sigma_RatioTo4Sigma.size(), &(vn8vn6Cent_4p5Sigma_RatioTo4Sigma[0]), &(vn8vn6_4p5Sigma_RatioTo4Sigma[0]), CERR, &(vn8vn6_4p5Sigma_RatioTo4Sigma_err[0]));

  formatGraph(grVn2_4p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         23, "grVn2_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn4_4p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 23, "grVn4_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn6_4p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         23, "grVn6_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8_4p5Sigma_RatioTo4Sigma,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 23, "grVn8_4p5Sigma_RatioTo4Sigma");
  formatGraph(grG1e_4p5Sigma_RatioTo4Sigma,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         23, "grG1e_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn6Vn4_4p5Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         23, "grVn6Vn4_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn4_4p5Sigma_RatioTo4Sigma, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  23, "grVn8Vn4_4p5Sigma_RatioTo4Sigma");
  formatGraph(grVn8Vn6_4p5Sigma_RatioTo4Sigma, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 23, "grVn8Vn6_4p5Sigma_RatioTo4Sigma");

  TLine * lone_vn2 = new TLine(grVn2_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn2_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn4 = new TLine(grVn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn6 = new TLine(grVn6_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn6_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn8 = new TLine(grVn8_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn8_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_g1e = new TLine(grG1e_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grG1e_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn6vn4 = new TLine(grVn6Vn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn6Vn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn8vn4 = new TLine(grVn8Vn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn8Vn4_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);
  TLine * lone_vn8vn6 = new TLine(grVn8Vn6_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmin(), 1., grVn8Vn6_3p5Sigma_RatioTo4Sigma->GetXaxis()->GetXmax(), 1.);

  lone_vn2->SetLineWidth(2);
  lone_vn4->SetLineWidth(2);
  lone_vn6->SetLineWidth(2);
  lone_vn8->SetLineWidth(2);
  lone_g1e->SetLineWidth(2);
  lone_vn6vn4->SetLineWidth(2);
  lone_vn8vn4->SetLineWidth(2);
  lone_vn8vn6->SetLineWidth(2);

  lone_vn2->SetLineStyle(2);
  lone_vn4->SetLineStyle(2);
  lone_vn6->SetLineStyle(2);
  lone_vn8->SetLineStyle(2);
  lone_g1e->SetLineStyle(2);
  lone_vn6vn4->SetLineStyle(2);
  lone_vn8vn4->SetLineStyle(2);
  lone_vn8vn6->SetLineStyle(2);

  TLegend * legVn2 = new TLegend(0.55, 0.2, 0.9, 0.55);
  legVn2->SetBorderSize(0);
  legVn2->SetFillStyle(0);
  legVn2->AddEntry(grVn2_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn2->AddEntry(grVn2_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn2->AddEntry(grVn2_4Sigma,    "4 #sigma Cut",   "lp");
  legVn2->AddEntry(grVn2_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn2->AddEntry(grVn2_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TLegend * legVn4 = new TLegend(0.55, 0.2, 0.9, 0.55);
  legVn4->SetBorderSize(0);
  legVn4->SetFillStyle(0);
  legVn4->AddEntry(grVn4_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn4->AddEntry(grVn4_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn4->AddEntry(grVn4_4Sigma,    "4 #sigma Cut",   "lp");
  legVn4->AddEntry(grVn4_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn4->AddEntry(grVn4_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TLegend * legVn6 = new TLegend(0.55, 0.2, 0.9, 0.55);
  legVn6->SetBorderSize(0);
  legVn6->SetFillStyle(0);
  legVn6->AddEntry(grVn6_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn6->AddEntry(grVn6_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn6->AddEntry(grVn6_4Sigma,    "4 #sigma Cut",   "lp");
  legVn6->AddEntry(grVn6_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn6->AddEntry(grVn6_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TLegend * legVn8 = new TLegend(0.55, 0.2, 0.9, 0.55);
  legVn8->SetBorderSize(0);
  legVn8->SetFillStyle(0);
  legVn8->AddEntry(grVn8_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn8->AddEntry(grVn8_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn8->AddEntry(grVn8_4Sigma,    "4 #sigma Cut",   "lp");
  legVn8->AddEntry(grVn8_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn8->AddEntry(grVn8_4p2Sigma,  "4.5 #sigma Cut", "lp");

  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 2000, 500);
  cCumu->Divide(4,1);

  cCumu->cd(1);
  grVn2Sys_4Sigma->Draw("apE2");
  grVn2_3p5Sigma->Draw("psame");
  grVn2_3p8Sigma->Draw("psame");
  grVn2_4Sigma->Draw("psame");
  grVn2_4p2Sigma->Draw("psame");
  grVn2_4p5Sigma->Draw("psame");
  legVn2->Draw("same");

  cCumu->cd(2);
  grVn4Sys_4Sigma->Draw("apE2");
  grVn4_3p5Sigma->Draw("psame");
  grVn4_3p8Sigma->Draw("psame");
  grVn4_4Sigma->Draw("psame");
  grVn4_4p2Sigma->Draw("psame");
  grVn4_4p5Sigma->Draw("psame");
  legVn4->Draw("same");

  cCumu->cd(3);
  grVn6Sys_4Sigma->Draw("apE2");
  grVn6_3p5Sigma->Draw("psame");
  grVn6_3p8Sigma->Draw("psame");
  grVn6_4Sigma->Draw("psame");
  grVn6_4p2Sigma->Draw("psame");
  grVn6_4p5Sigma->Draw("psame");
  legVn6->Draw("same");

  cCumu->cd(4);
  grVn8Sys_4Sigma->Draw("apE2");
  grVn8_3p5Sigma->Draw("psame");
  grVn8_3p8Sigma->Draw("psame");
  grVn8_4Sigma->Draw("psame");
  grVn8_4p2Sigma->Draw("psame");
  grVn8_4p5Sigma->Draw("psame");
  legVn8->Draw("same");

  cCumu->SaveAs("cCumu.pdf");

  //-- sysCumu
  TLegend * legVn2Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn2Sys->SetBorderSize(0);
  legVn2Sys->SetFillStyle(0);
  legVn2Sys->AddEntry(grVn2_3p8Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn2Sys->AddEntry(grVn2_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn2Sys->AddEntry(grVn2_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn2Sys->AddEntry(grVn2_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TLegend * legVn4Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn4Sys->SetBorderSize(0);
  legVn4Sys->SetFillStyle(0);
  legVn4Sys->AddEntry(grVn4_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn4Sys->AddEntry(grVn4_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn4Sys->AddEntry(grVn4_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn4Sys->AddEntry(grVn4_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TLegend * legVn6Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn6Sys->SetBorderSize(0);
  legVn6Sys->SetFillStyle(0);
  legVn6Sys->AddEntry(grVn6_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn6Sys->AddEntry(grVn6_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn6Sys->AddEntry(grVn6_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn6Sys->AddEntry(grVn6_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TLegend * legVn8Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn8Sys->SetBorderSize(0);
  legVn8Sys->SetFillStyle(0);
  legVn8Sys->AddEntry(grVn8_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn8Sys->AddEntry(grVn8_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn8Sys->AddEntry(grVn8_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn8Sys->AddEntry(grVn8_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TCanvas * cCumuSys = new TCanvas("cCumuSys", "cCumuSys", 2000, 500);
  cCumuSys->Divide(4,1);

  cCumuSys->cd(1);
  grVn2_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn2_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn2_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn2_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn2Sys->Draw("same");
  lone_vn2->Draw("same");

  cCumuSys->cd(2);
  grVn4_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn4_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn4_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn4_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn4Sys->Draw("same");
  lone_vn4->Draw("same");

  cCumuSys->cd(3);
  grVn6_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn6_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn6_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn6_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn6Sys->Draw("same");
  lone_vn6->Draw("same");

  cCumuSys->cd(4);
  grVn8_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn8_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn8_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn8_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn8Sys->Draw("same");
  lone_vn8->Draw("same");

  cCumuSys->SaveAs("cCumuSys.pdf");


  TLegend * legVn6Vn4 = new TLegend(0.3, 0.2, 0.65, 0.55);
  legVn6Vn4->SetBorderSize(0);
  legVn6Vn4->SetFillStyle(0);
  legVn6Vn4->AddEntry(grVn6Vn4_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_4Sigma,    "4 #sigma Cut",   "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TLegend * legVn8Vn4 = new TLegend(0.3, 0.2, 0.65, 0.55);
  legVn8Vn4->SetBorderSize(0);
  legVn8Vn4->SetFillStyle(0);
  legVn8Vn4->AddEntry(grVn8Vn4_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_4Sigma,    "4 #sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TLegend * legVn8Vn6 = new TLegend(0.2, 0.2, 0.55, 0.55);
  legVn8Vn6->SetBorderSize(0);
  legVn8Vn6->SetFillStyle(0);
  legVn8Vn6->AddEntry(grVn8Vn6_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_4Sigma,    "4 #sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);

  cCumuRatio->cd(1);
  grVn6Vn4Sys_4Sigma->Draw("apE2");
  grVn6Vn4_3p5Sigma->Draw("psame");
  grVn6Vn4_3p8Sigma->Draw("psame");
  grVn6Vn4_4Sigma->Draw("psame");
  grVn6Vn4_4p2Sigma->Draw("psame");
  grVn6Vn4_4p5Sigma->Draw("psame");
  legVn6Vn4->Draw("same");

  cCumuRatio->cd(2);
  grVn8Vn4Sys_4Sigma->Draw("apE2");
  grVn8Vn4_3p5Sigma->Draw("psame");
  grVn8Vn4_3p8Sigma->Draw("psame");
  grVn8Vn4_4Sigma->Draw("psame");
  grVn8Vn4_4p2Sigma->Draw("psame");
  grVn8Vn4_4p5Sigma->Draw("psame");
  legVn8Vn4->Draw("same");

  cCumuRatio->cd(3);
  grVn8Vn6Sys_4Sigma->Draw("apE2");
  grVn8Vn6_3p5Sigma->Draw("psame");
  grVn8Vn6_3p8Sigma->Draw("psame");
  grVn8Vn6_4Sigma->Draw("psame");
  grVn8Vn6_4p2Sigma->Draw("psame");
  grVn8Vn6_4p5Sigma->Draw("psame");
  legVn8Vn6->Draw("same");

  cCumuRatio->Update();
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- Ratio to default 
  TLegend * legVn6Vn4Sys = new TLegend(0.25, 0.2, 0.75, 0.45);
  legVn6Vn4Sys->SetBorderSize(0);
  legVn6Vn4Sys->SetFillStyle(0);
  legVn6Vn4Sys->AddEntry(grVn6Vn4_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn6Vn4Sys->AddEntry(grVn6Vn4_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn6Vn4Sys->AddEntry(grVn6Vn4_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn6Vn4Sys->AddEntry(grVn6Vn4_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TLegend * legVn8Vn4Sys = new TLegend(0.25, 0.2, 0.75, 0.45);
  legVn8Vn4Sys->SetBorderSize(0);
  legVn8Vn4Sys->SetFillStyle(0);
  legVn8Vn4Sys->AddEntry(grVn8Vn4_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn8Vn4Sys->AddEntry(grVn8Vn4_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn8Vn4Sys->AddEntry(grVn8Vn4_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn8Vn4Sys->AddEntry(grVn8Vn4_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TLegend * legVn8Vn6Sys = new TLegend(0.25, 0.2, 0.75, 0.45);
  legVn8Vn6Sys->SetBorderSize(0);
  legVn8Vn6Sys->SetFillStyle(0);
  legVn8Vn6Sys->AddEntry(grVn8Vn6_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legVn8Vn6Sys->AddEntry(grVn8Vn6_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legVn8Vn6Sys->AddEntry(grVn8Vn6_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legVn8Vn6Sys->AddEntry(grVn8Vn6_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);

  cCumuRatioSys->cd(1);
  grVn6Vn4_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn6Vn4_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn6Vn4_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn6Vn4_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn6Vn4Sys->Draw("same");
  lone_vn6vn4->Draw("same");

  cCumuRatioSys->cd(2);
  grVn8Vn4_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn8Vn4_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn8Vn4_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn8Vn4_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legVn8Vn4Sys->Draw("same");
  lone_vn8vn4->Draw("same");

  cCumuRatioSys->cd(3);
  grVn8Vn6_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grVn8Vn6_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grVn8Vn6_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grVn8Vn6_4p5Sigma_RatioTo4Sigma->Draw("psame");
  lone_vn8vn6->Draw("same");
  legVn8Vn6Sys->Draw("same");

  cCumuRatioSys->Update();
  cCumuRatioSys->SaveAs("cCumuRatioSys.pdf");


  TLegend * legG1e = new TLegend(0.2, 0.2, 0.55, 0.55);
  legG1e->SetBorderSize(0);
  legG1e->SetFillStyle(0);
  legG1e->AddEntry(grG1e_3p5Sigma,  "3.5 #sigma Cut", "lp");
  legG1e->AddEntry(grG1e_3p8Sigma,  "3.8 #sigma Cut", "lp");
  legG1e->AddEntry(grG1e_4Sigma,    "4 #sigma Cut",   "lp");
  legG1e->AddEntry(grG1e_4p2Sigma,  "4.2 #sigma Cut", "lp");
  legG1e->AddEntry(grG1e_4p5Sigma,  "4.5 #sigma Cut", "lp");

  TCanvas * cGamma1 = new TCanvas("cGamma1", "cGamma1", 500, 500);
  cGamma1->cd();
  grG1eSys_4Sigma->Draw("apE2");
  grG1e_3p5Sigma->Draw("psame");
  grG1e_3p8Sigma->Draw("psame");
  grG1e_4Sigma->Draw("psame");
  grG1e_4p2Sigma->Draw("psame");
  grG1e_4p5Sigma->Draw("psame");
  legG1e->Draw("same");
  cGamma1->SaveAs("cGamma1.pdf");

  //-- Ratio to Default
  TLegend * legG1eSys = new TLegend(0.5383, 0.6843, 0.9960, 0.8708);
  legG1eSys->SetBorderSize(0);
  legG1eSys->SetFillStyle(0);
  legG1eSys->AddEntry(grG1e_3p5Sigma_RatioTo4Sigma,  "3.5 #sigma Cut / Default", "lp");
  legG1eSys->AddEntry(grG1e_3p8Sigma_RatioTo4Sigma,  "3.8 #sigma Cut / Default", "lp");
  legG1eSys->AddEntry(grG1e_4p2Sigma_RatioTo4Sigma,  "4.2 #sigma Cut / Default", "lp");
  legG1eSys->AddEntry(grG1e_4p5Sigma_RatioTo4Sigma,  "4.5 #sigma Cut / Default", "lp");

  TCanvas * cGamma1Sys = new TCanvas("cGamma1Sys", "cGamma1Sys", 500, 500);
  cGamma1Sys->cd();
  grG1e_3p5Sigma_RatioTo4Sigma->Draw("ap");
  grG1e_3p8Sigma_RatioTo4Sigma->Draw("psame");
  grG1e_4p2Sigma_RatioTo4Sigma->Draw("psame");
  grG1e_4p5Sigma_RatioTo4Sigma->Draw("psame");
  legG1eSys->Draw("same");
  lone_g1e->Draw("same");
  cGamma1Sys->SaveAs("cGamma1Sys.pdf");

  //-- Big plot showing where the truncs are
  TLegend * legUnf = new TLegend(0.64, 0.2, 0.99, 0.55);
  legInit(legUnf);
  legUnf->AddEntry(l3p5s[0], "3.5#sigma", "l");
  legUnf->AddEntry(l3p8s[0], "3.8#sigma", "l");
  legUnf->AddEntry(l4s[0],   "4#sigma",   "l");
  legUnf->AddEntry(l4p2s[0], "4.2#sigma", "l");
  legUnf->AddEntry(l4p5s[0], "4.5#sigma", "l");

  TCanvas * cUnfBig = new TCanvas("cUnfBig", "cUnfBig", 2000, 1500);
  cUnfBig->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    cUnfBig->cd(icent+1);
    cUnfBig->cd(icent+1)->SetLogy();
    hFinalUnfold_Default[icent]->Draw();
    l3p5s[icent]->Draw("same");
    l3p8s[icent]->Draw("same");
    l4s[icent]->Draw("same");
    l4p2s[icent]->Draw("same");
    l4p5s[icent]->Draw("same");
    hFinalUnfold_Default[icent]->Draw("same");
    latex.DrawLatex(0.65, 0.88, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%") );
    if(icent == 0) legUnf->Draw("same");
  }

  cUnfBig->SaveAs("cUnfBig.pdf");

}

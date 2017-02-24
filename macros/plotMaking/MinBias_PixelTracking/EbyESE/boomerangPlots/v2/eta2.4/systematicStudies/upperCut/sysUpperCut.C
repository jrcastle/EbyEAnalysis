#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TExec.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"


#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

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
  double rMin6484 = 0.96;
  double rMax6484 = 1.02;
  double rMin86   = 0.998;
  double rMax86   = 1.012;
  double rMing1e  = -1.0;
  double rMaxg1e  = 2.0;

  TFile * fAna;
  TH1D * hObs[NCENT];
  TFile * fUnfold;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];

  TH1D * hFinalUnfold_Default[NCENT];
  TH1D * hFinalUnfold_Cut5Sigma[NCENT];
  TH1D * hFinalUnfold_Cut4Sigma[NCENT];
  TH1D * hFinalUnfold_Cut3Sigma[NCENT];

  TLatex latex;

  //-- Vectors galore!
  //-- Default
  std::vector<double> vn2_Default; 
  std::vector<double> vn2Cent_Default;
  std::vector<double> vn4_Default;
  std::vector<double> vn4Cent_Default;
  std::vector<double> vn6_Default;
  std::vector<double> vn6Cent_Default;
  std::vector<double> vn8_Default;
  std::vector<double> vn8Cent_Default;

  std::vector<double> g1e_Default;
  std::vector<double> g1eCent_Default;
  std::vector<double> vn6vn4_Default;
  std::vector<double> vn6vn4Cent_Default;
  std::vector<double> vn6vn4Npart_Default;
  std::vector<double> vn8vn4_Default;
  std::vector<double> vn8vn4Cent_Default;
  std::vector<double> vn8vn4Npart_Default;
  std::vector<double> vn8vn6_Default;
  std::vector<double> vn8vn6Cent_Default;

  //-- 5 Sigma
  std::vector<double> vn2_5Sigma;
  std::vector<double> vn2Cent_5Sigma;
  std::vector<double> vn4_5Sigma;
  std::vector<double> vn4Cent_5Sigma;
  std::vector<double> vn6_5Sigma;
  std::vector<double> vn6Cent_5Sigma;
  std::vector<double> vn8_5Sigma;
  std::vector<double> vn8Cent_5Sigma;

  std::vector<double> g1e_5Sigma;
  std::vector<double> g1eCent_5Sigma;
  std::vector<double> vn6vn4_5Sigma;
  std::vector<double> vn6vn4Cent_5Sigma;
  std::vector<double> vn6vn4Npart_5Sigma;
  std::vector<double> vn8vn4_5Sigma;
  std::vector<double> vn8vn4Cent_5Sigma;
  std::vector<double> vn8vn4Npart_5Sigma;
  std::vector<double> vn8vn6_5Sigma;
  std::vector<double> vn8vn6Cent_5Sigma;

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

  //-- 3 Sigma
  std::vector<double> vn2_3Sigma;
  std::vector<double> vn2Cent_3Sigma;
  std::vector<double> vn4_3Sigma;
  std::vector<double> vn4Cent_3Sigma;
  std::vector<double> vn6_3Sigma;
  std::vector<double> vn6Cent_3Sigma;
  std::vector<double> vn8_3Sigma;
  std::vector<double> vn8Cent_3Sigma;

  std::vector<double> g1e_3Sigma;
  std::vector<double> g1eCent_3Sigma;
  std::vector<double> vn6vn4_3Sigma;
  std::vector<double> vn6vn4Cent_3Sigma;
  std::vector<double> vn6vn4Npart_3Sigma;
  std::vector<double> vn8vn4_3Sigma;
  std::vector<double> vn8vn4Cent_3Sigma;
  std::vector<double> vn8vn4Npart_3Sigma;
  std::vector<double> vn8vn6_3Sigma;
  std::vector<double> vn8vn6Cent_3Sigma;

  TGraph * grVn2_Default;
  TGraph * grVn4_Default;
  TGraph * grVn6_Default;
  TGraph * grVn8_Default;
  TGraph * grG1e_Default;
  TGraph * grVn6Vn4_Default;
  TGraph * grVn8Vn4_Default;
  TGraph * grVn8Vn6_Default;

  TGraph * grVn2_5Sigma;
  TGraph * grVn4_5Sigma;
  TGraph * grVn6_5Sigma;
  TGraph * grVn8_5Sigma;
  TGraph * grG1e_5Sigma;
  TGraph * grVn6Vn4_5Sigma;
  TGraph * grVn8Vn4_5Sigma;
  TGraph * grVn8Vn6_5Sigma;

  TGraph * grVn2_4Sigma;
  TGraph * grVn4_4Sigma;
  TGraph * grVn6_4Sigma;
  TGraph * grVn8_4Sigma;
  TGraph * grG1e_4Sigma;
  TGraph * grVn6Vn4_4Sigma;
  TGraph * grVn8Vn4_4Sigma;
  TGraph * grVn8Vn6_4Sigma;

  TGraph * grVn2_3Sigma;
  TGraph * grVn4_3Sigma;
  TGraph * grVn6_3Sigma;
  TGraph * grVn8_3Sigma;
  TGraph * grG1e_3Sigma;
  TGraph * grVn6Vn4_3Sigma;
  TGraph * grVn8Vn4_3Sigma;
  TGraph * grVn8Vn6_3Sigma;

  //-- vs Npart
  TGraph * grVn6Vn4Npart_Default;
  TGraph * grVn8Vn4Npart_Default;
  TGraph * grVn6Vn4Npart_5Sigma;
  TGraph * grVn8Vn4Npart_5Sigma;
  TGraph * grVn6Vn4Npart_4Sigma;
  TGraph * grVn8Vn4Npart_4Sigma;
  TGraph * grVn6Vn4Npart_3Sigma;
  TGraph * grVn8Vn4Npart_3Sigma;

  //-- 5Sigma / Default
  std::vector<double> vn2_5Sigma_RatioToDefault;
  std::vector<double> vn2Cent_5Sigma_RatioToDefault;
  std::vector<double> vn4_5Sigma_RatioToDefault;
  std::vector<double> vn4Cent_5Sigma_RatioToDefault;
  std::vector<double> vn6_5Sigma_RatioToDefault;
  std::vector<double> vn6Cent_5Sigma_RatioToDefault;
  std::vector<double> vn8_5Sigma_RatioToDefault;
  std::vector<double> vn8Cent_5Sigma_RatioToDefault;

  std::vector<double> g1e_5Sigma_RatioToDefault;
  std::vector<double> g1eCent_5Sigma_RatioToDefault;
  std::vector<double> vn6vn4_5Sigma_RatioToDefault;
  std::vector<double> vn6vn4Cent_5Sigma_RatioToDefault;
  std::vector<double> vn8vn4_5Sigma_RatioToDefault;
  std::vector<double> vn8vn4Cent_5Sigma_RatioToDefault;
  std::vector<double> vn8vn6_5Sigma_RatioToDefault;
  std::vector<double> vn8vn6Cent_5Sigma_RatioToDefault;

  //-- 4Sigma / Default
  std::vector<double> vn2_4Sigma_RatioToDefault;
  std::vector<double> vn2Cent_4Sigma_RatioToDefault;
  std::vector<double> vn4_4Sigma_RatioToDefault;
  std::vector<double> vn4Cent_4Sigma_RatioToDefault;
  std::vector<double> vn6_4Sigma_RatioToDefault;
  std::vector<double> vn6Cent_4Sigma_RatioToDefault;
  std::vector<double> vn8_4Sigma_RatioToDefault;
  std::vector<double> vn8Cent_4Sigma_RatioToDefault;

  std::vector<double> g1e_4Sigma_RatioToDefault;
  std::vector<double> g1eCent_4Sigma_RatioToDefault;
  std::vector<double> vn6vn4_4Sigma_RatioToDefault;
  std::vector<double> vn6vn4Cent_4Sigma_RatioToDefault;
  std::vector<double> vn8vn4_4Sigma_RatioToDefault;
  std::vector<double> vn8vn4Cent_4Sigma_RatioToDefault;
  std::vector<double> vn8vn6_4Sigma_RatioToDefault;
  std::vector<double> vn8vn6Cent_4Sigma_RatioToDefault;

  //-- 3Sigma / Default
  std::vector<double> vn2_3Sigma_RatioToDefault;
  std::vector<double> vn2Cent_3Sigma_RatioToDefault;
  std::vector<double> vn4_3Sigma_RatioToDefault;
  std::vector<double> vn4Cent_3Sigma_RatioToDefault;
  std::vector<double> vn6_3Sigma_RatioToDefault;
  std::vector<double> vn6Cent_3Sigma_RatioToDefault;
  std::vector<double> vn8_3Sigma_RatioToDefault;
  std::vector<double> vn8Cent_3Sigma_RatioToDefault;

  std::vector<double> g1e_3Sigma_RatioToDefault;
  std::vector<double> g1eCent_3Sigma_RatioToDefault;
  std::vector<double> vn6vn4_3Sigma_RatioToDefault;
  std::vector<double> vn6vn4Cent_3Sigma_RatioToDefault;
  std::vector<double> vn8vn4_3Sigma_RatioToDefault;
  std::vector<double> vn8vn4Cent_3Sigma_RatioToDefault;
  std::vector<double> vn8vn6_3Sigma_RatioToDefault;
  std::vector<double> vn8vn6Cent_3Sigma_RatioToDefault;

  TGraph * grVn2_5Sigma_RatioToDefault;
  TGraph * grVn4_5Sigma_RatioToDefault;
  TGraph * grVn6_5Sigma_RatioToDefault;
  TGraph * grVn8_5Sigma_RatioToDefault;
  TGraph * grG1e_5Sigma_RatioToDefault;
  TGraph * grVn6Vn4_5Sigma_RatioToDefault;
  TGraph * grVn8Vn4_5Sigma_RatioToDefault;
  TGraph * grVn8Vn6_5Sigma_RatioToDefault;

  TGraph * grVn2_4Sigma_RatioToDefault;
  TGraph * grVn4_4Sigma_RatioToDefault;
  TGraph * grVn6_4Sigma_RatioToDefault;
  TGraph * grVn8_4Sigma_RatioToDefault;
  TGraph * grG1e_4Sigma_RatioToDefault;
  TGraph * grVn6Vn4_4Sigma_RatioToDefault;
  TGraph * grVn8Vn4_4Sigma_RatioToDefault;
  TGraph * grVn8Vn6_4Sigma_RatioToDefault;

  TGraph * grVn2_3Sigma_RatioToDefault;
  TGraph * grVn4_3Sigma_RatioToDefault;
  TGraph * grVn6_3Sigma_RatioToDefault;
  TGraph * grVn8_3Sigma_RatioToDefault;
  TGraph * grG1e_3Sigma_RatioToDefault;
  TGraph * grVn6Vn4_3Sigma_RatioToDefault;
  TGraph * grVn8Vn4_3Sigma_RatioToDefault;
  TGraph * grVn8Vn6_3Sigma_RatioToDefault;


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

  //-- Start fetching histograms....
  for(int icent = 0; icent < NCENT; icent++){

    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    for(int i = 0; i < NITER; i++){
      hUnfold[icent][i]     = (TH1D*) fUnfold->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i]     = (TH1D*) fUnfold->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");

      if( chi2 < 1.2 ){
        hFinalUnfold_Default[icent]   = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_c%i", icent) );
	hFinalUnfold_Cut5Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut5Sigma_c%i", icent) );
	hFinalUnfold_Cut4Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut4Sigma_c%i", icent) );
	hFinalUnfold_Cut3Sigma[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_Cut3Sigma_c%i", icent) );
	break;
      }

    } //-- End iter loop

  } //-- End cent loop

  //-- Now let's calculate some values!
  for(int icent = 0; icent < NCENT; icent++){

    //-- Default
    FixUnfold( hFinalUnfold_Default[icent] );
    EbyECumu cumuDef( hFinalUnfold_Default[icent] );
    double vn2_Def    = cumuDef.GetCumu_vn2();
    double vn4_Def    = cumuDef.GetCumu_vn4();
    double vn6_Def    = cumuDef.GetCumu_vn6();
    double vn8_Def    = cumuDef.GetCumu_vn8();
    double g1e_Def    = cumuDef.GetGamma1Exp();
    double vn6vn4_Def;
    double vn8vn4_Def;
    double vn8vn6_Def;

    if(vn2_Def != 0){
      vn2_Default.push_back(vn2_Def);
      vn2Cent_Default.push_back(centBinCenter[icent]);
    }
    if(vn4_Def != 0){
      vn4_Default.push_back(vn4_Def);
      vn4Cent_Default.push_back(centBinCenter[icent]);
    }
    if(vn6_Def != 0){
      vn6_Default.push_back(vn6_Def);
      vn6Cent_Default.push_back(centBinCenter[icent]);
    }
    if(vn8_Def != 0){
      vn8_Default.push_back(vn8_Def);
      vn8Cent_Default.push_back(centBinCenter[icent]);
    }
    if(g1e_Def != -10000){
      g1e_Default.push_back(g1e_Def);
      g1eCent_Default.push_back(centBinCenter[icent]);
    }
    if(vn4_Def != 0 && vn6_Def != 0){
      vn6vn4_Default.push_back( vn6_Def/vn4_Def );
      vn6vn4Cent_Default.push_back(centBinCenter[icent]);
      vn6vn4Npart_Default.push_back(Npart[icent]);
    }
    if(vn4_Def != 0 && vn8_Def!= 0){
      vn8vn4_Default.push_back( vn8_Def/vn4_Def );
      vn8vn4Cent_Default.push_back(centBinCenter[icent]);
      vn8vn4Npart_Default.push_back(Npart[icent]);
    }
    if(vn6_Def != 0 && vn8_Def!= 0){
      vn8vn6_Default.push_back( vn8_Def/vn6_Def );
      vn8vn6Cent_Default.push_back(centBinCenter[icent]);
    }

    //-- 5Sigma
    FIXUNFOLD( hFinalUnfold_Cut5Sigma[icent], 5 );
    EbyECumu cumu5Sig( hFinalUnfold_Cut5Sigma[icent] );
    double vn2_5Sig    = cumu5Sig.GetCumu_vn2();
    double vn4_5Sig    = cumu5Sig.GetCumu_vn4();
    double vn6_5Sig    = cumu5Sig.GetCumu_vn6();
    double vn8_5Sig    = cumu5Sig.GetCumu_vn8();
    double g1e_5Sig    = cumu5Sig.GetGamma1Exp();
    double vn6vn4_5Sig;
    double vn8vn4_5Sig;
    double vn8vn6_5Sig;

    if(vn2_5Sig != 0){
      vn2_5Sigma.push_back(vn2_5Sig);
      vn2Cent_5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_5Sig != 0){
      vn4_5Sigma.push_back(vn4_5Sig);
      vn4Cent_5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_5Sig != 0){
      vn6_5Sigma.push_back(vn6_5Sig);
      vn6Cent_5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_5Sig != 0){
      vn8_5Sigma.push_back(vn8_5Sig);
      vn8Cent_5Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_5Sig != -10000){
      g1e_5Sigma.push_back(g1e_5Sig);
      g1eCent_5Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_5Sig != 0 && vn6_5Sig != 0){
      vn6vn4_5Sigma.push_back( vn6_5Sig/vn4_5Sig );
      vn6vn4Cent_5Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_5Sigma.push_back(Npart[icent]);
    }
    if(vn4_5Sig != 0 && vn8_5Sig!= 0){
      vn8vn4_5Sigma.push_back( vn8_5Sig/vn4_5Sig );
      vn8vn4Cent_5Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_5Sigma.push_back(Npart[icent]);
    }
    if(vn6_5Sig != 0 && vn8_5Sig!= 0){
      vn8vn6_5Sigma.push_back( vn8_5Sig/vn6_5Sig );
      vn8vn6Cent_5Sigma.push_back(centBinCenter[icent]);
    }
    //-- 5Sigma/Default
    if(vn2_5Sig != 0 && vn2_Def != 0){
      vn2_5Sigma_RatioToDefault.push_back(vn2_5Sig / vn2_Def);
      vn2Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_5Sig != 0 && vn4_Def != 0){
      vn4_5Sigma_RatioToDefault.push_back(vn4_5Sig / vn4_Def);
      vn4Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_5Sig != 0 && vn6_Def != 0){
      vn6_5Sigma_RatioToDefault.push_back(vn6_5Sig / vn6_Def);
      vn6Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn8_5Sig != 0 && vn8_Def != 0){
      vn8_5Sigma_RatioToDefault.push_back(vn8_5Sig / vn8_Def);
      vn8Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(g1e_5Sig != -10000 && g1e_Def != 0){
      g1e_5Sigma_RatioToDefault.push_back(g1e_5Sig / g1e_Def);
      g1eCent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_5Sig != 0 && vn6_5Sig != 0 && vn4_Def != 0 && vn6_Def != 0){
      vn6vn4_5Sigma_RatioToDefault.push_back( (vn6_5Sig/vn4_5Sig) / (vn6_Def/vn4_Def) );
      vn6vn4Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_5Sig != 0 && vn8_5Sig!= 0 && vn4_Def != 0 && vn8_Def != 0){
      vn8vn4_5Sigma_RatioToDefault.push_back( (vn8_5Sig/vn4_5Sig) / (vn8_Def/vn4_Def) );
      vn8vn4Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_5Sig != 0 && vn8_5Sig!= 0 && vn6_Def != 0 && vn8_Def != 0){
      vn8vn6_5Sigma_RatioToDefault.push_back( (vn8_5Sig/vn6_5Sig) / (vn8_Def/vn6_Def) );
      vn8vn6Cent_5Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }

    //-- 4Sigma
    FIXUNFOLD( hFinalUnfold_Cut4Sigma[icent], 4 );
    EbyECumu cumu4Sig( hFinalUnfold_Cut4Sigma[icent] );
    double vn2_4Sig    = cumu4Sig.GetCumu_vn2();
    double vn4_4Sig    = cumu4Sig.GetCumu_vn4();
    double vn6_4Sig    = cumu4Sig.GetCumu_vn6();
    double vn8_4Sig    = cumu4Sig.GetCumu_vn8();
    double g1e_4Sig    = cumu4Sig.GetGamma1Exp();
    double vn6vn4_4Sig;
    double vn8vn4_4Sig;
    double vn8vn6_4Sig;

    if(vn2_4Sig != 0){
      vn2_4Sigma.push_back(vn2_4Sig);
      vn2Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig != 0){
      vn4_4Sigma.push_back(vn4_4Sig);
      vn4Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_4Sig != 0){
      vn6_4Sigma.push_back(vn6_4Sig);
      vn6Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_4Sig != 0){
      vn8_4Sigma.push_back(vn8_4Sig);
      vn8Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_4Sig != -10000){
      g1e_4Sigma.push_back(g1e_4Sig);
      g1eCent_4Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig != 0 && vn6_4Sig != 0){
      vn6vn4_4Sigma.push_back( vn6_4Sig/vn4_4Sig );
      vn6vn4Cent_4Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_4Sigma.push_back(Npart[icent]);
    }
    if(vn4_4Sig != 0 && vn8_4Sig!= 0){
      vn8vn4_4Sigma.push_back( vn8_4Sig/vn4_4Sig );
      vn8vn4Cent_4Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_4Sigma.push_back(Npart[icent]);
    }
    if(vn6_4Sig != 0 && vn8_4Sig!= 0){
      vn8vn6_4Sigma.push_back( vn8_4Sig/vn6_4Sig );
      vn8vn6Cent_4Sigma.push_back(centBinCenter[icent]);
    }
    //-- 4Sigma/Default
    if(vn2_4Sig != 0 && vn2_Def != 0){
      vn2_4Sigma_RatioToDefault.push_back(vn2_4Sig / vn2_Def);
      vn2Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig != 0 && vn4_Def != 0){
      vn4_4Sigma_RatioToDefault.push_back(vn4_4Sig / vn4_Def);
      vn4Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_4Sig != 0 && vn6_Def != 0){
      vn6_4Sigma_RatioToDefault.push_back(vn6_4Sig / vn6_Def);
      vn6Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn8_4Sig != 0 && vn8_Def != 0){
      vn8_4Sigma_RatioToDefault.push_back(vn8_4Sig / vn8_Def);
      vn8Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(g1e_4Sig != -10000 && g1e_Def != 0){
      g1e_4Sigma_RatioToDefault.push_back(g1e_4Sig / g1e_Def);
      g1eCent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig != 0 && vn6_4Sig != 0 && vn4_Def != 0 && vn6_Def != 0){
      vn6vn4_4Sigma_RatioToDefault.push_back( (vn6_4Sig/vn4_4Sig) / (vn6_Def/vn4_Def) );
      vn6vn4Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_4Sig != 0 && vn8_4Sig!= 0 && vn4_Def != 0 && vn8_Def != 0){
      vn8vn4_4Sigma_RatioToDefault.push_back( (vn8_4Sig/vn4_4Sig) / (vn8_Def/vn4_Def) );
      vn8vn4Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_4Sig != 0 && vn8_4Sig!= 0 && vn6_Def != 0 && vn8_Def != 0){
      vn8vn6_4Sigma_RatioToDefault.push_back( (vn8_4Sig/vn6_4Sig) / (vn8_Def/vn6_Def) );
      vn8vn6Cent_4Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }

    //-- 3Sigma
    FIXUNFOLD( hFinalUnfold_Cut3Sigma[icent], 3 );
    EbyECumu cumu3Sig( hFinalUnfold_Cut3Sigma[icent] );
    double vn2_3Sig    = cumu3Sig.GetCumu_vn2();
    double vn4_3Sig    = cumu3Sig.GetCumu_vn4();
    double vn6_3Sig    = cumu3Sig.GetCumu_vn6();
    double vn8_3Sig    = cumu3Sig.GetCumu_vn8();
    double g1e_3Sig    = cumu3Sig.GetGamma1Exp();
    double vn6vn4_3Sig;
    double vn8vn4_3Sig;
    double vn8vn6_3Sig;

    if(vn2_3Sig != 0){
      vn2_3Sigma.push_back(vn2_3Sig);
      vn2Cent_3Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3Sig != 0){
      vn4_3Sigma.push_back(vn4_3Sig);
      vn4Cent_3Sigma.push_back(centBinCenter[icent]);
    }
    if(vn6_3Sig != 0){
      vn6_3Sigma.push_back(vn6_3Sig);
      vn6Cent_3Sigma.push_back(centBinCenter[icent]);
    }
    if(vn8_3Sig != 0){
      vn8_3Sigma.push_back(vn8_3Sig);
      vn8Cent_3Sigma.push_back(centBinCenter[icent]);
    }
    if(g1e_3Sig != -10000){
      g1e_3Sigma.push_back(g1e_3Sig);
      g1eCent_3Sigma.push_back(centBinCenter[icent]);
    }
    if(vn4_3Sig != 0 && vn6_3Sig != 0){
      vn6vn4_3Sigma.push_back( vn6_3Sig/vn4_3Sig );
      vn6vn4Cent_3Sigma.push_back(centBinCenter[icent]);
      vn6vn4Npart_3Sigma.push_back(Npart[icent]);
    }
    if(vn4_3Sig != 0 && vn8_3Sig!= 0){
      vn8vn4_3Sigma.push_back( vn8_3Sig/vn4_3Sig );
      vn8vn4Cent_3Sigma.push_back(centBinCenter[icent]);
      vn8vn4Npart_3Sigma.push_back(Npart[icent]);
    }
    if(vn6_3Sig != 0 && vn8_3Sig!= 0){
      vn8vn6_3Sigma.push_back( vn8_3Sig/vn6_3Sig );
      vn8vn6Cent_3Sigma.push_back(centBinCenter[icent]);
    }
    //-- 3Sigma/Default
    if(vn2_3Sig != 0 && vn2_Def != 0){
      vn2_3Sigma_RatioToDefault.push_back(vn2_3Sig / vn2_Def);
      vn2Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_3Sig != 0 && vn4_Def != 0){
      vn4_3Sigma_RatioToDefault.push_back(vn4_3Sig / vn4_Def);
      vn4Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_3Sig != 0 && vn6_Def != 0){
      vn6_3Sigma_RatioToDefault.push_back(vn6_3Sig / vn6_Def);
      vn6Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn8_3Sig != 0 && vn8_Def != 0){
      vn8_3Sigma_RatioToDefault.push_back(vn8_3Sig / vn8_Def);
      vn8Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(g1e_3Sig != -10000 && g1e_Def != 0){
      g1e_3Sigma_RatioToDefault.push_back(g1e_3Sig / g1e_Def);
      g1eCent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_3Sig != 0 && vn6_3Sig != 0 && vn4_Def != 0 && vn6_Def != 0){
      vn6vn4_3Sigma_RatioToDefault.push_back( (vn6_3Sig/vn4_3Sig) / (vn6_Def/vn4_Def) );
      vn6vn4Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn4_3Sig != 0 && vn8_3Sig!= 0 && vn4_Def != 0 && vn8_Def != 0){
      vn8vn4_3Sigma_RatioToDefault.push_back( (vn8_3Sig/vn4_3Sig) / (vn8_Def/vn4_Def) );
      vn8vn4Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }
    if(vn6_3Sig != 0 && vn8_3Sig!= 0 && vn6_Def != 0 && vn8_Def != 0){
      vn8vn6_3Sigma_RatioToDefault.push_back( (vn8_3Sig/vn6_3Sig) / (vn8_Def/vn6_Def) );
      vn8vn6Cent_3Sigma_RatioToDefault.push_back(centBinCenter[icent]);
    }


  } //-- End Cent loop

  //-- Make sweet, sweet TGraphs
  //-- Default
  grVn2_Default    = new TGraph(vn2Cent_Default.size(), &(vn2Cent_Default[0]), &(vn2_Default[0]));
  grVn4_Default    = new TGraph(vn4Cent_Default.size(), &(vn4Cent_Default[0]), &(vn4_Default[0]));
  grVn6_Default    = new TGraph(vn6Cent_Default.size(), &(vn6Cent_Default[0]), &(vn6_Default[0]));
  grVn8_Default    = new TGraph(vn8Cent_Default.size(), &(vn8Cent_Default[0]), &(vn8_Default[0]));
  grG1e_Default    = new TGraph(g1eCent_Default.size(), &(g1eCent_Default[0]), &(g1e_Default[0]));
  grVn6Vn4_Default = new TGraph(vn6vn4Cent_Default.size(), &(vn6vn4Cent_Default[0]), &(vn6vn4_Default[0]));
  grVn8Vn4_Default = new TGraph(vn8vn4Cent_Default.size(), &(vn8vn4Cent_Default[0]), &(vn8vn4_Default[0]));
  grVn8Vn6_Default = new TGraph(vn8vn6Cent_Default.size(), &(vn8vn6Cent_Default[0]), &(vn8vn6_Default[0]));

  formatGraph(grVn2_Default,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         20, "grVn2_Default");
  formatGraph(grVn4_Default,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 20, "grVn4_Default");
  formatGraph(grVn6_Default,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         20, "grVn6_Default");
  formatGraph(grVn8_Default,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 20, "grVn8_Default");
  formatGraph(grG1e_Default,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         20, "grG1e_Default");
  formatGraph(grVn6Vn4_Default, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         20, "grVn6Vn4_Default");
  formatGraph(grVn8Vn4_Default, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  20, "grVn8Vn4_Default");
  formatGraph(grVn8Vn6_Default, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 20, "grVn8Vn6_Default");

  //-- 5Sigma
  grVn2_5Sigma    = new TGraph(vn2Cent_5Sigma.size(), &(vn2Cent_5Sigma[0]), &(vn2_5Sigma[0]));
  grVn4_5Sigma    = new TGraph(vn4Cent_5Sigma.size(), &(vn4Cent_5Sigma[0]), &(vn4_5Sigma[0]));
  grVn6_5Sigma    = new TGraph(vn6Cent_5Sigma.size(), &(vn6Cent_5Sigma[0]), &(vn6_5Sigma[0]));
  grVn8_5Sigma    = new TGraph(vn8Cent_5Sigma.size(), &(vn8Cent_5Sigma[0]), &(vn8_5Sigma[0]));
  grG1e_5Sigma    = new TGraph(g1eCent_5Sigma.size(), &(g1eCent_5Sigma[0]), &(g1e_5Sigma[0]));
  grVn6Vn4_5Sigma = new TGraph(vn6vn4Cent_5Sigma.size(), &(vn6vn4Cent_5Sigma[0]), &(vn6vn4_5Sigma[0]));
  grVn8Vn4_5Sigma = new TGraph(vn8vn4Cent_5Sigma.size(), &(vn8vn4Cent_5Sigma[0]), &(vn8vn4_5Sigma[0]));
  grVn8Vn6_5Sigma = new TGraph(vn8vn6Cent_5Sigma.size(), &(vn8vn6Cent_5Sigma[0]), &(vn8vn6_5Sigma[0]));

  formatGraph(grVn2_5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         25, "grVn2_5Sigma");
  formatGraph(grVn4_5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 25, "grVn4_5Sigma");
  formatGraph(grVn6_5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         25, "grVn6_5Sigma");
  formatGraph(grVn8_5Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 25, "grVn8_5Sigma");
  formatGraph(grG1e_5Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         25, "grG1e_5Sigma");
  formatGraph(grVn6Vn4_5Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         25, "grVn6Vn4_5Sigma");
  formatGraph(grVn8Vn4_5Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  25, "grVn8Vn4_5Sigma");
  formatGraph(grVn8Vn6_5Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 25, "grVn8Vn6_5Sigma");

  //-- 5Sigma/Default
  grVn2_5Sigma_RatioToDefault    = new TGraph(vn2Cent_5Sigma_RatioToDefault.size(), &(vn2Cent_5Sigma_RatioToDefault[0]), &(vn2_5Sigma_RatioToDefault[0]));
  grVn4_5Sigma_RatioToDefault    = new TGraph(vn4Cent_5Sigma_RatioToDefault.size(), &(vn4Cent_5Sigma_RatioToDefault[0]), &(vn4_5Sigma_RatioToDefault[0]));
  grVn6_5Sigma_RatioToDefault    = new TGraph(vn6Cent_5Sigma_RatioToDefault.size(), &(vn6Cent_5Sigma_RatioToDefault[0]), &(vn6_5Sigma_RatioToDefault[0]));
  grVn8_5Sigma_RatioToDefault    = new TGraph(vn8Cent_5Sigma_RatioToDefault.size(), &(vn8Cent_5Sigma_RatioToDefault[0]), &(vn8_5Sigma_RatioToDefault[0]));
  grG1e_5Sigma_RatioToDefault    = new TGraph(g1eCent_5Sigma_RatioToDefault.size(), &(g1eCent_5Sigma_RatioToDefault[0]), &(g1e_5Sigma_RatioToDefault[0]));
  grVn6Vn4_5Sigma_RatioToDefault = new TGraph(vn6vn4Cent_5Sigma_RatioToDefault.size(), &(vn6vn4Cent_5Sigma_RatioToDefault[0]), &(vn6vn4_5Sigma_RatioToDefault[0]));
  grVn8Vn4_5Sigma_RatioToDefault = new TGraph(vn8vn4Cent_5Sigma_RatioToDefault.size(), &(vn8vn4Cent_5Sigma_RatioToDefault[0]), &(vn8vn4_5Sigma_RatioToDefault[0]));
  grVn8Vn6_5Sigma_RatioToDefault = new TGraph(vn8vn6Cent_5Sigma_RatioToDefault.size(), &(vn8vn6Cent_5Sigma_RatioToDefault[0]), &(vn8vn6_5Sigma_RatioToDefault[0]));

  formatGraph(grVn2_5Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         25, "grVn2_5Sigma_RatioToDefault");
  formatGraph(grVn4_5Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 25, "grVn4_5Sigma_RatioToDefault");
  formatGraph(grVn6_5Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         25, "grVn6_5Sigma_RatioToDefault");
  formatGraph(grVn8_5Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 25, "grVn8_5Sigma_RatioToDefault");
  formatGraph(grG1e_5Sigma_RatioToDefault,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         25, "grG1e_5Sigma_RatioToDefault");
  formatGraph(grVn6Vn4_5Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         25, "grVn6Vn4_5Sigma_RatioToDefault");
  formatGraph(grVn8Vn4_5Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  25, "grVn8Vn4_5Sigma_RatioToDefault");
  formatGraph(grVn8Vn6_5Sigma_RatioToDefault, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 25, "grVn8Vn6_5Sigma_RatioToDefault");

  //-- 4Sigma
  grVn2_4Sigma    = new TGraph(vn2Cent_4Sigma.size(), &(vn2Cent_4Sigma[0]), &(vn2_4Sigma[0]));
  grVn4_4Sigma    = new TGraph(vn4Cent_4Sigma.size(), &(vn4Cent_4Sigma[0]), &(vn4_4Sigma[0]));
  grVn6_4Sigma    = new TGraph(vn6Cent_4Sigma.size(), &(vn6Cent_4Sigma[0]), &(vn6_4Sigma[0]));
  grVn8_4Sigma    = new TGraph(vn8Cent_4Sigma.size(), &(vn8Cent_4Sigma[0]), &(vn8_4Sigma[0]));
  grG1e_4Sigma    = new TGraph(g1eCent_4Sigma.size(), &(g1eCent_4Sigma[0]), &(g1e_4Sigma[0]));
  grVn6Vn4_4Sigma = new TGraph(vn6vn4Cent_4Sigma.size(), &(vn6vn4Cent_4Sigma[0]), &(vn6vn4_4Sigma[0]));
  grVn8Vn4_4Sigma = new TGraph(vn8vn4Cent_4Sigma.size(), &(vn8vn4Cent_4Sigma[0]), &(vn8vn4_4Sigma[0]));
  grVn8Vn6_4Sigma = new TGraph(vn8vn6Cent_4Sigma.size(), &(vn8vn6Cent_4Sigma[0]), &(vn8vn6_4Sigma[0]));

  formatGraph(grVn2_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         26, "grVn2_4Sigma");
  formatGraph(grVn4_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 26, "grVn4_4Sigma");
  formatGraph(grVn6_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         26, "grVn6_4Sigma");
  formatGraph(grVn8_4Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 26, "grVn8_4Sigma");
  formatGraph(grG1e_4Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         26, "grG1e_4Sigma");
  formatGraph(grVn6Vn4_4Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         26, "grVn6Vn4_4Sigma");
  formatGraph(grVn8Vn4_4Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  26, "grVn8Vn4_4Sigma");
  formatGraph(grVn8Vn6_4Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 26, "grVn8Vn6_4Sigma");

  //-- 4Sigma/Default
  grVn2_4Sigma_RatioToDefault    = new TGraph(vn2Cent_4Sigma_RatioToDefault.size(), &(vn2Cent_4Sigma_RatioToDefault[0]), &(vn2_4Sigma_RatioToDefault[0]));
  grVn4_4Sigma_RatioToDefault    = new TGraph(vn4Cent_4Sigma_RatioToDefault.size(), &(vn4Cent_4Sigma_RatioToDefault[0]), &(vn4_4Sigma_RatioToDefault[0]));
  grVn6_4Sigma_RatioToDefault    = new TGraph(vn6Cent_4Sigma_RatioToDefault.size(), &(vn6Cent_4Sigma_RatioToDefault[0]), &(vn6_4Sigma_RatioToDefault[0]));
  grVn8_4Sigma_RatioToDefault    = new TGraph(vn8Cent_4Sigma_RatioToDefault.size(), &(vn8Cent_4Sigma_RatioToDefault[0]), &(vn8_4Sigma_RatioToDefault[0]));
  grG1e_4Sigma_RatioToDefault    = new TGraph(g1eCent_4Sigma_RatioToDefault.size(), &(g1eCent_4Sigma_RatioToDefault[0]), &(g1e_4Sigma_RatioToDefault[0]));
  grVn6Vn4_4Sigma_RatioToDefault = new TGraph(vn6vn4Cent_4Sigma_RatioToDefault.size(), &(vn6vn4Cent_4Sigma_RatioToDefault[0]), &(vn6vn4_4Sigma_RatioToDefault[0]));
  grVn8Vn4_4Sigma_RatioToDefault = new TGraph(vn8vn4Cent_4Sigma_RatioToDefault.size(), &(vn8vn4Cent_4Sigma_RatioToDefault[0]), &(vn8vn4_4Sigma_RatioToDefault[0]));
  grVn8Vn6_4Sigma_RatioToDefault = new TGraph(vn8vn6Cent_4Sigma_RatioToDefault.size(), &(vn8vn6Cent_4Sigma_RatioToDefault[0]), &(vn8vn6_4Sigma_RatioToDefault[0]));

  formatGraph(grVn2_4Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         26, "grVn2_4Sigma_RatioToDefault");
  formatGraph(grVn4_4Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 26, "grVn4_4Sigma_RatioToDefault");
  formatGraph(grVn6_4Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         26, "grVn6_4Sigma_RatioToDefault");
  formatGraph(grVn8_4Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 26, "grVn8_4Sigma_RatioToDefault");
  formatGraph(grG1e_4Sigma_RatioToDefault,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         26, "grG1e_4Sigma_RatioToDefault");
  formatGraph(grVn6Vn4_4Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         26, "grVn6Vn4_4Sigma_RatioToDefault");
  formatGraph(grVn8Vn4_4Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  26, "grVn8Vn4_4Sigma_RatioToDefault");
  formatGraph(grVn8Vn6_4Sigma_RatioToDefault, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 26, "grVn8Vn6_4Sigma_RatioToDefault");


  //-- 3Sigma
  grVn2_3Sigma    = new TGraph(vn2Cent_3Sigma.size(), &(vn2Cent_3Sigma[0]), &(vn2_3Sigma[0]));
  grVn4_3Sigma    = new TGraph(vn4Cent_3Sigma.size(), &(vn4Cent_3Sigma[0]), &(vn4_3Sigma[0]));
  grVn6_3Sigma    = new TGraph(vn6Cent_3Sigma.size(), &(vn6Cent_3Sigma[0]), &(vn6_3Sigma[0]));
  grVn8_3Sigma    = new TGraph(vn8Cent_3Sigma.size(), &(vn8Cent_3Sigma[0]), &(vn8_3Sigma[0]));
  grG1e_3Sigma    = new TGraph(g1eCent_3Sigma.size(), &(g1eCent_3Sigma[0]), &(g1e_3Sigma[0]));
  grVn6Vn4_3Sigma = new TGraph(vn6vn4Cent_3Sigma.size(), &(vn6vn4Cent_3Sigma[0]), &(vn6vn4_3Sigma[0]));
  grVn8Vn4_3Sigma = new TGraph(vn8vn4Cent_3Sigma.size(), &(vn8vn4Cent_3Sigma[0]), &(vn8vn4_3Sigma[0]));
  grVn8Vn6_3Sigma = new TGraph(vn8vn6Cent_3Sigma.size(), &(vn8vn6Cent_3Sigma[0]), &(vn8vn6_3Sigma[0]));

  formatGraph(grVn2_3Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{2}", norder_),                    1,         28, "grVn2_3Sigma");
  formatGraph(grVn4_3Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{4}", norder_),                    kSpring+4, 28, "grVn4_3Sigma");
  formatGraph(grVn6_3Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{6}", norder_),                    6,         28, "grVn6_3Sigma");
  formatGraph(grVn8_3Sigma,    "Centrality %", vnCumuMin, vnCumuMax, Form("v_{%i}{8}", norder_),                    kOrange+7, 28, "grVn8_3Sigma");
  formatGraph(grG1e_3Sigma,    "Centrality %", g1eMin,    g1eMax,   "#gamma_{1}^{exp}",                             2,         28, "grG1e_3Sigma");
  formatGraph(grVn6Vn4_3Sigma, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 4,         28, "grVn6Vn4_3Sigma");
  formatGraph(grVn8Vn4_3Sigma, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), kGreen+2,  28, "grVn8Vn4_3Sigma");
  formatGraph(grVn8Vn6_3Sigma, "Centrality %", vn8vn6Min, vn8vn6Max, Form("v_{%i}{8}/v_{%i}{6}", norder_, norder_), kViolet-1, 28, "grVn8Vn6_3Sigma");

  //-- 3Sigma/Default
  grVn2_3Sigma_RatioToDefault    = new TGraph(vn2Cent_3Sigma_RatioToDefault.size(), &(vn2Cent_3Sigma_RatioToDefault[0]), &(vn2_3Sigma_RatioToDefault[0]));
  grVn4_3Sigma_RatioToDefault    = new TGraph(vn4Cent_3Sigma_RatioToDefault.size(), &(vn4Cent_3Sigma_RatioToDefault[0]), &(vn4_3Sigma_RatioToDefault[0]));
  grVn6_3Sigma_RatioToDefault    = new TGraph(vn6Cent_3Sigma_RatioToDefault.size(), &(vn6Cent_3Sigma_RatioToDefault[0]), &(vn6_3Sigma_RatioToDefault[0]));
  grVn8_3Sigma_RatioToDefault    = new TGraph(vn8Cent_3Sigma_RatioToDefault.size(), &(vn8Cent_3Sigma_RatioToDefault[0]), &(vn8_3Sigma_RatioToDefault[0]));
  grG1e_3Sigma_RatioToDefault    = new TGraph(g1eCent_3Sigma_RatioToDefault.size(), &(g1eCent_3Sigma_RatioToDefault[0]), &(g1e_3Sigma_RatioToDefault[0]));
  grVn6Vn4_3Sigma_RatioToDefault = new TGraph(vn6vn4Cent_3Sigma_RatioToDefault.size(), &(vn6vn4Cent_3Sigma_RatioToDefault[0]), &(vn6vn4_3Sigma_RatioToDefault[0]));
  grVn8Vn4_3Sigma_RatioToDefault = new TGraph(vn8vn4Cent_3Sigma_RatioToDefault.size(), &(vn8vn4Cent_3Sigma_RatioToDefault[0]), &(vn8vn4_3Sigma_RatioToDefault[0]));
  grVn8Vn6_3Sigma_RatioToDefault = new TGraph(vn8vn6Cent_3Sigma_RatioToDefault.size(), &(vn8vn6Cent_3Sigma_RatioToDefault[0]), &(vn8vn6_3Sigma_RatioToDefault[0]));

  formatGraph(grVn2_3Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 1,         28, "grVn2_3Sigma_RatioToDefault");
  formatGraph(grVn4_3Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kSpring+4, 28, "grVn4_3Sigma_RatioToDefault");
  formatGraph(grVn6_3Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", 6,         28, "grVn6_3Sigma_RatioToDefault");
  formatGraph(grVn8_3Sigma_RatioToDefault,    "Centrality %", rMinCumu, rMaxCumu, "Ratio to Default", kOrange+7, 28, "grVn8_3Sigma_RatioToDefault");
  formatGraph(grG1e_3Sigma_RatioToDefault,    "Centrality %", rMing1e,  rMaxg1e,  "Ratio to Default", 2,         28, "grG1e_3Sigma_RatioToDefault");
  formatGraph(grVn6Vn4_3Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", 4,         28, "grVn6Vn4_3Sigma_RatioToDefault");
  formatGraph(grVn8Vn4_3Sigma_RatioToDefault, "Centrality %", rMin6484, rMax6484, "Ratio to Default", kGreen+2,  28, "grVn8Vn4_3Sigma_RatioToDefault");
  formatGraph(grVn8Vn6_3Sigma_RatioToDefault, "Centrality %", rMin86,   rMax86,   "Ratio to Default", kViolet-1, 28, "grVn8Vn6_3Sigma_RatioToDefault");


  TLegend * legVn2 = new TLegend(0.2, 0.67, 0.5, 0.95);
  legVn2->SetBorderSize(0);
  legVn2->SetFillStyle(0);
  legVn2->AddEntry(grVn2_Default, "Default", "lp");
  legVn2->AddEntry(grVn2_5Sigma,  "5#sigma Cut", "lp");
  legVn2->AddEntry(grVn2_4Sigma,  "4#sigma Cut", "lp");
  legVn2->AddEntry(grVn2_3Sigma,  "3#sigma Cut", "lp");

  TLegend * legVn4 = new TLegend(0.2, 0.67, 0.5, 0.95);
  legVn4->SetBorderSize(0);
  legVn4->SetFillStyle(0);
  legVn4->AddEntry(grVn4_Default, "Default", "lp");
  legVn4->AddEntry(grVn4_5Sigma,  "5#sigma Cut", "lp");
  legVn4->AddEntry(grVn4_4Sigma,  "4#sigma Cut", "lp");
  legVn4->AddEntry(grVn4_3Sigma,  "3#sigma Cut", "lp");

  TLegend * legVn6 = new TLegend(0.2, 0.67, 0.5, 0.95);
  legVn6->SetBorderSize(0);
  legVn6->SetFillStyle(0);
  legVn6->AddEntry(grVn6_Default, "Default", "lp");
  legVn6->AddEntry(grVn6_5Sigma,  "5#sigma Cut", "lp");
  legVn6->AddEntry(grVn6_4Sigma,  "4#sigma Cut", "lp");
  legVn6->AddEntry(grVn6_3Sigma,  "3#sigma Cut", "lp");

  TLegend * legVn8 = new TLegend(0.2, 0.67, 0.5, 0.93);
  legVn8->SetBorderSize(0);
  legVn8->SetFillStyle(0);
  legVn8->AddEntry(grVn8_Default, "Default", "lp");
  legVn8->AddEntry(grVn8_5Sigma,  "5#sigma Cut", "lp");
  legVn8->AddEntry(grVn8_4Sigma,  "4#sigma Cut", "lp");
  legVn8->AddEntry(grVn8_3Sigma,  "3#sigma Cut", "lp");

  TCanvas * cCumu = new TCanvas("cCumu", "cCumu", 2000, 500);
  cCumu->Divide(4,1);
  cCumu->cd(1);
  grVn2_Default->Draw("ap");
  grVn2_5Sigma->Draw("psame");
  grVn2_4Sigma->Draw("psame");
  grVn2_3Sigma->Draw("psame");
  legVn2->Draw("same");
  cCumu->cd(2);
  grVn4_Default->Draw("ap");
  grVn4_5Sigma->Draw("psame");
  grVn4_4Sigma->Draw("psame");
  grVn4_3Sigma->Draw("psame");
  legVn4->Draw("same");
  cCumu->cd(3);
  grVn6_Default->Draw("ap");
  grVn6_5Sigma->Draw("psame");
  grVn6_4Sigma->Draw("psame");
  grVn6_3Sigma->Draw("psame");
  legVn6->Draw("same");
  cCumu->cd(4);
  grVn8_Default->Draw("ap");
  grVn8_5Sigma->Draw("psame");
  grVn8_4Sigma->Draw("psame");
  grVn8_3Sigma->Draw("psame");
  legVn8->Draw("same");
  cCumu->SaveAs("cCumu.pdf");

  //-- sysCumu
  TLegend * legVn2Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn2Sys->SetBorderSize(0);
  legVn2Sys->SetFillStyle(0);
  legVn2Sys->AddEntry(grVn2_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn2Sys->AddEntry(grVn2_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn2Sys->AddEntry(grVn2_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TLegend * legVn4Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn4Sys->SetBorderSize(0);
  legVn4Sys->SetFillStyle(0);
  legVn4Sys->AddEntry(grVn4_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn4Sys->AddEntry(grVn4_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn4Sys->AddEntry(grVn4_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TLegend * legVn6Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn6Sys->SetBorderSize(0);
  legVn6Sys->SetFillStyle(0);
  legVn6Sys->AddEntry(grVn6_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn6Sys->AddEntry(grVn6_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn6Sys->AddEntry(grVn6_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TLegend * legVn8Sys = new TLegend(0.6257, 0.1778, 0.9901, 0.4135);
  legVn8Sys->SetBorderSize(0);
  legVn8Sys->SetFillStyle(0);
  legVn8Sys->AddEntry(grVn8_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn8Sys->AddEntry(grVn8_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn8Sys->AddEntry(grVn8_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TCanvas * cCumuSys = new TCanvas("cCumuSys", "cCumuSys", 2000, 500);
  cCumuSys->Divide(4,1);
  cCumuSys->cd(1);
  grVn2_5Sigma_RatioToDefault->Draw("ap");
  grVn2_4Sigma_RatioToDefault->Draw("psame");
  grVn2_3Sigma_RatioToDefault->Draw("psame");
  legVn2Sys->Draw("same");
  cCumuSys->cd(2);
  grVn4_5Sigma_RatioToDefault->Draw("ap");
  grVn4_4Sigma_RatioToDefault->Draw("psame");
  grVn4_3Sigma_RatioToDefault->Draw("psame");
  legVn4Sys->Draw("same");
  cCumuSys->cd(3);
  grVn6_5Sigma_RatioToDefault->Draw("ap");
  grVn6_4Sigma_RatioToDefault->Draw("psame");
  grVn6_3Sigma_RatioToDefault->Draw("psame");
  legVn6Sys->Draw("same");
  cCumuSys->cd(4);
  grVn8_5Sigma_RatioToDefault->Draw("ap");
  grVn8_4Sigma_RatioToDefault->Draw("psame");
  grVn8_3Sigma_RatioToDefault->Draw("psame");
  legVn8Sys->Draw("same");
  cCumuSys->SaveAs("cCumuSys.pdf");


  TLegend * legVn6Vn4 = new TLegend(0.6, 0.67, 0.9, 0.93);
  legVn6Vn4->SetBorderSize(0);
  legVn6Vn4->SetFillStyle(0);
  legVn6Vn4->AddEntry(grVn6Vn4_Default, "Default", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_5Sigma,  "5#sigma Cut", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_4Sigma,  "4#sigma Cut", "lp");
  legVn6Vn4->AddEntry(grVn6Vn4_3Sigma,  "3#sigma Cut", "lp");

  TLegend * legVn8Vn4 = new TLegend(0.6, 0.67, 0.9, 0.93);
  legVn8Vn4->SetBorderSize(0);
  legVn8Vn4->SetFillStyle(0);
  legVn8Vn4->AddEntry(grVn8Vn4_Default, "Default", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_5Sigma,  "5#sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_4Sigma,  "4#sigma Cut", "lp");
  legVn8Vn4->AddEntry(grVn8Vn4_3Sigma,  "3#sigma Cut", "lp");

  TLegend * legVn8Vn6 = new TLegend(0.6, 0.67, 0.9, 0.93);
  legVn8Vn6->SetBorderSize(0);
  legVn8Vn6->SetFillStyle(0);
  legVn8Vn6->AddEntry(grVn8Vn6_Default, "Default", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_5Sigma,  "5#sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_4Sigma,  "4#sigma Cut", "lp");
  legVn8Vn6->AddEntry(grVn8Vn6_3Sigma,  "3#sigma Cut", "lp");


  TCanvas * cCumuRatio = new TCanvas("cCumuRatio", "cCumuRatio", 1500, 500);
  cCumuRatio->Divide(3,1);
  cCumuRatio->cd(1);
  grVn6Vn4_Default->Draw("ap");
  grVn6Vn4_5Sigma->Draw("psame");
  grVn6Vn4_4Sigma->Draw("psame");
  grVn6Vn4_3Sigma->Draw("psame");
  legVn6Vn4->Draw("same");
  cCumuRatio->cd(2);
  grVn8Vn4_Default->Draw("ap");
  grVn8Vn4_5Sigma->Draw("psame");
  grVn8Vn4_4Sigma->Draw("psame");
  grVn8Vn4_3Sigma->Draw("psame");
  legVn8Vn4->Draw("same");
  cCumuRatio->cd(3);
  grVn8Vn6_Default->Draw("ap");
  grVn8Vn6_5Sigma->Draw("psame");
  grVn8Vn6_4Sigma->Draw("psame");
  grVn8Vn6_3Sigma->Draw("psame");
  legVn8Vn6->Draw("same");
  cCumuRatio->Update();
  cCumuRatio->SaveAs("cCumuRatio.pdf");

  //-- Ratio to default 
  TLegend * legVn6Vn4Sys = new TLegend(0.5339, 0.7443, 0.9916, 0.9302);
  legVn6Vn4Sys->SetBorderSize(0);
  legVn6Vn4Sys->SetFillStyle(0);
  legVn6Vn4Sys->AddEntry(grVn6Vn4_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn6Vn4Sys->AddEntry(grVn6Vn4_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn6Vn4Sys->AddEntry(grVn6Vn4_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TLegend * legVn8Vn4Sys = new TLegend(0.5339, 0.7443, 0.9916, 0.9302);
  legVn8Vn4Sys->SetBorderSize(0);
  legVn8Vn4Sys->SetFillStyle(0);
  legVn8Vn4Sys->AddEntry(grVn8Vn4_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn8Vn4Sys->AddEntry(grVn8Vn4_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn8Vn4Sys->AddEntry(grVn8Vn4_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TLegend * legVn8Vn6Sys = new TLegend(0.5339, 0.7443, 0.9916, 0.9302);
  legVn8Vn6Sys->SetBorderSize(0);
  legVn8Vn6Sys->SetFillStyle(0);
  legVn8Vn6Sys->AddEntry(grVn8Vn6_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legVn8Vn6Sys->AddEntry(grVn8Vn6_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legVn8Vn6Sys->AddEntry(grVn8Vn6_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TCanvas * cCumuRatioSys = new TCanvas("cCumuRatioSys", "cCumuRatioSys", 1500, 500);
  cCumuRatioSys->Divide(3,1);
  cCumuRatioSys->cd(1);
  grVn6Vn4_5Sigma_RatioToDefault->Draw("ap");
  grVn6Vn4_4Sigma_RatioToDefault->Draw("psame");
  grVn6Vn4_3Sigma_RatioToDefault->Draw("psame");
  legVn6Vn4Sys->Draw("same");
  cCumuRatioSys->cd(2);
  grVn8Vn4_5Sigma_RatioToDefault->Draw("ap");
  grVn8Vn4_4Sigma_RatioToDefault->Draw("psame");
  grVn8Vn4_3Sigma_RatioToDefault->Draw("psame");
  legVn8Vn4Sys->Draw("same");
  cCumuRatioSys->cd(3);
  grVn8Vn6_5Sigma_RatioToDefault->Draw("ap");
  grVn8Vn6_4Sigma_RatioToDefault->Draw("psame");
  grVn8Vn6_3Sigma_RatioToDefault->Draw("psame");
  legVn8Vn6Sys->Draw("same");
  cCumuRatioSys->Update();
  cCumuRatioSys->SaveAs("cCumuRatioSys.pdf");


  TLegend * legG1e = new TLegend(0.6, 0.67, 0.9, 0.93);
  legG1e->SetBorderSize(0);
  legG1e->SetFillStyle(0);
  legG1e->AddEntry(grG1e_Default, "Default", "lp");
  legG1e->AddEntry(grG1e_5Sigma,  "5#sigma Cut", "lp");
  legG1e->AddEntry(grG1e_4Sigma,  "4#sigma Cut", "lp");
  legG1e->AddEntry(grG1e_3Sigma,  "3#sigma Cut", "lp");

  TCanvas * cGamma1 = new TCanvas("cGamma1", "cGamma1", 500, 500);
  cGamma1->cd();
  grG1e_Default->Draw("ap");
  grG1e_5Sigma->Draw("psame");
  grG1e_4Sigma->Draw("psame");
  grG1e_3Sigma->Draw("psame");
  legG1e->Draw("same");
  cGamma1->SaveAs("cGamma1.pdf");

  //-- Ratio to Default
  TLegend * legG1eSys = new TLegend(0.5383, 0.1843, 0.9960, 0.3708);
  legG1eSys->SetBorderSize(0);
  legG1eSys->SetFillStyle(0);
  legG1eSys->AddEntry(grG1e_5Sigma_RatioToDefault,  "5#sigma Cut / Default", "lp");
  legG1eSys->AddEntry(grG1e_4Sigma_RatioToDefault,  "4#sigma Cut / Default", "lp");
  legG1eSys->AddEntry(grG1e_3Sigma_RatioToDefault,  "3#sigma Cut / Default", "lp");

  TCanvas * cGamma1Sys = new TCanvas("cGamma1Sys", "cGamma1Sys", 500, 500);
  cGamma1Sys->cd();
  grG1e_5Sigma_RatioToDefault->Draw("ap");
  grG1e_4Sigma_RatioToDefault->Draw("psame");
  grG1e_3Sigma_RatioToDefault->Draw("psame");
  legG1eSys->Draw("same");
  cGamma1Sys->SaveAs("cGamma1Sys.pdf");


  //-------------------------------------------------------------------
  //-- ATLAS Comparison
  if( compATLAS ){

    double ratioMin = 0.9;
    double ratioMax = 1.05;

    double npartSysWidth[NCENT];
    for(int icent = 0; icent < NCENT; icent++) npartSysWidth[icent] = 5;

    //-- CMS 5.02 TeV, 0.3 < pT < 3.0 GeV, |eta| < 2.4
    //-- vn{6} / vn{4}
    grVn6Vn4Npart_Default = new TGraph(vn6vn4Npart_Default.size(), &(vn6vn4Npart_Default[0]), &(vn6vn4_Default[0]));
    grVn6Vn4Npart_5Sigma  = new TGraph(vn6vn4Npart_5Sigma.size(),  &(vn6vn4Npart_5Sigma[0]),  &(vn6vn4_5Sigma[0]));
    grVn6Vn4Npart_4Sigma  = new TGraph(vn6vn4Npart_4Sigma.size(),  &(vn6vn4Npart_4Sigma[0]),  &(vn6vn4_4Sigma[0]));
    grVn6Vn4Npart_3Sigma  = new TGraph(vn6vn4Npart_3Sigma.size(),  &(vn6vn4Npart_3Sigma[0]),  &(vn6vn4_3Sigma[0]));

    formatGraph(grVn6Vn4Npart_Default, "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 2, 20, "grVn6Vn4Npart_Default");
    formatGraph(grVn6Vn4Npart_5Sigma,  "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 2, 25, "grVn6Vn4Npart_5Sigma");
    formatGraph(grVn6Vn4Npart_4Sigma,  "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 2, 26, "grVn6Vn4Npart_4Sigma");
    formatGraph(grVn6Vn4Npart_3Sigma,  "Centrality %", vn6vn4Min, vn6vn4Max, Form("v_{%i}{6}/v_{%i}{4}", norder_, norder_), 2, 28, "grVn6Vn4Npart_3Sigma");

    //-- vn{8} / vn{4}
    grVn8Vn4Npart_Default = new TGraph(vn8vn4Npart_Default.size(), &(vn8vn4Npart_Default[0]), &(vn8vn4_Default[0]));
    grVn8Vn4Npart_5Sigma  = new TGraph(vn8vn4Npart_5Sigma.size(),  &(vn8vn4Npart_5Sigma[0]),  &(vn8vn4_5Sigma[0]));
    grVn8Vn4Npart_4Sigma  = new TGraph(vn8vn4Npart_4Sigma.size(),  &(vn8vn4Npart_4Sigma[0]),  &(vn8vn4_4Sigma[0]));
    grVn8Vn4Npart_3Sigma  = new TGraph(vn8vn4Npart_3Sigma.size(),  &(vn8vn4Npart_3Sigma[0]),  &(vn8vn4_3Sigma[0]));

    formatGraph(grVn8Vn4Npart_Default, "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), 2, 20, "grVn8Vn4Npart_Default");
    formatGraph(grVn8Vn4Npart_5Sigma,  "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), 2, 25, "grVn8Vn4Npart_5Sigma");
    formatGraph(grVn8Vn4Npart_4Sigma,  "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), 2, 26, "grVn8Vn4Npart_4Sigma");
    formatGraph(grVn8Vn4Npart_3Sigma,  "Centrality %", vn8vn4Min, vn8vn4Max, Form("v_{%i}{8}/v_{%i}{4}", norder_, norder_), 2, 28, "grVn8Vn4Npart_3Sigma");

    //-- ATLAS 2.76 TeV, 0.5 < pT < 20 GeV, |eta| < 2.5
    //-- vn{6} / vn{4}
    double ATLASVn6Vn4_xval[] = { 22.6, 46.1, 59.9, 76.1, 95.0, 117.0, 142.0, 170.0, 203.0, 239.5, 281.9, 330.3, 372.533 };
    double ATLASVn6Vn4_xerrminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ATLASVn6Vn4_xerrplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ATLASVn6Vn4_yval[] = { 0.99391, 0.97918, 0.98259, 0.98755, 0.9909, 0.99359, 0.9945, 0.99601, 0.99636, 0.99602, 0.99559, 0.99242, 0.98003 };
    double ATLASVn6Vn4_yerrminus[] = { 0.11599758618178226, 0.012629725095978931, 0.007083214595083224, 0.007002367885222826, 0.005813305084029222, 0.004373746791939377, 0.005512068667932213, 0.0036425772469502963, 0.006161371600544801, 0.006271722969009394, 0.0029355852908747176, 0.004781224111041021, 0.10590297446247673 };
    double ATLASVn6Vn4_yerrplus[] = { 0.11599758618178226, 0.012629725095978931, 0.007083214595083224, 0.007002367885222826, 0.005813305084029222, 0.004373746791939377, 0.005512068667932213, 0.0036425772469502963, 0.006161371600544801, 0.006271722969009394, 0.0029355852908747176, 0.004781224111041021, 0.10590297446247673 };
    double ATLASVn6Vn4_ystatminus[] = { 0.0262, 8.66E-4, 5.73E-4, 4.16E-4, 1.96E-4, 1.81E-4, 1.51E-4, 1.37E-4, 1.3E-4, 1.47E-4, 1.81E-4, 6.98E-4, 0.0138 };
    double ATLASVn6Vn4_ystatplus[] = { 0.0262, 8.66E-4, 5.73E-4, 4.16E-4, 1.96E-4, 1.81E-4, 1.51E-4, 1.37E-4, 1.3E-4, 1.47E-4, 1.81E-4, 6.98E-4, 0.0138 };
    int ATLASVn6Vn4_numpoints = 13;
    double ATLASVn6Vn4_ysysminus[ATLASVn6Vn4_numpoints];
    double ATLASVn6Vn4_ysysplus[ATLASVn6Vn4_numpoints];
    double ATLASVn6Vn4_xsyswidth[ATLASVn6Vn4_numpoints];
    for(int i = 0; i < ATLASVn6Vn4_numpoints; i++){
      ATLASVn6Vn4_ysysminus[i] = sqrt( pow(ATLASVn6Vn4_yerrminus[i], 2) - pow(ATLASVn6Vn4_ystatminus[i], 2) );
      ATLASVn6Vn4_ysysplus[i] = sqrt( pow(ATLASVn6Vn4_yerrplus[i], 2) - pow(ATLASVn6Vn4_ystatplus[i], 2) );
      ATLASVn6Vn4_xsyswidth[i] = 5.;
    }

    TGraphAsymmErrors * grvn6vn4Ratio_ATLASNpart = new TGraphAsymmErrors(ATLASVn6Vn4_numpoints, ATLASVn6Vn4_xval, ATLASVn6Vn4_yval, ATLASVn6Vn4_xerrminus, ATLASVn6Vn4_xerrplus, ATLASVn6Vn4_ystatminus, ATLASVn6Vn4_ystatplus);
    grvn6vn4Ratio_ATLASNpart->SetLineColor(1);
    grvn6vn4Ratio_ATLASNpart->SetMarkerColor(1);
    grvn6vn4Ratio_ATLASNpart->SetMarkerStyle(20);
    grvn6vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "N_{part}");
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
    grvn6vn4Ratio_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);

    TGraphAsymmErrors * grvn6vn4RatioSys_ATLASNpart = new TGraphAsymmErrors(ATLASVn6Vn4_numpoints, ATLASVn6Vn4_xval, ATLASVn6Vn4_yval, ATLASVn6Vn4_xsyswidth, ATLASVn6Vn4_xsyswidth, ATLASVn6Vn4_ysysminus, ATLASVn6Vn4_ysysplus);
    grvn6vn4RatioSys_ATLASNpart->SetLineColor(18);
    grvn6vn4RatioSys_ATLASNpart->SetMarkerColor(18);
    grvn6vn4RatioSys_ATLASNpart->SetMarkerStyle(21);
    grvn6vn4RatioSys_ATLASNpart->SetFillColor(18);
    grvn6vn4RatioSys_ATLASNpart->GetXaxis()->SetTitle( "N_{part}");
    grvn6vn4RatioSys_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{6} / v_{%i}{4}", norder_, norder_) );
    grvn6vn4RatioSys_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);

    //-- vn{8} / vn{4}
    double ATLASVn8Vn4_xval[] = { 22.6, 46.1, 59.9, 76.1, 95.0, 117.0, 142.0, 170.0, 203.0, 239.5, 281.9, 330.3, 372.533 };
    double ATLASVn8Vn4_xerrminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ATLASVn8Vn4_xerrplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ATLASVn8Vn4_yval[] = { 0.94991, 0.97382, 0.97875, 0.98429, 0.98796, 0.99096, 0.99222, 0.99354, 0.99451, 0.99202, 0.99349, 0.98923, 0.97844 };
    double ATLASVn8Vn4_yerrminus[] = { 0.0558578553114958, 0.021625993618791254, 0.02300806121775583, 0.02560347796687005, 0.026000784680466855, 0.024800667813589215, 0.02280050666103716, 0.024900399213667237, 0.019500481250471744, 0.03870026051591901, 0.020300879808520612, 0.027108182454749708, 0.08442866811693761 };
    double ATLASVn8Vn4_yerrplus[] = { 0.0558578553114958, 0.021625993618791254, 0.02300806121775583, 0.02560347796687005, 0.026000784680466855, 0.024800667813589215, 0.02280050666103716, 0.024900399213667237, 0.019500481250471744, 0.03870026051591901, 0.020300879808520612, 0.027108182454749708, 0.08442866811693761 };
    double ATLASVn8Vn4_ystatminus[] = { 0.0103, 0.00106, 6.09E-4, 4.22E-4, 2.02E-4, 1.82E-4, 1.52E-4, 1.41E-4, 1.37E-4, 1.42E-4, 1.89E-4, 6.66E-4, 0.0118 };
    double ATLASVn8Vn4_ystatplus[] = { 0.0103, 0.00106, 6.09E-4, 4.22E-4, 2.02E-4, 1.82E-4, 1.52E-4, 1.41E-4, 1.37E-4, 1.42E-4, 1.89E-4, 6.66E-4, 0.0118 };
    int ATLASVn8Vn4_numpoints = 13;
    double ATLASVn8Vn4_ysysminus[ATLASVn8Vn4_numpoints];
    double ATLASVn8Vn4_ysysplus[ATLASVn8Vn4_numpoints];
    double ATLASVn8Vn4_xsyswidth[ATLASVn8Vn4_numpoints];
    for(int i = 0; i < ATLASVn8Vn4_numpoints; i++){
      ATLASVn8Vn4_ysysminus[i] = sqrt( pow(ATLASVn8Vn4_yerrminus[i], 2) - pow(ATLASVn8Vn4_ystatminus[i], 2) );
      ATLASVn8Vn4_ysysplus[i] = sqrt( pow(ATLASVn8Vn4_yerrplus[i], 2) - pow(ATLASVn8Vn4_ystatplus[i], 2) );
      ATLASVn8Vn4_xsyswidth[i] = 5.;
    }

    TGraphAsymmErrors * grvn8vn4Ratio_ATLASNpart = new TGraphAsymmErrors(ATLASVn8Vn4_numpoints, ATLASVn8Vn4_xval, ATLASVn8Vn4_yval, ATLASVn8Vn4_xerrminus, ATLASVn8Vn4_xerrplus, ATLASVn8Vn4_ystatminus, ATLASVn8Vn4_ystatplus);
    grvn8vn4Ratio_ATLASNpart->SetLineColor(1);
    grvn8vn4Ratio_ATLASNpart->SetMarkerColor(1);
    grvn8vn4Ratio_ATLASNpart->SetMarkerStyle(21);
    grvn8vn4Ratio_ATLASNpart->GetXaxis()->SetTitle( "N_{part}");
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
    grvn8vn4Ratio_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);

    TGraphAsymmErrors * grvn8vn4RatioSys_ATLASNpart = new TGraphAsymmErrors(ATLASVn8Vn4_numpoints, ATLASVn8Vn4_xval, ATLASVn8Vn4_yval, ATLASVn8Vn4_xsyswidth, ATLASVn8Vn4_xsyswidth, ATLASVn8Vn4_ysysminus, ATLASVn8Vn4_ysysplus);
    grvn8vn4RatioSys_ATLASNpart->SetLineColor(18);
    grvn8vn4RatioSys_ATLASNpart->SetMarkerColor(18);
    grvn8vn4RatioSys_ATLASNpart->SetMarkerStyle(21);
    grvn8vn4RatioSys_ATLASNpart->SetFillColor(18);
    grvn8vn4RatioSys_ATLASNpart->GetXaxis()->SetTitle( "N_{part}");
    grvn8vn4RatioSys_ATLASNpart->GetYaxis()->SetTitle( Form("v_{%i}{8} / v_{%i}{4}", norder_, norder_) );
    grvn8vn4RatioSys_ATLASNpart->GetYaxis()->SetRangeUser(ratioMin, ratioMax);

    TLegend * leg64 = new TLegend(0.2051, 0.1973, 0.8430, 0.3897);
    leg64->SetFillStyle(0);
    leg64->SetBorderSize(0);
    leg64->AddEntry(grVn6Vn4Npart_Default,    "Default CMS", "lp");
    leg64->AddEntry(grVn6Vn4Npart_5Sigma,     "5#sigma CMS", "lp");
    leg64->AddEntry(grVn6Vn4Npart_4Sigma,     "4#sigma CMS", "lp");
    leg64->AddEntry(grVn6Vn4Npart_3Sigma,     "3#sigma CMS", "lp");
    leg64->AddEntry(grvn6vn4Ratio_ATLASNpart, "ATLAS, 0.5 < p_{T} < 20.0 GeV, |#eta| < 2.5", "lp");

    TLegend * leg84 = new TLegend(0.2051, 0.1973, 0.8430, 0.3897);
    leg84->SetFillStyle(0);
    leg84->SetBorderSize(0);
    leg84->AddEntry(grVn8Vn4Npart_Default,    "Default CMS", "lp");
    leg84->AddEntry(grVn8Vn4Npart_5Sigma,     "5#sigma CMS", "lp");
    leg84->AddEntry(grVn8Vn4Npart_4Sigma,     "4#sigma CMS", "lp");
    leg84->AddEntry(grVn8Vn4Npart_3Sigma,     "3#sigma CMS", "lp");
    leg84->AddEntry(grvn8vn4Ratio_ATLASNpart, "ATLAS, 0.5 < p_{T} < 20.0 GeV, |#eta| < 2.5", "lp");

    TCanvas * cATLASComp = new TCanvas("cATLASComp_Vn6Vn4", "cATLASComp_Vn6Vn4", 1000, 500);
    cATLASComp->Divide(2,1);
    cATLASComp->cd(1);
    grvn6vn4RatioSys_ATLASNpart->Draw("apE2");
    grvn6vn4Ratio_ATLASNpart->Draw("psame");
    grVn6Vn4Npart_Default->Draw("psame");
    grVn6Vn4Npart_5Sigma->Draw("psame");
    grVn6Vn4Npart_4Sigma->Draw("psame");
    grVn6Vn4Npart_3Sigma->Draw("psame");
    leg64->Draw("same");
    cATLASComp->cd(2);
    grvn8vn4RatioSys_ATLASNpart->Draw("apE2");
    grvn8vn4Ratio_ATLASNpart->Draw("psame");
    grVn8Vn4Npart_Default->Draw("psame");
    grVn8Vn4Npart_5Sigma->Draw("psame");
    grVn8Vn4Npart_4Sigma->Draw("psame");
    grVn8Vn4Npart_3Sigma->Draw("psame");
    leg84->Draw("same");
    cATLASComp->SaveAs("cATLASComp.pdf");

  }

}

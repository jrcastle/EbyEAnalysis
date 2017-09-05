#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace ebyese;

void SmoothSys(){

  const int norder_ = 2;
  bool propRespUncert = 1;

  TFile * fOut;
  TH1D * SmoothSysTotVn2;
  TH1D * SmoothSysTotVn4;
  TH1D * SmoothSysTotVn6;
  TH1D * SmoothSysTotVn8;
  TH1D * SmoothSysTotVn6Vn4;
  TH1D * SmoothSysTotVn8Vn4;
  TH1D * SmoothSysTotVn8Vn6;
  TH1D * SmoothSysTotVn46_Vn68;
  TH1D * SmoothSysTotG1E;
  TH1D * SmoothSysTotKn;
  TH1D * SmoothSysTotAlpha;
  TH1D * SmoothSysTotE0;

  TFile * fSysEllp;

  //-- SysReg
  TFile * fSysReg;
  TGraphErrors * grSysRegVn2;
  TGraphErrors * grSysRegVn4;
  TGraphErrors * grSysRegVn6;
  TGraphErrors * grSysRegVn8;
  TGraphErrors * grSysRegVn6Vn4;
  TGraphErrors * grSysRegVn8Vn4;
  TGraphErrors * grSysRegVn8Vn6;
  TGraphErrors * grSysRegVn46_Vn68;
  TGraphErrors * grSysRegG1E;
  TGraphErrors * grSysRegKn;
  TGraphErrors * grSysRegAlpha;
  TGraphErrors * grSysRegE0;

  double sysRegVn2[NCENT];
  double sysRegVn4[NCENT];
  double sysRegVn6[NCENT];
  double sysRegVn8[NCENT];
  double sysRegVn6Vn4[NCENT];
  double sysRegVn8Vn4[NCENT];
  double sysRegVn8Vn6[NCENT];
  double sysRegVn46_Vn68[NCENT];
  double sysRegG1E[NCENT];
  double sysRegKn[NCENT];
  double sysRegAlpha[NCENT];
  double sysRegE0[NCENT];

  double sysRegVn2e[NCENT];
  double sysRegVn4e[NCENT];
  double sysRegVn6e[NCENT];
  double sysRegVn8e[NCENT];
  double sysRegVn6Vn4e[NCENT];
  double sysRegVn8Vn4e[NCENT];
  double sysRegVn8Vn6e[NCENT];
  double sysRegVn46_Vn68e[NCENT];
  double sysRegG1Ee[NCENT];
  double sysRegKne[NCENT];
  double sysRegAlphae[NCENT];
  double sysRegE0e[NCENT];

  TGraphErrors * grTransSysRegVn2;
  TGraphErrors * grTransSysRegVn4;
  TGraphErrors * grTransSysRegVn6;
  TGraphErrors * grTransSysRegVn8;
  TGraphErrors * grTransSysRegVn6Vn4;
  TGraphErrors * grTransSysRegVn8Vn4;
  TGraphErrors * grTransSysRegVn8Vn6;
  TGraphErrors * grTransSysRegVn46_Vn68;
  TGraphErrors * grTransSysRegG1E;
  TGraphErrors * grTransSysRegKn;
  TGraphErrors * grTransSysRegAlpha;
  TGraphErrors * grTransSysRegE0;

  TF1 * fSmoothSysRegVn2;
  TF1 * fSmoothSysRegVn4;
  TF1 * fSmoothSysRegVn6;
  TF1 * fSmoothSysRegVn8;
  TF1 * fSmoothSysRegVn6Vn4;
  TF1 * fSmoothSysRegVn8Vn4;
  TF1 * fSmoothSysRegVn8Vn6;
  TF1 * fSmoothSysRegVn46_Vn68;
  TF1 * fSmoothSysRegG1E;
  TF1 * fSmoothSysRegKn;
  TF1 * fSmoothSysRegAlpha;
  TF1 * fSmoothSysRegE0;


  //-- SysResp
  TFile * fSysResp;
  TGraphErrors * grSysRespVn2;
  TGraphErrors * grSysRespVn4;
  TGraphErrors * grSysRespVn6;
  TGraphErrors * grSysRespVn8;
  TGraphErrors * grSysRespVn6Vn4;
  TGraphErrors * grSysRespVn8Vn4;
  TGraphErrors * grSysRespVn8Vn6;
  TGraphErrors * grSysRespVn46_Vn68;
  TGraphErrors * grSysRespG1E;
  TGraphErrors * grSysRespKn;
  TGraphErrors * grSysRespAlpha;
  TGraphErrors * grSysRespE0;

  double sysRespVn2[NCENT];
  double sysRespVn4[NCENT];
  double sysRespVn6[NCENT];
  double sysRespVn8[NCENT];
  double sysRespVn6Vn4[NCENT];
  double sysRespVn8Vn4[NCENT];
  double sysRespVn8Vn6[NCENT];
  double sysRespVn46_Vn68[NCENT];
  double sysRespG1E[NCENT];
  double sysRespKn[NCENT];
  double sysRespAlpha[NCENT];
  double sysRespE0[NCENT];

  double sysRespVn2e[NCENT];
  double sysRespVn4e[NCENT];
  double sysRespVn6e[NCENT];
  double sysRespVn8e[NCENT];
  double sysRespVn6Vn4e[NCENT];
  double sysRespVn8Vn4e[NCENT];
  double sysRespVn8Vn6e[NCENT];
  double sysRespVn46_Vn68e[NCENT];
  double sysRespG1Ee[NCENT];
  double sysRespKne[NCENT];
  double sysRespAlphae[NCENT];
  double sysRespE0e[NCENT];

  TGraphErrors * grTransSysRespVn2;
  TGraphErrors * grTransSysRespVn4;
  TGraphErrors * grTransSysRespVn6;
  TGraphErrors * grTransSysRespVn8;
  TGraphErrors * grTransSysRespVn6Vn4;
  TGraphErrors * grTransSysRespVn8Vn4;
  TGraphErrors * grTransSysRespVn8Vn6;
  TGraphErrors * grTransSysRespVn46_Vn68;
  TGraphErrors * grTransSysRespG1E;
  TGraphErrors * grTransSysRespKn;
  TGraphErrors * grTransSysRespAlpha;
  TGraphErrors * grTransSysRespE0;

  TF1 * fSmoothSysRespVn2;
  TF1 * fSmoothSysRespVn4;
  TF1 * fSmoothSysRespVn6;
  TF1 * fSmoothSysRespVn8;
  TF1 * fSmoothSysRespVn6Vn4;
  TF1 * fSmoothSysRespVn8Vn4;
  TF1 * fSmoothSysRespVn8Vn6;
  TF1 * fSmoothSysRespVn46_Vn68;
  TF1 * fSmoothSysRespG1E;
  TF1 * fSmoothSysRespKn;
  TF1 * fSmoothSysRespAlpha;
  TF1 * fSmoothSysRespE0;

  //-- SysNewCC
  TFile * fSysNewCC;
  TGraphErrors * grSysNewCCVn2;
  TGraphErrors * grSysNewCCVn4;
  TGraphErrors * grSysNewCCVn6;
  TGraphErrors * grSysNewCCVn8;
  TGraphErrors * grSysNewCCVn6Vn4;
  TGraphErrors * grSysNewCCVn8Vn4;
  TGraphErrors * grSysNewCCVn8Vn6;
  TGraphErrors * grSysNewCCVn46_Vn68;
  TGraphErrors * grSysNewCCG1E;
  TGraphErrors * grSysNewCCKn;
  TGraphErrors * grSysNewCCAlpha;
  TGraphErrors * grSysNewCCE0;

  double sysNewCCVn2[NCENT];
  double sysNewCCVn4[NCENT];
  double sysNewCCVn6[NCENT];
  double sysNewCCVn8[NCENT];
  double sysNewCCVn6Vn4[NCENT];
  double sysNewCCVn8Vn4[NCENT];
  double sysNewCCVn8Vn6[NCENT];
  double sysNewCCVn46_Vn68[NCENT];
  double sysNewCCG1E[NCENT];
  double sysNewCCKn[NCENT];
  double sysNewCCAlpha[NCENT];
  double sysNewCCE0[NCENT];

  double sysNewCCVn2e[NCENT];
  double sysNewCCVn4e[NCENT];
  double sysNewCCVn6e[NCENT];
  double sysNewCCVn8e[NCENT];
  double sysNewCCVn6Vn4e[NCENT];
  double sysNewCCVn8Vn4e[NCENT];
  double sysNewCCVn8Vn6e[NCENT];
  double sysNewCCVn46_Vn68e[NCENT];
  double sysNewCCG1Ee[NCENT];
  double sysNewCCKne[NCENT];
  double sysNewCCAlphae[NCENT];
  double sysNewCCE0e[NCENT];

  TGraphErrors * grTransSysNewCCVn2;
  TGraphErrors * grTransSysNewCCVn4;
  TGraphErrors * grTransSysNewCCVn6;
  TGraphErrors * grTransSysNewCCVn8;
  TGraphErrors * grTransSysNewCCVn6Vn4;
  TGraphErrors * grTransSysNewCCVn8Vn4;
  TGraphErrors * grTransSysNewCCVn8Vn6;
  TGraphErrors * grTransSysNewCCVn46_Vn68;
  TGraphErrors * grTransSysNewCCG1E;
  TGraphErrors * grTransSysNewCCKn;
  TGraphErrors * grTransSysNewCCAlpha;
  TGraphErrors * grTransSysNewCCE0;

  TF1 * fSmoothSysNewCCVn2;
  TF1 * fSmoothSysNewCCVn4;
  TF1 * fSmoothSysNewCCVn6;
  TF1 * fSmoothSysNewCCVn8;
  TF1 * fSmoothSysNewCCVn6Vn4;
  TF1 * fSmoothSysNewCCVn8Vn4;
  TF1 * fSmoothSysNewCCVn8Vn6;
  TF1 * fSmoothSysNewCCVn46_Vn68;
  TF1 * fSmoothSysNewCCG1E;
  TF1 * fSmoothSysNewCCKn;
  TF1 * fSmoothSysNewCCAlpha;
  TF1 * fSmoothSysNewCCE0;


  //-- SysTkQ
  TFile * fSysTkQ;

  //-- Loose
  TGraphErrors * grSysTkQLVn2;
  TGraphErrors * grSysTkQLVn4;
  TGraphErrors * grSysTkQLVn6;
  TGraphErrors * grSysTkQLVn8;
  TGraphErrors * grSysTkQLVn6Vn4;
  TGraphErrors * grSysTkQLVn8Vn4;
  TGraphErrors * grSysTkQLVn8Vn6;
  TGraphErrors * grSysTkQLVn46_Vn68;
  TGraphErrors * grSysTkQLG1E;
  TGraphErrors * grSysTkQLKn;
  TGraphErrors * grSysTkQLAlpha;
  TGraphErrors * grSysTkQLE0;

  //-- Tight
  TGraphErrors * grSysTkQTVn2;
  TGraphErrors * grSysTkQTVn4;
  TGraphErrors * grSysTkQTVn6;
  TGraphErrors * grSysTkQTVn8;
  TGraphErrors * grSysTkQTVn6Vn4;
  TGraphErrors * grSysTkQTVn8Vn4;
  TGraphErrors * grSysTkQTVn8Vn6;
  TGraphErrors * grSysTkQTVn46_Vn68;
  TGraphErrors * grSysTkQTG1E;
  TGraphErrors * grSysTkQTKn;
  TGraphErrors * grSysTkQTAlpha;
  TGraphErrors * grSysTkQTE0;

  //-- Merge
  double sysTkQVn2[NCENT];
  double sysTkQVn4[NCENT];
  double sysTkQVn6[NCENT];
  double sysTkQVn8[NCENT];
  double sysTkQVn6Vn4[NCENT];
  double sysTkQVn8Vn4[NCENT];
  double sysTkQVn8Vn6[NCENT];
  double sysTkQVn46_Vn68[NCENT];
  double sysTkQG1E[NCENT];
  double sysTkQKn[NCENT];
  double sysTkQAlpha[NCENT];
  double sysTkQE0[NCENT];

  double sysTkQVn2e[NCENT];
  double sysTkQVn4e[NCENT];
  double sysTkQVn6e[NCENT];
  double sysTkQVn8e[NCENT];
  double sysTkQVn6Vn4e[NCENT];
  double sysTkQVn8Vn4e[NCENT];
  double sysTkQVn8Vn6e[NCENT];
  double sysTkQVn46_Vn68e[NCENT];
  double sysTkQG1Ee[NCENT];
  double sysTkQKne[NCENT];
  double sysTkQAlphae[NCENT];
  double sysTkQE0e[NCENT];

  TGraphErrors * grTransSysTkQVn2;
  TGraphErrors * grTransSysTkQVn4;
  TGraphErrors * grTransSysTkQVn6;
  TGraphErrors * grTransSysTkQVn8;
  TGraphErrors * grTransSysTkQVn6Vn4;
  TGraphErrors * grTransSysTkQVn8Vn4;
  TGraphErrors * grTransSysTkQVn8Vn6;
  TGraphErrors * grTransSysTkQVn46_Vn68;
  TGraphErrors * grTransSysTkQG1E;
  TGraphErrors * grTransSysTkQKn;
  TGraphErrors * grTransSysTkQAlpha;
  TGraphErrors * grTransSysTkQE0;

  TF1 * fSmoothSysTkQVn2;
  TF1 * fSmoothSysTkQVn4;
  TF1 * fSmoothSysTkQVn6;
  TF1 * fSmoothSysTkQVn8;
  TF1 * fSmoothSysTkQVn6Vn4;
  TF1 * fSmoothSysTkQVn8Vn4;
  TF1 * fSmoothSysTkQVn8Vn6;
  TF1 * fSmoothSysTkQVn46_Vn68;
  TF1 * fSmoothSysTkQG1E;
  TF1 * fSmoothSysTkQKn;
  TF1 * fSmoothSysTkQAlpha;
  TF1 * fSmoothSysTkQE0;

  //-- SysVtx
  TFile * fSysVtx;

  //-- Vtx3
  TGraphErrors * grSysVtx3Vn2;
  TGraphErrors * grSysVtx3Vn4;
  TGraphErrors * grSysVtx3Vn6;
  TGraphErrors * grSysVtx3Vn8;
  TGraphErrors * grSysVtx3Vn6Vn4;
  TGraphErrors * grSysVtx3Vn8Vn4;
  TGraphErrors * grSysVtx3Vn8Vn6;
  TGraphErrors * grSysVtx3Vn46_Vn68;
  TGraphErrors * grSysVtx3G1E;
  TGraphErrors * grSysVtx3Kn;
  TGraphErrors * grSysVtx3Alpha;
  TGraphErrors * grSysVtx3E0;

  //-- Vtx3_15
  TGraphErrors * grSysVtx3_15Vn2;
  TGraphErrors * grSysVtx3_15Vn4;
  TGraphErrors * grSysVtx3_15Vn6;
  TGraphErrors * grSysVtx3_15Vn8;
  TGraphErrors * grSysVtx3_15Vn6Vn4;
  TGraphErrors * grSysVtx3_15Vn8Vn4;
  TGraphErrors * grSysVtx3_15Vn8Vn6;
  TGraphErrors * grSysVtx3_15Vn46_Vn68;
  TGraphErrors * grSysVtx3_15G1E;
  TGraphErrors * grSysVtx3_15Kn;
  TGraphErrors * grSysVtx3_15Alpha;
  TGraphErrors * grSysVtx3_15E0;

  //-- Merge
  double sysVtxVn2[NCENT];
  double sysVtxVn4[NCENT];
  double sysVtxVn6[NCENT];
  double sysVtxVn8[NCENT];
  double sysVtxVn6Vn4[NCENT];
  double sysVtxVn8Vn4[NCENT];
  double sysVtxVn8Vn6[NCENT];
  double sysVtxVn46_Vn68[NCENT];
  double sysVtxG1E[NCENT];
  double sysVtxKn[NCENT];
  double sysVtxAlpha[NCENT];
  double sysVtxE0[NCENT];

  double sysVtxVn2e[NCENT];
  double sysVtxVn4e[NCENT];
  double sysVtxVn6e[NCENT];
  double sysVtxVn8e[NCENT];
  double sysVtxVn6Vn4e[NCENT];
  double sysVtxVn8Vn4e[NCENT];
  double sysVtxVn8Vn6e[NCENT];
  double sysVtxVn46_Vn68e[NCENT];
  double sysVtxG1Ee[NCENT];
  double sysVtxKne[NCENT];
  double sysVtxAlphae[NCENT];
  double sysVtxE0e[NCENT];

  TGraphErrors * grTransSysVtxVn2;
  TGraphErrors * grTransSysVtxVn4;
  TGraphErrors * grTransSysVtxVn6;
  TGraphErrors * grTransSysVtxVn8;
  TGraphErrors * grTransSysVtxVn6Vn4;
  TGraphErrors * grTransSysVtxVn8Vn4;
  TGraphErrors * grTransSysVtxVn8Vn6;
  TGraphErrors * grTransSysVtxVn46_Vn68;
  TGraphErrors * grTransSysVtxG1E;
  TGraphErrors * grTransSysVtxKn;
  TGraphErrors * grTransSysVtxAlpha;
  TGraphErrors * grTransSysVtxE0;

  TF1 * fSmoothSysVtxVn2;
  TF1 * fSmoothSysVtxVn4;
  TF1 * fSmoothSysVtxVn6;
  TF1 * fSmoothSysVtxVn8;
  TF1 * fSmoothSysVtxVn6Vn4;
  TF1 * fSmoothSysVtxVn8Vn4;
  TF1 * fSmoothSysVtxVn8Vn6;
  TF1 * fSmoothSysVtxVn46_Vn68;
  TF1 * fSmoothSysVtxG1E;
  TF1 * fSmoothSysVtxKn;
  TF1 * fSmoothSysVtxAlpha;
  TF1 * fSmoothSysVtxE0;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Initialize output file
  fOut = new TFile("SmoothSysTot.root", "recreate");
  fOut->cd();
  SmoothSysTotVn2       = new TH1D("SmoothSysTotVn2",       "SmoothSysTotVn2",       NCENT, centbinsDefault);
  SmoothSysTotVn4       = new TH1D("SmoothSysTotVn4",       "SmoothSysTotVn4",       NCENT, centbinsDefault);
  SmoothSysTotVn6       = new TH1D("SmoothSysTotVn6",       "SmoothSysTotVn6",       NCENT, centbinsDefault);
  SmoothSysTotVn8       = new TH1D("SmoothSysTotVn8",       "SmoothSysTotVn8",       NCENT, centbinsDefault);
  SmoothSysTotVn6Vn4    = new TH1D("SmoothSysTotVn6Vn4",    "SmoothSysTotVn6Vn4",    NCENT, centbinsDefault);
  SmoothSysTotVn8Vn4    = new TH1D("SmoothSysTotVn8Vn4",    "SmoothSysTotVn8Vn4",    NCENT, centbinsDefault);
  SmoothSysTotVn8Vn6    = new TH1D("SmoothSysTotVn8Vn6",    "SmoothSysTotVn8Vn6",    NCENT, centbinsDefault);
  SmoothSysTotVn46_Vn68 = new TH1D("SmoothSysTotVn46_Vn68", "SmoothSysTotVn46_Vn68", NCENT, centbinsDefault);
  SmoothSysTotG1E       = new TH1D("SmoothSysTotG1E",       "SmoothSysTotG1E",       NCENT, centbinsDefault);
  SmoothSysTotKn        = new TH1D("SmoothSysTotKn",        "SmoothSysTotKn",        NCENT, centbinsDefault);
  SmoothSysTotAlpha     = new TH1D("SmoothSysTotAlpha",     "SmoothSysTotAlpha",     NCENT, centbinsDefault);
  SmoothSysTotE0        = new TH1D("SmoothSysTotE0",        "SmoothSysTotE0",        NCENT, centbinsDefault);

  //-- Get Graphs
  //-- SysReg
  fSysReg = new TFile("chi2Cutoff/SysReg.root");
  grSysRegVn2       = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn2");
  grSysRegVn4       = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn4");
  grSysRegVn6       = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn6");
  grSysRegVn8       = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn8");
  grSysRegVn6Vn4    = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn6vn4");
  grSysRegVn8Vn4    = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn8vn4");
  grSysRegVn8Vn6    = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn8vn6");
  grSysRegVn46_Vn68 = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_vn46_vn68");
  grSysRegG1E       = (TGraphErrors*) fSysReg->Get("grRatioChi2_21_gamma1exp");

  //-- SysResp
  fSysResp = new TFile("responseElements/SysResp.root");
  grSysRespVn2       = (TGraphErrors*) fSysResp->Get("grVn2DoSys_RatioToNominal");
  grSysRespVn4       = (TGraphErrors*) fSysResp->Get("grVn4DoSys_RatioToNominal");
  grSysRespVn6       = (TGraphErrors*) fSysResp->Get("grVn6DoSys_RatioToNominal");
  grSysRespVn8       = (TGraphErrors*) fSysResp->Get("grVn8DoSys_RatioToNominal");
  grSysRespVn6Vn4    = (TGraphErrors*) fSysResp->Get("grVn6Vn4DoSys_RatioToNominal");
  grSysRespVn8Vn4    = (TGraphErrors*) fSysResp->Get("grVn8Vn4DoSys_RatioToNominal");
  grSysRespVn8Vn6    = (TGraphErrors*) fSysResp->Get("grVn8Vn6DoSys_RatioToNominal");
  grSysRespVn46_Vn68 = (TGraphErrors*) fSysResp->Get("grVn46_Vn68DoSys_RatioToNominal");
  grSysRespG1E       = (TGraphErrors*) fSysResp->Get("grGamma1ExpDoSys_RatioToNominal");

  //-- SysNewCC
  fSysNewCC = new TFile("clusCompatTune/SysNewCC.root");
  grSysNewCCVn2       = (TGraphErrors*) fSysNewCC->Get("grVn22pct_RatioToDefault");
  grSysNewCCVn4       = (TGraphErrors*) fSysNewCC->Get("grVn42pct_RatioToDefault");
  grSysNewCCVn6       = (TGraphErrors*) fSysNewCC->Get("grVn62pct_RatioToDefault");
  grSysNewCCVn8       = (TGraphErrors*) fSysNewCC->Get("grVn82pct_RatioToDefault");
  grSysNewCCVn6Vn4    = (TGraphErrors*) fSysNewCC->Get("grVn6Vn42pct_RatioToDefault");
  grSysNewCCVn8Vn4    = (TGraphErrors*) fSysNewCC->Get("grVn8Vn42pct_RatioToDefault");
  grSysNewCCVn8Vn6    = (TGraphErrors*) fSysNewCC->Get("grVn8Vn62pct_RatioToDefault");
  grSysNewCCVn46_Vn68 = (TGraphErrors*) fSysNewCC->Get("grVn46_Vn682pct_RatioToDefault");
  grSysNewCCG1E       = (TGraphErrors*) fSysNewCC->Get("grG1E2pct_RatioToDefault");

  //-- SysTkQ
  fSysTkQ = new TFile("tkQuality/SysTkQuality.root");
  grSysTkQLVn2       = (TGraphErrors*) fSysTkQ->Get("grVn2Loose_RatiotoNominal");
  grSysTkQLVn4       = (TGraphErrors*) fSysTkQ->Get("grVn4Loose_RatiotoNominal");
  grSysTkQLVn6       = (TGraphErrors*) fSysTkQ->Get("grVn6Loose_RatiotoNominal");
  grSysTkQLVn8       = (TGraphErrors*) fSysTkQ->Get("grVn8Loose_RatiotoNominal");
  grSysTkQLVn6Vn4    = (TGraphErrors*) fSysTkQ->Get("grVn6Vn4Loose_RatiotoNominal");
  grSysTkQLVn8Vn4    = (TGraphErrors*) fSysTkQ->Get("grVn8Vn4Loose_RatiotoNominal");
  grSysTkQLVn8Vn6    = (TGraphErrors*) fSysTkQ->Get("grVn8Vn6Loose_RatiotoNominal");
  grSysTkQLVn46_Vn68 = (TGraphErrors*) fSysTkQ->Get("grVn46_Vn68Loose_RatiotoNominal");
  grSysTkQLG1E       = (TGraphErrors*) fSysTkQ->Get("grGamma1ExpLoose_RatiotoNominal");

  grSysTkQTVn2       = (TGraphErrors*) fSysTkQ->Get("grVn2Tight_RatiotoNominal");
  grSysTkQTVn4       = (TGraphErrors*) fSysTkQ->Get("grVn4Tight_RatiotoNominal");
  grSysTkQTVn6       = (TGraphErrors*) fSysTkQ->Get("grVn6Tight_RatiotoNominal");
  grSysTkQTVn8       = (TGraphErrors*) fSysTkQ->Get("grVn8Tight_RatiotoNominal");
  grSysTkQTVn6Vn4    = (TGraphErrors*) fSysTkQ->Get("grVn6Vn4Tight_RatiotoNominal");
  grSysTkQTVn8Vn4    = (TGraphErrors*) fSysTkQ->Get("grVn8Vn4Tight_RatiotoNominal");
  grSysTkQTVn8Vn6    = (TGraphErrors*) fSysTkQ->Get("grVn8Vn6Tight_RatiotoNominal");
  grSysTkQTVn46_Vn68 = (TGraphErrors*) fSysTkQ->Get("grVn46_Vn68Tight_RatiotoNominal");
  grSysTkQTG1E       = (TGraphErrors*) fSysTkQ->Get("grGamma1ExpTight_RatiotoNominal");

  //-- SysVtx
  fSysVtx = new TFile("vtxCut/SysVtx.root");
  grSysVtx3Vn2       = (TGraphErrors*) fSysVtx->Get("grVn2Vtx3_RatiotoVtx15");
  grSysVtx3Vn4       = (TGraphErrors*) fSysVtx->Get("grVn4Vtx3_RatiotoVtx15");
  grSysVtx3Vn6       = (TGraphErrors*) fSysVtx->Get("grVn6Vtx3_RatiotoVtx15");
  grSysVtx3Vn8       = (TGraphErrors*) fSysVtx->Get("grVn8Vtx3_RatiotoVtx15");
  grSysVtx3Vn6Vn4    = (TGraphErrors*) fSysVtx->Get("grVn6Vn4Vtx3_RatiotoVtx15");
  grSysVtx3Vn8Vn4    = (TGraphErrors*) fSysVtx->Get("grVn8Vn4Vtx3_RatiotoVtx15");
  grSysVtx3Vn8Vn6    = (TGraphErrors*) fSysVtx->Get("grVn8Vn6Vtx3_RatiotoVtx15");
  grSysVtx3Vn46_Vn68 = (TGraphErrors*) fSysVtx->Get("grVn46_Vn68Vtx3_RatiotoVtx15");
  grSysVtx3G1E       = (TGraphErrors*) fSysVtx->Get("grGamma1ExpVtx3_RatiotoVtx15");

  grSysVtx3_15Vn2       = (TGraphErrors*) fSysVtx->Get("grVn2Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn4       = (TGraphErrors*) fSysVtx->Get("grVn4Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn6       = (TGraphErrors*) fSysVtx->Get("grVn6Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn8       = (TGraphErrors*) fSysVtx->Get("grVn8Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn6Vn4    = (TGraphErrors*) fSysVtx->Get("grVn6Vn4Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn8Vn4    = (TGraphErrors*) fSysVtx->Get("grVn8Vn4Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn8Vn6    = (TGraphErrors*) fSysVtx->Get("grVn8Vn6Vtx3_15_RatiotoVtx15");
  grSysVtx3_15Vn46_Vn68 = (TGraphErrors*) fSysVtx->Get("grVn46_Vn68Vtx3_15_RatiotoVtx15");
  grSysVtx3_15G1E       = (TGraphErrors*) fSysVtx->Get("grGamma1ExpVtx3_15_RatiotoVtx15");

  //-- Ellp
  fSysEllp = new TFile("SysEllpParms.root");
  grSysRegKn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_SysReg");
  grSysRegAlpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_SysReg");
  grSysRegE0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_SysReg");

  grSysRespKn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_SysResp");
  grSysRespAlpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_SysResp");
  grSysRespE0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_SysResp");

  grSysNewCCKn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_SysNewCC");
  grSysNewCCAlpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_SysNewCC");
  grSysNewCCE0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_SysNewCC");

  grSysTkQLKn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_LooseSysTkQ");
  grSysTkQLAlpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_LooseSysTkQ");
  grSysTkQLE0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_LooseSysTkQ");

  grSysTkQTKn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_TightSysTkQ");
  grSysTkQTAlpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_TightSysTkQ");
  grSysTkQTE0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_TightSysTkQ");

  grSysVtx3Kn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_Vtx3SysVtx");
  grSysVtx3Alpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_Vtx3SysVtx");
  grSysVtx3E0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_Vtx3SysVtx");

  grSysVtx3_15Kn    = (TGraphErrors*) fSysEllp->Get("grFitKn_RatioToNom_Vtx3_15SysVtx");
  grSysVtx3_15Alpha = (TGraphErrors*) fSysEllp->Get("grFitAlpha_RatioToNom_Vtx3_15SysVtx");
  grSysVtx3_15E0    = (TGraphErrors*) fSysEllp->Get("grFitE0_RatioToNom_Vtx3_15SysVtx");

  //-- Ave Sys
  double aveSys_00_10_Vn2 = 0.;
  double aveSys_10_20_Vn2 = 0.;
  double aveSys_20_30_Vn2 = 0.;
  double aveSys_30_50_Vn2 = 0.;
  double aveSys_50_60_Vn2 = 0.;

  double aveSys_00_10_Vn4 = 0.;
  double aveSys_10_20_Vn4 = 0.;
  double aveSys_20_30_Vn4 = 0.;
  double aveSys_30_50_Vn4 = 0.;
  double aveSys_50_60_Vn4 = 0.;

  double aveSys_00_10_Vn6 = 0.;
  double aveSys_10_20_Vn6 = 0.;
  double aveSys_20_30_Vn6 = 0.;
  double aveSys_30_50_Vn6 = 0.;
  double aveSys_50_60_Vn6 = 0.;

  double aveSys_00_10_Vn8 = 0.;
  double aveSys_10_20_Vn8 = 0.;
  double aveSys_20_30_Vn8 = 0.;
  double aveSys_30_50_Vn8 = 0.;
  double aveSys_50_60_Vn8 = 0.;

  double aveSys_00_10_Vn6Vn4 = 0.;
  double aveSys_10_20_Vn6Vn4 = 0.;
  double aveSys_20_30_Vn6Vn4 = 0.;
  double aveSys_30_50_Vn6Vn4 = 0.;
  double aveSys_50_60_Vn6Vn4 = 0.;

  double aveSys_00_10_Vn8Vn4 = 0.;
  double aveSys_10_20_Vn8Vn4 = 0.;
  double aveSys_20_30_Vn8Vn4 = 0.;
  double aveSys_30_50_Vn8Vn4 = 0.;
  double aveSys_50_60_Vn8Vn4 = 0.;

  double aveSys_00_10_Vn8Vn6 = 0.;
  double aveSys_10_20_Vn8Vn6 = 0.;
  double aveSys_20_30_Vn8Vn6 = 0.;
  double aveSys_30_50_Vn8Vn6 = 0.;
  double aveSys_50_60_Vn8Vn6 = 0.;

  double aveSys_00_10_Vn46_Vn68 = 0.;
  double aveSys_10_20_Vn46_Vn68 = 0.;
  double aveSys_20_30_Vn46_Vn68 = 0.;
  double aveSys_30_50_Vn46_Vn68 = 0.;
  double aveSys_50_60_Vn46_Vn68 = 0.;

  double aveSys_00_10_G1E = 0.;
  double aveSys_10_20_G1E = 0.;
  double aveSys_20_30_G1E = 0.;
  double aveSys_30_50_G1E = 0.;
  double aveSys_50_60_G1E = 0.;

  double aveSys_15_25_Kn = 0.;
  double aveSys_25_50_Kn = 0.;
  double aveSys_50_60_Kn = 0.;

  double aveSys_15_25_Alpha = 0.;
  double aveSys_25_50_Alpha = 0.;
  double aveSys_50_60_Alpha = 0.;

  double aveSys_15_25_E0 = 0.;
  double aveSys_25_50_E0 = 0.;
  double aveSys_50_60_E0 = 0.;

  //-- Start grabbing data points and translating them...
  for(int icent = 0; icent < NCENT; icent++){

    //-- SysReg
    sysRegVn2[icent]       = fabs( 1. - grSysRegVn2->GetY()[icent] );
    sysRegVn4[icent]       = fabs( 1. - grSysRegVn4->GetY()[icent] );
    sysRegVn6[icent]       = fabs( 1. - grSysRegVn6->GetY()[icent] );
    sysRegVn8[icent]       = fabs( 1. - grSysRegVn8->GetY()[icent] );
    sysRegVn6Vn4[icent]    = fabs( 1. - grSysRegVn6Vn4->GetY()[icent] );
    sysRegVn8Vn4[icent]    = fabs( 1. - grSysRegVn8Vn4->GetY()[icent] );
    sysRegVn8Vn6[icent]    = fabs( 1. - grSysRegVn8Vn6->GetY()[icent] );
    sysRegVn46_Vn68[icent] = fabs( 1. - grSysRegVn46_Vn68->GetY()[icent] );
    sysRegG1E[icent]       = fabs( 1. - grSysRegG1E->GetY()[icent] );
    sysRegKn[icent]        = fabs( 1. - grSysRegKn->GetY()[icent] );
    sysRegAlpha[icent]     = fabs( 1. - grSysRegAlpha->GetY()[icent] );
    sysRegE0[icent]        = fabs( 1. - grSysRegE0->GetY()[icent] );

    sysRegVn2e[icent]       = grSysRegVn2->GetErrorY(icent);
    sysRegVn4e[icent]       = grSysRegVn4->GetErrorY(icent);
    sysRegVn6e[icent]       = grSysRegVn6->GetErrorY(icent);
    sysRegVn8e[icent]       = grSysRegVn8->GetErrorY(icent);
    sysRegVn6Vn4e[icent]    = grSysRegVn6Vn4->GetErrorY(icent);
    sysRegVn8Vn4e[icent]    = grSysRegVn8Vn4->GetErrorY(icent);
    sysRegVn8Vn6e[icent]    = grSysRegVn8Vn6->GetErrorY(icent);
    sysRegVn46_Vn68e[icent] = grSysRegVn46_Vn68->GetErrorY(icent);
    sysRegG1Ee[icent]       = grSysRegG1E->GetErrorY(icent);
    sysRegKne[icent]        = grSysRegKn->GetErrorY(icent);
    sysRegAlphae[icent]     = grSysRegAlpha->GetErrorY(icent);
    sysRegE0e[icent]        = grSysRegE0->GetErrorY(icent);

    //-- SysResp
    sysRespVn2[icent]       = fabs( 1. - grSysRespVn2->GetY()[icent] );
    sysRespVn4[icent]       = fabs( 1. - grSysRespVn4->GetY()[icent] );
    sysRespVn6[icent]       = fabs( 1. - grSysRespVn6->GetY()[icent] );
    sysRespVn8[icent]       = fabs( 1. - grSysRespVn8->GetY()[icent] );
    sysRespVn6Vn4[icent]    = fabs( 1. - grSysRespVn6Vn4->GetY()[icent] );
    sysRespVn8Vn4[icent]    = fabs( 1. - grSysRespVn8Vn4->GetY()[icent] );
    sysRespVn8Vn6[icent]    = fabs( 1. - grSysRespVn8Vn6->GetY()[icent] );
    sysRespVn46_Vn68[icent] = fabs( 1. - grSysRespVn46_Vn68->GetY()[icent] );
    sysRespG1E[icent]       = fabs( 1. - grSysRespG1E->GetY()[icent] );
    sysRespKn[icent]        = fabs( 1. - grSysRespKn->GetY()[icent] );
    sysRespAlpha[icent]     = fabs( 1. - grSysRespAlpha->GetY()[icent] );
    sysRespE0[icent]        = fabs( 1. - grSysRespE0->GetY()[icent] );

    sysRespVn2e[icent]       = grSysRespVn2->GetErrorY(icent);
    sysRespVn4e[icent]       = grSysRespVn4->GetErrorY(icent);
    sysRespVn6e[icent]       = grSysRespVn6->GetErrorY(icent);
    sysRespVn8e[icent]       = grSysRespVn8->GetErrorY(icent);
    sysRespVn6Vn4e[icent]    = grSysRespVn6Vn4->GetErrorY(icent);
    sysRespVn8Vn4e[icent]    = grSysRespVn8Vn4->GetErrorY(icent);
    sysRespVn8Vn6e[icent]    = grSysRespVn8Vn6->GetErrorY(icent);
    sysRespVn46_Vn68e[icent] = grSysRespVn46_Vn68->GetErrorY(icent);
    sysRespG1Ee[icent]       = grSysRespG1E->GetErrorY(icent);
    sysRespKne[icent]        = grSysRespKn->GetErrorY(icent);
    sysRespAlphae[icent]     = grSysRespAlpha->GetErrorY(icent);
    sysRespE0e[icent]        = grSysRespE0->GetErrorY(icent);

    //-- SysNewCC
    sysNewCCVn2[icent]       = fabs( 1. - grSysNewCCVn2->GetY()[icent] );
    sysNewCCVn4[icent]       = fabs( 1. - grSysNewCCVn4->GetY()[icent] );
    sysNewCCVn6[icent]       = fabs( 1. - grSysNewCCVn6->GetY()[icent] );
    sysNewCCVn8[icent]       = fabs( 1. - grSysNewCCVn8->GetY()[icent] );
    sysNewCCVn6Vn4[icent]    = fabs( 1. - grSysNewCCVn6Vn4->GetY()[icent] );
    sysNewCCVn8Vn4[icent]    = fabs( 1. - grSysNewCCVn8Vn4->GetY()[icent] );
    sysNewCCVn8Vn6[icent]    = fabs( 1. - grSysNewCCVn8Vn6->GetY()[icent] );
    sysNewCCVn46_Vn68[icent] = fabs( 1. - grSysNewCCVn46_Vn68->GetY()[icent] );
    sysNewCCG1E[icent]       = fabs( 1. - grSysNewCCG1E->GetY()[icent] );
    sysNewCCKn[icent]        = fabs( 1. - grSysNewCCKn->GetY()[icent] );
    sysNewCCAlpha[icent]     = fabs( 1. - grSysNewCCAlpha->GetY()[icent] );
    sysNewCCE0[icent]        = fabs( 1. - grSysNewCCE0->GetY()[icent] );

    sysNewCCVn2e[icent]       = grSysNewCCVn2->GetErrorY(icent);
    sysNewCCVn4e[icent]       = grSysNewCCVn4->GetErrorY(icent);
    sysNewCCVn6e[icent]       = grSysNewCCVn6->GetErrorY(icent);
    sysNewCCVn8e[icent]       = grSysNewCCVn8->GetErrorY(icent);
    sysNewCCVn6Vn4e[icent]    = grSysNewCCVn6Vn4->GetErrorY(icent);
    sysNewCCVn8Vn4e[icent]    = grSysNewCCVn8Vn4->GetErrorY(icent);
    sysNewCCVn8Vn6e[icent]    = grSysNewCCVn8Vn6->GetErrorY(icent);
    sysNewCCVn46_Vn68e[icent] = grSysNewCCVn46_Vn68->GetErrorY(icent);
    sysNewCCG1Ee[icent]       = grSysNewCCG1E->GetErrorY(icent);
    sysNewCCKne[icent]        = grSysNewCCKn->GetErrorY(icent);
    sysNewCCAlphae[icent]     = grSysNewCCAlpha->GetErrorY(icent);
    sysNewCCE0e[icent]        = grSysNewCCE0->GetErrorY(icent);

    //-- TkQL
    double tkQLVn2       = fabs( 1. - grSysTkQLVn2->GetY()[icent] );
    double tkQLVn4       = fabs( 1. - grSysTkQLVn4->GetY()[icent] );
    double tkQLVn6       = fabs( 1. - grSysTkQLVn6->GetY()[icent] );
    double tkQLVn8       = fabs( 1. - grSysTkQLVn8->GetY()[icent] );
    double tkQLVn6Vn4    = fabs( 1. - grSysTkQLVn6Vn4->GetY()[icent] );
    double tkQLVn8Vn4    = fabs( 1. - grSysTkQLVn8Vn4->GetY()[icent] );
    double tkQLVn8Vn6    = fabs( 1. - grSysTkQLVn8Vn6->GetY()[icent] );
    double tkQLVn46_Vn68 = fabs( 1. - grSysTkQLVn46_Vn68->GetY()[icent] );
    double tkQLG1E       = fabs( 1. - grSysTkQLG1E->GetY()[icent] );
    double tkQLKn        = fabs( 1. - grSysTkQLKn->GetY()[icent] );
    double tkQLAlpha     = fabs( 1. - grSysTkQLAlpha->GetY()[icent] );
    double tkQLE0        = fabs( 1. - grSysTkQLE0->GetY()[icent] );

    double tkQLVn2e       = grSysTkQLVn2->GetErrorY(icent);
    double tkQLVn4e       = grSysTkQLVn4->GetErrorY(icent);
    double tkQLVn6e       = grSysTkQLVn6->GetErrorY(icent);
    double tkQLVn8e       = grSysTkQLVn8->GetErrorY(icent);
    double tkQLVn6Vn4e    = grSysTkQLVn6Vn4->GetErrorY(icent);
    double tkQLVn8Vn4e    = grSysTkQLVn8Vn4->GetErrorY(icent);
    double tkQLVn8Vn6e    = grSysTkQLVn8Vn6->GetErrorY(icent);
    double tkQLVn46_Vn68e = grSysTkQLVn46_Vn68->GetErrorY(icent);
    double tkQLG1Ee       = grSysTkQLG1E->GetErrorY(icent);
    double tkQLKne        = grSysTkQLKn->GetErrorY(icent);
    double tkQLAlphae     = grSysTkQLAlpha->GetErrorY(icent);
    double tkQLE0e        = grSysTkQLE0->GetErrorY(icent);

    //-- TkQT
    double tkQTVn2       = fabs( 1. - grSysTkQTVn2->GetY()[icent] );
    double tkQTVn4       = fabs( 1. - grSysTkQTVn4->GetY()[icent] );
    double tkQTVn6       = fabs( 1. - grSysTkQTVn6->GetY()[icent] );
    double tkQTVn8       = fabs( 1. - grSysTkQTVn8->GetY()[icent] );
    double tkQTVn6Vn4    = fabs( 1. - grSysTkQTVn6Vn4->GetY()[icent] );
    double tkQTVn8Vn4    = fabs( 1. - grSysTkQTVn8Vn4->GetY()[icent] );
    double tkQTVn8Vn6    = fabs( 1. - grSysTkQTVn8Vn6->GetY()[icent] );
    double tkQTVn46_Vn68 = fabs( 1. - grSysTkQTVn46_Vn68->GetY()[icent] );
    double tkQTG1E       = fabs( 1. - grSysTkQTG1E->GetY()[icent] );
    double tkQTKn        = fabs( 1. - grSysTkQTKn->GetY()[icent] );
    double tkQTAlpha     = fabs( 1. - grSysTkQTAlpha->GetY()[icent] );
    double tkQTE0        = fabs( 1. - grSysTkQTE0->GetY()[icent] );

    double tkQTVn2e       = grSysTkQTVn2->GetErrorY(icent);
    double tkQTVn4e       = grSysTkQTVn4->GetErrorY(icent);
    double tkQTVn6e       = grSysTkQTVn6->GetErrorY(icent);
    double tkQTVn8e       = grSysTkQTVn8->GetErrorY(icent);
    double tkQTVn6Vn4e    = grSysTkQTVn6Vn4->GetErrorY(icent);
    double tkQTVn8Vn4e    = grSysTkQTVn8Vn4->GetErrorY(icent);
    double tkQTVn8Vn6e    = grSysTkQTVn8Vn6->GetErrorY(icent);
    double tkQTVn46_Vn68e = grSysTkQTVn46_Vn68->GetErrorY(icent);
    double tkQTG1Ee       = grSysTkQTG1E->GetErrorY(icent);
    double tkQTKne        = grSysTkQTKn->GetErrorY(icent);
    double tkQTAlphae     = grSysTkQTAlpha->GetErrorY(icent);
    double tkQTE0e        = grSysTkQTE0->GetErrorY(icent);

    //-- Merge
    if(tkQLVn2 > tkQTVn2){
      sysTkQVn2[icent]  = tkQLVn2;
      sysTkQVn2e[icent] = tkQLVn2e;
    }
    else{
      sysTkQVn2[icent]  = tkQTVn2;
      sysTkQVn2e[icent] = tkQTVn2e;
    }

    if(tkQLVn4 > tkQTVn4){
      sysTkQVn4[icent]  = tkQLVn4;
      sysTkQVn4e[icent] = tkQLVn4e;
    }
    else{
      sysTkQVn4[icent]  = tkQTVn4;
      sysTkQVn4e[icent] = tkQTVn4e;
    }

    if(tkQLVn6 > tkQTVn6){
      sysTkQVn6[icent]  = tkQLVn6;
      sysTkQVn6e[icent] = tkQLVn6e;
    }
    else{
      sysTkQVn6[icent]  = tkQTVn6;
      sysTkQVn6e[icent] = tkQTVn6e;
    }

    if(tkQLVn8 > tkQTVn8){
      sysTkQVn8[icent]  = tkQLVn8;
      sysTkQVn8e[icent] = tkQLVn8e;
    }
    else{
      sysTkQVn8[icent]  = tkQTVn8;
      sysTkQVn8e[icent] = tkQTVn8e;
    }

    if(tkQLVn6Vn4 > tkQTVn6Vn4){
      sysTkQVn6Vn4[icent]  = tkQLVn6Vn4;
      sysTkQVn6Vn4e[icent] = tkQLVn6Vn4e;
    }
    else{
      sysTkQVn6Vn4[icent]  = tkQTVn6Vn4;
      sysTkQVn6Vn4e[icent] = tkQTVn6Vn4e;
    }

    if(tkQLVn8Vn4 > tkQTVn8Vn4){
      sysTkQVn8Vn4[icent]  = tkQLVn8Vn4;
      sysTkQVn8Vn4e[icent] = tkQLVn8Vn4e;
    }
    else{
      sysTkQVn8Vn4[icent]  = tkQTVn8Vn4;
      sysTkQVn8Vn4e[icent] = tkQTVn8Vn4e;
    }

    if(tkQLVn8Vn6 > tkQTVn8Vn6){
      sysTkQVn8Vn6[icent]  = tkQLVn8Vn6;
      sysTkQVn8Vn6e[icent] = tkQLVn8Vn6e;
    }
    else{
      sysTkQVn8Vn6[icent]  = tkQTVn8Vn6;
      sysTkQVn8Vn6e[icent] = tkQTVn8Vn6e;
    }

    if(tkQLVn46_Vn68 > tkQTVn46_Vn68){
      sysTkQVn46_Vn68[icent]  = tkQLVn46_Vn68;
      sysTkQVn46_Vn68e[icent] = tkQLVn46_Vn68e;
    }
    else{
      sysTkQVn46_Vn68[icent]  = tkQTVn46_Vn68;
      sysTkQVn46_Vn68e[icent] = tkQTVn46_Vn68e;
    }

    if(tkQLG1E > tkQTG1E){
      sysTkQG1E[icent]  = tkQLG1E;
      sysTkQG1Ee[icent] = tkQLG1Ee;
    }
    else{
      sysTkQG1E[icent]  = tkQTG1E;
      sysTkQG1Ee[icent] = tkQTG1Ee;
    }

    if(tkQLKn > tkQTKn){
      sysTkQKn[icent]  = tkQLKn;
      sysTkQKne[icent] = tkQLKne;
    }
    else{
      sysTkQKn[icent]  = tkQTKn;
      sysTkQKne[icent] = tkQTKne;
    }

    if(tkQLAlpha > tkQTAlpha){
      sysTkQAlpha[icent]  = tkQLAlpha;
      sysTkQAlphae[icent] = tkQLAlphae;
    }
    else{
      sysTkQAlpha[icent]  = tkQTAlpha;
      sysTkQAlphae[icent] = tkQTAlphae;
    }

    if(tkQLE0 > tkQTE0){
      sysTkQE0[icent]  = tkQLE0;
      sysTkQE0e[icent] = tkQLE0e;
    }
    else{
      sysTkQE0[icent]  = tkQTE0;
      sysTkQE0e[icent] = tkQTE0e;
    } 

    //-- Vtx3
    double Vtx3Vn2       = fabs( 1. - grSysVtx3Vn2->GetY()[icent] );
    double Vtx3Vn4       = fabs( 1. - grSysVtx3Vn4->GetY()[icent] );
    double Vtx3Vn6       = fabs( 1. - grSysVtx3Vn6->GetY()[icent] );
    double Vtx3Vn8       = fabs( 1. - grSysVtx3Vn8->GetY()[icent] );
    double Vtx3Vn6Vn4    = fabs( 1. - grSysVtx3Vn6Vn4->GetY()[icent] );
    double Vtx3Vn8Vn4    = fabs( 1. - grSysVtx3Vn8Vn4->GetY()[icent] );
    double Vtx3Vn8Vn6    = fabs( 1. - grSysVtx3Vn8Vn6->GetY()[icent] );
    double Vtx3Vn46_Vn68 = fabs( 1. - grSysVtx3Vn46_Vn68->GetY()[icent] );
    double Vtx3G1E       = fabs( 1. - grSysVtx3G1E->GetY()[icent] );
    double Vtx3Kn        = fabs( 1. - grSysVtx3Kn->GetY()[icent] );
    double Vtx3Alpha     = fabs( 1. - grSysVtx3Alpha->GetY()[icent] );
    double Vtx3E0        = fabs( 1. - grSysVtx3E0->GetY()[icent] );

    double Vtx3Vn2e       = grSysVtx3Vn2->GetErrorY(icent);
    double Vtx3Vn4e       = grSysVtx3Vn4->GetErrorY(icent);
    double Vtx3Vn6e       = grSysVtx3Vn6->GetErrorY(icent);
    double Vtx3Vn8e       = grSysVtx3Vn8->GetErrorY(icent);
    double Vtx3Vn6Vn4e    = grSysVtx3Vn6Vn4->GetErrorY(icent);
    double Vtx3Vn8Vn4e    = grSysVtx3Vn8Vn4->GetErrorY(icent);
    double Vtx3Vn8Vn6e    = grSysVtx3Vn8Vn6->GetErrorY(icent);
    double Vtx3Vn46_Vn68e = grSysVtx3Vn46_Vn68->GetErrorY(icent);
    double Vtx3G1Ee       = grSysVtx3G1E->GetErrorY(icent);
    double Vtx3Kne        = grSysVtx3Kn->GetErrorY(icent);
    double Vtx3Alphae     = grSysVtx3Alpha->GetErrorY(icent);
    double Vtx3E0e        = grSysVtx3E0->GetErrorY(icent);

    //-- Vtx3_15
    double Vtx3_15Vn2       = fabs( 1. - grSysVtx3_15Vn2->GetY()[icent] );
    double Vtx3_15Vn4       = fabs( 1. - grSysVtx3_15Vn4->GetY()[icent] );
    double Vtx3_15Vn6       = fabs( 1. - grSysVtx3_15Vn6->GetY()[icent] );
    double Vtx3_15Vn8       = fabs( 1. - grSysVtx3_15Vn8->GetY()[icent] );
    double Vtx3_15Vn6Vn4    = fabs( 1. - grSysVtx3_15Vn6Vn4->GetY()[icent] );
    double Vtx3_15Vn8Vn4    = fabs( 1. - grSysVtx3_15Vn8Vn4->GetY()[icent] );
    double Vtx3_15Vn8Vn6    = fabs( 1. - grSysVtx3_15Vn8Vn6->GetY()[icent] );
    double Vtx3_15Vn46_Vn68 = fabs( 1. - grSysVtx3_15Vn46_Vn68->GetY()[icent] );
    double Vtx3_15G1E       = fabs( 1. - grSysVtx3_15G1E->GetY()[icent] );
    double Vtx3_15Kn        = fabs( 1. - grSysVtx3_15Kn->GetY()[icent] );
    double Vtx3_15Alpha     = fabs( 1. - grSysVtx3_15Alpha->GetY()[icent] );
    double Vtx3_15E0        = fabs( 1. - grSysVtx3_15E0->GetY()[icent] );

    double Vtx3_15Vn2e       = grSysVtx3_15Vn2->GetErrorY(icent);
    double Vtx3_15Vn4e       = grSysVtx3_15Vn4->GetErrorY(icent);
    double Vtx3_15Vn6e       = grSysVtx3_15Vn6->GetErrorY(icent);
    double Vtx3_15Vn8e       = grSysVtx3_15Vn8->GetErrorY(icent);
    double Vtx3_15Vn6Vn4e    = grSysVtx3_15Vn6Vn4->GetErrorY(icent);
    double Vtx3_15Vn8Vn4e    = grSysVtx3_15Vn8Vn4->GetErrorY(icent);
    double Vtx3_15Vn8Vn6e    = grSysVtx3_15Vn8Vn6->GetErrorY(icent);
    double Vtx3_15Vn46_Vn68e = grSysVtx3_15Vn46_Vn68->GetErrorY(icent);
    double Vtx3_15G1Ee       = grSysVtx3_15G1E->GetErrorY(icent);
    double Vtx3_15Kne        = grSysVtx3_15Kn->GetErrorY(icent);
    double Vtx3_15Alphae     = grSysVtx3_15Alpha->GetErrorY(icent);
    double Vtx3_15E0e        = grSysVtx3_15E0->GetErrorY(icent);

    //-- Merge
    if(Vtx3Vn2 > Vtx3_15Vn2){
      sysVtxVn2[icent]  = Vtx3Vn2;
      sysVtxVn2e[icent] = Vtx3Vn2e;
    }
    else{
      sysVtxVn2[icent]  = Vtx3_15Vn2;
      sysVtxVn2e[icent] = Vtx3_15Vn2e;
    }

    if(Vtx3Vn4 > Vtx3_15Vn4){
      sysVtxVn4[icent]  = Vtx3Vn4;
      sysVtxVn4e[icent] = Vtx3Vn4e;
    }
    else{
      sysVtxVn4[icent]  = Vtx3_15Vn4;
      sysVtxVn4e[icent] = Vtx3_15Vn4e;
    }

    if(Vtx3Vn6 > Vtx3_15Vn6){
      sysVtxVn6[icent]  = Vtx3Vn6;
      sysVtxVn6e[icent] = Vtx3Vn6e;
    }
    else{
      sysVtxVn6[icent]  = Vtx3_15Vn6;
      sysVtxVn6e[icent] = Vtx3_15Vn6e;
    }

    if(Vtx3Vn8 > Vtx3_15Vn8){
      sysVtxVn8[icent]  = Vtx3Vn8;
      sysVtxVn8e[icent] = Vtx3Vn8e;
    }
    else{
      sysVtxVn8[icent]  = Vtx3_15Vn8;
      sysVtxVn8e[icent] = Vtx3_15Vn8e;
    }

    if(Vtx3Vn6Vn4 > Vtx3_15Vn6Vn4){
      sysVtxVn6Vn4[icent]  = Vtx3Vn6Vn4;
      sysVtxVn6Vn4e[icent] = Vtx3Vn6Vn4e;
    }
    else{
      sysVtxVn6Vn4[icent]  = Vtx3_15Vn6Vn4;
      sysVtxVn6Vn4e[icent] = Vtx3_15Vn6Vn4e;
    }

    if(Vtx3Vn8Vn4 > Vtx3_15Vn8Vn4){
      sysVtxVn8Vn4[icent]  = Vtx3Vn8Vn4;
      sysVtxVn8Vn4e[icent] = Vtx3Vn8Vn4e;
    }
    else{
      sysVtxVn8Vn4[icent]  = Vtx3_15Vn8Vn4;
      sysVtxVn8Vn4e[icent] = Vtx3_15Vn8Vn4e;
    }

    if(Vtx3Vn8Vn6 > Vtx3_15Vn8Vn6){
      sysVtxVn8Vn6[icent]  = Vtx3Vn8Vn6;
      sysVtxVn8Vn6e[icent] = Vtx3Vn8Vn6e;
    }
    else{
      sysVtxVn8Vn6[icent]  = Vtx3_15Vn8Vn6;
      sysVtxVn8Vn6e[icent] = Vtx3_15Vn8Vn6e;
    }

    if(Vtx3Vn46_Vn68 > Vtx3_15Vn46_Vn68){
      sysVtxVn46_Vn68[icent]  = Vtx3Vn46_Vn68;
      sysVtxVn46_Vn68e[icent] = Vtx3Vn46_Vn68e;
    }
    else{
      sysVtxVn46_Vn68[icent]  = Vtx3_15Vn46_Vn68;
      sysVtxVn46_Vn68e[icent] = Vtx3_15Vn46_Vn68e;
    }

    if(Vtx3G1E > Vtx3_15G1E){
      sysVtxG1E[icent]  = Vtx3G1E;
      sysVtxG1Ee[icent] = Vtx3G1Ee;
    }
    else{
      sysVtxG1E[icent]  = Vtx3_15G1E;
      sysVtxG1Ee[icent] = Vtx3_15G1Ee;
    }

    if(Vtx3Kn > Vtx3_15Kn){
      sysVtxKn[icent]  = Vtx3Kn;
      sysVtxKne[icent] = Vtx3Kne;
    }
    else{
      sysVtxKn[icent]  = Vtx3_15Kn;
      sysVtxKne[icent] = Vtx3_15Kne;
    }

    if(Vtx3Alpha > Vtx3_15Alpha){
      sysVtxAlpha[icent]  = Vtx3Alpha;
      sysVtxAlphae[icent] = Vtx3Alphae;
    }
    else{
      sysVtxAlpha[icent]  = Vtx3_15Alpha;
      sysVtxAlphae[icent] = Vtx3_15Alphae;
    }

    if(Vtx3E0 > Vtx3_15E0){
      sysVtxE0[icent]  = Vtx3E0;
      sysVtxE0e[icent] = Vtx3E0e;
    }
    else{
      sysVtxE0[icent]  = Vtx3_15E0;
      sysVtxE0e[icent] = Vtx3_15E0e;
    }

  } //-- End Cent loop

  //-- Make translated TGraphs
  //-- SysReg
  grTransSysRegVn2       = new TGraphErrors(NCENT, centBinCenter, sysRegVn2,       CERR, sysRegVn2e);
  grTransSysRegVn4       = new TGraphErrors(NCENT, centBinCenter, sysRegVn4,       CERR, sysRegVn4e);
  grTransSysRegVn6       = new TGraphErrors(NCENT, centBinCenter, sysRegVn6,       CERR, sysRegVn6e);
  grTransSysRegVn8       = new TGraphErrors(NCENT, centBinCenter, sysRegVn8,       CERR, sysRegVn8e);
  grTransSysRegVn6Vn4    = new TGraphErrors(NCENT, centBinCenter, sysRegVn6Vn4,    CERR, sysRegVn6Vn4e);
  grTransSysRegVn8Vn4    = new TGraphErrors(NCENT, centBinCenter, sysRegVn8Vn4,    CERR, sysRegVn8Vn4e);
  grTransSysRegVn8Vn6    = new TGraphErrors(NCENT, centBinCenter, sysRegVn8Vn6,    CERR, sysRegVn8Vn6e);
  grTransSysRegVn46_Vn68 = new TGraphErrors(NCENT, centBinCenter, sysRegVn46_Vn68, CERR, sysRegVn46_Vn68e);
  grTransSysRegG1E       = new TGraphErrors(NCENT, centBinCenter, sysRegG1E,       CERR, sysRegG1Ee);
  grTransSysRegKn        = new TGraphErrors(NCENT, centBinCenter, sysRegKn,        CERR, sysRegKne);
  grTransSysRegAlpha     = new TGraphErrors(NCENT, centBinCenter, sysRegAlpha,     CERR, sysRegAlphae);
  grTransSysRegE0        = new TGraphErrors(NCENT, centBinCenter, sysRegE0,        CERR, sysRegE0e);

  formatGraph(grTransSysRegVn2,       "Centrality %", 0., 1., Form("v_{%i}{2} Sys Uncert %s", norder_, "%"),                                                                  9,         20, "grTransSysRegVn2");
  formatGraph(grTransSysRegVn4,       "Centrality %", 0., 1., Form("v_{%i}{4} Sys Uncert %s", norder_, "%"),                                                                  kSpring+4, 20, "grTransSysRegVn4");
  formatGraph(grTransSysRegVn6,       "Centrality %", 0., 1., Form("v_{%i}{6} Sys Uncert %s", norder_, "%"),                                                                  6,         20, "grTransSysRegVn6");
  formatGraph(grTransSysRegVn8,       "Centrality %", 0., 1., Form("v_{%i}{8} Sys Uncert %s", norder_, "%"),                                                                  kOrange+7, 20, "grTransSysRegVn8");
  formatGraph(grTransSysRegVn6Vn4,    "Centrality %", 0., 1., Form("v_{%i}{6}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               4,         20, "grTransSysRegVn6Vn4");
  formatGraph(grTransSysRegVn8Vn4,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               kGreen+2,  20, "grTransSysRegVn8Vn4");
  formatGraph(grTransSysRegVn8Vn6,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{6} Sys Uncert %s", norder_, norder_, "%"),                                               kViolet-1, 20, "grTransSysRegVn8Vn6");
  formatGraph(grTransSysRegVn46_Vn68, "Centrality %", 0., 1., Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) Sys Uncert %s", norder_, norder_, norder_, norder_, "%"), kGray+2,   20, "grTransSysRegVn46_Vn68");
  formatGraph(grTransSysRegG1E,       "Centrality %", 0., 1., "#gamma_{1}^{exp} Sys Uncert %",                                                                                2,         20, "grTransSysRegG1E");
  formatGraph(grTransSysRegKn,        "Centrality %", 0., 1., "k_{n} Sys Uncert %",                                                                                           2,         20, "grTransSysRegKn");
  formatGraph(grTransSysRegAlpha,     "Centrality %", 0., 1., "#alpha Sys Uncert %",                                                                                          2,         20, "grTransSysRegAlpha");
  formatGraph(grTransSysRegE0,        "Centrality %", 0., 1., "#epsilon_{0} Sys Uncert %",                                                                                    2,         20, "grTransSysRegE0");


  //-- SysResp
  grTransSysRespVn2       = new TGraphErrors(NCENT, centBinCenter, sysRespVn2,       CERR, sysRespVn2e);
  grTransSysRespVn4       = new TGraphErrors(NCENT, centBinCenter, sysRespVn4,       CERR, sysRespVn4e);
  grTransSysRespVn6       = new TGraphErrors(NCENT, centBinCenter, sysRespVn6,       CERR, sysRespVn6e);
  grTransSysRespVn8       = new TGraphErrors(NCENT, centBinCenter, sysRespVn8,       CERR, sysRespVn8e);
  grTransSysRespVn6Vn4    = new TGraphErrors(NCENT, centBinCenter, sysRespVn6Vn4,    CERR, sysRespVn6Vn4e);
  grTransSysRespVn8Vn4    = new TGraphErrors(NCENT, centBinCenter, sysRespVn8Vn4,    CERR, sysRespVn8Vn4e);
  grTransSysRespVn8Vn6    = new TGraphErrors(NCENT, centBinCenter, sysRespVn8Vn6,    CERR, sysRespVn8Vn6e);
  grTransSysRespVn46_Vn68 = new TGraphErrors(NCENT, centBinCenter, sysRespVn46_Vn68, CERR, sysRespVn46_Vn68e);
  grTransSysRespG1E       = new TGraphErrors(NCENT, centBinCenter, sysRespG1E,       CERR, sysRespG1Ee);
  grTransSysRespKn        = new TGraphErrors(NCENT, centBinCenter, sysRespKn,        CERR, sysRespKne);
  grTransSysRespAlpha     = new TGraphErrors(NCENT, centBinCenter, sysRespAlpha,     CERR, sysRespAlphae);
  grTransSysRespE0        = new TGraphErrors(NCENT, centBinCenter, sysRespE0,        CERR, sysRespE0e);

  formatGraph(grTransSysRespVn2,       "Centrality %", 0., 1., Form("v_{%i}{2} Sys Uncert %s", norder_, "%"),                                                                  9,         20, "grTransSysRespVn2");
  formatGraph(grTransSysRespVn4,       "Centrality %", 0., 1., Form("v_{%i}{4} Sys Uncert %s", norder_, "%"),                                                                  kSpring+4, 20, "grTransSysRespVn4");
  formatGraph(grTransSysRespVn6,       "Centrality %", 0., 1., Form("v_{%i}{6} Sys Uncert %s", norder_, "%"),                                                                  6,         20, "grTransSysRespVn6");
  formatGraph(grTransSysRespVn8,       "Centrality %", 0., 1., Form("v_{%i}{8} Sys Uncert %s", norder_, "%"),                                                                  kOrange+7, 20, "grTransSysRespVn8");
  formatGraph(grTransSysRespVn6Vn4,    "Centrality %", 0., 1., Form("v_{%i}{6}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               4,         20, "grTransSysRespVn6Vn4");
  formatGraph(grTransSysRespVn8Vn4,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               kGreen+2,  20, "grTransSysRespVn8Vn4");
  formatGraph(grTransSysRespVn8Vn6,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{6} Sys Uncert %s", norder_, norder_, "%"),                                               kViolet-1, 20, "grTransSysRespVn8Vn6");
  formatGraph(grTransSysRespVn46_Vn68, "Centrality %", 0., 1., Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) Sys Uncert %s", norder_, norder_, norder_, norder_, "%"), kGray+2,   20, "grTransSysRespVn46_Vn68");
  formatGraph(grTransSysRespG1E,       "Centrality %", 0., 1., "#gamma_{1}^{exp} Sys Uncert %",                                                                                2,         20, "grTransSysRespG1E");
  formatGraph(grTransSysRespKn,        "Centrality %", 0., 1., "k_{n} Sys Uncert %",                                                                                           2,         20, "grTransSysRespKn");
  formatGraph(grTransSysRespAlpha,     "Centrality %", 0., 1., "#alpha Sys Uncert %",                                                                                          2,         20, "grTransSysRespAlpha");
  formatGraph(grTransSysRespE0,        "Centrality %", 0., 1., "#epsilon_{0} Sys Uncert %",                                                                                    2,         20, "grTransSysRespE0");

  //-- SysNewCC
  grTransSysNewCCVn2       = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn2,       CERR, sysNewCCVn2e);
  grTransSysNewCCVn4       = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn4,       CERR, sysNewCCVn4e);
  grTransSysNewCCVn6       = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn6,       CERR, sysNewCCVn6e);
  grTransSysNewCCVn8       = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn8,       CERR, sysNewCCVn8e);
  grTransSysNewCCVn6Vn4    = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn6Vn4,    CERR, sysNewCCVn6Vn4e);
  grTransSysNewCCVn8Vn4    = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn8Vn4,    CERR, sysNewCCVn8Vn4e);
  grTransSysNewCCVn8Vn6    = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn8Vn6,    CERR, sysNewCCVn8Vn6e);
  grTransSysNewCCVn46_Vn68 = new TGraphErrors(NCENT, centBinCenter, sysNewCCVn46_Vn68, CERR, sysNewCCVn46_Vn68e);
  grTransSysNewCCG1E       = new TGraphErrors(NCENT, centBinCenter, sysNewCCG1E,       CERR, sysNewCCG1Ee);
  grTransSysNewCCKn        = new TGraphErrors(NCENT, centBinCenter, sysNewCCKn,        CERR, sysNewCCKne);
  grTransSysNewCCAlpha     = new TGraphErrors(NCENT, centBinCenter, sysNewCCAlpha,     CERR, sysNewCCAlphae);
  grTransSysNewCCE0        = new TGraphErrors(NCENT, centBinCenter, sysNewCCE0,        CERR, sysNewCCE0e);

  formatGraph(grTransSysNewCCVn2,       "Centrality %", 0., 1., Form("v_{%i}{2} Sys Uncert %s", norder_, "%"),                                                                  9,         20, "grTransSysNewCCVn2");
  formatGraph(grTransSysNewCCVn4,       "Centrality %", 0., 1., Form("v_{%i}{4} Sys Uncert %s", norder_, "%"),                                                                  kSpring+4, 20, "grTransSysNewCCVn4");
  formatGraph(grTransSysNewCCVn6,       "Centrality %", 0., 1., Form("v_{%i}{6} Sys Uncert %s", norder_, "%"),                                                                  6,         20, "grTransSysNewCCVn6");
  formatGraph(grTransSysNewCCVn8,       "Centrality %", 0., 1., Form("v_{%i}{8} Sys Uncert %s", norder_, "%"),                                                                  kOrange+7, 20, "grTransSysNewCCVn8");
  formatGraph(grTransSysNewCCVn6Vn4,    "Centrality %", 0., 1., Form("v_{%i}{6}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               4,         20, "grTransSysNewCCVn6Vn4");
  formatGraph(grTransSysNewCCVn8Vn4,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               kGreen+2,  20, "grTransSysNewCCVn8Vn4");
  formatGraph(grTransSysNewCCVn8Vn6,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{6} Sys Uncert %s", norder_, norder_, "%"),                                               kViolet-1, 20, "grTransSysNewCCVn8Vn6");
  formatGraph(grTransSysNewCCVn46_Vn68, "Centrality %", 0., 1., Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) Sys Uncert %s", norder_, norder_, norder_, norder_, "%"), kGray+2,   20, "grTransSysNewCCVn46_Vn68");
  formatGraph(grTransSysNewCCG1E,       "Centrality %", 0., 1., "#gamma_{1}^{exp} Sys Uncert %",                                                                                2,         20, "grTransSysNewCCG1E");
  formatGraph(grTransSysNewCCKn,        "Centrality %", 0., 1., "k_{n} Sys Uncert %",                                                                                           2,         20, "grTransSysNewCCKn");
  formatGraph(grTransSysNewCCAlpha,     "Centrality %", 0., 1., "#alpha Sys Uncert %",                                                                                          2,         20, "grTransSysNewCCAlpha");
  formatGraph(grTransSysNewCCE0,        "Centrality %", 0., 1., "#epsilon_{0} Sys Uncert %",                                                                                    2,         20, "grTransSysNewCCE0");

  //-- SysTkQ
  grTransSysTkQVn2       = new TGraphErrors(NCENT, centBinCenter, sysTkQVn2,       CERR, sysTkQVn2e);
  grTransSysTkQVn4       = new TGraphErrors(NCENT, centBinCenter, sysTkQVn4,       CERR, sysTkQVn4e);
  grTransSysTkQVn6       = new TGraphErrors(NCENT, centBinCenter, sysTkQVn6,       CERR, sysTkQVn6e);
  grTransSysTkQVn8       = new TGraphErrors(NCENT, centBinCenter, sysTkQVn8,       CERR, sysTkQVn8e);
  grTransSysTkQVn6Vn4    = new TGraphErrors(NCENT, centBinCenter, sysTkQVn6Vn4,    CERR, sysTkQVn6Vn4e);
  grTransSysTkQVn8Vn4    = new TGraphErrors(NCENT, centBinCenter, sysTkQVn8Vn4,    CERR, sysTkQVn8Vn4e);
  grTransSysTkQVn8Vn6    = new TGraphErrors(NCENT, centBinCenter, sysTkQVn8Vn6,    CERR, sysTkQVn8Vn6e);
  grTransSysTkQVn46_Vn68 = new TGraphErrors(NCENT, centBinCenter, sysTkQVn46_Vn68, CERR, sysTkQVn46_Vn68e);
  grTransSysTkQG1E       = new TGraphErrors(NCENT, centBinCenter, sysTkQG1E,       CERR, sysTkQG1Ee);
  grTransSysTkQKn        = new TGraphErrors(NCENT, centBinCenter, sysTkQKn,        CERR, sysTkQKne);
  grTransSysTkQAlpha     = new TGraphErrors(NCENT, centBinCenter, sysTkQAlpha,     CERR, sysTkQAlphae);
  grTransSysTkQE0        = new TGraphErrors(NCENT, centBinCenter, sysTkQE0,        CERR, sysTkQE0e);

  formatGraph(grTransSysTkQVn2,       "Centrality %", 0., 1., Form("v_{%i}{2} Sys Uncert %s", norder_, "%"),                                                                  9,         20, "grTransSysTkQVn2");
  formatGraph(grTransSysTkQVn4,       "Centrality %", 0., 1., Form("v_{%i}{4} Sys Uncert %s", norder_, "%"),                                                                  kSpring+4, 20, "grTransSysTkQVn4");
  formatGraph(grTransSysTkQVn6,       "Centrality %", 0., 1., Form("v_{%i}{6} Sys Uncert %s", norder_, "%"),                                                                  6,         20, "grTransSysTkQVn6");
  formatGraph(grTransSysTkQVn8,       "Centrality %", 0., 1., Form("v_{%i}{8} Sys Uncert %s", norder_, "%"),                                                                  kOrange+7, 20, "grTransSysTkQVn8");
  formatGraph(grTransSysTkQVn6Vn4,    "Centrality %", 0., 1., Form("v_{%i}{6}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               4,         20, "grTransSysTkQVn6Vn4");
  formatGraph(grTransSysTkQVn8Vn4,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               kGreen+2,  20, "grTransSysTkQVn8Vn4");
  formatGraph(grTransSysTkQVn8Vn6,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{6} Sys Uncert %s", norder_, norder_, "%"),                                               kViolet-1, 20, "grTransSysTkQVn8Vn6");
  formatGraph(grTransSysTkQVn46_Vn68, "Centrality %", 0., 1., Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) Sys Uncert %s", norder_, norder_, norder_, norder_, "%"), kGray+2,   20, "grTransSysTkQVn46_Vn68");
  formatGraph(grTransSysTkQG1E,       "Centrality %", 0., 1., "#gamma_{1}^{exp} Sys Uncert %",                                                                                2,         20, "grTransSysTkQG1E");
  formatGraph(grTransSysTkQKn,        "Centrality %", 0., 1., "k_{n} Sys Uncert %",                                                                                           2,         20, "grTransSysTkQKn");
  formatGraph(grTransSysTkQAlpha,     "Centrality %", 0., 1., "#alpha Sys Uncert %",                                                                                          2,         20, "grTransSysTkQAlpha");
  formatGraph(grTransSysTkQE0,        "Centrality %", 0., 1., "#epsilon_{0} Sys Uncert %",                                                                                    2,         20, "grTransSysTkQE0");

  //-- SysVtx
  grTransSysVtxVn2       = new TGraphErrors(NCENT, centBinCenter, sysVtxVn2,       CERR, sysVtxVn2e);
  grTransSysVtxVn4       = new TGraphErrors(NCENT, centBinCenter, sysVtxVn4,       CERR, sysVtxVn4e);
  grTransSysVtxVn6       = new TGraphErrors(NCENT, centBinCenter, sysVtxVn6,       CERR, sysVtxVn6e);
  grTransSysVtxVn8       = new TGraphErrors(NCENT, centBinCenter, sysVtxVn8,       CERR, sysVtxVn8e);
  grTransSysVtxVn6Vn4    = new TGraphErrors(NCENT, centBinCenter, sysVtxVn6Vn4,    CERR, sysVtxVn6Vn4e);
  grTransSysVtxVn8Vn4    = new TGraphErrors(NCENT, centBinCenter, sysVtxVn8Vn4,    CERR, sysVtxVn8Vn4e);
  grTransSysVtxVn8Vn6    = new TGraphErrors(NCENT, centBinCenter, sysVtxVn8Vn6,    CERR, sysVtxVn8Vn6e);
  grTransSysVtxVn46_Vn68 = new TGraphErrors(NCENT, centBinCenter, sysVtxVn46_Vn68, CERR, sysVtxVn46_Vn68e);
  grTransSysVtxG1E       = new TGraphErrors(NCENT, centBinCenter, sysVtxG1E,       CERR, sysVtxG1Ee);
  grTransSysVtxKn        = new TGraphErrors(NCENT, centBinCenter, sysVtxKn,        CERR, sysVtxKne);
  grTransSysVtxAlpha     = new TGraphErrors(NCENT, centBinCenter, sysVtxAlpha,     CERR, sysVtxAlphae);
  grTransSysVtxE0        = new TGraphErrors(NCENT, centBinCenter, sysVtxE0,        CERR, sysVtxE0e);

  formatGraph(grTransSysVtxVn2,       "Centrality %", 0., 1., Form("v_{%i}{2} Sys Uncert %s", norder_, "%"),                                                                  9,         20, "grTransSysVtxVn2");
  formatGraph(grTransSysVtxVn4,       "Centrality %", 0., 1., Form("v_{%i}{4} Sys Uncert %s", norder_, "%"),                                                                  kSpring+4, 20, "grTransSysVtxVn4");
  formatGraph(grTransSysVtxVn6,       "Centrality %", 0., 1., Form("v_{%i}{6} Sys Uncert %s", norder_, "%"),                                                                  6,         20, "grTransSysVtxVn6");
  formatGraph(grTransSysVtxVn8,       "Centrality %", 0., 1., Form("v_{%i}{8} Sys Uncert %s", norder_, "%"),                                                                  kOrange+7, 20, "grTransSysVtxVn8");
  formatGraph(grTransSysVtxVn6Vn4,    "Centrality %", 0., 1., Form("v_{%i}{6}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               4,         20, "grTransSysVtxVn6Vn4");
  formatGraph(grTransSysVtxVn8Vn4,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{4} Sys Uncert %s", norder_, norder_, "%"),                                               kGreen+2,  20, "grTransSysVtxVn8Vn4");
  formatGraph(grTransSysVtxVn8Vn6,    "Centrality %", 0., 1., Form("v_{%i}{8}/v_{%i}{6} Sys Uncert %s", norder_, norder_, "%"),                                               kViolet-1, 20, "grTransSysVtxVn8Vn6");
  formatGraph(grTransSysVtxVn46_Vn68, "Centrality %", 0., 1., Form("(v_{%i}{4} - v_{%i}{6})/(v_{%i}{6} - v_{%i}{8}) Sys Uncert %s", norder_, norder_, norder_, norder_, "%"), kGray+2,   20, "grTransSysVtxVn46_Vn68");
  formatGraph(grTransSysVtxG1E,       "Centrality %", 0., 1., "#gamma_{1}^{exp} Sys Uncert %",                                                                                2,         20, "grTransSysVtxG1E");
  formatGraph(grTransSysVtxKn,        "Centrality %", 0., 1., "k_{n} Sys Uncert %",                                                                                           2,         20, "grTransSysVtxKn");
  formatGraph(grTransSysVtxAlpha,     "Centrality %", 0., 1., "#alpha Sys Uncert %",                                                                                          2,         20, "grTransSysVtxAlpha");
  formatGraph(grTransSysVtxE0,        "Centrality %", 0., 1., "#epsilon_{0} Sys Uncert %",                                                                                    2,         20, "grTransSysVtxE0");

  //-- Smooth Fits
  //-- SysReg
  fSmoothSysRegVn2       = new TF1("fSmoothSysRegVn2",       "pol5", 0,  60);
  fSmoothSysRegVn4       = new TF1("fSmoothSysRegVn4",       "pol5", 0,  60);
  fSmoothSysRegVn6       = new TF1("fSmoothSysRegVn6",       "pol5", 0,  60);
  fSmoothSysRegVn8       = new TF1("fSmoothSysRegVn8",       "pol5", 0,  60);
  fSmoothSysRegVn6Vn4    = new TF1("fSmoothSysRegVn6Vn4",    "pol5", 0,  60);
  fSmoothSysRegVn8Vn4    = new TF1("fSmoothSysRegVn8Vn4",    "pol5", 0,  60);
  fSmoothSysRegVn8Vn6    = new TF1("fSmoothSysRegVn8Vn6",    "pol5", 0,  60);
  fSmoothSysRegVn46_Vn68 = new TF1("fSmoothSysRegVn46_Vn68", "pol5", 0,  60);
  fSmoothSysRegG1E       = new TF1("fSmoothSysRegG1E",       "pol5", 0,  60);
  fSmoothSysRegKn        = new TF1("fSmoothSysRegKn",        "pol5", 15, 60);
  fSmoothSysRegAlpha     = new TF1("fSmoothSysRegAlpha",     "pol5", 15, 60);
  fSmoothSysRegE0        = new TF1("fSmoothSysRegE0",        "pol5", 15, 60);

  grTransSysRegVn2       -> Fit("fSmoothSysRegVn2",       "Q0", "", 5,  60);
  grTransSysRegVn4       -> Fit("fSmoothSysRegVn4",       "Q0", "", 5,  60);
  grTransSysRegVn6       -> Fit("fSmoothSysRegVn6",       "Q0", "", 5,  60);
  grTransSysRegVn8       -> Fit("fSmoothSysRegVn8",       "Q0", "", 5,  60);
  grTransSysRegVn6Vn4    -> Fit("fSmoothSysRegVn6Vn4",    "Q0", "", 5,  60);
  grTransSysRegVn8Vn4    -> Fit("fSmoothSysRegVn8Vn4",    "Q0", "", 5,  60);
  grTransSysRegVn8Vn6    -> Fit("fSmoothSysRegVn8Vn6",    "Q0", "", 5,  60);
  grTransSysRegVn46_Vn68 -> Fit("fSmoothSysRegVn46_Vn68", "Q0", "", 5,  60);
  grTransSysRegG1E       -> Fit("fSmoothSysRegG1E",       "Q0", "", 5,  60);
  grTransSysRegKn        -> Fit("fSmoothSysRegKn",        "Q0", "", 15, 60);
  grTransSysRegAlpha     -> Fit("fSmoothSysRegAlpha",     "Q0", "", 15, 60);
  grTransSysRegE0        -> Fit("fSmoothSysRegE0",        "Q0", "", 15, 60);

  //-- SysResp
  fSmoothSysRespVn2       = new TF1("fSmoothSysRespVn2",       "pol5", 0,  60);
  fSmoothSysRespVn4       = new TF1("fSmoothSysRespVn4",       "pol5", 0,  60);
  fSmoothSysRespVn6       = new TF1("fSmoothSysRespVn6",       "pol5", 0,  60);
  fSmoothSysRespVn8       = new TF1("fSmoothSysRespVn8",       "pol5", 0,  60);
  fSmoothSysRespVn6Vn4    = new TF1("fSmoothSysRespVn6Vn4",    "pol5", 0,  60);
  fSmoothSysRespVn8Vn4    = new TF1("fSmoothSysRespVn8Vn4",    "pol5", 0,  60);
  fSmoothSysRespVn8Vn6    = new TF1("fSmoothSysRespVn8Vn6",    "pol5", 0,  60);
  fSmoothSysRespVn46_Vn68 = new TF1("fSmoothSysRespVn46_Vn68", "pol5", 0,  60);
  fSmoothSysRespG1E       = new TF1("fSmoothSysRespG1E",       "pol5", 0,  60);
  fSmoothSysRespKn        = new TF1("fSmoothSysRespKn",        "pol5", 15, 60);
  fSmoothSysRespAlpha     = new TF1("fSmoothSysRespAlpha",     "pol5", 15, 60);
  fSmoothSysRespE0        = new TF1("fSmoothSysRespE0",        "pol5", 15, 60);

  grTransSysRespVn2       -> Fit("fSmoothSysRespVn2",       "Q0", "", 5,  60);
  grTransSysRespVn4       -> Fit("fSmoothSysRespVn4",       "Q0", "", 5,  60);
  grTransSysRespVn6       -> Fit("fSmoothSysRespVn6",       "Q0", "", 5,  60);
  grTransSysRespVn8       -> Fit("fSmoothSysRespVn8",       "Q0", "", 5,  60);
  grTransSysRespVn6Vn4    -> Fit("fSmoothSysRespVn6Vn4",    "Q0", "", 5,  60);
  grTransSysRespVn8Vn4    -> Fit("fSmoothSysRespVn8Vn4",    "Q0", "", 5,  60);
  grTransSysRespVn8Vn6    -> Fit("fSmoothSysRespVn8Vn6",    "Q0", "", 5,  60);
  grTransSysRespVn46_Vn68 -> Fit("fSmoothSysRespVn46_Vn68", "Q0", "", 5,  60);
  grTransSysRespG1E       -> Fit("fSmoothSysRespG1E",       "Q0", "", 5,  60);
  grTransSysRespKn        -> Fit("fSmoothSysRespKn",        "Q0", "", 15, 60);
  grTransSysRespAlpha     -> Fit("fSmoothSysRespAlpha",     "Q0", "", 15, 60);
  grTransSysRespE0        -> Fit("fSmoothSysRespE0",        "Q0", "", 15, 60);

  //-- SysNewCC
  fSmoothSysNewCCVn2       = new TF1("fSmoothSysNewCCVn2",       "pol5", 0,  60);
  fSmoothSysNewCCVn4       = new TF1("fSmoothSysNewCCVn4",       "pol5", 0,  60);
  fSmoothSysNewCCVn6       = new TF1("fSmoothSysNewCCVn6",       "pol5", 0,  60);
  fSmoothSysNewCCVn8       = new TF1("fSmoothSysNewCCVn8",       "pol5", 0,  60);
  fSmoothSysNewCCVn6Vn4    = new TF1("fSmoothSysNewCCVn6Vn4",    "pol5", 0,  60);
  fSmoothSysNewCCVn8Vn4    = new TF1("fSmoothSysNewCCVn8Vn4",    "pol5", 0,  60);
  fSmoothSysNewCCVn8Vn6    = new TF1("fSmoothSysNewCCVn8Vn6",    "pol5", 0,  60);
  fSmoothSysNewCCVn46_Vn68 = new TF1("fSmoothSysNewCCVn46_Vn68", "pol5", 0,  60);
  fSmoothSysNewCCG1E       = new TF1("fSmoothSysNewCCG1E",       "pol5", 0,  60);
  fSmoothSysNewCCKn        = new TF1("fSmoothSysNewCCKn",        "pol5", 15, 60);
  fSmoothSysNewCCAlpha     = new TF1("fSmoothSysNewCCAlpha",     "pol5", 15, 60);
  fSmoothSysNewCCE0        = new TF1("fSmoothSysNewCCE0",        "pol5", 15, 60);

  grTransSysNewCCVn2       -> Fit("fSmoothSysNewCCVn2",       "Q0", "", 5,  60);
  grTransSysNewCCVn4       -> Fit("fSmoothSysNewCCVn4",       "Q0", "", 5,  60);
  grTransSysNewCCVn6       -> Fit("fSmoothSysNewCCVn6",       "Q0", "", 5,  60);
  grTransSysNewCCVn8       -> Fit("fSmoothSysNewCCVn8",       "Q0", "", 5,  60);
  grTransSysNewCCVn6Vn4    -> Fit("fSmoothSysNewCCVn6Vn4",    "Q0", "", 5,  60);
  grTransSysNewCCVn8Vn4    -> Fit("fSmoothSysNewCCVn8Vn4",    "Q0", "", 5,  60);
  grTransSysNewCCVn8Vn6    -> Fit("fSmoothSysNewCCVn8Vn6",    "Q0", "", 5,  60);
  grTransSysNewCCVn46_Vn68 -> Fit("fSmoothSysNewCCVn46_Vn68", "Q0", "", 5,  60);
  grTransSysNewCCG1E       -> Fit("fSmoothSysNewCCG1E",       "Q0", "", 5,  60);
  grTransSysNewCCKn        -> Fit("fSmoothSysNewCCKn",        "Q0", "", 15, 60);
  grTransSysNewCCAlpha     -> Fit("fSmoothSysNewCCAlpha",     "Q0", "", 15, 60);
  grTransSysNewCCE0        -> Fit("fSmoothSysNewCCE0",        "Q0", "", 15, 60);

  //-- SysTkQ
  fSmoothSysTkQVn2       = new TF1("fSmoothSysTkQVn2",       "pol5", 0,  60);
  fSmoothSysTkQVn4       = new TF1("fSmoothSysTkQVn4",       "pol5", 0,  60);
  fSmoothSysTkQVn6       = new TF1("fSmoothSysTkQVn6",       "pol5", 0,  60);
  fSmoothSysTkQVn8       = new TF1("fSmoothSysTkQVn8",       "pol5", 0,  60);
  fSmoothSysTkQVn6Vn4    = new TF1("fSmoothSysTkQVn6Vn4",    "pol5", 0,  60);
  fSmoothSysTkQVn8Vn4    = new TF1("fSmoothSysTkQVn8Vn4",    "pol5", 0,  60);
  fSmoothSysTkQVn8Vn6    = new TF1("fSmoothSysTkQVn8Vn6",    "pol5", 0,  60);
  fSmoothSysTkQVn46_Vn68 = new TF1("fSmoothSysTkQVn46_Vn68", "pol5", 0,  60);
  fSmoothSysTkQG1E       = new TF1("fSmoothSysTkQG1E",       "pol5", 0,  60);
  fSmoothSysTkQKn        = new TF1("fSmoothSysTkQKn",        "pol5", 15, 60);
  fSmoothSysTkQAlpha     = new TF1("fSmoothSysTkQAlpha",     "pol5", 15, 60);
  fSmoothSysTkQE0        = new TF1("fSmoothSysTkQE0",        "pol5", 15, 60);

  grTransSysTkQVn2       -> Fit("fSmoothSysTkQVn2",       "Q0", "", 5,  60);
  grTransSysTkQVn4       -> Fit("fSmoothSysTkQVn4",       "Q0", "", 5,  60);
  grTransSysTkQVn6       -> Fit("fSmoothSysTkQVn6",       "Q0", "", 5,  60);
  grTransSysTkQVn8       -> Fit("fSmoothSysTkQVn8",       "Q0", "", 5,  60);
  grTransSysTkQVn6Vn4    -> Fit("fSmoothSysTkQVn6Vn4",    "Q0", "", 5,  60);
  grTransSysTkQVn8Vn4    -> Fit("fSmoothSysTkQVn8Vn4",    "Q0", "", 5,  60);
  grTransSysTkQVn8Vn6    -> Fit("fSmoothSysTkQVn8Vn6",    "Q0", "", 5,  60);
  grTransSysTkQVn46_Vn68 -> Fit("fSmoothSysTkQVn46_Vn68", "Q0", "", 5,  60);
  grTransSysTkQG1E       -> Fit("fSmoothSysTkQG1E",       "Q0", "", 5,  60);
  grTransSysTkQKn        -> Fit("fSmoothSysTkQKn",        "Q0", "", 15, 60);
  grTransSysTkQAlpha     -> Fit("fSmoothSysTkQAlpha",     "Q0", "", 15, 60);
  grTransSysTkQE0        -> Fit("fSmoothSysTkQE0",        "Q0", "", 15, 60);

  //-- SysVtx
  fSmoothSysVtxVn2       = new TF1("fSmoothSysVtxVn2",       "pol5", 0,  60);
  fSmoothSysVtxVn4       = new TF1("fSmoothSysVtxVn4",       "pol5", 0,  60);
  fSmoothSysVtxVn6       = new TF1("fSmoothSysVtxVn6",       "pol5", 0,  60);
  fSmoothSysVtxVn8       = new TF1("fSmoothSysVtxVn8",       "pol5", 0,  60);
  fSmoothSysVtxVn6Vn4    = new TF1("fSmoothSysVtxVn6Vn4",    "pol5", 0,  60);
  fSmoothSysVtxVn8Vn4    = new TF1("fSmoothSysVtxVn8Vn4",    "pol5", 0,  60);
  fSmoothSysVtxVn8Vn6    = new TF1("fSmoothSysVtxVn8Vn6",    "pol5", 0,  60);
  fSmoothSysVtxVn46_Vn68 = new TF1("fSmoothSysVtxVn46_Vn68", "pol5", 0,  60);
  fSmoothSysVtxG1E       = new TF1("fSmoothSysVtxG1E",       "pol5", 0,  60);
  fSmoothSysVtxKn        = new TF1("fSmoothSysVtxKn",        "pol5", 15, 60);
  fSmoothSysVtxAlpha     = new TF1("fSmoothSysVtxAlpha",     "pol5", 15, 60);
  fSmoothSysVtxE0        = new TF1("fSmoothSysVtxE0",        "pol5", 15, 60);

  grTransSysVtxVn2       -> Fit("fSmoothSysVtxVn2",       "Q0", "", 5,  60);
  grTransSysVtxVn4       -> Fit("fSmoothSysVtxVn4",       "Q0", "", 5,  60);
  grTransSysVtxVn6       -> Fit("fSmoothSysVtxVn6",       "Q0", "", 5,  60);
  grTransSysVtxVn8       -> Fit("fSmoothSysVtxVn8",       "Q0", "", 5,  60);
  grTransSysVtxVn6Vn4    -> Fit("fSmoothSysVtxVn6Vn4",    "Q0", "", 5,  60);
  grTransSysVtxVn8Vn4    -> Fit("fSmoothSysVtxVn8Vn4",    "Q0", "", 5,  60);
  grTransSysVtxVn8Vn6    -> Fit("fSmoothSysVtxVn8Vn6",    "Q0", "", 5,  60);
  grTransSysVtxVn46_Vn68 -> Fit("fSmoothSysVtxVn46_Vn68", "Q0", "", 5,  60);
  grTransSysVtxG1E       -> Fit("fSmoothSysVtxG1E",       "Q0", "", 5,  60);
  grTransSysVtxKn        -> Fit("fSmoothSysVtxKn",        "Q0", "", 15, 60);
  grTransSysVtxAlpha     -> Fit("fSmoothSysVtxAlpha",     "Q0", "", 15, 60);
  grTransSysVtxE0        -> Fit("fSmoothSysVtxE0",        "Q0", "", 15, 60);

  for(int icent = 0; icent < NCENT; icent++){

    double sysRegVn2       = fSmoothSysRegVn2->Eval( centBinCenter[icent] );
    double sysRegVn4       = fSmoothSysRegVn4->Eval( centBinCenter[icent] );
    double sysRegVn6       = fSmoothSysRegVn6->Eval( centBinCenter[icent] );
    double sysRegVn8       = fSmoothSysRegVn8->Eval( centBinCenter[icent] );
    double sysRegVn6Vn4    = fSmoothSysRegVn6Vn4->Eval( centBinCenter[icent] );
    double sysRegVn8Vn4    = fSmoothSysRegVn8Vn4->Eval( centBinCenter[icent] );
    double sysRegVn8Vn6    = fSmoothSysRegVn8Vn6->Eval( centBinCenter[icent] );
    double sysRegVn46_Vn68 = fSmoothSysRegVn46_Vn68->Eval( centBinCenter[icent] );
    double sysRegG1E       = fSmoothSysRegG1E->Eval( centBinCenter[icent] );

    double sysRespVn2       = fSmoothSysRespVn2->Eval( centBinCenter[icent] );
    double sysRespVn4       = fSmoothSysRespVn4->Eval( centBinCenter[icent] );
    double sysRespVn6       = fSmoothSysRespVn6->Eval( centBinCenter[icent] );
    double sysRespVn8       = fSmoothSysRespVn8->Eval( centBinCenter[icent] );
    double sysRespVn6Vn4    = fSmoothSysRespVn6Vn4->Eval( centBinCenter[icent] );
    double sysRespVn8Vn4    = fSmoothSysRespVn8Vn4->Eval( centBinCenter[icent] );
    double sysRespVn8Vn6    = fSmoothSysRespVn8Vn6->Eval( centBinCenter[icent] );
    double sysRespVn46_Vn68 = fSmoothSysRespVn46_Vn68->Eval( centBinCenter[icent] );
    double sysRespG1E       = fSmoothSysRespG1E->Eval( centBinCenter[icent] );

    if(!propRespUncert){
      sysRespVn2       = 0;
      sysRespVn4       = 0;
      sysRespVn6       = 0;
      sysRespVn8       = 0;
      sysRespVn6Vn4    = 0;
      sysRespVn8Vn4    = 0;
      sysRespVn8Vn6    = 0;
      sysRespVn46_Vn68 = 0;
      sysRespG1E       = 0;
    }

    double sysNewCCVn2       = fSmoothSysNewCCVn2->Eval( centBinCenter[icent] );
    double sysNewCCVn4       = fSmoothSysNewCCVn4->Eval( centBinCenter[icent] );
    double sysNewCCVn6       = fSmoothSysNewCCVn6->Eval( centBinCenter[icent] );
    double sysNewCCVn8       = fSmoothSysNewCCVn8->Eval( centBinCenter[icent] );
    double sysNewCCVn6Vn4    = fSmoothSysNewCCVn6Vn4->Eval( centBinCenter[icent] );
    double sysNewCCVn8Vn4    = fSmoothSysNewCCVn8Vn4->Eval( centBinCenter[icent] );
    double sysNewCCVn8Vn6    = fSmoothSysNewCCVn8Vn6->Eval( centBinCenter[icent] );
    double sysNewCCVn46_Vn68 = fSmoothSysNewCCVn46_Vn68->Eval( centBinCenter[icent] );
    double sysNewCCG1E       = fSmoothSysNewCCG1E->Eval( centBinCenter[icent] );

    double sysTkQVn2       = fSmoothSysTkQVn2->Eval( centBinCenter[icent] );
    double sysTkQVn4       = fSmoothSysTkQVn4->Eval( centBinCenter[icent] );
    double sysTkQVn6       = fSmoothSysTkQVn6->Eval( centBinCenter[icent] );
    double sysTkQVn8       = fSmoothSysTkQVn8->Eval( centBinCenter[icent] );
    double sysTkQVn6Vn4    = fSmoothSysTkQVn6Vn4->Eval( centBinCenter[icent] );
    double sysTkQVn8Vn4    = fSmoothSysTkQVn8Vn4->Eval( centBinCenter[icent] );
    double sysTkQVn8Vn6    = fSmoothSysTkQVn8Vn6->Eval( centBinCenter[icent] );
    double sysTkQVn46_Vn68 = fSmoothSysTkQVn46_Vn68->Eval( centBinCenter[icent] );
    double sysTkQG1E       = fSmoothSysTkQG1E->Eval( centBinCenter[icent] );

    double sysVtxVn2       = fSmoothSysVtxVn2->Eval( centBinCenter[icent] );
    double sysVtxVn4       = fSmoothSysVtxVn4->Eval( centBinCenter[icent] );
    double sysVtxVn6       = fSmoothSysVtxVn6->Eval( centBinCenter[icent] );
    double sysVtxVn8       = fSmoothSysVtxVn8->Eval( centBinCenter[icent] );
    double sysVtxVn6Vn4    = fSmoothSysVtxVn6Vn4->Eval( centBinCenter[icent] );
    double sysVtxVn8Vn4    = fSmoothSysVtxVn8Vn4->Eval( centBinCenter[icent] );
    double sysVtxVn8Vn6    = fSmoothSysVtxVn8Vn6->Eval( centBinCenter[icent] );
    double sysVtxVn46_Vn68 = fSmoothSysVtxVn46_Vn68->Eval( centBinCenter[icent] );
    double sysVtxG1E       = fSmoothSysVtxG1E->Eval( centBinCenter[icent] );

    //-- Total
    double sysTotVn2       = sqrt( pow(sysRegVn2,2) + pow(sysRespVn2,2) + pow(sysNewCCVn2,2) + pow(sysTkQVn2,2) + pow(sysVtxVn2,2) + pow(0.005,2) );
    double sysTotVn4       = sqrt( pow(sysRegVn4,2) + pow(sysRespVn4,2) + pow(sysNewCCVn4,2) + pow(sysTkQVn4,2) + pow(sysVtxVn4,2) + pow(0.005,2) );
    double sysTotVn6       = sqrt( pow(sysRegVn6,2) + pow(sysRespVn6,2) + pow(sysNewCCVn6,2) + pow(sysTkQVn6,2) + pow(sysVtxVn6,2) + pow(0.005,2) );
    double sysTotVn8       = sqrt( pow(sysRegVn8,2) + pow(sysRespVn8,2) + pow(sysNewCCVn8,2) + pow(sysTkQVn8,2) + pow(sysVtxVn8,2) + pow(0.005,2) );
    double sysTotVn6Vn4    = sqrt( pow(sysRegVn6Vn4,2) + pow(sysRespVn6Vn4,2) + pow(sysNewCCVn6Vn4,2) + pow(sysTkQVn6Vn4,2) + pow(sysVtxVn6Vn4,2) );
    double sysTotVn8Vn4    = sqrt( pow(sysRegVn8Vn4,2) + pow(sysRespVn8Vn4,2) + pow(sysNewCCVn8Vn4,2) + pow(sysTkQVn8Vn4,2) + pow(sysVtxVn8Vn4,2) );
    double sysTotVn8Vn6    = sqrt( pow(sysRegVn8Vn6,2) + pow(sysRespVn8Vn6,2) + pow(sysNewCCVn8Vn6,2) + pow(sysTkQVn8Vn6,2) + pow(sysVtxVn8Vn6,2) );
    double sysTotVn46_Vn68 = sqrt( pow(sysRegVn46_Vn68,2) + pow(sysRespVn46_Vn68,2) + pow(sysNewCCVn46_Vn68,2) + pow(sysTkQVn46_Vn68,2) + pow(sysVtxVn46_Vn68,2) );
    double sysTotG1E       = sqrt( pow(sysRegG1E,2) + pow(sysRespG1E,2) + pow(sysNewCCG1E,2) + pow(sysTkQG1E,2) + pow(sysVtxG1E,2) );

    SmoothSysTotVn2       -> SetBinContent(icent+1, sysTotVn2);
    SmoothSysTotVn4       -> SetBinContent(icent+1, sysTotVn4);
    SmoothSysTotVn6       -> SetBinContent(icent+1, sysTotVn6);
    SmoothSysTotVn8       -> SetBinContent(icent+1, sysTotVn8);
    SmoothSysTotVn6Vn4    -> SetBinContent(icent+1, sysTotVn6Vn4);
    SmoothSysTotVn8Vn4    -> SetBinContent(icent+1, sysTotVn8Vn4);
    SmoothSysTotVn8Vn6    -> SetBinContent(icent+1, sysTotVn8Vn6);
    SmoothSysTotVn46_Vn68 -> SetBinContent(icent+1, sysTotVn46_Vn68);
    SmoothSysTotG1E       -> SetBinContent(icent+1, sysTotG1E);

    if(centBinCenter[icent] < 10. && centBinCenter[icent] > 5.){
      aveSys_00_10_Vn2       += sysTotVn2;
      aveSys_00_10_Vn4       += sysTotVn4;
      aveSys_00_10_Vn6       += sysTotVn6;
      aveSys_00_10_Vn8       += sysTotVn8;
      aveSys_00_10_Vn6Vn4    += sysTotVn6Vn4;
      aveSys_00_10_Vn8Vn4    += sysTotVn8Vn4;
      aveSys_00_10_Vn8Vn6    += sysTotVn8Vn6;
      aveSys_00_10_Vn46_Vn68 += sysTotVn46_Vn68;  
      aveSys_00_10_G1E       += sysTotG1E;
    }
    if(centBinCenter[icent] > 10. && centBinCenter[icent] < 20.){
      aveSys_10_20_Vn2       += sysTotVn2;
      aveSys_10_20_Vn4       += sysTotVn4;
      aveSys_10_20_Vn6       += sysTotVn6;
      aveSys_10_20_Vn8       += sysTotVn8;
      aveSys_10_20_Vn6Vn4    += sysTotVn6Vn4;
      aveSys_10_20_Vn8Vn4    += sysTotVn8Vn4;
      aveSys_10_20_Vn8Vn6    += sysTotVn8Vn6;
      aveSys_10_20_Vn46_Vn68 += sysTotVn46_Vn68;
      aveSys_10_20_G1E       += sysTotG1E;
    }
    if(centBinCenter[icent] > 20. && centBinCenter[icent] < 30.){
      aveSys_20_30_Vn2       += sysTotVn2;
      aveSys_20_30_Vn4       += sysTotVn4;
      aveSys_20_30_Vn6       += sysTotVn6;
      aveSys_20_30_Vn8       += sysTotVn8;
      aveSys_20_30_Vn6Vn4    += sysTotVn6Vn4;
      aveSys_20_30_Vn8Vn4    += sysTotVn8Vn4;
      aveSys_20_30_Vn8Vn6    += sysTotVn8Vn6;
      aveSys_20_30_Vn46_Vn68 += sysTotVn46_Vn68;
      aveSys_20_30_G1E       += sysTotG1E;
    }
    if(centBinCenter[icent] > 30. && centBinCenter[icent] < 50.){
      aveSys_30_50_Vn2       += sysTotVn2;
      aveSys_30_50_Vn4       += sysTotVn4;
      aveSys_30_50_Vn6       += sysTotVn6;
      aveSys_30_50_Vn8       += sysTotVn8;
      aveSys_30_50_Vn6Vn4    += sysTotVn6Vn4;
      aveSys_30_50_Vn8Vn4    += sysTotVn8Vn4;
      aveSys_30_50_Vn8Vn6    += sysTotVn8Vn6;
      aveSys_30_50_Vn46_Vn68 += sysTotVn46_Vn68;
      aveSys_30_50_G1E       += sysTotG1E;
    }
    if(centBinCenter[icent] > 50.){
      aveSys_50_60_Vn2       += sysTotVn2;
      aveSys_50_60_Vn4       += sysTotVn4;
      aveSys_50_60_Vn6       += sysTotVn6;
      aveSys_50_60_Vn8       += sysTotVn8;
      aveSys_50_60_Vn6Vn4    += sysTotVn6Vn4;
      aveSys_50_60_Vn8Vn4    += sysTotVn8Vn4;
      aveSys_50_60_Vn8Vn6    += sysTotVn8Vn6;
      aveSys_50_60_Vn46_Vn68 += sysTotVn46_Vn68;
      aveSys_50_60_G1E       += sysTotG1E;
    }

    //-- Ellp 
    if(icent >= 3 ){
      double sysRegKn     = fSmoothSysRegKn->Eval( centBinCenter[icent] );
      double sysRegAlpha  = fSmoothSysRegAlpha->Eval( centBinCenter[icent] );
      double sysRegE0     = fSmoothSysRegE0->Eval( centBinCenter[icent] );

      double sysRespKn     = fSmoothSysRespKn->Eval( centBinCenter[icent] );
      double sysRespAlpha  = fSmoothSysRespAlpha->Eval( centBinCenter[icent] );
      double sysRespE0     = fSmoothSysRespE0->Eval( centBinCenter[icent] );

      if(!propRespUncert){
	sysRespKn     = 0;
	sysRespAlpha  = 0;
	sysRespE0     = 0;
      }

      double sysNewCCKn     = fSmoothSysNewCCKn->Eval( centBinCenter[icent] );
      double sysNewCCAlpha  = fSmoothSysNewCCAlpha->Eval( centBinCenter[icent] );
      double sysNewCCE0     = fSmoothSysNewCCE0->Eval( centBinCenter[icent] );

      double sysTkQKn     = fSmoothSysTkQKn->Eval( centBinCenter[icent] );
      double sysTkQAlpha  = fSmoothSysTkQAlpha->Eval( centBinCenter[icent] );
      double sysTkQE0     = fSmoothSysTkQE0->Eval( centBinCenter[icent] );

      double sysVtxKn     = fSmoothSysVtxKn->Eval( centBinCenter[icent] );
      double sysVtxAlpha  = fSmoothSysVtxAlpha->Eval( centBinCenter[icent] );
      double sysVtxE0     = fSmoothSysVtxE0->Eval( centBinCenter[icent] );

      double sysTotKn    = sqrt( pow(sysRegKn,2) + pow(sysRespKn,2) + pow(sysNewCCKn,2) + pow(sysTkQKn,2) + pow(sysVtxKn,2) );
      double sysTotAlpha = sqrt( pow(sysRegAlpha,2) + pow(sysRespAlpha,2) + pow(sysNewCCAlpha,2) + pow(sysTkQAlpha,2) + pow(sysVtxAlpha,2) );
      double sysTotE0    = sqrt( pow(sysRegE0,2) + pow(sysRespE0,2) + pow(sysNewCCE0,2) + pow(sysTkQE0,2) + pow(sysVtxE0,2) );

      SmoothSysTotKn    -> SetBinContent(icent+1, sysTotKn);
      SmoothSysTotAlpha -> SetBinContent(icent+1, sysTotAlpha);
      SmoothSysTotE0    -> SetBinContent(icent+1, sysTotE0);

      if(centBinCenter[icent] > 15. && centBinCenter[icent] < 25.){
        aveSys_15_25_Kn     += sysTotKn;
        aveSys_15_25_Alpha  += sysTotAlpha;
        aveSys_15_25_E0     += sysTotE0;
      }
      if(centBinCenter[icent] > 25. && centBinCenter[icent] < 50.){
	aveSys_25_50_Kn     += sysTotKn;
	aveSys_25_50_Alpha  += sysTotAlpha;
	aveSys_25_50_E0     += sysTotE0;
      }
      if(centBinCenter[icent] > 50. && centBinCenter[icent] < 60.){
	aveSys_50_60_Kn     += sysTotKn;
	aveSys_50_60_Alpha  += sysTotAlpha;
	aveSys_50_60_E0     += sysTotE0;
      }

    } //-- End if


  }

  aveSys_00_10_Vn2  /= 1.;
  aveSys_10_20_Vn2  /= 2.;
  aveSys_20_30_Vn2  /= 2.;
  aveSys_30_50_Vn2  /= 4.;
  aveSys_50_60_Vn2  /= 2.;

  aveSys_00_10_Vn4  /= 1.;
  aveSys_10_20_Vn4  /= 2.;
  aveSys_20_30_Vn4  /= 2.;
  aveSys_30_50_Vn4  /= 4.;
  aveSys_50_60_Vn4  /= 2.;

  aveSys_00_10_Vn6  /= 1.;
  aveSys_10_20_Vn6  /= 2.;
  aveSys_20_30_Vn6  /= 2.;
  aveSys_30_50_Vn6  /= 4.;
  aveSys_50_60_Vn6  /= 2.;

  aveSys_00_10_Vn8  /= 1.;
  aveSys_10_20_Vn8  /= 2.;
  aveSys_20_30_Vn8  /= 2.;
  aveSys_30_50_Vn8  /= 4.;
  aveSys_50_60_Vn8  /= 2.;

  aveSys_00_10_Vn6Vn4  /= 1.;
  aveSys_10_20_Vn6Vn4  /= 2.;
  aveSys_20_30_Vn6Vn4  /= 2.;
  aveSys_30_50_Vn6Vn4  /= 4.;
  aveSys_50_60_Vn6Vn4  /= 2.;

  aveSys_00_10_Vn8Vn4  /= 1.;
  aveSys_10_20_Vn8Vn4  /= 2.;
  aveSys_20_30_Vn8Vn4  /= 2.;
  aveSys_30_50_Vn8Vn4  /= 4.;
  aveSys_50_60_Vn8Vn4  /= 2.;

  aveSys_00_10_Vn8Vn6  /= 1.;
  aveSys_10_20_Vn8Vn6  /= 2.;
  aveSys_20_30_Vn8Vn6  /= 2.;
  aveSys_30_50_Vn8Vn6  /= 4.;
  aveSys_50_60_Vn8Vn6  /= 2.;

  aveSys_00_10_Vn46_Vn68  /= 1.;
  aveSys_10_20_Vn46_Vn68  /= 2.;
  aveSys_20_30_Vn46_Vn68  /= 2.;
  aveSys_30_50_Vn46_Vn68  /= 4.;
  aveSys_50_60_Vn46_Vn68  /= 2.;

  aveSys_00_10_G1E  /= 1.;
  aveSys_10_20_G1E  /= 2.;
  aveSys_20_30_G1E  /= 2.;
  aveSys_30_50_G1E  /= 4.;
  aveSys_50_60_G1E  /= 2.;

  aveSys_15_25_Kn  /= 2.;
  aveSys_25_50_Kn  /= 5.;
  aveSys_50_60_Kn  /= 2.;

  aveSys_15_25_Alpha  /= 2.;
  aveSys_25_50_Alpha  /= 5.;
  aveSys_50_60_Alpha  /= 2.;

  aveSys_15_25_E0  /= 2.;
  aveSys_25_50_E0  /= 5.;
  aveSys_50_60_E0  /= 2.;

  //-- make percents
  aveSys_00_10_Vn2 *= 100.;
  aveSys_10_20_Vn2 *= 100.;
  aveSys_20_30_Vn2 *= 100.;
  aveSys_30_50_Vn2 *= 100.;
  aveSys_50_60_Vn2 *= 100.;

  aveSys_00_10_Vn4 *= 100.;
  aveSys_10_20_Vn4 *= 100.;
  aveSys_20_30_Vn4 *= 100.;
  aveSys_30_50_Vn4 *= 100.;
  aveSys_50_60_Vn4 *= 100.;

  aveSys_00_10_Vn6 *= 100.;
  aveSys_10_20_Vn6 *= 100.;
  aveSys_20_30_Vn6 *= 100.;
  aveSys_30_50_Vn6 *= 100.;
  aveSys_50_60_Vn6 *= 100.;

  aveSys_00_10_Vn8 *= 100.;
  aveSys_10_20_Vn8 *= 100.;
  aveSys_20_30_Vn8 *= 100.;
  aveSys_30_50_Vn8 *= 100.;
  aveSys_50_60_Vn8 *= 100.;

  aveSys_00_10_Vn6Vn4 *= 100.;
  aveSys_10_20_Vn6Vn4 *= 100.;
  aveSys_20_30_Vn6Vn4 *= 100.;
  aveSys_30_50_Vn6Vn4 *= 100.;
  aveSys_50_60_Vn6Vn4 *= 100.;

  aveSys_00_10_Vn8Vn4 *= 100.;
  aveSys_10_20_Vn8Vn4 *= 100.;
  aveSys_20_30_Vn8Vn4 *= 100.;
  aveSys_30_50_Vn8Vn4 *= 100.;
  aveSys_50_60_Vn8Vn4 *= 100.;

  aveSys_00_10_Vn8Vn6 *= 100.;
  aveSys_10_20_Vn8Vn6 *= 100.;
  aveSys_20_30_Vn8Vn6 *= 100.;
  aveSys_30_50_Vn8Vn6 *= 100.;
  aveSys_50_60_Vn8Vn6 *= 100.;

  aveSys_00_10_Vn46_Vn68 *= 100.;
  aveSys_10_20_Vn46_Vn68 *= 100.;
  aveSys_20_30_Vn46_Vn68 *= 100.;
  aveSys_30_50_Vn46_Vn68 *= 100.;
  aveSys_50_60_Vn46_Vn68 *= 100.;

  aveSys_00_10_G1E *= 100.;
  aveSys_10_20_G1E *= 100.;
  aveSys_20_30_G1E *= 100.;
  aveSys_30_50_G1E *= 100.;
  aveSys_50_60_G1E *= 100.;

  aveSys_15_25_Kn *= 100.;
  aveSys_25_50_Kn *= 100.;
  aveSys_50_60_Kn *= 100.;

  aveSys_15_25_Alpha *= 100.;
  aveSys_25_50_Alpha *= 100.;
  aveSys_50_60_Alpha *= 100.;

  aveSys_15_25_E0 *= 100.;
  aveSys_25_50_E0 *= 100.;
  aveSys_50_60_E0 *= 100.;

  std::cout << "\n\n======================================================" << std::endl;
  std::cout << "Cent"  << "\tVn2      " << "\tVn4      " << "\tVn6      " << "\tVn8      " << "\tG1E      " << "\tVn6Vn4   " << "\tVn8Vn4   "<< "\tVn8Vn6   " << std::endl;
  std::cout << "0--10\\%  & "  << Form("%.1f", aveSys_00_10_Vn2)    << " & " << Form("%.1f", aveSys_00_10_Vn4)    << " & " << Form("%.1f", aveSys_00_10_Vn6) 
	    << " & "           << Form("%.1f", aveSys_00_10_Vn8)    << " & " << Form("%.1f", aveSys_00_10_G1E)    << " & " << Form("%.1f", aveSys_00_10_Vn6Vn4) 
	    << " & "           << Form("%.1f", aveSys_00_10_Vn8Vn4) << " & " << Form("%.1f", aveSys_00_10_Vn8Vn6) << " \\\\" << std::endl;
  std::cout << "10--20\\% & "  << Form("%.1f", aveSys_10_20_Vn2)    << " & " << Form("%.1f", aveSys_10_20_Vn4)    << " & " << Form("%.1f", aveSys_10_20_Vn6)
            << " & "           << Form("%.1f", aveSys_10_20_Vn8)    << " & " << Form("%.1f", aveSys_10_20_G1E)    << " & " << Form("%.1f", aveSys_10_20_Vn6Vn4)
            << " & "           << Form("%.1f", aveSys_10_20_Vn8Vn4) << " & " << Form("%.1f", aveSys_10_20_Vn8Vn6) << " \\\\" << std::endl;
  std::cout << "20--30\\% & "  << Form("%.1f", aveSys_20_30_Vn2)    << " & " << Form("%.1f", aveSys_20_30_Vn4)    << " & " << Form("%.1f", aveSys_20_30_Vn6)
            << " & "           << Form("%.1f", aveSys_20_30_Vn8)    << " & " << Form("%.1f", aveSys_20_30_G1E)    << " & " << Form("%.1f", aveSys_20_30_Vn6Vn4)
            << " & "           << Form("%.1f", aveSys_20_30_Vn8Vn4) << " & " << Form("%.1f", aveSys_20_30_Vn8Vn6) << " \\\\" << std::endl;
  std::cout << "30--50\\% & "  << Form("%.1f", aveSys_30_50_Vn2)    << " & " << Form("%.1f", aveSys_30_50_Vn4)    << " & " << Form("%.1f", aveSys_30_50_Vn6)
            << " & "           << Form("%.1f", aveSys_30_50_Vn8)    << " & " << Form("%.1f", aveSys_30_50_G1E)    << " & " << Form("%.1f", aveSys_30_50_Vn6Vn4)
            << " & "           << Form("%.1f", aveSys_30_50_Vn8Vn4) << " & " << Form("%.1f", aveSys_30_50_Vn8Vn6) << " \\\\" << std::endl;
  std::cout << "50--60\\% & "  << Form("%.1f", aveSys_50_60_Vn2)    << " & " << Form("%.1f", aveSys_50_60_Vn4)    << " & " << Form("%.1f", aveSys_50_60_Vn6)
            << " & "           << Form("%.1f", aveSys_50_60_Vn8)    << " & " << Form("%.1f", aveSys_50_60_G1E)    << " & " << Form("%.1f", aveSys_50_60_Vn6Vn4)
            << " & "           << Form("%.1f", aveSys_50_60_Vn8Vn4) << " & " << Form("%.1f", aveSys_50_60_Vn8Vn6) << " \\\\" << std::endl;


  std::cout << "\n\n======================================================" << std::endl;
  std::cout << "Cent"  << "\tKn   " << "\tAlpha" << "\tE0   " << std::endl;
  std::cout << "15-25%\t" << Form("%.1f", aveSys_15_25_Kn)    << "\t" << Form("%.1f", aveSys_15_25_Alpha)    << "\t" << Form("%.1f", aveSys_15_25_E0) << std::endl;
  std::cout << "25-50%\t" << Form("%.1f", aveSys_25_50_Kn)    << "\t" << Form("%.1f", aveSys_25_50_Alpha)    << "\t" << Form("%.1f", aveSys_25_50_E0) << std::endl;
  std::cout << "50-60%\t" << Form("%.1f", aveSys_50_60_Kn)    << "\t" << Form("%.1f", aveSys_50_60_Alpha)    << "\t" << Form("%.1f", aveSys_50_60_E0) << std::endl;

  fOut->Write();

}

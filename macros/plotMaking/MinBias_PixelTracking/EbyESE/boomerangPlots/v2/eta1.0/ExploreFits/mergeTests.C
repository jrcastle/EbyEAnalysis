#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

#include <iostream>

using namespace hi;
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


}

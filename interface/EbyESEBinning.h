#ifndef __EbyESEBinning__
#define __EbyESEBinning__

#include "TGraphErrors.h"
#include "TMatrixD.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TString.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

namespace ebyese{

  //-- Analyzer output tree
  const TString fAnaTreeName = "/rfs/jcastle/PbPb2015/PixelTracking_MB2/EbyETree_pixel.root";
  const TString fileSplit    = "/rfs/jcastle/PbPb2015/PixelTracking_MB2/SplitTree.root";

  //-- Analyzer pt binning
  //static const int    nptbinsDefault = 7;
  //const double ptbinsDefault[] = {
  //  0.30,    0.50,    1.00,    1.25,    1.50,
  //  2.00,    2.50,    3.00
  //};
  static const int    nptbinsDefault = 2;
  const double ptbinsDefault[] = {
    0.30,    1.00,    3.00
  };

  //-- Analyzer eta binning
  //static const int    netabinsDefault = 14;
  //const double etabinsDefault[] = {
  //  -2.4,     -2.0,     -1.6,     -1.2,     -1.0,
  //  -0.8,     -0.4,      0.0,      0.4,      0.8,
  //   1.0,      1.2,      1.6,      2.0,      2.4
  //};
  static const int    netabinsDefault = 4;
  const double etabinsDefault[] = {
    -2.4,    -1.0,    0.0,    1.0,    2.4
  };

  //-- Unfolding pt binning 
  //static const int NPT  = 7;
  //const double pt_min[NPT] = {
  //  0.30,    0.50,    1.00,    1.25,    1.50,
  //  2.00,    2.50
  //};
  //const double pt_max[NPT] ={
  //  0.50,    1.00,    1.25,    1.50,    2.00,
  //  2.50,    3.00
  //};
  static const int NPT  = 2;
  const double pt_min[NPT] = {
    0.30,    1.00
  };
  const double pt_max[NPT] ={
    1.00,    3.00
  };

  //-- Centrality binning, colors and markers
  static const int NCENT    = 12;
  const double centBinWidth = 5.;
  const int cent_min[NCENT] = {
    0,      5,    10,     15,     20,
    25,    30,    35,     40,     45,
    50,    55
  };
  const int cent_max[NCENT] = {
    5,     10,    15,    20,    25,
    30,    35,    40,    45,    50,
    55,    60
  };
  const double centbinsDefault[] = {
    0,      5,    10,     15,     20,
    25,    30,    35,     40,     45,
    50,    55,    60
  };
  const double centBinCenter[] = {
     2.5,     7.5,    12.5,    17.5,    22.5,
    27.5,    32.5,    37.5,    42.5,    47.5,
    52.5,    57.5
  };
  const double centBinErr[] = {
    1.0,    1.0,    1.0,    1.0,    1.0,
    1.0,    1.0,    1.0,    1.0,    1.0,
    1.0,    1.0
  };
  const double CERR[NCENT] = {
    0., 0., 0., 0., 0., 
    0., 0., 0., 0., 0., 
    0., 0.
  };
  const double Npart[NCENT] = {
    384.4,    333.4,    285.4,    242.9,    205.7,
    172.7,    144.1,    118.7,    96.46,    77.35,
    60.81,    46.89
  };
  /*
  static const int NCENT    = 6;
  const double centBinWidth = 10.;
  const int cent_min[NCENT] = {
    0,     10,    20,    30,    40,
    50
  };
  const int cent_max[NCENT] = {
    10,    20,    30,    40,    50,
    60
  };
  const double centbinsDefault[] = {
    0,     10,    20,    30,    40,
    50,    60
  };
  const double centBinCenter[] = {
    5,     15,    25,    35,    45,
    55
  };
  const double centBinErr[] = {
    1.0,    1.0,    1.0,    1.0,    1.0,
    1.0    
  };
  const double CERR[NCENT] = {
    0., 0., 0., 0., 0., 
    0.
  };
  */
  TH1D hCentBins("hCentBins", "hCentBins", NCENT, centbinsDefault);
  const int centCol[]  = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3, kOrange-2};
  const int centMark[] = {24, 29, 27, 34, 25, 20, 26, 33, 28, 23, 32, 21};

  //-- File Splitting for getting a handle on statistical errors
  const int NSPLIT = 10;

  //-- Number of QN bins
  static const int NQN = 10;
  const int qnCol[] = {kBlack, kBlue, kViolet-1, kMagenta+1, kRed, kOrange-3, kGreen+1, kCyan+2, kRed+3, kGreen+3};
  const int qnMrk[] = {24, 29, 27, 34, 25, 20, 26, 33, 28, 23};

  //-- Harmonics for cross-correlations
  static const int NVN = 3;
  const int vn_[] = {2,3,4};

  //-- VN binning
  const int    NBins    = 150;
  const double vnMax[]  = {0.0, 0.0, 0.6, 0.4, 0.4}; //-- First two elements are 0.0 because this is accessed in code by vnMax[norder_] where norder_ = 2, 3, 4

  const double v2Max[NCENT] = {
    0.168,    0.200,    0.224,    0.240,    0.272,
    0.296,    0.320,    0.336,    0.360,    0.392,
    0.424,    0.448
  };
  double binw = 0.008;
  int NBinsV2[NCENT] = {
    static_cast<int>( v2Max[0]/binw),    static_cast<int>( v2Max[1]/binw),    static_cast<int>( v2Max[2]/binw),    static_cast<int>( v2Max[3]/binw),    static_cast<int>( v2Max[4]/binw),
    static_cast<int>( v2Max[5]/binw),    static_cast<int>( v2Max[6]/binw),    static_cast<int>( v2Max[7]/binw),    static_cast<int>( v2Max[8]/binw),    static_cast<int>( v2Max[9]/binw),
    static_cast<int>( v2Max[10]/binw),    static_cast<int>( v2Max[11]/binw)
  };

  /*
  const double v2Max[NCENT] = {
    0.213,    0.214,    0.231,    0.261,    0.273,
    0.300,    0.333,    0.375,    0.386,    0.450,
    0.480,    0.533
  };
  double binw[NCENT] = {
    0.019,    0.021,    0.023,    0.025,    0.027,
    0.030,    0.033,    0.037,    0.042,    0.048,
    0.055,    0.064
  };
  int NBinsV2[NCENT] = {
    static_cast<int>( v2Max[0]/binw[0]),    static_cast<int>( v2Max[1]/binw[1]),    static_cast<int>( v2Max[2]/binw[2]),    static_cast<int>( v2Max[3]/binw[3]),    static_cast<int>( v2Max[4]/binw[4]),
    static_cast<int>( v2Max[5]/binw[5]),    static_cast<int>( v2Max[6]/binw[6]),    static_cast<int>( v2Max[7]/binw[7]),    static_cast<int>( v2Max[8]/binw[8]),    static_cast<int>( v2Max[9]/binw[9]),
    static_cast<int>( v2Max[10]/binw[10]),  static_cast<int>( v2Max[11]/binw[11])
  };
  */

  //-- Iterations and colors for unfolding
  const int NITER        = 11;
  const int iter[]       = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
  const double diter[]   = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
  const double iterErr[] = {0, 0, 0, 0,  0,  0,  0,   0,   0,   0,    0};
  const int col[]        = {kOrange-2, 38, 46, kGreen+3, kCyan, kMagenta, 30, kViolet-1, 28, kBlue, kRed};

  //-- SVD kregs
  const int NKREG = 12;
  //const int NKREG = 8;

  //-- Chi2 cutoff scenarios
  //-- [0] = loose
  //-- [1] = nominal
  //-- [2] = tight
  const int NCHICUT        = 3;
  const double chi2Cut[]   = {1.5, 1.2, 1.0};
  string chi2CutScenario[] = {"loose", "nominal", "tight"};

  const int NCHI2 = 12;
  const double chi2Cutoff[] = {
    1.0,    1.2,    1.5,    2.0,    3.0,
    4.0,    5.0,    6.0,    7.0,    8.0,
    9.0,   10.0
  };
  const double chi2Cutoffe[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  //-- ------------------------- Helpful Functions -------------------------
  void formatGraph(TGraphErrors * g, string xT, double yMin, double yMax, string yT, int Color, int mrkS, string name){
    g->GetXaxis()->SetTitle( xT.data() );
    g->GetYaxis()->SetRangeUser(yMin, yMax);
    g->GetYaxis()->SetTitle( yT.data() );
    g->SetLineColor( Color );
    g->SetMarkerColor( Color );
    g->SetMarkerStyle( mrkS );
    g->GetXaxis()->SetNdivisions(509);
    g->SetName( name.data() );
  }
  void formatGraph(TGraph * g, string xT, double yMin, double yMax, string yT, int Color, int mrkS, string name){
    g->GetXaxis()->SetTitle( xT.data() );
    g->GetYaxis()->SetRangeUser(yMin, yMax);
    g->GetYaxis()->SetTitle( yT.data() );
    g->SetLineColor( Color );
    g->SetMarkerColor( Color );
    g->SetMarkerStyle( mrkS );
    g->GetXaxis()->SetNdivisions(509);
    g->SetName( name.data() );
  }

  void crap(int i){
    std::cout << "DEBUG " << i << std::endl;
  }
  void crap(string s){
    std::cout << "DEBUG: " << s.data() << std::endl;
  }

  bool compatibleWithOne(double ratio, double err){
    double high = ratio+err;
    double low  = ratio-err;
    if( ratio >= 1.0 && low < 1.0 )  return true;
    if( ratio < 1.0 && high >= 1.0 ) return true;
    return false;
  }

  bool compatibleWithZero(double ratio, double err){
    double low  = ratio-err;
    if( ratio >= 0. && low < 0. )  return true;
    return false;
  }
  /*
  void FixUnfold(TH1D * h){
    bool threshold = 0;
    int nb = h->GetNbinsX();
    for(int ibin = 1; ibin <= nb; ibin++){
      if( ibin < nb/10 ) continue;
      double bc = h->GetBinContent(ibin);
      if( !threshold && bc < 1. ) threshold = 1;
      if( threshold ){
        h->SetBinContent(ibin, 0);
        h->SetBinError(ibin, 0);
      }
    }
  }
  */
  void FixUnfold(TH1D * h){
    bool threshold = 0;
    int nb         = h->GetNbinsX();
    double sig4    = h->GetMean() + 4.*h->GetRMS();
    for(int ibin = 1; ibin <= nb; ibin++){
      double bc   = h->GetBinCenter(ibin);
      double binc = h->GetBinContent(ibin);
      if( !threshold && bc >= sig4 ) threshold = 1;
      if( threshold ){
        h->SetBinContent(ibin, 0);
        h->SetBinError(ibin, 0);
      }
    }
  }

  void  M2H(TMatrixD M, TH2D * H){
    int nx = M.GetNrows();
    int ny = M.GetNcols();
    for(int x = 0; x < nx; x++){
      for(int y = 0; y < ny; y++){
	H->SetBinContent(x+1, y+1, M(y,x));
      }
    }
  }

  TMatrixD H2M(TH2D * H){
    int nx = H->GetXaxis()->GetNbins();
    int ny = H->GetYaxis()->GetNbins();
    TMatrixD M(ny, nx);
    for(int x = 0; x < nx; x++){
      for(int y = 0; y < ny; y++){
	M(y,x) = H->GetBinContent(x+1, y+1);
      }
    }
    return M;
  }

  void legInit( TLegend * l ){
    l->SetBorderSize(0);
    l->SetFillStyle(0);
  }



}
#endif

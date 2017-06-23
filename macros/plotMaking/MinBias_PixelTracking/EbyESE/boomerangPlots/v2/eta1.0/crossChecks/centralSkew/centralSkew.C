#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/SysTablesEbyESE.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;
using namespace sysebyese;

double pEllP(double * x, double * par){

  int nbin = 180;
  //int nbin = 360;

  //-- [0] = e0
  //-- [1] = alpha
  //-- [2] = kn
  //-- [3] = Scale 
  //-- xx  = vn

  double e0    = par[0];
  double alpha = par[1];
  double kn    = par[2];
  double scale = par[3];
  double eccn  = x[0] / kn;


  double pi = TMath::Pi();
  double p1 = (scale * 2. * alpha * eccn / pi / kn) * pow( 1 - e0*e0, alpha + 0.5);
  //double p1 = (2. * alpha * eccn / pi / kn) * pow( 1 - e0*e0, alpha + 0.5);

  //-- [0] = xx
  //-- [1] = e0
  //-- [2] = alpha
  //-- [3] = kn
  //TF1 f("f", "pow(1 - [0]*[0]/[3]/[3], [2] - 1) / pow(1 - [1]*[0]*TMath::Cos(x)/[3]/[3], 2*[2]+1)", 0., 2*pi);
  //f.SetParameters(xx, e0, alpha, kn);
  //double integ = f.Integral(0, pi);
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 0; i <= nbin; i++){
    double phi = 0. + dphi*(double)(i);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  //integ += 0.5 * dphi * ( pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(0), -2.*alpha-1 ) + pow(1-eccn*eccn, alpha-1) *pow( 1-e0*eccn*cos(pi), -2.*alpha-1 ) );
  double ellp = p1 * integ;

  return ellp;

}

void centralSkew(){

  int N = 1000000;
  TH1D::SetDefaultSumw2();

  TF1 * f;
  TH1D * hBG;

  TLatex latex;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  //f = new TF1("f", "[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))", 0.0, 0.3);
  //f->SetParameters(1.0, 0.0, 0.02);

  f = new TF1("f", pEllP, 0.0, 0.22, 4);
  //f->SetParameters(0.05, 200, 4.33292e-01, 3.92495e-03);
  f->SetParameters(2.93223e-01, 1.01156e+01, 2.70041e-01, 3.94381e-03);

  hBG = new TH1D("hBG", "hBG", 150, 0.0, 0.3);
  hBG->GetXaxis()->SetTitle("v_{2}");

  for(int i = 0; i < N; i++){
    double v2 = f->GetRandom();
    hBG->Fill(v2);
  }

  EbyECumu cumu( hBG );
  double gamma1exp = cumu.GetGamma1Exp();

  std::cout << gamma1exp << "\t" << hBG->GetSkewness() << std::endl;

  TCanvas * c = new TCanvas("c", "c", 500, 500);
  c->cd();
  c->SetLogy();
  hBG->Draw();



}

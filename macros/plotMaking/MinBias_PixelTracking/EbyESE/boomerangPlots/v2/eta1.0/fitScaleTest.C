#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

double fk(double k, double alpha, double ecc0){
  double p1 = alpha / ( alpha + k );
  double p2 = pow( 1 - pow(ecc0,2), k );
  double p3 = ROOT::Math::hyperg( k+0.5, k, alpha+k+1, pow(ecc0,2) );
  double f  = p1*p2*p3;
  return f;
}

void fitScaleTest(){

  double ecc0[NCENT]  = {0.0,   0.0,   0.0,   0.177, 0.195, 0.216, 0.226, 0.240, 0.265, 0.259, 0.271, 0.279};
  double alpha[NCENT] = {2.8e5, 2.8e5, 2.8e5, 69.3,  62.4,  50.3,  44.0,  33.4,  23.0,  21.0,  15.6,  12.3};
  double kn[NCENT]    = {16.2,  16.2,  16.2,  0.366, 0.383, 0.379, 0.385, 0.371, 0.338, 0.341, 0.314, 0.289};

  TFile * fRes;
  TGraphErrors * grVn2;
  TGraphErrors * grVn4;
  TGraphErrors * grVn6;


  TGraph * grFitEn6En4;
  TGraph * grFitEn4En2;
  double fitEn6En4[NCENT];
  double fitEn4En2[NCENT];

  TGraph * grVn6Vn4;
  TGraph * grVn4Vn2;
  double vn6vn4[NCENT];
  double vn4vn2[NCENT];

  //
  // MAIN
  //
  setTDRStyle();

  fRes = new TFile( "systematicStudies/PhysicsResults.root" );
  grVn2 = (TGraphErrors*) fRes->Get("grVn2Raw");
  grVn4 = (TGraphErrors*) fRes->Get("grVn4Raw");
  grVn6 = (TGraphErrors*) fRes->Get("grVn6Raw");

  for(int icent = 3; icent < NCENT; icent++){


    vn6vn4[icent] = grVn6->GetY()[icent] / grVn4->GetY()[icent];
    vn4vn2[icent] = grVn4->GetY()[icent] / grVn2->GetY()[icent];

    double a = alpha[icent];
    double e = ecc0[icent];

    double en2 = pow( 1 - fk(1,a,e), 1./2.);
    double en4 = pow( 1 - 2.*fk(1,a,e) + 2.*pow(fk(1,a,e), 2) - fk(2,a,e), 1./4.);
    double en6 = pow( 1 + (9./2.)*pow(fk(1,a,e),2) - 3.*pow(fk(1,a,e),3) + 3.*fk(1,a,e)*( (3./4.)*fk(2,a,e)-1 ) - (3./2.)*fk(2,a,e) - (1./4.)*fk(3,a,e), 1./6.);

    fitEn6En4[icent] = en6 / en4;
    fitEn4En2[icent] = en4 / en2; 

    std::cout << en4 / en2 <<std::endl;

  } //-- End Cent loop

  grFitEn6En4 = new TGraph(NCENT, centBinCenter, fitEn6En4);
  grFitEn4En2 = new TGraph(NCENT, centBinCenter, fitEn4En2);
  grVn6Vn4 = new TGraph(NCENT, centBinCenter, vn6vn4);
  grVn4Vn2 = new TGraph(NCENT, centBinCenter, vn4vn2);

  formatGraph(grFitEn6En4, "Centrality %", 0.95, 1.02, "Cumu Ratio", 1, 24, "grFitEn6En4");
  formatGraph(grFitEn4En2, "Centrality %", 0.80, 0.90, "Cumu Ratio", 1, 25, "grFitEn4En2");
  formatGraph(grVn6Vn4,    "Centrality %", 0.95, 1.02, "Cumu Ratio", 2, 20, "grVn6Vn4");
  formatGraph(grVn4Vn2,    "Centrality %", 0.80, 0.90, "Cumu Ratio", 2, 21, "grVn4Vn2");

  TLegend * l1 = new TLegend(0.18, 0.77, 0.49, 0.93);
  legInit( l1 );
  l1->AddEntry(grVn4Vn2,    "v_{2}{4} / v_{2}{2}",               "p");
  l1->AddEntry(grFitEn4En2, "#epsilon_{2}{4} / #epsilon_{2}{2}", "p");

  TLegend * l2 = new TLegend(0.18, 0.77, 0.49, 0.93);
  legInit( l2 );
  l2->AddEntry(grVn6Vn4,    "v_{2}{6} / v_{2}{4}",              "p");
  l2->AddEntry(grFitEn6En4, "#epsilon_{2}{6} / #epsilon_{2}{4}", "p");



  TCanvas * c = new TCanvas("c", "c", 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  grVn4Vn2->Draw("ap");
  grFitEn4En2->Draw("psame");
  l1->Draw("same");
  c->cd(2);
  grVn6Vn4->Draw("ap");
  grFitEn6En4->Draw("psame");
  l2->Draw("same");
  c->SaveAs("plots/fitScaleTest.pdf");

}

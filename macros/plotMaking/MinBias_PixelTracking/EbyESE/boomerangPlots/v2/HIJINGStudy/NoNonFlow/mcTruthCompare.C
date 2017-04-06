#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

double pEllP(double * x, double * par){
  int nbin = 50;
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
  
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 1; i <= nbin; i++){
    double phi = (2*(double)i-1.)*pi/(2.*(double)nbin);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }
  double ellp = p1 * integ;
  return ellp;
}


void mcTruthCompare(){

  const int c = 0;
  const int norder_ = 2;

  TFile * fAna;
  TH1D * hObs[NCENT];
  TH1D * hV2True[NCENT];

  TFile * fUnf;
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hRefold[NCENT][NITER];
  TH1D * hFinalUnfold[NCENT];

  TFile * fStat;
  TH1D * hVarianceOfMean_Vn4;
  TH1D * hVarianceOfMean_Vn6;
  TH1D * hVarianceOfMean_Vn8;
  TH1D * hVarianceOfMean_Vn4True;
  TH1D * hVarianceOfMean_Vn6True;
  TH1D * hVarianceOfMean_Vn8True;

  TF1 * fMCTruth;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(20);

  //-- Grab analyzer and unfold files
  fAna = new TFile("AnalyzerResults/CastleEbyE.root");
  fUnf = new TFile( Form("UnfoldResults/dataResp/data%i.root", norder_) );

  //-- Get Stat Uncerts
  fStat = new TFile(Form("StatErr/StatisticalUncertainties_v%i.root", norder_));
  hVarianceOfMean_Vn4     = (TH1D*) fStat->Get("hVarianceOfMean_Vn4");
  hVarianceOfMean_Vn6     = (TH1D*) fStat->Get("hVarianceOfMean_Vn6");
  hVarianceOfMean_Vn8     = (TH1D*) fStat->Get("hVarianceOfMean_Vn8");
  hVarianceOfMean_Vn4True = (TH1D*) fStat->Get("hVarianceOfMean_Vn4True");
  hVarianceOfMean_Vn6True = (TH1D*) fStat->Get("hVarianceOfMean_Vn6True");
  hVarianceOfMean_Vn8True = (TH1D*) fStat->Get("hVarianceOfMean_Vn8True");


  //-- Set up the MC Truth
  double ecc0_8pct_v2  = 2.16286e-01;
  double alpha_8pct_v2 = 5.02657e+01;
  double kn_8pct_v2    = 3.79222e-01;
  fMCTruth = new TF1("fMCTruth", pEllP, 0.0, kn_8pct_v2, 4); //-- There's an asymptote at vn = kn, don't let the TF1 range exceed that value
  fMCTruth->SetParNames("ecc0", "alpha", "kn", "Scale");
  fMCTruth->SetParameters(ecc0_8pct_v2, alpha_8pct_v2, kn_8pct_v2, 1.0);
  fMCTruth->SetLineColor(2);

  //-- Loop over centbins
  for(int icent = 0; icent < NCENT; icent++){

    if(icent != c) continue;

    //-- Grab the observed and truth  hist
    hObs[icent]    = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
    hV2True[icent] = (TH1D*) fAna->Get( Form("qwebye/hV2True_c%i", icent) );
    hV2True[icent]->SetLineColor(2);
    hV2True[icent]->SetMarkerColor(2);

    //-- Loop over iterations and get the final unfolding iteration
    for(int i = 0; i < NITER; i++){
      hUnfold[icent][i] = (TH1D*) fUnf->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold[icent][i]->GetXaxis()->SetRange( 1, hUnfold[icent][i]->FindBin(0.3) );
      hRefold[icent][i] = (TH1D*) fUnf->Get( Form("hrefold%i_c%i", iter[i], icent) );

      double chi2 = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");
      std::cout<<"Grabbing iter = " << iter[i] << "\t Refold Chi2 = " << chi2 <<std::endl;
      if(chi2 < 1.2){
	hFinalUnfold[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form("hFinalUnfold_c%i", icent) );
	break;
      }
      if(i == NITER-1){
	hFinalUnfold[icent] = (TH1D*) hUnfold[icent][i]->Clone(Form("hFinalUnfold_c%i", icent) );
        break;
      } 

    } //-- End iter loop

  } //-- End cent loop

  //-- Quan Cumu compare
  FixUnfold(hFinalUnfold[c]);
  EbyECumu cumu( hFinalUnfold[c] );
  double v24  = cumu.GetCumu_vn4();
  double v26  = cumu.GetCumu_vn6();
  double v28  = cumu.GetCumu_vn8();
  double v24e = sqrt(hVarianceOfMean_Vn4->GetBinContent(c+1));
  double v26e = sqrt(hVarianceOfMean_Vn6->GetBinContent(c+1));
  double v28e = sqrt(hVarianceOfMean_Vn8->GetBinContent(c+1));

  double ebyeCumu[3]  = {v24, v26, v28};
  double ebyeCumue[3] = {v24e, v26e, v28e};

  EbyECumu cumuT(hV2True[c]);
  double v24t  = cumuT.GetCumu_vn4();
  double v26t  = cumuT.GetCumu_vn6();
  double v28t  = cumuT.GetCumu_vn8();
  double v24te = sqrt(hVarianceOfMean_Vn4True->GetBinContent(c+1));
  double v26te = sqrt(hVarianceOfMean_Vn6True->GetBinContent(c+1));
  double v28te = sqrt(hVarianceOfMean_Vn8True->GetBinContent(c+1));

  double trueCumu[3]  = {v24t, v26t, v28t};
  double trueCumue[3] = {v24te, v26te, v28te};

  double quanV24 = 0.0827373390644;
  double quanV26 = 0.0824165143178;
  double quanV28 = 0.0824136316097;

  double quanV24e = 6.35622389118e-05;
  double quanV26e = 5.69858415971e-05;
  double quanV28e = 5.38852799655e-05;

  double quanCumu[3]  = {quanV24,  quanV26,  quanV28};
  double quanCumue[3] = {quanV24e, quanV26e, quanV28e};

  TH1D * grEbyE = new TH1D("grEbyE", "grEbyE", 3, 1, 3);
  TH1D * grQuan = new TH1D("grQuan", "grQuan", 3, 1, 3);
  TH1D * grTrue = new TH1D("grTrue", "grTrue", 3, 1, 3);

  string st[3] = {"v_{2}{4}", "v_{2}{6}", "v_{2}{8}"};

  for(int i = 0; i<3; i++){
    grEbyE->SetBinContent(i+1, ebyeCumu[i]);
    grEbyE->SetBinError(i+1, ebyeCumue[i]);
 
    grTrue->SetBinContent(i+1, trueCumu[i]);
    grTrue->SetBinError(i+1, trueCumue[i]);

    grQuan->SetBinContent(i+1, quanCumu[i]);
    grQuan->SetBinError(i+1, quanCumue[i]);

    grEbyE->GetXaxis()->SetBinLabel(i+1, st[i].data());
    grTrue->GetXaxis()->SetBinLabel(i+1, st[i].data());
    grQuan->GetXaxis()->SetBinLabel(i+1, st[i].data());

  }
  grEbyE->SetMinimum(0.082);
  grEbyE->SetMaximum(0.084);
  grTrue->SetMinimum(0.082);
  grTrue->SetMaximum(0.084);
  grQuan->SetMinimum(0.082);
  grQuan->SetMaximum(0.084);

  grEbyE->GetXaxis()->SetLabelSize(0.08);
  grTrue->GetXaxis()->SetLabelSize(0.08);
  grQuan->GetXaxis()->SetLabelSize(0.08);

  grEbyE->GetYaxis()->SetNdivisions(508);
  grTrue->GetYaxis()->SetNdivisions(508);
  grQuan->GetYaxis()->SetNdivisions(508);

  grEbyE->SetLineColor(4);
  grEbyE->SetMarkerColor(4);

  grQuan->SetLineColor(2);
  grQuan->SetMarkerColor(2);

  //-- Draw MC compare 
  TLegend * leg = new TLegend(0.67, 0.78, 0.92, 0.93);
  legInit( leg );
  leg->AddEntry(hFinalUnfold[c], "Unfolding", "lp");
  leg->AddEntry(hV2True[c],      "MC Truth",  "lp");

  hFinalUnfold[c]->GetXaxis()->SetTitle("v_{2}");
  hFinalUnfold[c]->GetYaxis()->SetTitle("Events");

  TCanvas * cRawDistComp = new TCanvas("cRawDistComp", "cRawDistComp", 500, 500);
  cRawDistComp->cd();
  cRawDistComp->SetLogy();
  hFinalUnfold[c]->Draw();
  hV2True[c]->Draw("same");
  leg->Draw("same");
  cRawDistComp->SaveAs("cRawDistComp.pdf");

  //-- Draw Quan Compare
  TLegend * legQ = new TLegend(0.67, 0.78, 0.92, 0.93);
  legInit(legQ);
  legQ->AddEntry(grTrue, "Truth", "lp");
  legQ->AddEntry(grEbyE, "EbyE",  "lp");
  legQ->AddEntry(grQuan, "QCumu", "lp");

  TCanvas * cQuan = new TCanvas("cQuan", "cQuan", 500, 500);
  cQuan->cd();
  grTrue->Draw();
  grEbyE->Draw("same");
  grQuan->Draw("same");
  legQ->Draw("same");
  latex.DrawLatex(0.2, 0.89, "#bf{Flow}");
  cQuan->SaveAs("QCumuComp.pdf");

}

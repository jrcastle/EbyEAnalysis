#include "TDecompSVD.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TCanvas.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/DAgostiniUnfold/DAgostiniUnfold.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void BottomLineTest(){

  string fin = "EllPFits_Cubic2_Test1.root";

  bool cubic = 1;

  const int norder_ = 2;

  TFile * fAna;
  TH1D * hObs[NCENT];

  TFile * fFit;
  TF1 * fEllP[NCENT];
  TF1 * fBG[NCENT];
  TH1D * hEllP[NCENT];
  TH1D * hBG[NCENT];
  TH1D * hSmearEllP[NCENT];
  TH1D * hSmearBG[NCENT];

  TFile * fUnf;
  TH1D * hUnf[NCENT];
  TH2D * hCov[NCENT];
  TH2D * hResp[NCENT];
  TH1D * hNormFactor;

  double chi2SmearEllP[NCENT];
  double chi2SmearBG[NCENT];
  double chi2UnfEllP[NCENT];
  double chi2UnfBG[NCENT];

  TGraph * grEllPChi2;
  TGraph * grBGChi2;

  TH1D * hRatio_EllPData[NCENT];
  TH1D * hRatio_BGData[NCENT];

  TLatex latex;
  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );
  fUnf = new TFile( Form("systematicStudies/SysUnfoldDistns_v%i.root", norder_) );
  fFit = 0;

  fFit = new TFile( fin.data() );


  hNormFactor = (TH1D*) fUnf->Get("hNormFactor");

  for(int icent = 3; icent < NCENT; icent++){

    //-- Get the obs hist
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    //-- Get Fits
    fEllP[icent] = (TF1*) fFit->Get( Form("fEllP_c%i", icent) );
    fBG[icent]   = (TF1*) fFit->Get( Form("fBG_c%i", icent) );

    //-- Get final unfold dist, cov, and resp matrices
    hUnf[icent]  = (TH1D*) fUnf->Get( Form("hFinalUnfoldStat_c%i", icent) );
    hUnf[icent]->SetLineColor(1);
    hUnf[icent]->SetMarkerColor(1);
    hCov[icent]  = (TH2D*) fUnf->Get( Form("finalCovMat_c%i", icent) );
    hResp[icent] = (TH2D*) fUnf->Get( Form("hresp_c%i", icent) );

    //-- Set up fit hists
    hEllP[icent] = (TH1D*) hUnf[icent]->Clone( Form("hEllP_c%i", icent) );
    hEllP[icent]->Reset();
    hEllP[icent]->SetLineColor(2);
    hBG[icent] = (TH1D*) hUnf[icent]->Clone( Form("hEllP_c%i", icent) );
    hBG[icent]->Reset();
    hBG[icent]->SetLineColor(7);

    //-- Unnormalize the unfolded hist
    hUnf[icent]->Scale(1./hNormFactor->GetBinContent(icent+1));

    //-- Convert fits to hists
    for(int i = 1; i < NBins; i++){
      double vn    = hUnf[icent]->GetBinCenter(i);
      double pEllP = fEllP[icent]->Eval( vn );
      double pBG   = fBG[icent]->Eval( vn );

      hEllP[icent]->SetBinContent(i, pEllP);
      hBG[icent]->SetBinContent(i, pBG);
    } //-- end bin loop

    //-- Ratio fit/Data:
    hRatio_EllPData[icent] = (TH1D*) hEllP[icent]->Clone( Form("hRatio_EllPData_c%i", icent) );
    hRatio_BGData[icent]   = (TH1D*) hBG[icent]->Clone( Form("hRatio_BGData_c%i", icent) );

    hRatio_EllPData[icent]->SetLineColor(2);
    hRatio_EllPData[icent]->SetMarkerColor(2);
    hRatio_EllPData[icent]->SetMarkerStyle(20);
    hRatio_EllPData[icent]->SetMinimum(0.8);
    hRatio_EllPData[icent]->SetMaximum(1.2);
    hRatio_EllPData[icent]->GetXaxis()->SetRange(1, hRatio_EllPData[icent]->FindBin(0.28));

    hRatio_BGData[icent]->SetLineColor(4);
    hRatio_BGData[icent]->SetMarkerColor(4);
    hRatio_BGData[icent]->SetMarkerStyle(21);
    hRatio_BGData[icent]->SetMinimum(0.8);
    hRatio_BGData[icent]->SetMaximum(1.2);
    hRatio_BGData[icent]->GetXaxis()->SetRange(1, hRatio_BGData[icent]->FindBin(0.28));

    hRatio_EllPData[icent]->Divide( hUnf[icent] );
    hRatio_BGData[icent]->Divide( hUnf[icent]);



    //-- Refold the fit histograms
    DAgostiniUnfold unfolder( hResp[icent] );
    hSmearEllP[icent] = (TH1D*) unfolder.refold( hEllP[icent], Form("hrefoldellp_c%i", icent) );
    hSmearBG[icent]   = (TH1D*) unfolder.refold( hBG[icent], Form("hrefoldbg_c%i", icent) );

    //-- Smeared space chi2
    chi2SmearEllP[icent] = hSmearEllP[icent]->Chi2Test(hObs[icent], "CHI2/NDF");
    chi2SmearBG[icent]   = hSmearBG[icent]->Chi2Test(hObs[icent], "CHI2/NDF");

    //-- Unfold space chi2
    /*
    TMatrixD MCov = H2M( hCov[icent] );
    TDecompSVD svd( MCov );
    TMatrixD MCovInv = svd.Invert();
    TVectorD VresEllP(NBins);
    TVectorD VresBG(NBins);
    for(int i = 0; i < NBins; i++){
      if( hUnf[icent]->GetBinContent(i+1) == 0 ) continue;
      VresEllP(i) = hUnf[icent]->GetBinContent(i+1) - hEllP[icent]->GetBinContent(i+1);
      VresBG(i) = hUnf[icent]->GetBinContent(i+1) - hBG[icent]->GetBinContent(i+1);
    }


    //-- Ellp
    TVectorD VresEllP2 = VresEllP;
    VresEllP2 *= MCovInv;
    chi2UnfEllP[icent] = VresEllP2 * VresEllP;

    //-- BG
    TVectorD VresBG2 = VresBG;
    VresBG2 *= MCovInv;
    chi2UnfBG[icent] = VresBG2 * VresBG;
    */

  } //-- End cent loop

  for(int icent = 0; icent < 3; icent++){
    chi2SmearEllP[icent] = -1.;
    chi2SmearBG[icent] = -1.;
  }

  grEllPChi2 = new TGraph(NCENT, centBinCenter, chi2SmearEllP);
  grBGChi2   = new TGraph(NCENT, centBinCenter, chi2SmearBG);

  formatGraph(grEllPChi2, "Centrality %", 0., 10., "#chi^{2}/NDF", 2, 20, "grEllPChi2");
  formatGraph(grBGChi2,   "Centrality %", 0., 10., "#chi^{2}/NDF", 7, 21, "grBGChi2");

  TFile * fOut = 0;
  if( cubic ) fOut = new TFile("SmearSpaceChi2_Cubic.root", "recreate");
  else        fOut = new TFile("SmearSpaceChi2.root", "recreate");  

  fOut->cd();
  grEllPChi2->Write();
  grBGChi2->Write();

  TCanvas * c = new TCanvas("c", "c", 500, 500);
  c->cd();
  grEllPChi2->Draw("ap");
  grBGChi2->Draw("psame");

  TCanvas * cc = new TCanvas("cc", "cc", 2000, 1500);
  cc->Divide(3,2);

  TLine * lone = new TLine(hRatio_EllPData[0]->GetXaxis()->GetXmin(), 1., hRatio_EllPData[0]->GetXaxis()->GetXmax(), 1.);
  lone->SetLineStyle(2);
  lone->SetLineWidth(2);
  std::cout<<hRatio_EllPData[0]->GetXaxis()->GetXmin()<< "\t" << hRatio_EllPData[0]->GetXaxis()->GetXmax()<<std::endl;

  int panel = 1;
  for(int icent = 3; icent < NCENT; icent++){
    if(icent == 4 || icent == 6 || icent == 8 || icent == 10) continue;
    if(panel == 3) panel++;
    cc->cd(panel);
    hRatio_EllPData[icent]->Draw();
    hRatio_BGData[icent]->Draw("same");
    lone->Draw("same");
    latex.DrawLatex(0.2, 0.8, Form("Cent. %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    panel++;
  }
  cc->SaveAs("FitRatio.pdf");

}

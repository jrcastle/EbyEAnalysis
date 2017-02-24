#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void ptResDATAComp(){

  const int cbin = 1;

  //-- DATA RECO 
  TFile * fAna_DATA_RECO;
  TH1D * hObs_DATA_RECO[NCENT];

  TFile * fUnf_DATA_RECO;
  TH1D * hUnfold_DATA_RECO[NCENT][NITER];
  TH1D * hRefold_DATA_RECO[NCENT][NITER];
  int finalUnfoldIter_DATA_RECO[NCENT];

  //-- DATA RECO_SMEAR
  TFile * fAna_DATA_RECO_SMEAR;
  TH1D * hObs_DATA_RECO_SMEAR[NCENT];

  TFile * fUnf_DATA_RECO_SMEAR;
  TH1D * hUnfold_DATA_RECO_SMEAR[NCENT][NITER];
  TH1D * hRefold_DATA_RECO_SMEAR[NCENT][NITER];
  int finalUnfoldIter_DATA_RECO_SMEAR[NCENT];

  //-- RECO/RECO_SMEAR Comparison
  TH1D * hRatioObsRECO_SMEAR_RECO_DATA[NCENT];
  TH1D * hRatioUnfRECO_SMEAR_RECO_DATA[NCENT];
  TH1D * hRatio[NCENT];

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Get Ana Files
  fAna_DATA_RECO  = new TFile("DATA_RECO/AnalyzerResults/CastleEbyE.root");
  fAna_DATA_RECO_SMEAR = new TFile("DATA_RECO_SMEAR/AnalyzerResults/CastleEbyE.root");

  //-- Get Unf files
  fUnf_DATA_RECO  = new TFile("DATA_RECO/UnfoldResults/dataResp/data2.root");
  fUnf_DATA_RECO_SMEAR = new TFile("DATA_RECO_SMEAR/UnfoldResults/dataResp/data2.root");

  //-- Cent loop
  for(int icent = 0; icent < NCENT; icent++){

    if(icent != cbin) continue;

    //-- Get Obs hists
    hObs_DATA_RECO[icent]    = (TH1D*) fAna_DATA_RECO->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs_DATA_RECO[icent]->SetLineColor(2);
    hObs_DATA_RECO[icent]->SetMarkerColor(2);
    hObs_DATA_RECO[icent]->SetMarkerStyle(20);

    hObs_DATA_RECO_SMEAR[icent]   = (TH1D*) fAna_DATA_RECO_SMEAR->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs_DATA_RECO_SMEAR[icent]->SetLineColor(2);
    hObs_DATA_RECO_SMEAR[icent]->SetMarkerColor(2);
    hObs_DATA_RECO_SMEAR[icent]->SetMarkerStyle(24);

    bool iterStop_DATA_RECO  = false;
    bool iterStop_DATA_RECO_SMEAR = false;

    //-- iter loop
    for(int i = 0; i < NITER; i++){

      double chi2;

      //-- Get unfold iteration
      hUnfold_DATA_RECO[icent][i]    = (TH1D*) fUnf_DATA_RECO->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold_DATA_RECO[icent][i]->SetLineColor(4);
      hUnfold_DATA_RECO[icent][i]->SetMarkerColor(4);
      hUnfold_DATA_RECO[icent][i]->SetMarkerStyle(21);
      hUnfold_DATA_RECO[icent][i]->Scale(1./hUnfold_DATA_RECO[icent][i]->Integral());
      hUnfold_DATA_RECO[icent][i]->SetMinimum(10e-5);
      hUnfold_DATA_RECO[icent][i]->SetMaximum(1);

      hUnfold_DATA_RECO_SMEAR[icent][i]   = (TH1D*) fUnf_DATA_RECO_SMEAR->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfold_DATA_RECO_SMEAR[icent][i]->SetLineColor(4);
      hUnfold_DATA_RECO_SMEAR[icent][i]->SetMarkerColor(4);
      hUnfold_DATA_RECO_SMEAR[icent][i]->SetMarkerStyle(25);
      hUnfold_DATA_RECO_SMEAR[icent][i]->Scale(1./hUnfold_DATA_RECO_SMEAR[icent][i]->Integral());
      hUnfold_DATA_RECO_SMEAR[icent][i]->SetMinimum(10e-5);
      hUnfold_DATA_RECO_SMEAR[icent][i]->SetMaximum(1);


      //-- Get refold iteration 
      hRefold_DATA_RECO[icent][i]  = (TH1D*) fUnf_DATA_RECO->Get( Form("hrefold%i_c%i", iter[i],icent) );
      hRefold_DATA_RECO_SMEAR[icent][i] = (TH1D*) fUnf_DATA_RECO_SMEAR->Get( Form("hrefold%i_c%i", iter[i],icent) );

      //-- DATA_RECO
      chi2 = hRefold_DATA_RECO[icent][i]->Chi2Test(hObs_DATA_RECO[icent], "CHI2/NDF");
      if( chi2 < 1.2 && !iterStop_DATA_RECO ){
        iterStop_DATA_RECO               = true;
        finalUnfoldIter_DATA_RECO[icent] = i;
      }
      if( i == NITER-1 && !iterStop_DATA_RECO ){
        iterStop_DATA_RECO               = true;
	finalUnfoldIter_DATA_RECO[icent] = i;
      }

      //-- DATA_RECO_SMEAR
      chi2 = hRefold_DATA_RECO_SMEAR[icent][i]->Chi2Test(hObs_DATA_RECO_SMEAR[icent], "CHI2/NDF");
      if( chi2 < 1.2 && !iterStop_DATA_RECO_SMEAR ){
	iterStop_DATA_RECO_SMEAR               = true;
	finalUnfoldIter_DATA_RECO_SMEAR[icent] = i;
      }
      if( i == NITER-1 && !iterStop_DATA_RECO_SMEAR ){
	iterStop_DATA_RECO_SMEAR               = true;
        finalUnfoldIter_DATA_RECO_SMEAR[icent] = i;
      }

    } //-- End iter loop


  } //-- End cent loop

  //-- =========================================================

  //-- Calculate Raios between gen and reco

  TCanvas * c[NCENT];
  TCanvas * cR[NCENT];

  TLegend * leg_Observed = new TLegend(0.182, 0.198, 0.923, 0.356);
  legInit( leg_Observed );
  leg_Observed->AddEntry(hObs_DATA_RECO[cbin],  "DATA RECO:  Observed", "lp");
  leg_Observed->AddEntry(hObs_DATA_RECO_SMEAR[cbin], "DATA RECO Smeared: Observed", "lp");

  TLegend * leg_Unfold = new TLegend(0.182, 0.198, 0.923, 0.356);
  legInit( leg_Unfold );
  leg_Unfold->AddEntry(hUnfold_DATA_RECO[cbin][0],  "DATA RECO:  Unfold", "lp");
  leg_Unfold->AddEntry(hUnfold_DATA_RECO_SMEAR[cbin][0], "DATA RECO Smeared: Unfold", "lp");


  for(int icent = 0; icent < NCENT; icent++){

    if(icent != cbin) continue;

    int i_DATA_RECO = finalUnfoldIter_DATA_RECO[icent];
    int i_DATA_RECO_SMEAR = finalUnfoldIter_DATA_RECO_SMEAR[icent];

    //-- Obs Ratio
    hObs_DATA_RECO[icent]->Scale(1./hObs_DATA_RECO[icent]->Integral());
    hObs_DATA_RECO_SMEAR[icent]->Scale(1./hObs_DATA_RECO_SMEAR[icent]->Integral());

    hRatioObsRECO_SMEAR_RECO_DATA[icent] = (TH1D*) hObs_DATA_RECO_SMEAR[icent]->Clone( Form("hRatioObsRECO_SMEAR_RECO_DATA_c%i", icent) );
    hRatioObsRECO_SMEAR_RECO_DATA[icent]->Divide( hObs_DATA_RECO[icent] );    
    hRatioObsRECO_SMEAR_RECO_DATA[icent]->GetYaxis()->SetTitle("Ratio: RECO Smeared/RECO");
    hRatioObsRECO_SMEAR_RECO_DATA[icent]->SetMinimum(0.1);
    hRatioObsRECO_SMEAR_RECO_DATA[icent]->SetMaximum(1.9);

    //-- Unfold Ratio
    hRatioUnfRECO_SMEAR_RECO_DATA[icent] = (TH1D*) hUnfold_DATA_RECO_SMEAR[icent][i_DATA_RECO_SMEAR]->Clone( Form("hRatioUnfRECO_RECO_DATA_c%i", icent) );
    hRatioUnfRECO_SMEAR_RECO_DATA[icent]->Divide( hUnfold_DATA_RECO[icent][i_DATA_RECO] );
    hRatioUnfRECO_SMEAR_RECO_DATA[icent]->GetYaxis()->SetTitle("Ratio: RECO Smeared/RECO");
    hRatioUnfRECO_SMEAR_RECO_DATA[icent]->SetMinimum(0.1);
    hRatioUnfRECO_SMEAR_RECO_DATA[icent]->SetMaximum(1.9);

    //-- Extracted parameter ratios:
    EbyECumu cumuDATA_RECO( hUnfold_DATA_RECO[icent][i_DATA_RECO] );
    EbyECumu cumuDATA_RECO_SMEAR( hUnfold_DATA_RECO_SMEAR[icent][i_DATA_RECO_SMEAR] );

    double vn2r = cumuDATA_RECO_SMEAR.GetCumu_vn2() / cumuDATA_RECO.GetCumu_vn2();
    double vn4r = cumuDATA_RECO_SMEAR.GetCumu_vn4() / cumuDATA_RECO.GetCumu_vn4();
    double vn6r = cumuDATA_RECO_SMEAR.GetCumu_vn6() / cumuDATA_RECO.GetCumu_vn6();
    double vn8r = cumuDATA_RECO_SMEAR.GetCumu_vn8() / cumuDATA_RECO.GetCumu_vn8();
    double vn6vn4r = ( cumuDATA_RECO_SMEAR.GetCumu_vn6() / cumuDATA_RECO_SMEAR.GetCumu_vn4() ) / ( cumuDATA_RECO.GetCumu_vn6() / cumuDATA_RECO.GetCumu_vn4() );
    double vn8vn4r = ( cumuDATA_RECO_SMEAR.GetCumu_vn8() / cumuDATA_RECO_SMEAR.GetCumu_vn4() ) / ( cumuDATA_RECO.GetCumu_vn8() / cumuDATA_RECO.GetCumu_vn4() );
    double vn8vn6r = ( cumuDATA_RECO_SMEAR.GetCumu_vn8() / cumuDATA_RECO_SMEAR.GetCumu_vn6() ) / ( cumuDATA_RECO.GetCumu_vn8() / cumuDATA_RECO.GetCumu_vn6() );
    double g1er = cumuDATA_RECO_SMEAR.GetGamma1Exp() / cumuDATA_RECO.GetGamma1Exp();

    hRatio[icent] = new TH1D(Form("hRatio_c%i", icent), Form("hRatio_c%i", icent), 8, 1, 8);
    hRatio[icent]->SetMinimum(0.95);
    hRatio[icent]->SetMaximum(1.02);
    hRatio[icent]->GetYaxis()->SetTitle("Ratio: RECO Smeared/RECO");

    hRatio[icent]->GetXaxis()->SetBinLabel(1, "v_{2}{2}");
    hRatio[icent]->GetXaxis()->SetBinLabel(2, "v_{2}{4}");
    hRatio[icent]->GetXaxis()->SetBinLabel(3, "v_{2}{6}");
    hRatio[icent]->GetXaxis()->SetBinLabel(4, "v_{2}{8}");
    hRatio[icent]->GetXaxis()->SetBinLabel(5, "#gamma_{1}^{exp}");
    hRatio[icent]->GetXaxis()->SetBinLabel(6, "v_{2}{6}/v_{2}{4}");
    hRatio[icent]->GetXaxis()->SetBinLabel(7, "v_{2}{8}/v_{2}{4}");
    hRatio[icent]->GetXaxis()->SetBinLabel(8, "v_{2}{8}/v_{2}{6}");

    hRatio[icent]->SetBinContent(1, vn2r);
    hRatio[icent]->SetBinContent(2, vn4r);
    hRatio[icent]->SetBinContent(3, vn6r);
    hRatio[icent]->SetBinContent(4, vn8r);
    hRatio[icent]->SetBinContent(5, g1er);
    hRatio[icent]->SetBinContent(6, vn6vn4r);
    hRatio[icent]->SetBinContent(7, vn8vn4r);
    hRatio[icent]->SetBinContent(8, vn8vn6r);

    std::cout << "v22r = " << vn2r <<std::endl;
    std::cout << "v24r = " << vn4r <<std::endl;
    std::cout << "v26r = " << vn6r <<std::endl;
    std::cout << "v28r = " << vn8r <<std::endl;
    std::cout << "v26vn4r = " << vn6vn4r <<std::endl;
    std::cout << "v28vn4r = " << vn8vn4r <<std::endl;
    std::cout << "v28vn6r = " << vn8vn6r <<std::endl;
    std::cout << "g1er = " << g1er <<std::endl;
    std::cout<<"Ratio Mean: RECO_SMEAR/RECO = " << hUnfold_DATA_RECO_SMEAR[icent][i_DATA_RECO_SMEAR]->GetMean() / hUnfold_DATA_RECO[icent][i_DATA_RECO]->GetMean() << std::endl;

    //-- DRAW
    c[icent] = new TCanvas( Form("c%i", icent), Form("c%i", icent), 1000, 1000 );
    c[icent]->Divide(2,2);

    c[icent]->cd(1);
    c[icent]->cd(1)->SetLogy();
    hObs_DATA_RECO[icent]->Draw();
    hObs_DATA_RECO_SMEAR[icent]->Draw("same");
    leg_Observed->Draw("same");
    latex.DrawLatex(0.65, 0.88, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));


    c[icent]->cd(2);
    c[icent]->cd(2)->SetLogy();
    hUnfold_DATA_RECO[icent][i_DATA_RECO]->Draw();
    hUnfold_DATA_RECO_SMEAR[icent][i_DATA_RECO_SMEAR]->Draw("same");
    leg_Unfold->Draw("same");
    latex.DrawLatex(0.65, 0.88, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));

    c[icent]->cd(3);
    hRatioObsRECO_SMEAR_RECO_DATA[icent]->Draw();
    latex.DrawLatex(0.46, 0.88, "Observed Distributions");
    latex.DrawLatex(0.65, 0.8, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));

    c[icent]->cd(4);
    hRatioUnfRECO_SMEAR_RECO_DATA[icent]->Draw();
    latex.DrawLatex(0.52, 0.88, "Unfold Distributions");
    latex.DrawLatex(0.65, 0.8, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));

    c[icent]->SaveAs( Form("PtResDATA_c%i.pdf", icent) );

    cR[icent] = new TCanvas(Form("cR_c%i", icent), Form("cR_c%i", icent), 500, 500);
    cR[icent]->cd();
    cR[icent]->SetRightMargin(0.09);
    hRatio[icent]->Draw();
    cR[icent]->SaveAs( Form("dataParmRatio_c%i.pdf", icent) );


  } //-- End cent loop


}

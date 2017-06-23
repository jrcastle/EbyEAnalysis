#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/ATLAS_PV2.h"

using namespace hi;
using namespace ebyese;
using namespace atlas_pv2;

void compareAtlasPV2(){

  const int norder_ = 2;

  //-- ATLAS
  TH1D * hATLAS[NCENT];

  //-- CMS
  TFile * fFinalUnf;
  TH1D * hFinalUnfold[NCENT];
  TH1D * hFinalUnfoldSys[NCENT];

  TH1D * hFinalUnfoldShift[NCENT];

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  fFinalUnf = new TFile( Form("systematicStudies/SysUnfoldDistns_v%i.root", norder_) );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get CMS hists
    hFinalUnfold[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldStatAndSys_c%i", icent) );
    hFinalUnfold[icent]->SetLineColor(2);
    hFinalUnfold[icent]->SetMarkerColor(2);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->SetFillColor(15);

    hFinalUnfoldSys[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldSys_c%i", icent) );
    hFinalUnfoldSys[icent]->SetMarkerColor(2);
    hFinalUnfoldSys[icent]->SetLineColor(2);
    hFinalUnfoldSys[icent]->SetMarkerStyle(20);
    hFinalUnfoldSys[icent]->SetFillColor(15);
    hFinalUnfoldSys[icent]->SetStats(0);
    hFinalUnfoldSys[icent]->GetXaxis()->SetTitleSize(0.06);
    hFinalUnfoldSys[icent]->GetXaxis()->SetLabelSize(0.06);
    hFinalUnfoldSys[icent]->GetXaxis()->SetNdivisions(507);

    hFinalUnfoldShift[icent] = (TH1D*) hFinalUnfold[icent]->Clone( Form("hFinalUnfoldShift_c%i", icent) );
    hFinalUnfoldShift[icent]->Reset();

    ATLAS_PV2[icent]->SetLineColor(1);
    ATLAS_PV2[icent]->SetMarkerColor(1);
    ATLAS_PV2[icent]->SetMarkerStyle(20);

    double atlasmax = TMath::MaxElement(ATLAS_PV2[icent]->GetN(), ATLAS_PV2[icent]->GetY());
    hFinalUnfold[icent]->Scale( atlasmax / hFinalUnfold[icent]->GetMaximum() );

    //-- Shift CMS to have the same mean as ATLAS
    double meanATLAS = ATLAS_PV2[icent]->GetMean(1);
    double meanCMS = hFinalUnfold[icent]->GetMean();
    int nb = hFinalUnfold[icent]->GetNbinsX();
    
    for(int i = 1; i <= nb; i++){
      double vn   = hFinalUnfold[icent]->GetBinCenter(i);
      double pvn  = hFinalUnfold[icent]->GetBinContent(i);
      double pvne = hFinalUnfold[icent]->GetBinError(i);
      /*
      double shiftedVn = vn + (meanATLAS - meanCMS);
      int shiftedVnBin = hFinalUnfold[icent]->FindBin( shiftedVn );

      hFinalUnfoldShift[icent]->SetBinContent(shiftedVnBin, pvn);
      hFinalUnfoldShift[icent]->SetBinError(shiftedVnBin, pvne);
      */
      double shiftedVn = vn * (meanATLAS / meanCMS);
      int shiftedVnBin = hFinalUnfold[icent]->FindBin( shiftedVn );

      hFinalUnfoldShift[icent]->Fill(shiftedVn, pvn);
      //hFinalUnfoldShift[icent]->SetBinError(shiftedVnBin, pvne);

    }

  } //-- End cent loop

  TLegend * leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  legInit( leg );
  leg->AddEntry(hFinalUnfold[0], "CMS",   "lp");
  leg->AddEntry(ATLAS_PV2[0],    "ATLAS", "lp");

  TCanvas * c = new TCanvas("c", "c", 2000, 1500);
  c->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    c->cd(icent+1);
    c->cd(icent+1)->SetLogy();
    hFinalUnfold[icent]->Draw();
    ATLAS_PV2[icent]->Draw("psame");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.65, 0.88, Form("Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );

    std::cout << Form("=========== Cent %i ===========", icent) << std::endl;
    std::cout << "MEAN CMS   = " << hFinalUnfoldShift[icent]->GetMean() << std::endl;
    std::cout << "MEAN ATLAS = " << ATLAS_PV2[icent]->GetMean(1) << std::endl;

  }

  c->SaveAs("plots/compareATLAS_pv2.pdf");


}

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/ATLAS_PV2.h"

using namespace hi;
using namespace ebyese;
using namespace atlas_pv2;

void compareAtlasPV2(){

  bool ATLASreg = 0;

  const int norder_ = 2;

  //-- CMS
  TFile * fAna;
  TH1D * hObs[NCENT];

  TFile * fUnf;
  TH1D * hRefold[NCENT][NITER];
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hFinalUnfold[NCENT];
  TH1D * hUnfold_ATLASreg[NCENT];
  TFile * fOut;

  TLatex latex3;

  //
  // MAIN
  //
  setTDRStyle();
  latex3.SetNDC();

  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );
  fUnf = new TFile( "UnfoldResults/dataResp/data2.root" );
  fOut = new TFile( "FinalUnfoldDistns.root", "recreate" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get Obs
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    //-- Get ATLAS reg hist
    hUnfold_ATLASreg[icent] = (TH1D*) fUnf->Get( Form("hreco128_c%i", icent) );

    //-- Get Unf
    for(int i = 0; i < NITER; i++){

      hUnfold[icent][i] = (TH1D*) fUnf->Get( Form("hreco%i_c%i", iter[i], icent) );
      hRefold[icent][i] = (TH1D*) fUnf->Get( Form("hrefold%i_c%i", iter[i], icent) );
      double chi2 = hRefold[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");

      if( chi2 <= 1.2 ){
	hFinalUnfold[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form( "hFinalUnfold_c%i", icent) );
	break;
      }
      if(i == NITER-1){
	hFinalUnfold[icent] = (TH1D*) hUnfold[icent][i]->Clone( Form( "hFinalUnfold_c%i", icent) );
        break;
      }

    } //-- End iter loop

    FixUnfold( hFinalUnfold[icent] );
    FixUnfold( hUnfold_ATLASreg[icent] );

    fOut->cd();
    hFinalUnfold[icent]->Write( Form("hFinalUnfoldStat_c%i", icent) );
    hUnfold_ATLASreg[icent]->Write( Form("hFinalUnfoldATLASreg_c%i", icent) );

    //-- Cosmetics
    hFinalUnfold[icent]->SetLineColor(2);
    hFinalUnfold[icent]->SetMarkerColor(2);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->GetXaxis()->SetRange(1, hFinalUnfold[icent]->FindBin(0.3));

    hUnfold_ATLASreg[icent]->SetLineColor(4);
    hUnfold_ATLASreg[icent]->SetMarkerColor(4);
    hUnfold_ATLASreg[icent]->SetMarkerStyle(20);
    hUnfold_ATLASreg[icent]->GetXaxis()->SetRange(1, hUnfold_ATLASreg[icent]->FindBin(0.3));

    ATLASPV2_Stat[icent]->SetLineColor(1);
    ATLASPV2_Stat[icent]->SetMarkerColor(1);
    ATLASPV2_Stat[icent]->SetMarkerStyle(20);
    ATLASPV2_Stat[icent]->GetXaxis()->SetNdivisions(507);
    double atlasmax = TMath::MaxElement(ATLASPV2_Stat[icent]->GetN(), ATLASPV2_Stat[icent]->GetY());
    hFinalUnfold[icent]->Scale( atlasmax / hFinalUnfold[icent]->GetMaximum() );
    hUnfold_ATLASreg[icent]->Scale( atlasmax / hUnfold_ATLASreg[icent]->GetMaximum() );

  } //-- End cent loop


  //-- CMS Refold Reg VS ATLAS Data
  TLegend * leg = new TLegend(0.42, 0.71, 0.71, 0.90);
  legInit( leg );
  leg->AddEntry(hFinalUnfold[0],  "CMS 5.02 TeV",   "lp");
  leg->AddEntry(ATLASPV2_Stat[0], "ATLAS 2.76 TeV", "lp");

  TCanvas * c = new TCanvas("c", "c", 2000, 1500);
  c->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    c->cd(icent+1);
    c->cd(icent+1)->SetLogy();
    if( ATLASreg ) hUnfold_ATLASreg[icent]->Draw();
    else           hFinalUnfold[icent]->Draw();
    ATLASPV2_Stat[icent]->Draw("psame");
    if(icent == 0) leg->Draw("same");
    latex3.DrawLatex(0.2, 0.22, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );

    ATLASPV2_Stat[icent]->GetYaxis()->SetRangeUser(ATLASPV2_Stat[icent]->GetYaxis()->GetXmin(), 2.*ATLASPV2_Stat[icent]->GetYaxis()->GetXmax() );
  }

  leg->SetTextFont(43);
  leg->SetTextSize(23);
  c->Update();
  c->SaveAs("ATLASComp.pdf");

  //-- CMS Refold Reg VS CMS ATLAS Reg
  TLegend * leg2 = new TLegend(0.42, 0.71, 0.71, 0.90);
  legInit( leg2 );
  leg2->AddEntry(hFinalUnfold[0],     "CMS Refold Reg", "lp");
  leg2->AddEntry(hUnfold_ATLASreg[0], "CMS ATLAS Reg",  "lp");

  TCanvas * c2 = new TCanvas("c2", "c2", 2000, 1500);
  c2->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    c2->cd(icent+1);
    c2->cd(icent+1)->SetLogy();
    hUnfold_ATLASreg[icent]->Draw();
    hFinalUnfold[icent]->Draw("same");
    if(icent == 0) leg2->Draw("same");
    latex3.DrawLatex(0.2, 0.22, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
  }

  leg2->SetTextFont(43);
  leg2->SetTextSize(23);
  c2->Update();
  c2->SaveAs("CMSRegComp.pdf");

}

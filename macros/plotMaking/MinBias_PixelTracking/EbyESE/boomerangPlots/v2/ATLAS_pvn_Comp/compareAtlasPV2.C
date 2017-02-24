#include "TPaveStats.h"
#include "TMath.h"
#include "TExec.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/ATLAS_PV2.h"

#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;
using namespace atlas_pv2;

void compareAtlasPV2(){

  const int norder_ = 2;

  //-- CMS
  TFile * fAna;
  TH1D * hObs[NCENT];

  TFile * fUnf;
  TH1D * hRefold[NCENT][NITER];
  TH1D * hUnfold[NCENT][NITER];
  TH1D * hFinalUnfold[NCENT];



  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  fAna = new TFile( "AnalyzerResults/CastleEbyE.root" );
  fUnf = new TFile( "UnfoldResults/dataResp/data2.root" );

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get Obs
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

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

    //-- Cosmetics
    hFinalUnfold[icent]->SetLineColor(2);
    hFinalUnfold[icent]->SetMarkerColor(2);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->SetFillColor(15);

    ATLAS_PV2[icent]->SetLineColor(1);
    ATLAS_PV2[icent]->SetMarkerColor(1);
    ATLAS_PV2[icent]->SetMarkerStyle(20);

    double atlasmax = TMath::MaxElement(ATLAS_PV2[icent]->GetN(), ATLAS_PV2[icent]->GetY());
    hFinalUnfold[icent]->Scale( atlasmax / hFinalUnfold[icent]->GetMaximum() );

    //-- Shift CMS to have the same mean as ATLAS
    /*
    double meanATLAS = ATLAS_PV2[icent]->GetMean(1);
    double meanCMS = hFinalUnfold[icent]->GetMean();
    int nb = hFinalUnfold[icent]->GetNbinsX();
    
    for(int i = 1; i <= nb; i++){
      double vn   = hFinalUnfold[icent]->GetBinCenter(i);
      double pvn  = hFinalUnfold[icent]->GetBinContent(i);
      double pvne = hFinalUnfold[icent]->GetBinError(i);
      double shiftedVn = vn * (meanATLAS / meanCMS);
      int shiftedVnBin = hFinalUnfold[icent]->FindBin( shiftedVn );
      hFinalUnfoldShift[icent]->Fill(shiftedVn, pvn);
    }
    */
  } //-- End cent loop

  //TLegend * leg = new T

  TCanvas * c = new TCanvas("c", "c", 2000, 1500);
  c->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    c->cd(icent+1);
    c->cd(icent+1)->SetLogy();
    hFinalUnfold[icent]->Draw();
    ATLAS_PV2[icent]->Draw("psame");

    std::cout << Form("=========== Cent %i ===========", icent) << std::endl;
    std::cout << "MEAN CMS   = " << hFinalUnfoldShift[icent]->GetMean() << std::endl;
    std::cout << "MEAN ATLAS = " << ATLAS_PV2[icent]->GetMean(1) << std::endl;

  }




}

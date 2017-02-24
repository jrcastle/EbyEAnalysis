#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void ellpCorrelations(){

  TFile * fFit;

  const int Nsig = 3;
  TGraph * grKnE0[Nsig][NCENT];
  TGraph * grKnAlpha[Nsig][NCENT];
  TGraph * grAlphaE0[Nsig][NCENT];
  int contCol[Nsig] = {1, 2, 4};

  TCanvas * cEllPCor[NCENT];
  TLegend * legCont;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();


  fFit = new TFile("EllPFits.root");

  for(int icent = 3; icent < NCENT; icent++){

    cEllPCor[icent] = new TCanvas( Form("cEllPCor_c%i", icent), Form("cEllPCor_c%i", icent), 1500, 500);
    cEllPCor[icent]->Divide(3,1);

    if(icent == 3){
      legCont = new TLegend(0.2, 0.6, 0.5, 0.9);
      legInit( legCont );
    }

    for(int isig = 0; isig < Nsig; isig++){
      //grKnE0[isig][icent]    = (TGraph*) fFit->Get( Form("grKnE0_%is_c%i", isig+1, icent) );
      //grKnAlpha[isig][icent] = (TGraph*) fFit->Get( Form("grKnAlpha_%is_c%i", isig+1, icent) );
      grAlphaE0[isig][icent] = (TGraph*) fFit->Get( Form("grAlphaE0_%is_c%i", isig+1, icent) );
      //if(icent == 3) legCont->AddEntry(grKnE0[isig][icent], Form("%i#sigma", isig+1), "l");
    } //-- end sig loop

  } // end cent loop

  //-- DRAW
  for(int icent = 3; icent < NCENT; icent++){

    for(int isig = Nsig-1; isig >= 0; isig--){
      cEllPCor[icent]->cd(1);
      /*
      if(isig == Nsig-1) grKnE0[isig][icent]->Draw("alf");
      else               grKnE0[isig][icent]->Draw("lfsame");

      cEllPCor[icent]->cd(2);
      if(isig == Nsig-1) grKnAlpha[isig][icent]->Draw("alf");
      else               grKnAlpha[isig][icent]->Draw("lfsame");
      */
      cEllPCor[icent]->cd(3);
      if(isig == Nsig-1) grAlphaE0[isig][icent]->Draw("alf");
      else               grAlphaE0[isig][icent]->Draw("lfsame");
    }

    cEllPCor[icent]->cd(1);
    latex.DrawLatex(0.65, 0.88, Form("Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );

    cEllPCor[icent]->cd(2);
    legCont->Draw("same");

    cEllPCor[icent]->SaveAs( Form("plots/fitCorr/cEllPCor_c%i.pdf", icent ) );

  }
}

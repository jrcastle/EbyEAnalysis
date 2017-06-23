#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace hi;
using namespace ebyese;

void ellpCorrelations(){

  bool CUBIC   = 0;
  int SCENARIO = 2;

  double knMin     = 0.;
  double knMax     = 0.75;
  double knMin2    = 0.27;
  double knMax2    = 0.4;
  double knPrMin   = -0.01;
  double knPrMax   = 0.15;
  double alphaMin  = 0.;
  double alphaMax  = 130;
  double alphaMin2 = 0.;
  double alphaMax2 = 80;
  double E0Min     = 0.;
  double E0Max     = 0.5;
  double E0Min2    = 0.15;
  double E0Max2    = 0.2999;

  TFile * fFit;

  const int Nsig = 3;
  TGraph * grE0Kn[Nsig][NCENT];
  TGraph * grE0KnPr[Nsig][NCENT];
  TGraph * grE0Alpha[Nsig][NCENT];
  TGraph * grAlphaKn[Nsig][NCENT];
  TGraph * grAlphaKnPr[Nsig][NCENT];
  TGraph * grKnKnPr[Nsig][NCENT];

  TGraph * grE0Kn_Dummy;
  TGraph * grE0Kn_Dummy2;
  TGraph * grE0KnPr_Dummy;
  TGraph * grE0Alpha_Dummy;
  TGraph * grE0Alpha_Dummy2;
  TGraph * grAlphaKn_Dummy;
  TGraph * grAlphaKn_Dummy2;
  TGraph * grAlphaKnPr_Dummy;
  TGraph * grKnKnPr_Dummy;

  int contCol[Nsig] = {1, 2, 4};
  TCanvas * cEllPCor[NCENT];
  TLegend * legCont;

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  fFit = 0;
  if( CUBIC && SCENARIO == 1)      fFit = new TFile("EllPFits_Cubic1.root");
  else if( CUBIC && SCENARIO == 2) fFit = new TFile("EllPFits_Cubic2.root");
  else                             fFit = new TFile("EllPFits.root");

  for(int icent = 3; icent < NCENT; icent++){

    cEllPCor[icent] = 0;
    if( CUBIC ){
      cEllPCor[icent] = new TCanvas( Form("cEllPCor_c%i", icent), Form("cEllPCor_c%i", icent), 1500, 1000);
      cEllPCor[icent]->Divide(3,2);
    }
    else{
      cEllPCor[icent] = new TCanvas( Form("cEllPCor_c%i", icent), Form("cEllPCor_c%i", icent), 1500, 500);
      cEllPCor[icent]->Divide(3,1);
    }
    if(icent == 3){
      legCont = new TLegend(0.2, 0.6, 0.5, 0.9);
      legInit( legCont );
    }

    for(int isig = 0; isig < Nsig; isig++){
      grE0Kn[isig][icent]    = (TGraph*) fFit->Get( Form("grE0Kn_%is_c%i", isig+1, icent) );
      grAlphaKn[isig][icent] = (TGraph*) fFit->Get( Form("grAlphaKn_%is_c%i", isig+1, icent) );
      grE0Alpha[isig][icent] = (TGraph*) fFit->Get( Form("grE0Alpha_%is_c%i", isig+1, icent) );

      if( CUBIC ){
	grE0KnPr[isig][icent]    = (TGraph*) fFit->Get( Form("grE0KnPr_%is_c%i", isig+1, icent) );
	grAlphaKnPr[isig][icent] = (TGraph*) fFit->Get( Form("grAlphaKnPr_%is_c%i", isig+1, icent) );
	grKnKnPr[isig][icent]    = (TGraph*) fFit->Get( Form("grKnKnPr_%is_c%i", isig+1, icent) );
      }

      if(icent == 3) legCont->AddEntry(grE0Kn[isig][icent], Form("%i#sigma", isig+1), "l");
    } //-- end sig loop

  } // end cent loop

  //-- Make dummies for cases where the contours cannot be calculated
  double x[1] = {0.1};
  double y[1] = {-1.};
  grE0Kn_Dummy      = new TGraph(1,x,y);
  grE0Kn_Dummy2     = new TGraph(1,x,y);
  grE0KnPr_Dummy    = new TGraph(1,x,y);
  grE0Alpha_Dummy   = new TGraph(1,x,y);
  grE0Alpha_Dummy2  = new TGraph(1,x,y);
  grAlphaKn_Dummy   = new TGraph(1,x,y);
  grAlphaKn_Dummy2  = new TGraph(1,x,y);
  grAlphaKnPr_Dummy = new TGraph(1,x,y);
  grKnKnPr_Dummy    = new TGraph(1,x,y);

  formatGraph(grE0Kn_Dummy,      "k_{2}",    E0Min,    E0Max,    "#epsilon_{0}", 1, 20, "grE0Kn_Dummy");
  formatGraph(grE0Kn_Dummy2,     "k_{2}",    E0Min2,   E0Max2,   "#epsilon_{0}", 1, 20, "grE0Kn_Dummy2");
  formatGraph(grE0KnPr_Dummy,    "#kappa\'", E0Min,    E0Max,    "#epsilon_{0}", 1, 20, "grE0KnPr_Dummy");
  formatGraph(grE0Alpha_Dummy,   "#alpha",   E0Min,    E0Max,    "#epsilon_{0}", 1, 20, "grE0Alpha_Dummy");
  formatGraph(grE0Alpha_Dummy2,  "#alpha",   E0Min2,   E0Max2,   "#epsilon_{0}", 1, 20, "grE0Alpha_Dummy2");
  formatGraph(grAlphaKn_Dummy,   "k_{2}",    alphaMin, alphaMax, "#alpha",       1, 20, "grAlphaKn_Dummy");
  formatGraph(grAlphaKn_Dummy2,  "k_{2}",    alphaMin2, alphaMax2, "#alpha",       1, 20, "grAlphaKn_Dummy2");
  formatGraph(grAlphaKnPr_Dummy, "#kappa\'", alphaMin, alphaMax, "#alpha",       1, 20, "grAlphaKnPr_Dummy");
  formatGraph(grKnKnPr_Dummy,    "#kappa\'", knMin,    knMax,    "k_{2}",        1, 20, "grKnKnPr_Dummy");

  grE0Kn_Dummy      -> GetXaxis() -> SetLimits(knMin,     knMax);
  grE0Kn_Dummy2     -> GetXaxis() -> SetLimits(knMin2,    knMax2);
  grE0KnPr_Dummy    -> GetXaxis() -> SetLimits(knPrMin,   knPrMax);
  grE0Alpha_Dummy   -> GetXaxis() -> SetLimits(alphaMin,  alphaMax);
  grE0Alpha_Dummy2  -> GetXaxis() -> SetLimits(alphaMin2, alphaMax2);
  grAlphaKn_Dummy   -> GetXaxis() -> SetLimits(knMin,     knMax);
  grAlphaKn_Dummy2  -> GetXaxis() -> SetLimits(knMin2,    knMax2);
  grAlphaKnPr_Dummy -> GetXaxis() -> SetLimits(knPrMin,   knPrMax);
  grKnKnPr_Dummy    -> GetXaxis() -> SetLimits(knPrMin,   knPrMax);

  //-- DRAW
  for(int icent = 3; icent < NCENT; icent++){

    bool firsSigE0Kn      = false;
    bool firsSigE0KnPr    = false;
    bool firsSigE0Alpha   = false;
    bool firsSigAlphaKn   = false;
    bool firsSigAlphaKnPr = false;
    bool firsSigKnKnPr    = false;

    for(int isig = Nsig-1; isig >= 0; isig--){

      if( CUBIC ){
	cEllPCor[icent]->cd(1);
	if( grE0Kn[isig][icent] ){
	  if( !firsSigE0Kn ){
	    grE0Kn[isig][icent]->Draw("alf");
	    firsSigE0Kn = true;
	  }
	  else grE0Kn[isig][icent]->Draw("lfsame");
	}

	cEllPCor[icent]->cd(2);
	if( grE0KnPr[isig][icent] ){
	  if( !firsSigE0KnPr ){
	    grE0KnPr[isig][icent]->Draw("alf");
	    firsSigE0KnPr = true;
	  }
	  else grE0KnPr[isig][icent]->Draw("lfsame");
	}

	cEllPCor[icent]->cd(3);
	if( grE0Alpha[isig][icent] ){
	  if( !firsSigE0Alpha ){
	    grE0Alpha[isig][icent]->Draw("alf");
	    firsSigE0Alpha = true;
	  }
	  else grE0Alpha[isig][icent]->Draw("lfsame");
	}

        cEllPCor[icent]->cd(4);
	if( grAlphaKn[isig][icent] ){
	  if( !firsSigAlphaKn ){
	    grAlphaKn[isig][icent]->Draw("alf");
	    firsSigAlphaKn = true;
	  }
	  else grAlphaKn[isig][icent]->Draw("lfsame");
	}

	cEllPCor[icent]->cd(5);
	if( grAlphaKnPr[isig][icent] ){
	  if( !firsSigAlphaKnPr ){
	    grAlphaKnPr[isig][icent]->Draw("alf");
	    firsSigAlphaKnPr = true;
	  }
	  else grAlphaKnPr[isig][icent]->Draw("lfsame");
	}

        cEllPCor[icent]->cd(6);
	if( grKnKnPr[isig][icent] ){
	  if( !firsSigKnKnPr ){
	    grKnKnPr[isig][icent]->Draw("alf");
	    firsSigKnKnPr = true;
	  }
	  else grKnKnPr[isig][icent]->Draw("lfsame");
	}
      } //-- End if( CUBIC )
      else{
	cEllPCor[icent]->cd(1);
	if(isig == Nsig-1) grE0Kn[isig][icent]->Draw("alf");
	else               grE0Kn[isig][icent]->Draw("lfsame");

	cEllPCor[icent]->cd(2);
	if(isig == Nsig-1) grAlphaKn[isig][icent]->Draw("alf");
	else               grAlphaKn[isig][icent]->Draw("lfsame");
	
	cEllPCor[icent]->cd(3);
	if(isig == Nsig-1) grE0Alpha[isig][icent]->Draw("alf");
	else               grE0Alpha[isig][icent]->Draw("lfsame");
      } //-- End isig loop

      if( CUBIC ){

	if( !grE0Kn[0][icent] && !grE0Kn[1][icent] && !grE0Kn[2][icent] ){
	  cEllPCor[icent]->cd(1);
	  grE0Kn_Dummy->Draw("ap");
	  latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
	if( !grE0KnPr[0][icent] && !grE0KnPr[1][icent] && !grE0KnPr[2][icent] ){
          cEllPCor[icent]->cd(2);
          grE0KnPr_Dummy->Draw("ap");
          latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
	if( !grE0Alpha[0][icent] && !grE0Alpha[1][icent] && !grE0Alpha[2][icent] ){
          cEllPCor[icent]->cd(3);
          grE0Alpha_Dummy->Draw("ap");
          latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
	if( !grAlphaKn[0][icent] && !grAlphaKn[1][icent] && !grAlphaKn[2][icent] ){
          cEllPCor[icent]->cd(4);
          grAlphaKn_Dummy->Draw("ap");
          latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
	if( !grAlphaKnPr[0][icent] && !grAlphaKnPr[1][icent] && !grAlphaKnPr[2][icent] ){
          cEllPCor[icent]->cd(5);
          grAlphaKnPr_Dummy->Draw("ap");
          latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
	if( !grKnKnPr[0][icent] && !grKnKnPr[1][icent] && !grKnKnPr[2][icent] ){
          cEllPCor[icent]->cd(6);
          grKnKnPr_Dummy->Draw("ap");
          latex.DrawLatex(0.25, 0.5, "Unable to determine contours");
	}
      } //-- End if( CUBIC )

    } //-- End cent loop

    cEllPCor[icent]->cd(1);
    latex.DrawLatex(0.65, 0.88, Form("Cent %i-%i%s", cent_min[icent], cent_max[icent], "%") );

    cEllPCor[icent]->cd(2);
    legCont->Draw("same");

    if( CUBIC && SCENARIO == 1) cEllPCor[icent]->SaveAs( Form("Contours/Cubic/cEllPCor_Cubic1_c%i.pdf", icent ) );
    if( CUBIC && SCENARIO == 2) cEllPCor[icent]->SaveAs( Form("Contours/Cubic/cEllPCor_Cubic2_c%i.pdf", icent ) );
    if( !CUBIC )                cEllPCor[icent]->SaveAs( Form("Contours/Linear/cEllPCor_c%i.pdf", icent ) );

  }

  //-- Big Corr plot E0 Vs Alpha
  if( !CUBIC ){
    TCanvas * cContourMerge = new TCanvas("cContourMerge", "cContourMerge", 1500, 500);
    cContourMerge->Divide(3,1);

    cContourMerge->cd(1);
    grE0Alpha_Dummy2->Draw("ap");
    for(int icent = 3; icent < NCENT; icent++){
      for(int isig = Nsig-1; isig >= 0; isig--){
	grE0Alpha[isig][icent]->Draw("lfsame");
      }
    }

    cContourMerge->cd(2);
    grE0Kn_Dummy2->Draw("ap");
    for(int icent = 3; icent < NCENT; icent++){
      for(int isig = Nsig-1; isig >= 0; isig--){
        grE0Kn[isig][icent]->Draw("lfsame");
      }
    }

    cContourMerge->cd(3);
    grAlphaKn_Dummy2->Draw("ap");
    for(int icent = 3; icent < NCENT; icent++){
      for(int isig = Nsig-1; isig >= 0; isig--){
        grAlphaKn[isig][icent]->Draw("lfsame");
      }
    }

    cContourMerge->SaveAs("Contours/Linear/cContourMerge.pdf");

  }


}

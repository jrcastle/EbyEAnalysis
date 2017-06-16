#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void compareObs(){

  bool full = 1;
  bool sub0 = 0;
  bool sub1 = 0;

  TFile * vtx15;
  TH1D * hObs_vtx15[NCENT];
  TH1D * hObs_0_vtx15[NCENT];
  TH1D * hObs_1_vtx15[NCENT];
  TH1D * hRat01_vtx15[NCENT];
  TH2D * hSmearVec_vtx15[NCENT];

  TFile * vtx3;
  TH1D * hObs_vtx3[NCENT];
  TH1D * hObs_0_vtx3[NCENT];
  TH1D * hObs_1_vtx3[NCENT];
  TH1D * hObs_vtx3_Rto15[NCENT];
  TH1D * hObs_0_vtx3_Rto15[NCENT];
  TH1D * hObs_1_vtx3_Rto15[NCENT];
  TH1D * hRat01_vtx3[NCENT];
  TH2D * hSmearVec_vtx3[NCENT];

  TFile * vtx3_15;
  TH1D * hObs_vtx3_15[NCENT];
  TH1D * hObs_0_vtx3_15[NCENT];
  TH1D * hObs_1_vtx3_15[NCENT];
  TH1D * hObs_vtx3_15_Rto15[NCENT];
  TH1D * hObs_0_vtx3_15_Rto15[NCENT];
  TH1D * hObs_1_vtx3_15_Rto15[NCENT];
  TH1D * hRat01_vtx3_15[NCENT];
  TH2D * hSmearVec_vtx3_15[NCENT];

  TH1D * hObs_vtx3_Rto3_15[NCENT];
  TH1D * hObs_0_vtx3_Rto3_15[NCENT];
  TH1D * hObs_1_vtx3_Rto3_15[NCENT];

  TLatex latex;
  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  vtx15   = new TFile("AnalyzerResults/CastleEbyE.root");
  vtx3    = new TFile("AnalyzerResults_vtx_leq_3/CastleEbyE.root");
  vtx3_15 = new TFile("AnalyzerResults_vtx3_15/CastleEbyE.root");

  TCanvas* c[NCENT];
  TLegend * l = new TLegend(0.7, 0.7, 0.9, 0.9);
  legInit(l);
  TLine * lin = new TLine(0, 1, 0.6, 1);
  lin->SetLineStyle(2);

  for(int icent = 0; icent < NCENT; icent++){
    if(icent > 0) break;
    // 15
    hObs_vtx15[icent]   = (TH1D*) vtx15->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs_0_vtx15[icent] = (TH1D*) vtx15->Get( Form("qwebye/hVnSub0_c%i", icent) );
    hObs_1_vtx15[icent] = (TH1D*) vtx15->Get( Form("qwebye/hVnSub1_c%i", icent) );

    hObs_vtx15[icent]->Scale(1./hObs_vtx15[icent]->Integral());
    hObs_0_vtx15[icent]->Scale(1./hObs_0_vtx15[icent]->Integral());
    hObs_1_vtx15[icent]->Scale(1./hObs_1_vtx15[icent]->Integral());

    hObs_vtx15[icent]->SetLineColor(1);
    hObs_0_vtx15[icent]->SetLineColor(1);
    hObs_1_vtx15[icent]->SetLineColor(1);

    hObs_vtx15[icent]->SetMarkerColor(1);
    hObs_0_vtx15[icent]->SetMarkerColor(1);
    hObs_1_vtx15[icent]->SetMarkerColor(1);

    hRat01_vtx15[icent] = (TH1D*) hObs_0_vtx15[icent]->Clone( Form("hRat01_vtx15_c%i", icent) );
    hRat01_vtx15[icent]->Divide(hObs_1_vtx15[icent]);
    hRat01_vtx15[icent]->SetMinimum(0);  hRat01_vtx15[icent]->SetMaximum(2);

    // 3
    hObs_vtx3[icent]   = (TH1D*) vtx3->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs_0_vtx3[icent] = (TH1D*) vtx3->Get( Form("qwebye/hVnSub0_c%i", icent) );
    hObs_1_vtx3[icent] = (TH1D*) vtx3->Get( Form("qwebye/hVnSub1_c%i", icent) );

    hObs_vtx3[icent]->Scale(1./hObs_vtx3[icent]->Integral());
    hObs_0_vtx3[icent]->Scale(1./hObs_0_vtx3[icent]->Integral());
    hObs_1_vtx3[icent]->Scale(1./hObs_1_vtx3[icent]->Integral());

    hObs_vtx3[icent]->SetLineColor(2);
    hObs_0_vtx3[icent]->SetLineColor(2);
    hObs_1_vtx3[icent]->SetLineColor(2);

    hObs_vtx3[icent]->SetMarkerColor(2);
    hObs_0_vtx3[icent]->SetMarkerColor(2);
    hObs_1_vtx3[icent]->SetMarkerColor(2);

    hObs_vtx3_Rto15[icent]   = (TH1D*) hObs_vtx3[icent]->Clone(Form("hObs_vtx3_Rto15_c%i", icent));
    hObs_0_vtx3_Rto15[icent] = (TH1D*) hObs_0_vtx3[icent]->Clone(Form("hObs_0_vtx3_Rto15_c%i", icent));
    hObs_1_vtx3_Rto15[icent] = (TH1D*) hObs_1_vtx3[icent]->Clone(Form("hObs_1_vtx3_Rto15_c%i", icent));

    hObs_vtx3_Rto15[icent]->Divide( hObs_vtx15[icent] );
    hObs_0_vtx3_Rto15[icent]->Divide( hObs_0_vtx15[icent] );
    hObs_1_vtx3_Rto15[icent]->Divide( hObs_1_vtx15[icent] );

    hObs_vtx3_Rto15[icent]->SetMinimum(0); hObs_vtx3_Rto15[icent]->SetMaximum(2);
    hObs_0_vtx3_Rto15[icent]->SetMinimum(0); hObs_0_vtx3_Rto15[icent]->SetMaximum(2);
    hObs_1_vtx3_Rto15[icent]->SetMinimum(0); hObs_1_vtx3_Rto15[icent]->SetMaximum(2);

    hRat01_vtx3[icent] = (TH1D*) hObs_0_vtx3[icent]->Clone( Form("hRat01_vtx3_c%i", icent) );
    hRat01_vtx3[icent]->Divide(hObs_1_vtx3[icent]);
    hRat01_vtx3[icent]->SetMinimum(0); hRat01_vtx3[icent]->SetMaximum(2);

    // 3-15
    hObs_vtx3_15[icent]   = (TH1D*) vtx3_15->Get( Form("qwebye/hVnFull_c%i", icent) );
    hObs_0_vtx3_15[icent] = (TH1D*) vtx3_15->Get( Form("qwebye/hVnSub0_c%i", icent) );
    hObs_1_vtx3_15[icent] = (TH1D*) vtx3_15->Get( Form("qwebye/hVnSub1_c%i", icent) );

    hObs_vtx3_15[icent]->Scale(1./hObs_vtx3_15[icent]->Integral());
    hObs_0_vtx3_15[icent]->Scale(1./hObs_0_vtx3_15[icent]->Integral());
    hObs_1_vtx3_15[icent]->Scale(1./hObs_1_vtx3_15[icent]->Integral());

    hObs_vtx3_15[icent]->SetLineColor(4);
    hObs_0_vtx3_15[icent]->SetLineColor(4);
    hObs_1_vtx3_15[icent]->SetLineColor(4);

    hObs_vtx3_15[icent]->SetMarkerColor(4);
    hObs_0_vtx3_15[icent]->SetMarkerColor(4);
    hObs_1_vtx3_15[icent]->SetMarkerColor(4);

    hObs_vtx3_15_Rto15[icent]   = (TH1D*) hObs_vtx3_15[icent]->Clone(Form("hObs_vtx3_15_Rto15_c%i", icent));
    hObs_0_vtx3_15_Rto15[icent] = (TH1D*) hObs_0_vtx3_15[icent]->Clone(Form("hObs_0_vtx3_15_Rto15_c%i", icent));
    hObs_1_vtx3_15_Rto15[icent] = (TH1D*) hObs_1_vtx3_15[icent]->Clone(Form("hObs_1_vtx3_15_Rto15_c%i", icent));

    hObs_vtx3_15_Rto15[icent]->Divide( hObs_vtx15[icent] );
    hObs_0_vtx3_15_Rto15[icent]->Divide( hObs_0_vtx15[icent] );
    hObs_1_vtx3_15_Rto15[icent]->Divide( hObs_1_vtx15[icent] );

    hObs_vtx3_15_Rto15[icent]->SetMinimum(0); hObs_vtx3_15_Rto15[icent]->SetMaximum(2);
    hObs_0_vtx3_15_Rto15[icent]->SetMinimum(0); hObs_0_vtx3_15_Rto15[icent]->SetMaximum(2);
    hObs_1_vtx3_15_Rto15[icent]->SetMinimum(0); hObs_1_vtx3_15_Rto15[icent]->SetMaximum(2);

    hRat01_vtx3_15[icent] = (TH1D*) hObs_0_vtx3_15[icent]->Clone( Form("hRat01_vtx3_15_c%i", icent) );
    hRat01_vtx3_15[icent]->Divide(hObs_1_vtx3_15[icent]);
    hRat01_vtx3_15[icent]->SetMinimum(0); hRat01_vtx3_15[icent]->SetMaximum(2);

    // 3 / 3-15
    hObs_vtx3_Rto3_15[icent]   = (TH1D*) hObs_vtx3[icent]->Clone( Form("hObs_vtx3_Rto3_15_c%i", icent) );
    hObs_0_vtx3_Rto3_15[icent] = (TH1D*) hObs_0_vtx3[icent]->Clone( Form("hObs_0_vtx3_Rto3_15_c%i", icent) );
    hObs_1_vtx3_Rto3_15[icent] = (TH1D*) hObs_1_vtx3[icent]->Clone( Form("hObs_1_vtx3_Rto3_15_c%i", icent) );

    hObs_vtx3_Rto3_15[icent]->Divide( hObs_vtx3_15[icent] );
    hObs_0_vtx3_Rto3_15[icent]->Divide( hObs_0_vtx3_15[icent] );
    hObs_1_vtx3_Rto3_15[icent]->Divide( hObs_1_vtx3_15[icent] );

    hObs_vtx3_Rto3_15[icent]->SetMinimum(0); hObs_vtx3_Rto3_15[icent]->SetMaximum(2);
    hObs_0_vtx3_Rto3_15[icent]->SetMinimum(0); hObs_0_vtx3_Rto3_15[icent]->SetMaximum(2);
    hObs_1_vtx3_Rto3_15[icent]->SetMinimum(0); hObs_1_vtx3_Rto3_15[icent]->SetMaximum(2);


    std::cout << Form("Vz < 15 NEvents     = %.1f", hObs_vtx15[icent]->GetEntries()) << std::endl; 
    std::cout << Form("Vz < 3 NEvents      = %.1f", hObs_vtx3[icent]->GetEntries()) << std::endl;
    std::cout << Form("3 < Vz < 15 NEvents = %.1f", hObs_vtx3_15[icent]->GetEntries()) << std::endl;

    //-- Smear vec
    hSmearVec_vtx15[icent]   = (TH2D*) vtx15->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );
    hSmearVec_vtx3[icent]    = (TH2D*) vtx3->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );
    hSmearVec_vtx3_15[icent] = (TH2D*) vtx3_15->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );

    double smearW15 = ( hSmearVec_vtx15[icent]->ProjectionX()->GetRMS() + hSmearVec_vtx15[icent]->ProjectionY()->GetRMS() ) / 2.;
    double smearW3= ( hSmearVec_vtx3[icent]->ProjectionX()->GetRMS() +hSmearVec_vtx3[icent]->ProjectionY()->GetRMS()) / 2.;
    double smearW3_15= ( hSmearVec_vtx3_15[icent]->ProjectionX()->GetRMS() +hSmearVec_vtx3_15[icent]->ProjectionY()->GetRMS()) / 2.;

    std::cout << Form("Vz < 15 SmearW    = %.4f", smearW15) << std::endl;
    std::cout << Form("Vz < 3 SmearW     = %.4f", smearW3) << std::endl;
    std::cout << Form("3< Vz < 15 SmearW = %.4f", smearW3_15) << std::endl;

    if(icent == 0){
      l->AddEntry(hObs_vtx15[icent],   "|v_{z}| < 15 cm");
      l->AddEntry(hObs_vtx3[icent],    "|v_{z}| < 3 cm");
      l->AddEntry(hObs_vtx3_15[icent], "3 < |v_{z}| < 15 cm");
    }

    if( full ){
      c[icent] = new TCanvas(Form("c%i", icent),Form("c%i", icent),1000,1000);
      c[icent]->Divide(2,2);
      c[icent]->cd(1);
      c[icent]->SetLogy();
      hObs_vtx15[icent]->Draw();
      hObs_vtx3[icent]->Draw("same");
      hObs_vtx3_15[icent]->Draw("same");
      l->Draw("same");
      latex.DrawLatex(0.2, 0.8, "p(v_{n}^{obs}) Full Event");


      c[icent]->cd(3);
      hObs_0_vtx3_Rto15[icent]->Draw();
      hObs_0_vtx3_15_Rto15[icent]->Draw("same");
      latex.DrawLatex(0.2, 0.8,"Ratio to |v_{z}| < 15 Full Event");
      lin->Draw("same");

      c[icent]->cd(2);
      hObs_vtx3_Rto3_15[icent]->Draw();
      latex.DrawLatex(0.2, 0.8,"#frac{|v_{z}| < 3}{3 < |v_{z}| < 15} Full Event");
      lin->Draw("same");

      c[icent]->cd(4);
      hRat01_vtx15[icent]->Draw();
      hRat01_vtx3[icent]->Draw("same");
      hRat01_vtx3_15[icent]->Draw("same");
      lin->Draw("same");
      latex.DrawLatex(0.2, 0.8,"sub0 / sub1");

    }
    if( sub0 ){
      c[icent] = new TCanvas(Form("c%i", icent),Form("c%i", icent),1000,1000);
      c[icent]->Divide(2,2);
      c[icent]->cd(1);
      c[icent]->SetLogy();
      hObs_0_vtx15[icent]->Draw();
      hObs_0_vtx3[icent]->Draw("same");
      hObs_0_vtx3_15[icent]->Draw("same");
      l->Draw("same");
      latex.DrawLatex(0.2, 0.8, "p(v_{n}^{obs}) Subevent0");


      c[icent]->cd(3);
      hObs_0_vtx3_Rto15[icent]->Draw();
      hObs_0_vtx3_15_Rto15[icent]->Draw("same");
      latex.DrawLatex(0.2, 0.8,"Ratio to |v_{z}| < 15 Subevent0");
      lin->Draw("same");

      c[icent]->cd(2);
      hObs_0_vtx3_Rto3_15[icent]->Draw();
      latex.DrawLatex(0.2, 0.8,"#frac{|v_{z}| < 3}{3 < |v_{z}| < 15} Subevent0");
      lin->Draw("same");
    }
    if( sub1 ){
      c[icent] = new TCanvas(Form("c%i", icent),Form("c%i", icent),1000,1000);
      c[icent]->Divide(2,2);
      c[icent]->cd(1);
      c[icent]->SetLogy();
      hObs_1_vtx15[icent]->Draw();
      hObs_1_vtx3[icent]->Draw("same");
      hObs_1_vtx3_15[icent]->Draw("same");
      l->Draw("same");
      latex.DrawLatex(0.2, 0.8, "p(v_{n}^{obs}) Subevent1");


      c[icent]->cd(3);
      hObs_1_vtx3_Rto15[icent]->Draw();
      hObs_1_vtx3_15_Rto15[icent]->Draw("same");
      latex.DrawLatex(0.2, 0.8,"Ratio to |v_{z}| < 15 Subevent1");
      lin->Draw("same");

      c[icent]->cd(2);
      hObs_1_vtx3_Rto3_15[icent]->Draw();
      latex.DrawLatex(0.2, 0.8,"#frac{|v_{z}| < 3}{3 < |v_{z}| < 15} Subevent1");
      lin->Draw("same");
    }



  }








}

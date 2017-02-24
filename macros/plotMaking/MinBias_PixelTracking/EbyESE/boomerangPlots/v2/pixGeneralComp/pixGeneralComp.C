#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"

#include "/home/j550c590/tdrstyle.C"

#include <iostream>


void pixGeneralComp(){

  string clusTrunc = "2p0";

  const int N = 11;
  const int cmin[N] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int cmax[N] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const int ccol[N] = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};


  TFile * fMult;
  TH2D * hMultCent_eta24_pt03_3_hiPixelAndGeneral;
  TH2D * hMultCent_eta10_pt03_3_hiPixelAndGeneral;
  TH2D * hMultCent_eta24_pt03_3_hiPixel;
  TH2D * hMultCent_eta10_pt03_3_hiPixel;
  TH2D * hMultCent_eta24_pt03_3_hiGeneral;
  TH2D * hMultCent_eta10_pt03_3_hiGeneral;
  
  TH1D * hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[N];
  TH1D * hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[N];
  TH1D * hMultCent_eta24_pt03_3_hiPixel_proj[N];
  TH1D * hMultCent_eta10_pt03_3_hiPixel_proj[N];
  TH1D * hMultCent_eta24_pt03_3_hiGeneral_proj[N];
  TH1D * hMultCent_eta10_pt03_3_hiGeneral_proj[N];

  TLatex latex;

  //
  // MAIN
  //
  TH1D::SetDefaultSumw2();
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  fMult = new TFile( Form("MultCent_ClusVtx%s.root", clusTrunc.data()) );
  hMultCent_eta24_pt03_3_hiPixelAndGeneral = (TH2D*) fMult->Get( "multcentana/hMultCent_eta24_pt03_3_hiPixelAndGeneral" );
  hMultCent_eta10_pt03_3_hiPixelAndGeneral = (TH2D*) fMult->Get( "multcentana/hMultCent_eta10_pt03_3_hiPixelAndGeneral" );
  hMultCent_eta24_pt03_3_hiPixel           = (TH2D*) fMult->Get( "multcentana/hMultCent_eta24_pt03_3_hiPixel" );
  hMultCent_eta10_pt03_3_hiPixel           = (TH2D*) fMult->Get( "multcentana/hMultCent_eta10_pt03_3_hiPixel" );
  hMultCent_eta24_pt03_3_hiGeneral         = (TH2D*) fMult->Get( "multcentana/hMultCent_eta24_pt03_3_hiGeneral" );
  hMultCent_eta10_pt03_3_hiGeneral         = (TH2D*) fMult->Get( "multcentana/hMultCent_eta24_pt03_3_hiGeneral" );

  TCanvas * cMultCent_pt03_3_hiPixelAndGeneral = new TCanvas("cMultCent_pt03_3_hiPixelAndGeneral","cMultCent_pt03_3_hiPixelAndGeneral",1000,500);
  cMultCent_pt03_3_hiPixelAndGeneral->Divide(2,1);
  cMultCent_pt03_3_hiPixelAndGeneral->cd(1)->SetLogy();
  cMultCent_pt03_3_hiPixelAndGeneral->cd(2)->SetLogy();

  TCanvas * cMultCent_pt03_3_hiPixel = new TCanvas("cMultCent_pt03_3_hiPixel","cMultCent_pt03_3_hiPixel",1000,500);
  cMultCent_pt03_3_hiPixel->Divide(2,1);
  cMultCent_pt03_3_hiPixel->cd(1)->SetLogy();
  cMultCent_pt03_3_hiPixel->cd(2)->SetLogy();

  TCanvas * cMultCent_pt03_3_hiGeneral = new TCanvas("cMultCent_pt03_3_hiGeneral","cMultCent_pt03_3_hiGeneral",1000,500);
  cMultCent_pt03_3_hiGeneral->Divide(2,1);
  cMultCent_pt03_3_hiGeneral->cd(1)->SetLogy();
  cMultCent_pt03_3_hiGeneral->cd(2)->SetLogy();

  for(int icent = 0; icent < N; icent++){

    hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent] = (TH1D*) hMultCent_eta24_pt03_3_hiPixelAndGeneral->ProjectionY( Form("hMultCent_eta24_pt03_3_hiPixelAndGeneralc%i", icent), 2*cmin[icent], 2*cmax[icent] );
    hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent] = (TH1D*) hMultCent_eta10_pt03_3_hiPixelAndGeneral->ProjectionY( Form("hMultCent_eta10_pt03_3_hiPixelAndGeneral_c%i", icent), 2*cmin[icent], 2*cmax[icent] );
    hMultCent_eta24_pt03_3_hiPixel_proj[icent]           = (TH1D*) hMultCent_eta24_pt03_3_hiPixel->ProjectionY( Form("hMultCent_eta24_pt03_3_hiPixel_c%i", icent), 2*cmin[icent], 2*cmax[icent] );
    hMultCent_eta10_pt03_3_hiPixel_proj[icent]           = (TH1D*) hMultCent_eta10_pt03_3_hiPixel->ProjectionY( Form("hMultCent_eta10_pt03_3_hiPixel_c%i", icent), 2*cmin[icent], 2*cmax[icent] );
    hMultCent_eta24_pt03_3_hiGeneral_proj[icent]         = (TH1D*) hMultCent_eta24_pt03_3_hiGeneral->ProjectionY( Form("hMultCent_eta24_pt03_3_hiGeneral_c%i", icent), 2*cmin[icent], 2*cmax[icent] );
    hMultCent_eta10_pt03_3_hiGeneral_proj[icent]         = (TH1D*) hMultCent_eta10_pt03_3_hiGeneral->ProjectionY( Form("hMultCent_eta10_pt03_3_hiGeneral_c%i", icent), 2*cmin[icent], 2*cmax[icent] );

    /*
    hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->Scale(1./hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->Integral());
    hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->Scale(1./hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->Integral());
    hMultCent_eta24_pt03_3_hiPixel_proj[icent]->Scale(1./hMultCent_eta24_pt03_3_hiPixel_proj[icent]->Integral());
    hMultCent_eta10_pt03_3_hiPixel_proj[icent]->Scale(1./hMultCent_eta10_pt03_3_hiPixel_proj[icent]->Integral());
    hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->Scale(1./hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->Integral());
    hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->Scale(1./hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->Integral());

    std::cout<< Form("Cent %i-%i", cmin[icent], cmax[icent]) << "\t Integral = " << hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->Integral() << std::endl;
    */
    hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->SetLineColor(ccol[icent]);
    hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->SetLineColor(ccol[icent]);
    hMultCent_eta24_pt03_3_hiPixel_proj[icent]->SetLineColor(ccol[icent]);
    hMultCent_eta10_pt03_3_hiPixel_proj[icent]->SetLineColor(ccol[icent]);
    hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->SetLineColor(ccol[icent]);
    hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->SetLineColor(ccol[icent]);

    hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->SetMarkerColor(1);
    hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->SetMarkerColor(1);
    hMultCent_eta24_pt03_3_hiPixel_proj[icent]->SetMarkerColor(1);
    hMultCent_eta10_pt03_3_hiPixel_proj[icent]->SetMarkerColor(1);
    hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->SetMarkerColor(1);
    hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->SetMarkerColor(1);

    hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->SetMaximum( 20.*hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[0]->GetMaximum() );
    hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->SetMaximum( 20.*hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[0]->GetMaximum() );
    hMultCent_eta24_pt03_3_hiPixel_proj[icent]->SetMaximum( 20.*hMultCent_eta24_pt03_3_hiPixel_proj[0]->GetMaximum() );
    hMultCent_eta10_pt03_3_hiPixel_proj[icent]->SetMaximum( 20.*hMultCent_eta10_pt03_3_hiPixel_proj[0]->GetMaximum() );
    hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->SetMaximum( 20.*hMultCent_eta24_pt03_3_hiGeneral_proj[0]->GetMaximum() );
    hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->SetMaximum( 20.*hMultCent_eta10_pt03_3_hiGeneral_proj[0]->GetMaximum() );

    if(icent == 0){
      cMultCent_pt03_3_hiPixelAndGeneral->cd(1);
      hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->Draw();
      cMultCent_pt03_3_hiPixelAndGeneral->cd(2);
      hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->Draw();

      cMultCent_pt03_3_hiPixel->cd(1);
      hMultCent_eta24_pt03_3_hiPixel_proj[icent]->Draw();
      cMultCent_pt03_3_hiPixel->cd(2);
      hMultCent_eta10_pt03_3_hiPixel_proj[icent]->Draw();

      cMultCent_pt03_3_hiGeneral->cd(1);
      hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->Draw();
      cMultCent_pt03_3_hiGeneral->cd(2);
      hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->Draw();
    }
    else{
      cMultCent_pt03_3_hiPixelAndGeneral->cd(1);
      hMultCent_eta24_pt03_3_hiPixelAndGeneral_proj[icent]->Draw("same");
      cMultCent_pt03_3_hiPixelAndGeneral->cd(2);
      hMultCent_eta10_pt03_3_hiPixelAndGeneral_proj[icent]->Draw("same");

      cMultCent_pt03_3_hiPixel->cd(1);
      hMultCent_eta24_pt03_3_hiPixel_proj[icent]->Draw("same");
      cMultCent_pt03_3_hiPixel->cd(2);
      hMultCent_eta10_pt03_3_hiPixel_proj[icent]->Draw("same");

      cMultCent_pt03_3_hiGeneral->cd(1);
      hMultCent_eta24_pt03_3_hiGeneral_proj[icent]->Draw("same");
      cMultCent_pt03_3_hiGeneral->cd(2);
      hMultCent_eta10_pt03_3_hiGeneral_proj[icent]->Draw("same");
    }

  } //-- End cent loop

  cMultCent_pt03_3_hiPixelAndGeneral->cd(1);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.57, 0.70, "hiPixelAndGeneral");

  cMultCent_pt03_3_hiPixelAndGeneral->cd(2);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 1.0");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.57, 0.70, "hiPixelAndGeneral");

  cMultCent_pt03_3_hiPixel->cd(1);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.8, 0.70, "hiPixel");

  cMultCent_pt03_3_hiPixel->cd(2);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 1.0");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.8, 0.70, "hiPixel");

  cMultCent_pt03_3_hiGeneral->cd(1);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.75, 0.70, "hiGeneral");

  cMultCent_pt03_3_hiGeneral->cd(2);
  latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
  latex.DrawLatex(0.78, 0.82, "|#eta| < 1.0");
  latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
  latex.DrawLatex(0.75, 0.70, "hiGeneral");


  cMultCent_pt03_3_hiPixelAndGeneral->SaveAs( Form("cMultCent_pt03_3_hiPixelAndGeneral_%s.pdf", clusTrunc.data()) );
  //cMultCent_pt03_3_hiPixel->SaveAs("cMultCent_pt03_3_hiPixel.pdf");
  //cMultCent_pt03_3_hiGeneral->SaveAs("cMultCent_pt03_3_hiGeneral.pdf");

}

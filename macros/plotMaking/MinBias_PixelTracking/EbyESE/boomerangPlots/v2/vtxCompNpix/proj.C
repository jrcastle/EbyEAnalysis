#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TLegend.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void proj(){

  const int ll = 2;
  const int lh = 3;

  const int hl = 101;
  const int hh = 110;

  //TH1D::SetDefaultSumw2();

  TFile * fPix;
  TFile * fMC;

  TH2D * hVtxCompNpix_Pix;
  TH1D * hVtxComp_Pix_Low;
  TH1D * hVtxComp_Pix_High;
  TH1D * hVtxComp_Pix_Low_RebinReshift;
  TH1D * hVtxComp_Pix_High_RebinReshift;

  TH2D * hVtxCompNpix_MC;
  TH1D * hVtxComp_MC_Low;
  TH1D * hVtxComp_MC_High;
  TH1D * hVtxComp_MC_Low_RebinReshift;
  TH1D * hVtxComp_MC_High_RebinReshift;


  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  fPix = new TFile("MultCent_Pixel.root");
  hVtxCompNpix_Pix  = (TH2D*) fPix->Get("multcentana/hClusVtxCompInteg");
  hVtxComp_Pix_Low  = (TH1D*) hVtxCompNpix_Pix->ProjectionY("hVtxComp_Pix_Low", ll, lh);
  hVtxComp_Pix_High = (TH1D*) hVtxCompNpix_Pix->ProjectionY("hVtxComp_Pix_High", hl, hh);

  fMC  = new TFile("MultCent_MC.root");
  hVtxCompNpix_MC  = (TH2D*) fMC->Get("multcentana/hClusVtxCompInteg");
  hVtxComp_MC_Low  = (TH1D*) hVtxCompNpix_MC->ProjectionY("hVtxComp_MC_Low", ll, lh);
  hVtxComp_MC_High = (TH1D*) hVtxCompNpix_MC->ProjectionY("hVtxComp_MC_High", hl, hh);

  int bin_ll = hVtxCompNpix_MC->GetXaxis()->GetBinLowEdge(ll);
  int bin_lh = hVtxCompNpix_MC->GetXaxis()->GetBinLowEdge(lh+1);
  int bin_hl = hVtxCompNpix_MC->GetXaxis()->GetBinLowEdge(hl);
  int bin_hh = hVtxCompNpix_MC->GetXaxis()->GetBinLowEdge(hh+1);

  hVtxComp_Pix_Low->GetYaxis()->SetTitle("Events");
  hVtxComp_Pix_Low->SetLineColor(2);
  hVtxComp_Pix_Low->SetMarkerColor(2);
  hVtxComp_Pix_Low->SetMarkerStyle(20);
  //hVtxComp_Pix_Low->Scale(1./hVtxComp_Pix_Low->Integral());

  hVtxComp_MC_Low->GetYaxis()->SetTitle("Events");
  hVtxComp_MC_Low->SetLineColor(4);
  hVtxComp_MC_Low->SetMarkerColor(4);
  hVtxComp_MC_Low->SetMarkerStyle(24);
  //hVtxComp_MC_Low->Scale(1./hVtxComp_MC_Low->Integral());

  hVtxComp_Pix_High->GetYaxis()->SetTitle("Events");
  hVtxComp_Pix_High->SetLineColor(2);
  hVtxComp_Pix_High->SetMarkerColor(2);
  hVtxComp_Pix_High->SetMarkerStyle(21);
  //hVtxComp_Pix_High->Scale(1./hVtxComp_Pix_High->Integral());

  hVtxComp_MC_High->GetYaxis()->SetTitle("Events");
  hVtxComp_MC_High->SetLineColor(4);
  hVtxComp_MC_High->SetMarkerColor(4);
  hVtxComp_MC_High->SetMarkerStyle(25);
  //hVtxComp_MC_High->Scale(1./hVtxComp_MC_High->Integral());

  hVtxComp_Pix_Low_RebinReshift  = new TH1D("hVtxComp_Pix_Low_RebinReshift",  "hVtxComp_Pix_Low_RebinReshift",  100, 0, 20);
  hVtxComp_Pix_Low_RebinReshift->GetXaxis()->SetTitle("clus-vtx compatibility");
  hVtxComp_Pix_Low_RebinReshift->GetYaxis()->SetTitle("Events");
  hVtxComp_Pix_Low_RebinReshift->SetLineColor(2);
  hVtxComp_Pix_Low_RebinReshift->SetMarkerColor(2);
  hVtxComp_Pix_Low_RebinReshift->SetMarkerStyle(20);

  hVtxComp_Pix_High_RebinReshift = new TH1D("hVtxComp_Pix_High_RebinReshift", "hVtxComp_Pix_High_RebinReshift", 300, 0, 20);
  hVtxComp_Pix_High_RebinReshift->GetXaxis()->SetTitle("clus-vtx compatibility");
  hVtxComp_Pix_High_RebinReshift->GetYaxis()->SetTitle("Events");
  hVtxComp_Pix_High_RebinReshift->SetLineColor(2);
  hVtxComp_Pix_High_RebinReshift->SetMarkerColor(2);
  hVtxComp_Pix_High_RebinReshift->SetMarkerStyle(21);
  int highULBin = hVtxComp_Pix_High_RebinReshift->FindBin(6);
  hVtxComp_Pix_High_RebinReshift->GetXaxis()->SetRange(1,highULBin);

  hVtxComp_MC_Low_RebinReshift   = new TH1D("hVtxComp_MC_Low_RebinReshift",   "hVtxComp_MC_Low_RebinReshift",   100, 0, 20);
  hVtxComp_MC_Low_RebinReshift->GetXaxis()->SetTitle("clus-vtx compatibility");
  hVtxComp_MC_Low_RebinReshift->GetYaxis()->SetTitle("Events");
  hVtxComp_MC_Low_RebinReshift->SetLineColor(4);
  hVtxComp_MC_Low_RebinReshift->SetMarkerColor(4);
  hVtxComp_MC_Low_RebinReshift->SetMarkerStyle(24);

  hVtxComp_MC_High_RebinReshift  = new TH1D("hVtxComp_MC_High_RebinReshift",  "hVtxComp_MC_High_RebinReshift",  300, 0, 20);
  hVtxComp_MC_High_RebinReshift->GetXaxis()->SetTitle("clus-vtx compatibility");
  hVtxComp_MC_High_RebinReshift->GetYaxis()->SetTitle("Events");
  hVtxComp_MC_High_RebinReshift->SetLineColor(4);
  hVtxComp_MC_High_RebinReshift->SetMarkerColor(4);
  hVtxComp_MC_High_RebinReshift->SetMarkerStyle(25);
  hVtxComp_MC_High_RebinReshift->GetXaxis()->SetRange(1,highULBin);

  int N = hVtxComp_MC_High->GetNbinsX();
  double meanpixlow  = hVtxComp_Pix_Low->GetMean();
  double meanpixhigh = hVtxComp_Pix_High->GetMean();
  double meanmclow   = hVtxComp_MC_Low->GetMean();
  double meanmchigh  = hVtxComp_MC_High->GetMean();

  for(int i = 1; i <= N; i++){

    double clusvtxcomp = hVtxComp_Pix_Low->GetBinCenter(i);
    double pixlow  = hVtxComp_Pix_Low->GetBinContent(i);
    double pixhigh = hVtxComp_Pix_High->GetBinContent(i);
    double mclow   = hVtxComp_MC_Low->GetBinContent(i);
    double mchigh  = hVtxComp_MC_High->GetBinContent(i);

    double shiftClusVtxLow  = clusvtxcomp - fabs(meanmclow - meanpixlow);
    double shiftClusVtxHigh = clusvtxcomp - fabs(meanmchigh - meanpixhigh);

    hVtxComp_Pix_Low_RebinReshift->Fill(clusvtxcomp, pixlow);
    hVtxComp_Pix_High_RebinReshift->Fill(clusvtxcomp, pixhigh);
    hVtxComp_MC_Low_RebinReshift->Fill(shiftClusVtxLow, mclow);
    hVtxComp_MC_High_RebinReshift->Fill(shiftClusVtxHigh, mchigh);

  }

  hVtxComp_Pix_Low_RebinReshift->Scale(1./hVtxComp_Pix_Low_RebinReshift->Integral());
  hVtxComp_Pix_High_RebinReshift->Scale(1./hVtxComp_Pix_High_RebinReshift->Integral());
  hVtxComp_MC_Low_RebinReshift->Scale(1./hVtxComp_MC_Low_RebinReshift->Integral());
  hVtxComp_MC_High_RebinReshift->Scale(1./hVtxComp_MC_High_RebinReshift->Integral());

  std::cout << "Pix Low Mean = " << hVtxComp_Pix_Low_RebinReshift->GetMean() << "\t MC Low Mean = " << hVtxComp_MC_Low_RebinReshift->GetMean() << std::endl;
  std::cout << "Pix High Mean = " << hVtxComp_Pix_High_RebinReshift->GetMean() << "\t MC High Mean = " << hVtxComp_MC_High_RebinReshift->GetMean() << std::endl;

  TLegend * leg = new TLegend(0.2, 0.5, 0.5, 0.7);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hVtxComp_MC_Low, "MC", "l");
  leg->AddEntry(hVtxComp_Pix_Low, "DATA", "l");

  TCanvas * c = new TCanvas("c", "c", 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  c->cd(1)->SetLogy();
  //hVtxComp_MC_Low->Draw();
  //hVtxComp_Pix_Low->Draw("same");
  hVtxComp_MC_Low_RebinReshift->Draw();
  hVtxComp_Pix_Low_RebinReshift->Draw("same");
  latex.DrawLatex(0.5, 0.8,  Form("%i < N_{pixel} < %i", bin_ll, bin_lh) );
  c->cd(2);
  c->cd(2)->SetLogy();
  //hVtxComp_MC_High->Draw();
  //hVtxComp_Pix_High->Draw("same");
  hVtxComp_MC_High_RebinReshift->Draw();
  hVtxComp_Pix_High_RebinReshift->Draw("same");
  latex.DrawLatex(0.2, 0.85,  Form("%i < N_{pixel} < %i", bin_hl, bin_hh) );
  leg->Draw("same");
  c->SaveAs("Clus_Vtx_proj_DataVSMC.pdf");

}




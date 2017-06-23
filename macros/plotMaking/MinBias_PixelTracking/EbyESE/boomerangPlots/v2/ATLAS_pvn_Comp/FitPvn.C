#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/ATLAS_PV2.h"

using namespace hi;
using namespace ebyese;
using namespace atlas_pv2;

int BIN = 1;
string fname = "EllPFits.root";

double pEllP(double * x, double * par){

  int nbin = 50;
  //-- [0] = e0
  //-- [1] = alpha
  //-- [2] = kn
  //-- [3] = Scale 
  //-- xx  = vn
  double e0    = par[0];
  double alpha = par[1];
  double kn    = par[2];
  double scale = par[3];
  double eccn  = x[0] / kn;


  double pi = TMath::Pi();
  double p1 = (scale * 2. * alpha * eccn / pi / kn) * pow( 1 - e0*e0, alpha + 0.5);
  double integ = 0.;
  double dphi = pi/(double)nbin;
  for(int i = 1; i <= nbin; i++){
    double phi = (2*(double)i-1.)*pi/(2.*(double)nbin);
    integ += dphi * pow(1-eccn*eccn, alpha-1) * pow( 1-e0*eccn*cos(phi), -2.*alpha-1 );
  }

  double ellp = p1 * integ;

  return ellp;

}

void FitPvn(){

  bool dosys      = 0;
  bool ATLASreg   = 1;
  bool fixKn      = 0;
  bool fixAlpha   = 0;
  // kn 0.40
  // al 71.1
  // e0 .17
  //-- Free kn
  //double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.441, 0.394, 0.392, 0.364, 0.352, 0.344, 0.315, 0.280, 0.260};
  //double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 91.1,  56.3,  45.0,  30.2,  23.3,  18.6,  12.7,  8.1,   6.1};
  //double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.179, 0.227, 0.251, 0.287, 0.307, 0.314, 0.333, 0.0,   0.0};
  //double vnmax[NCENT]   = {0.3,   0.3,   0.3,   0.3,   0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.28,  0.26};
  //-- ATLAS KN
  double knGuess[NCENT] = {16.2,  16.2,  16.2,  0.441, 0.396, 0.394, 0.366, 0.355, 0.346, 0.319, 0.292, 0.275};
  double alGuess[NCENT] = {2.8e5, 2.8e5, 2.8e5, 91.1,  56.3,  45.0,  30.2,  23.3,  18.6,  12.7,  8.1,   6.1};
  double e0Guess[NCENT] = {0.0,   0.0,   0.0,   0.179, 0.227, 0.251, 0.287, 0.307, 0.314, 0.333, 0.333, 0.333};
  double vnmax[NCENT]   = {0.3,   0.3,   0.3,   0.3,   0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.28,  0.26};

  int norder_  = 2;
  double tkEta = 2.4;

  TLatex latex;
  TLatex latex2;
  TLatex latex3;

  //-- Unfold histos
  TFile * fFinalUnf;
  TH1D * hFinalUnfold[NCENT];
  TH1D * hDummy[NCENT];

  TF1 * fBG[NCENT];
  TF1 * fEllP[NCENT];
  TF1 * fEllPATLAS[NCENT];

  TFile * fBottomLine;
  TGraph * grChi2EllP;
  TGraph * grChi2BG;

  //-- Parms vs npart
  double fitKn[NCENT];
  double fitAlpha[NCENT];
  double fitE0[NCENT];

  double fitKn_err[NCENT];
  double fitAlpha_err[NCENT];
  double fitE0_err[NCENT];

  TGraphErrors * grFitKn;
  TGraphErrors * grFitAlpha;
  TGraphErrors * grFitE0;

  //-- ATLAS
  double fitATLASKn[NCENT];
  double fitATLASAlpha[NCENT];
  double fitATLASE0[NCENT];

  double fitATLASKn_err[NCENT];
  double fitATLASAlpha_err[NCENT];
  double fitATLASE0_err[NCENT];

  TGraphErrors * grFitATLASKn;
  TGraphErrors * grFitATLASAlpha;
  TGraphErrors * grFitATLASE0;

  //-- Ollitrault ATLAS
  TGraphErrors * grOlliATLASKn;
  TGraphErrors * grOlliATLASAlpha;
  TGraphErrors * grOlliATLASE0;

  TFile * fOut;

  //
  // MAIN
  //
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(20);
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(23);
  latex3.SetNDC();
  latex3.SetTextFont(43);
  latex3.SetTextSize(26);

  setTDRStyle();
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");

  fFinalUnf = new TFile( "FinalUnfoldDistns.root" );

  fOut = new TFile(fname.data(), "recreate");
  //fBottomLine = new TFile("SmearSpaceChi2.root");
  //grChi2EllP = (TGraph*) fBottomLine->Get("grEllPChi2");
  //grChi2BG   = (TGraph*) fBottomLine->Get("grBGChi2");

  for(int icent = 0; icent < NCENT; icent++){

    //if(icent != BIN) continue;

    fitKn[icent]    = -1;
    fitAlpha[icent] = -1;
    fitE0[icent]    = -1;

    fitATLASKn[icent]    = -1;
    fitATLASAlpha[icent] = -1;
    fitATLASE0[icent]    = -1;

    hDummy[icent] = new TH1D(Form("hDummy_c%i", icent), Form("hDummy_c%i",icent), 152, 0-binw, vnMax[norder_]);
    hDummy[icent]->GetXaxis()->SetTitle("v_{2}");
    hDummy[icent]->GetYaxis()->SetTitle("p(v_{2})");

    if(ATLASreg) hFinalUnfold[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldATLASreg_c%i", icent) );
    else         hFinalUnfold[icent] = (TH1D*) fFinalUnf->Get( Form("hFinalUnfoldStat_c%i", icent) );
    hFinalUnfold[icent]->SetLineColor(2);
    hFinalUnfold[icent]->SetMarkerColor(2);
    hFinalUnfold[icent]->SetMarkerStyle(20);
    hFinalUnfold[icent]->GetXaxis()->SetTitle("v_{2}");
    hFinalUnfold[icent]->GetYaxis()->SetTitle("p(v_{2})");
    hFinalUnfold[icent]->GetXaxis()->SetNdivisions(507);
    hFinalUnfold[icent]->GetXaxis()->SetRange(1, hFinalUnfold[icent]->FindBin(vnmax[icent]));
    hDummy[icent]->GetXaxis()->SetRange(1, hFinalUnfold[icent]->FindBin(0.3));


    //if( icent != 3 && icent != 5 && icent !=7 ) continue;
    if( icent < 3 ) continue;
    std::cout<<Form("-------- Centbin %i --------", icent)<<std::endl;

    fBG[icent] = new TF1(Form("fBG_c%i", icent), "[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))", 0.0, 0.4);
    fBG[icent]->SetParLimits(1, 0., 1.);
    fBG[icent]->SetParLimits(2, 0., 1.);
    fBG[icent]->SetLineColor(4);
    fBG[icent]->SetLineWidth(2);

    //-- BG
    EbyECumu cumu(hFinalUnfold[icent]);
    double vn2 = cumu.GetCumu_vn2();
    double vn4 = cumu.GetCumu_vn4();
    double BGmu    = vn4;
    double BGDelta = sqrt( (vn2*vn2 - vn4*vn4)/2. ) ;
    if(icent != 6){
      BGmu = hFinalUnfold[icent]->GetMean();
      BGDelta = hFinalUnfold[icent]->GetRMS();
    }
    fBG[icent]->SetParameters(hFinalUnfold[icent]->GetMaximum(), BGmu, BGDelta);

    //-- EllP
    fEllP[icent] = new TF1(Form("fEllP_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllP[icent]->SetLineColor(4);
    fEllP[icent]->SetLineWidth(2);

    fEllPATLAS[icent] = new TF1(Form("fEllPATLAS_c%i", icent), pEllP, 0.0, vnmax[icent], 4);
    fEllPATLAS[icent]->SetLineColor(4);
    fEllPATLAS[icent]->SetLineWidth(2);

    //-- [0] = ecc0; [1] = Alpha; [2] = kn; [3] = Scale;
    fEllP[icent]->SetParLimits(0, 0., 1.);
    fEllP[icent]->SetParLimits(2, 0., 1.);
    if(fixKn) fEllP[icent]->FixParameter(2, knGuess[icent]);
    fEllP[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], hFinalUnfold[icent]->GetMaximum());
    fEllP[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");

    fEllPATLAS[icent]->SetParLimits(0, 0., 1.);
    fEllPATLAS[icent]->SetParLimits(2, 0., 1.);
    fEllPATLAS[icent]->SetParameters(e0Guess[icent], alGuess[icent], knGuess[icent], TMath::MaxElement(ATLASPV2_Stat[icent]->GetN(), ATLASPV2_Stat[icent]->GetY()) );
    fEllPATLAS[icent]->SetParNames("ecc0", "alpha", "kn", "Scale");

    //-- Fit!
    std::cout << "CMS" << std::endl;
    hFinalUnfold[icent]->Fit( Form("fEllP_c%i", icent), "0", "", 0.0, vnmax[icent]);

    std::cout << "ATLAS" << std::endl;
    ATLASPV2_Stat[icent]->Fit( Form("fEllPATLAS_c%i", icent), "0", "", 0.0, vnmax[icent]);

    //-- Increase hist maximums post-fit.
    hFinalUnfold[icent]->SetMaximum( 5.*hFinalUnfold[icent]->GetMaximum() );
    ATLASPV2_Stat[icent]->SetMarkerStyle(20);
    ATLASPV2_Stat[icent]->SetMarkerColor(1);
    ATLASPV2_Stat[icent]->SetLineColor(1);
    ATLASPV2_Stat[icent]->GetXaxis()->SetNdivisions(507);

    fOut->cd();
    fEllP[icent]->Write();
    fEllPATLAS[icent]->Write();

    fitKn[icent]    = fEllP[icent]->GetParameter(2);
    fitAlpha[icent] = fEllP[icent]->GetParameter(1);
    fitE0[icent]    = fEllP[icent]->GetParameter(0);

    fitKn_err[icent]    = fEllP[icent]->GetParError(2);
    fitAlpha_err[icent] = fEllP[icent]->GetParError(1);
    fitE0_err[icent]    = fEllP[icent]->GetParError(0);

    fitATLASKn[icent]    = fEllPATLAS[icent]->GetParameter(2);
    fitATLASAlpha[icent] = fEllPATLAS[icent]->GetParameter(1);
    fitATLASE0[icent]    = fEllPATLAS[icent]->GetParameter(0);

    fitATLASKn_err[icent]    = fEllPATLAS[icent]->GetParError(2);
    fitATLASAlpha_err[icent] = fEllPATLAS[icent]->GetParError(1);
    fitATLASE0_err[icent]    = fEllPATLAS[icent]->GetParError(0);

  } // End cent loop

  //-- Make Graphs
  double c_err[NCENT];
  for(int icent = 0; icent < NCENT; icent++) c_err[icent] = 0;


  //-- CMS Parms
  grFitKn    = new TGraphErrors(NCENT, centBinCenter, fitKn,    c_err, fitKn_err);
  grFitAlpha = new TGraphErrors(NCENT, centBinCenter, fitAlpha, c_err, fitAlpha_err);
  grFitE0    = new TGraphErrors(NCENT, centBinCenter, fitE0,    c_err, fitE0_err);

  formatGraph(grFitKn,    "Centrality %", 0., 0.8,  "k_{n}",        1, 20, "grFitKn");
  formatGraph(grFitAlpha, "Centrality %", 0., 130,  "#alpha",       2, 21, "grFitAlpha");
  formatGraph(grFitE0,    "Centrality %", 0., 0.5,  "#epsilon_{0}", 4, 22, "grFitE0");

  //-- ATLAS Parms
  grFitATLASKn    = new TGraphErrors(NCENT, centBinCenter, fitATLASKn,    c_err, fitATLASKn_err);
  grFitATLASAlpha = new TGraphErrors(NCENT, centBinCenter, fitATLASAlpha, c_err, fitATLASAlpha_err);
  grFitATLASE0    = new TGraphErrors(NCENT, centBinCenter, fitATLASE0,    c_err, fitATLASE0_err);

  formatGraph(grFitATLASKn,    "Centrality %", 0., 0.8,  "k_{n}",         1, 24, "grFitATLASKn");
  formatGraph(grFitATLASAlpha, "Centrality %", 0., 130,  "#alpha",        2, 25, "grFitATLASAlpha");
  formatGraph(grFitATLASE0,    "Centrality %", 0., 0.5,  "#epsilon_{0}",  4, 26, "grFitATLASE0");

  //-- Ollitrault ATLAS
  grOlliATLASKn    = new TGraphErrors("k_n_ATLAS_green.txt",   "%lg %lg %lg");
  grOlliATLASAlpha = new TGraphErrors("alpha_ATLAS_green.txt", "%lg %lg %lg");
  grOlliATLASE0    = new TGraphErrors("ecc0_ATLAS_green.txt",  "%lg %lg %lg");

  formatGraph(grOlliATLASKn,    "Centrality %", 0.1, 1.0,  "k_{n}",        1, 24, "grOlliATLASKn");
  formatGraph(grOlliATLASAlpha, "Centrality %", 0.1, 120,  "#alpha",       2, 25, "grOlliATLASAlpha");
  formatGraph(grOlliATLASE0,    "Centrality %", 0.1, 0.5,  "#epsilon_{0}", 4, 26, "grOlliATLASE0");

  fOut->cd();
  grFitKn->Write();
  grFitAlpha->Write();
  grFitE0->Write();


  TLegend * legKn = new TLegend(0.2, 0.2, 0.65, 0.35);
  legKn->SetBorderSize(0);
  legKn->SetFillStyle(0);
  legKn->AddEntry(grFitKn,       "CMS ",                 "lp");
  legKn->AddEntry(grFitATLASKn,  "ATLAS 2.76 TeV",       "lp");
  legKn->AddEntry(grOlliATLASKn, "Ollitrault ATLAS Fit", "l");

  TLegend * legAlpha = new TLegend(0.18, 0.2, 0.63, 0.35);
  //TLegend * legAlpha = new TLegend(0.5, 0.75, 0.95, 0.9);
  legAlpha->SetBorderSize(0);
  legAlpha->SetFillStyle(0);
  legAlpha->AddEntry(grFitAlpha,       "CMS ",                 "lp");
  legAlpha->AddEntry(grFitATLASAlpha,  "ATLAS 2.76 TeV",       "lp");
  legAlpha->AddEntry(grOlliATLASAlpha, "Ollitrault ATLAS Fit", "l");

  TLegend * legEcc0 = new TLegend(0.2, 0.2, 0.65, 0.35);
  legEcc0->SetBorderSize(0);
  legEcc0->SetFillStyle(0);
  legEcc0->AddEntry(grFitE0,       "CMS ",                 "lp");
  legEcc0->AddEntry(grFitATLASE0,  "ATLAS 2.76 TeV",       "lp");
  legEcc0->AddEntry(grOlliATLASE0, "Ollitrault ATLAS Fit", "l");

  TCanvas * cParmSummary = new TCanvas("cParmSummary", "cParmSummary", 1500, 500);
  cParmSummary->Divide(3,1);

  cParmSummary->cd(1);
  grFitKn->Draw("ap");
  grOlliATLASKn->Draw("lsame");
  grFitKn->Draw("psame");
  grFitATLASKn->Draw("psame");
  legKn->Draw("same");
  latex.DrawLatex(0.2, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.2, 0.82, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex.DrawLatex(0.2, 0.76, Form("|#eta| < %.1f", tkEta));
  latex.DrawLatex(0.2, 0.70, Form("%.1f < p_{T} < %.1f GeV/c", 0.5, 20.0) );

  cParmSummary->cd(2);
  grFitE0->Draw("ap");
  grOlliATLASE0->Draw("lsame");
  grFitE0->Draw("psame");
  grFitATLASE0->Draw("psame");
  legEcc0->Draw("same");
    
  cParmSummary->cd(3);
  grFitAlpha->Draw("ap");
  grOlliATLASAlpha->Draw("lsame");
  grFitAlpha->Draw("psame");
  grFitATLASAlpha->Draw("psame");
  legAlpha->Draw("same");
  
  legKn->SetTextFont(43);
  legKn->SetTextSize(23);
  legAlpha->SetTextFont(43);
  legAlpha->SetTextSize(23);
  legEcc0->SetTextFont(43);
  legEcc0->SetTextSize(23);
  cParmSummary->Update();
  cParmSummary->SaveAs(Form("EllPParmSummary_v%i.pdf",norder_));

  //-- Fit Summaries

  //-- CMS
  TCanvas * cFitSummaryCMS = new TCanvas("cFitSummaryCMS", "cFitSummaryCMS", 1500, 1500);
  cFitSummaryCMS->Divide(3,3);
  int icanv = 1;
  for(int icent = 3; icent < NCENT; icent++){
    cFitSummaryCMS->cd(icanv);
    cFitSummaryCMS->cd(icanv)->SetLogy();
    hFinalUnfold[icent]->Draw();
    fEllP[icent]->Draw("same");  
    if(icent == 3) latex3.DrawLatex(0.2, 0.86, "#bf{CMS} #sqrt{s_{NN}} = 5.02 TeV" );
    latex3.DrawLatex(0.2, 0.22, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    icanv++;
  }
  cFitSummaryCMS->SaveAs("cFitSummaryCMS.pdf");

  //-- ATLAS
  TCanvas * cFitSummaryATLAS = new TCanvas("cFitSummaryATLAS", "cFitSummaryATLAS", 1500, 1500);
  cFitSummaryATLAS->Divide(3,3);
  icanv = 1;
  for(int icent = 3; icent < NCENT; icent++){
    cFitSummaryATLAS->cd(icanv);
    cFitSummaryATLAS->cd(icanv)->SetLogy();
    ATLASPV2_Stat[icent]->Draw("ap");
    fEllPATLAS[icent]->Draw("same");
    ATLASPV2_Stat[icent]->GetYaxis()->SetRangeUser(ATLASPV2_Stat[icent]->GetYaxis()->GetXmin(), 5.*ATLASPV2_Stat[icent]->GetYaxis()->GetXmax() );
    ATLASPV2_Stat[icent]->GetXaxis()->SetLimits(0., vnmax[icent]);
    if(icent ==3) latex3.DrawLatex(0.2, 0.86, "#bf{ATLAS} #sqrt{s_{NN}} = 2.76 TeV" );
    latex3.DrawLatex(0.2, 0.22, Form("#bf{Cent. %i - %i%s}", cent_min[icent], cent_max[icent], "%") );
    icanv++;
  }
  cFitSummaryATLAS->Update();
  cFitSummaryATLAS->SaveAs("cFitSummaryATLAS.pdf");

}

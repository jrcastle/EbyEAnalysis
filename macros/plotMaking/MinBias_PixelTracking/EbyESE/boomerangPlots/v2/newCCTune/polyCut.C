#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void polyCut(){

  bool makeProj = 0;

  const int NCLUS = 8;
  const double clusTruncVal[NCLUS] = {2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1};
  const double clusterPar_         = 0.006;
  const int col[]                  = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};

  const int N = 29;
  const double NPixMin[N] = {150, 300, 450, 600, 750, 900, 1050, 1650, 1800, 1950, 2100, 2250, 12025, 21800, 31575, 41350, 51125, 53080, 55035, 56990, 58945, 60900, 62855, 64810, 66765, 68720, 70675, 80450, 90225};
  const double NPixMax[N] = {300, 450, 600, 750, 900, 1050, 1650, 1800, 1950, 2100, 2250, 12025, 21800, 31575, 41350, 51125, 53080, 55035, 56990, 58945, 60900, 62855, 64810, 66765, 68720, 70675, 80450, 90225, 100000};
  const double CCMax      = 8.;


  TFile * f;
  TH2D * hCCvsNpix;
  TH1D * hCCproj[N];

  double CCcut_0p5pct[N];
  double CCcut_1pct[N];
  double CCcut_2pct[N];
  double CCcut_3pct[N];
  double NPix[N];

  TGraph * grCCcut_0p5pct;
  TGraph * grCCcut_1pct;
  TGraph * grCCcut_2pct;
  TGraph * grCCcut_3pct;

  TF1 * polyCut_0p5pct;
  TF1 * polyCut_1pct;
  TF1 * polyCut_2pct;
  TF1 * polyCut_3pct;

  TF1 * lowPixLine_0p5pct;
  TF1 * lowPixLine_1pct;
  TF1 * lowPixLine_2pct;
  TF1 * lowPixLine_3pct;

  TF1 * highPixLine_0p5pct;
  TF1 * highPixLine_1pct;
  TF1 * highPixLine_2pct;
  TF1 * highPixLine_3pct;

  TLine * l_0p5pct[N];
  TLine * l_1pct[N];
  TLine * l_2pct[N];
  TLine * l_3pct[N];

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  f = new TFile("MultCent_NoCCFilter.root");
  hCCvsNpix = (TH2D*) f->Get("multcentana/hClusVtxCompInteg");
  hCCvsNpix->GetYaxis()->SetTitle("Clus-Vtx Compatibility");
  hCCvsNpix->GetYaxis()->SetTitleSize(0.06);
  hCCvsNpix->GetYaxis()->SetLabelSize(0.05);
  hCCvsNpix->GetXaxis()->SetTitle("N_{Pixel}");
  hCCvsNpix->GetXaxis()->SetTitleSize(0.06);
  hCCvsNpix->GetXaxis()->SetLabelSize(0.05);
  hCCvsNpix->GetXaxis()->SetNdivisions(509);
  int bccmax = hCCvsNpix->GetYaxis()->FindBin(CCMax);
  hCCvsNpix->GetYaxis()->SetRange(1, bccmax);

  for(int i = 0; i < N; i++){

    NPix[i] = (NPixMin[i] + NPixMax[i])/2.;
    int bmin = hCCvsNpix->GetXaxis()->FindBin( NPixMin[i] );
    int bmax = hCCvsNpix->GetXaxis()->FindBin( NPixMax[i] );
    hCCproj[i] = (TH1D*) hCCvsNpix->ProjectionY( Form("hCCproj_NPix_%i_%i", (int)NPixMin[i], (int)NPixMax[i]), bmin, bmax );
    //int bccmax = hCCproj[i]->GetXaxis()->FindBin(CCMax);
    //hCCproj[i]->GetXaxis()->SetRange(1, bccmax);
    int nbins = hCCproj[i]->GetNbinsX();

    double fullIntegral = hCCproj[i]->Integral();

    //-- 0.5%
    for(int j = 1; j <= nbins; j++){
      double cc = hCCproj[i]->GetBinCenter(j);
      double subIntegral = hCCproj[i]->Integral(1, j);
      if( subIntegral / fullIntegral > 0.005 ){
	CCcut_0p5pct[i] = cc;
	break;
      }
    }

    //-- 1.0%
    for(int j = 1; j <= nbins; j++){
      double cc = hCCproj[i]->GetBinCenter(j);
      double subIntegral = hCCproj[i]->Integral(1, j);
      if( subIntegral / fullIntegral > 0.01 ){
        CCcut_1pct[i] = cc;
        break;
      }
    }

    //-- 2.0%
    for(int j = 1; j <= nbins; j++){
      double cc = hCCproj[i]->GetBinCenter(j);
      double subIntegral = hCCproj[i]->Integral(1, j);
      if( subIntegral / fullIntegral > 0.02 ){
        CCcut_2pct[i] = cc;
        break;
      }
    }

    //-- 3.0%
    for(int j = 1; j <= nbins; j++){
      double cc = hCCproj[i]->GetBinCenter(j);
      double subIntegral = hCCproj[i]->Integral(1, j);
      if( subIntegral / fullIntegral > 0.03 ){
        CCcut_3pct[i] = cc;
        break;
      }
    }

  }

  grCCcut_0p5pct = new TGraph(N, NPix, CCcut_0p5pct);
  grCCcut_0p5pct->GetXaxis()->SetTitle("N_{Pixel}");
  grCCcut_0p5pct->GetYaxis()->SetTitle("Clus-Vtx Compatibility");

  grCCcut_1pct = new TGraph(N, NPix, CCcut_1pct);
  grCCcut_1pct->GetXaxis()->SetTitle("N_{Pixel}");
  grCCcut_1pct->GetYaxis()->SetTitle("Clus-Vtx Compatibility");

  grCCcut_2pct = new TGraph(N, NPix, CCcut_2pct);
  grCCcut_2pct->GetXaxis()->SetTitle("N_{Pixel}");
  grCCcut_2pct->GetYaxis()->SetTitle("Clus-Vtx Compatibility");

  grCCcut_3pct = new TGraph(N, NPix, CCcut_3pct);
  grCCcut_3pct->GetXaxis()->SetTitle("N_{Pixel}");
  grCCcut_3pct->GetYaxis()->SetTitle("Clus-Vtx Compatibility");

  /*
  polyCut_0p5pct = new TF1("polyCut_0p5pct", "[0] + [1]*TMath::Log(x)", 1000, 60000);  
  polyCut_0p5pct->SetLineColor(2);
  polyCut_0p5pct->SetLineWidth(2);
  grCCcut_0p5pct->Fit("polyCut_0p5pct", "0", "", 1000, 60000);

  polyCut_1pct = new TF1("polyCut_1pct", "[0] + [1]*TMath::Log(x)", 1000, 60000);
  polyCut_1pct->SetLineColor(7);
  polyCut_1pct->SetLineWidth(2);
  grCCcut_1pct->Fit("polyCut_1pct", "0", "", 1000, 60000);
  */

  polyCut_0p5pct = new TF1("polyCut_0p5pct", "pol5", 1000, 60900);
  polyCut_0p5pct->SetLineColor(2);
  polyCut_0p5pct->SetLineWidth(2);
  grCCcut_0p5pct->Fit("polyCut_0p5pct", "0", "", 1000, 60900);

  lowPixLine_0p5pct = new TF1("lowPixLine_0p5pct", "pol1", 150, 1000);
  double lowPixSlope     = (polyCut_0p5pct->Eval(1000.) - CCcut_0p5pct[0]) / (1000. - NPix[0] );
  double lowPixIntercept = CCcut_0p5pct[0] - lowPixSlope * NPix[0];
  lowPixLine_0p5pct->SetParameters(lowPixIntercept, lowPixSlope);
  lowPixLine_0p5pct->SetLineColor(2);
  lowPixLine_0p5pct->SetLineWidth(2);

  highPixLine_0p5pct = new TF1("highPixLine_0p5pct", "pol1", 60900, 100000);
  highPixLine_0p5pct->SetParameters(polyCut_0p5pct->Eval(60900), 0);
  highPixLine_0p5pct->SetLineColor(2);
  highPixLine_0p5pct->SetLineWidth(2);

  std::cout << "0.5% Line p0 = " << lowPixIntercept << std::endl;
  std::cout << "0.5% Line p1 = " << lowPixSlope << std::endl;
  std::cout << "0.5% ClusterTrunc = " << polyCut_0p5pct->Eval(60900.) << std::endl;
  std::cout << "0.5% nHitsLineTruc = " << 1000 <<std::endl;

  polyCut_1pct = new TF1("polyCut_1pct", "pol5", 2100, 60900);
  polyCut_1pct->SetLineColor(7);
  polyCut_1pct->SetLineWidth(2);
  grCCcut_1pct->Fit("polyCut_1pct", "0", "", 2100, 60900);

  lowPixLine_1pct = new TF1("lowPixLine_1pct", "pol1", 150, 2100);
  lowPixSlope     = (polyCut_1pct->Eval(2100.) - CCcut_1pct[0]) / (2100. - NPix[0] );
  lowPixIntercept = CCcut_1pct[0] - lowPixSlope * NPix[0];
  lowPixLine_1pct->SetParameters(lowPixIntercept, lowPixSlope);
  lowPixLine_1pct->SetLineColor(7);
  lowPixLine_1pct->SetLineWidth(2);

  highPixLine_1pct = new TF1("highPixLine_1pct", "pol1", 60900, 100000);
  highPixLine_1pct->SetParameters(polyCut_1pct->Eval(60900), 0);
  highPixLine_1pct->SetLineColor(7);
  highPixLine_1pct->SetLineWidth(2);

  std::cout << "1% Line p0 = " << lowPixIntercept << std::endl;
  std::cout << "1% Line p1 = " << lowPixSlope << std::endl;
  std::cout << "1% ClusterTrunc = " << polyCut_1pct->Eval(60900.) << std::endl;
  std::cout << "1% nHitsLineTruc = " << 2100 << std::endl;
 
  polyCut_2pct = new TF1("polyCut_2pct", "pol5", 1000, 60000);
  polyCut_2pct->SetLineColor(8);
  polyCut_2pct->SetLineWidth(2);
  grCCcut_2pct->Fit("polyCut_2pct", "0", "", 1000, 60000);

  lowPixLine_2pct = new TF1("lowPixLine_2pct", "pol1", 150, 1000);
  lowPixSlope     = (polyCut_2pct->Eval(1000.) - CCcut_2pct[0]) / (1000. - NPix[0] );
  lowPixIntercept = CCcut_2pct[0] - lowPixSlope * NPix[0];
  lowPixLine_2pct->SetParameters(lowPixIntercept, lowPixSlope);
  lowPixLine_2pct->SetLineColor(8);
  lowPixLine_2pct->SetLineWidth(2);

  highPixLine_2pct = new TF1("highPixLine_2pct", "pol1", 60000, 100000);
  highPixLine_2pct->SetParameters(polyCut_2pct->Eval(60000), 0);
  highPixLine_2pct->SetLineColor(8);
  highPixLine_2pct->SetLineWidth(2);

  std::cout << "2% Line p0 = " << lowPixIntercept << std::endl;
  std::cout << "2% Line p1 = " << lowPixSlope << std::endl;
  std::cout << "2% ClusterTrunc = " << polyCut_2pct->Eval(60000.) << std::endl;
  std::cout << "2% nHitsLineTruc = " << 1000 <<std::endl;

  polyCut_3pct = new TF1("polyCut_3pct", "pol5", 1000, 58945); //60000
  polyCut_3pct->SetLineColor(6);
  polyCut_3pct->SetLineWidth(2);
  grCCcut_3pct->Fit("polyCut_3pct", "0", "", 1000, 58945);

  lowPixLine_3pct = new TF1("lowPixLine_3pct", "pol1", 150, 1000);
  lowPixSlope     = (polyCut_3pct->Eval(1000.) - CCcut_3pct[0]) / (1000. - NPix[0] );
  lowPixIntercept = CCcut_3pct[0] - lowPixSlope * NPix[0];
  lowPixLine_3pct->SetParameters(lowPixIntercept, lowPixSlope);
  lowPixLine_3pct->SetLineColor(6);
  lowPixLine_3pct->SetLineWidth(2);

  highPixLine_3pct = new TF1("highPixLine_3pct", "pol1", 58945, 100000);
  highPixLine_3pct->SetParameters(polyCut_3pct->Eval(58945), 0);
  highPixLine_3pct->SetLineColor(6);
  highPixLine_3pct->SetLineWidth(2);

  std::cout << "3% Line p0 = " << lowPixIntercept << std::endl;
  std::cout << "3% Line p1 = " << lowPixSlope << std::endl;
  std::cout << "3% ClusterTrunc = " << polyCut_3pct->Eval(58945.) << std::endl;
  std::cout << "3% nHitsLineTruc = " << 1000 <<std::endl;

  TCanvas * cFit_0p5pct = new TCanvas("cFit_0p5pct", "cFit_0p5pct", 500, 500);
  cFit_0p5pct->cd();
  grCCcut_0p5pct->Draw("ap");
  polyCut_0p5pct->Draw("same");
  lowPixLine_0p5pct->Draw("same");
  highPixLine_0p5pct->Draw("same");
  cFit_0p5pct->SaveAs("cFit_0p5pct.pdf");

  TCanvas * cFit_1pct = new TCanvas("cFit_1pct", "cFit_1pct", 500, 500);
  cFit_1pct->cd();
  grCCcut_1pct->Draw("ap");
  polyCut_1pct->Draw("same");
  lowPixLine_1pct->Draw("same");
  highPixLine_1pct->Draw("same");
  cFit_1pct->SaveAs("cFit_1pct.pdf");

  TCanvas * cFit_2pct = new TCanvas("cFit_2pct", "cFit_2pct", 500, 500);
  cFit_2pct->cd()->SetRightMargin(0.1);
  grCCcut_2pct->Draw("ap");
  polyCut_2pct->Draw("same");
  lowPixLine_2pct->Draw("same");
  highPixLine_2pct->Draw("same");
  cFit_2pct->SaveAs("cFit_2pct.pdf");

  TCanvas * cFit_3pct = new TCanvas("cFit_3pct", "cFit_3pct", 500, 500);
  cFit_3pct->cd();
  grCCcut_3pct->Draw("ap");
  polyCut_3pct->Draw("same");
  lowPixLine_3pct->Draw("same");
  highPixLine_3pct->Draw("same");
  cFit_3pct->SaveAs("cFit_3pct.pdf");

  TLegend * leg = new TLegend(0.8233, 0.7110, 0.9719, 0.8882);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(polyCut_0p5pct, "0.5%", "l");
  leg->AddEntry(polyCut_1pct, "1.0%","l");
  leg->AddEntry(polyCut_2pct, "2.0%","l");
  leg->AddEntry(polyCut_3pct, "3.0%","l");

  TCanvas * c2DHist = new TCanvas("c2DHist", "c2DHist", 500, 500);
  c2DHist->cd()->SetRightMargin(0.1);
  hCCvsNpix->Draw();
  c2DHist->Update();
  TPaletteAxis* palette = (TPaletteAxis*) hCCvsNpix->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.75);
  palette->SetY1NDC(0.6);
  palette->SetX2NDC(0.8);
  palette->SetY2NDC(0.9);
  /*
  polyCut_0p5pct->Draw("same");
  polyCut_1pct->Draw("same");
  polyCut_2pct->Draw("same");
  polyCut_3pct->Draw("same");
  lowPixLine_0p5pct->Draw("same");
  lowPixLine_1pct->Draw("same");
  lowPixLine_2pct->Draw("same");
  lowPixLine_3pct->Draw("same");
  highPixLine_0p5pct->Draw("same");
  highPixLine_1pct->Draw("same");
  highPixLine_2pct->Draw("same");
  highPixLine_3pct->Draw("same");
  leg->Draw("same");
  */
  c2DHist->SaveAs("c2DHist.pdf");


  TCanvas * cProj[N];
  if( makeProj ){
    for(int i = 0; i < N; i++){

      cProj[i] = new TCanvas(Form("cProj_%i", i), Form("cProj_%i", i), 500, 500);
      cProj[i]->cd();
      cProj[i]->SetLogy();
      hCCproj[i]->Draw(); 

      l_0p5pct[i] = new TLine(CCcut_0p5pct[i], 0, CCcut_0p5pct[i], hCCproj[i]->GetMaximum() );
      l_0p5pct[i]->SetLineColor(2);
      l_0p5pct[i]->SetLineStyle(2);
      l_0p5pct[i]->SetLineWidth(2);
      l_0p5pct[i]->Draw("same");

      l_1pct[i] = new TLine(CCcut_1pct[i], 0, CCcut_1pct[i], hCCproj[i]->GetMaximum() );
      l_1pct[i]->SetLineColor(7);
      l_1pct[i]->SetLineStyle(2);
      l_1pct[i]->SetLineWidth(2);
      l_1pct[i]->Draw("same");

      l_2pct[i] = new TLine(CCcut_2pct[i], 0, CCcut_2pct[i], hCCproj[i]->GetMaximum() );
      l_2pct[i]->SetLineColor(8);
      l_2pct[i]->SetLineStyle(2);
      l_2pct[i]->SetLineWidth(2);
      l_2pct[i]->Draw("same");

      l_3pct[i] = new TLine(CCcut_3pct[i], 0, CCcut_3pct[i], hCCproj[i]->GetMaximum() );
      l_3pct[i]->SetLineColor(6);
      l_3pct[i]->SetLineStyle(2);
      l_3pct[i]->SetLineWidth(2);
      l_3pct[i]->Draw("same");

      latex.DrawLatex(0.55, 0.9, Form("%i < N_{Pixel} < %i", (int)NPixMin[i], (int)NPixMax[i]) );
      leg->Draw("same");
      cProj[i]->SaveAs( Form("cProj%i.pdf", i) );
    
    }
  }
}

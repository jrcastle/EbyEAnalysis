#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void clusTree(){

  const int NCENT = 11; 
  const int cmin[NCENT]  = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int cmax[NCENT]  = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const double cbins[]   = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const int ccol[NCENT]  = {kRed+2, kGreen+1, kBlue+2, kMagenta, kOrange-3, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};
  const int Nbins[NCENT] = {300, 250, 250, 250, 175, 125, 100, 100, 100, 75, 50};
  const int mmax[NCENT]  = {6000, 5000, 4500, 3500, 2750, 2250, 1000, 500, 250, 150, 100};
  TH1D hCent("hCent", "hCent", NCENT, cbins);

  const int NCLUS = 8;
  string clusTrunc[NCLUS]    = {"2p0", "2p3", "2p6", "2p9", "3p2", "3p5", "3p8", "4p1"};
  double clusTruncVal[NCLUS] = {2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1};

  TFile * f;
  TTree * tree;
  double centval;
  unsigned int Mult_pt03_3_eta24;
  unsigned int Mult_pt03_3_eta10;
  double clusVtxQual;

  int NEvt[NCLUS][NCENT];
  int NEvtTot[NCLUS];

  TH1D * hMult24[NCLUS][NCENT];
  TH1D * hMult10[NCLUS][NCENT];

  TCanvas * cMultCent[NCLUS];
  TCanvas * cMultClus[NCENT];

  TLatex latex;

  //
  // MAIN
  //
  TH1D::SetDefaultSumw2();
  setTDRStyle();
  latex.SetNDC();
  latex.SetTextFont(43);
  latex.SetTextSize(23);

  f = new TFile("/rfs/jcastle/PbPb2015/PixelTracking_MB2/correction/clusTree.root");
  tree = (TTree*) f->Get("multcentana/tree");
  tree->SetBranchAddress("Cent",              &centval);
  tree->SetBranchAddress("Mult_pt03_3_eta24", &Mult_pt03_3_eta24);
  tree->SetBranchAddress("Mult_pt03_3_eta10", &Mult_pt03_3_eta10);
  tree->SetBranchAddress("clusVtxQual",       &clusVtxQual);

  //-- Initialize Event Counters
  for(int i = 0; i < NCLUS; i++){
    NEvtTot[i] = 0;
    for(int c = 0; c < NCENT; c++){
      NEvt[i][c] = 0;

      hMult24[i][c] = new TH1D(Form("hMult24_c%i_%i_%s", cmin[c], cmax[c], clusTrunc[i].data()), Form("hMult24_c%i_%i_%s", cmin[c], cmax[c], clusTrunc[i].data()), Nbins[c], 0, mmax[c]);
      hMult24[i][c]->GetXaxis()->SetNdivisions(509);
      hMult24[i][c]->GetXaxis()->SetTitle("Event Multiplicity");
      hMult24[i][c]->GetYaxis()->SetTitle("N_{Events}");
      hMult24[i][c]->SetMarkerColor(1);
      hMult24[i][c]->SetMarkerStyle(6);

      hMult10[i][c] = new TH1D(Form("hMult10_c%i_%i_%s", cmin[c], cmax[c], clusTrunc[i].data()), Form("hMult10_c%i_%i_%s", cmin[c], cmax[c], clusTrunc[i].data()), Nbins[c], 0, mmax[c]);
      hMult10[i][c]->GetXaxis()->SetNdivisions(509);
      hMult10[i][c]->GetXaxis()->SetTitle("Event Multiplicity");
      hMult10[i][c]->GetYaxis()->SetTitle("N_{Events}");
      hMult10[i][c]->SetMarkerColor(1);
      hMult10[i][c]->SetMarkerStyle(6);
    }
  }


  int N = tree->GetEntries();
  for(int ievent = 0; ievent < N; ievent++) {
        
    if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<100.*(ievent+1)/(double)N<<"% Completed"<<endl;
    tree->GetEntry(ievent);
    int icent = hCent.FindBin(centval)-1;

    for(int i = 0; i < NCLUS; i++){
      if(clusVtxQual > clusTruncVal[i]){
	NEvtTot[i]++;
	NEvt[i][icent]++;
	hMult24[i][icent]->Fill(Mult_pt03_3_eta24);
	hMult10[i][icent]->Fill(Mult_pt03_3_eta10);
      } //-- End if(clusVtxQual > clusTruncVal[i])
    } //-- End clus loop

  } //-- End event loop


  //-- Make plots

  //-- Group by Centrality
  TLegend * legCent = new TLegend(0.75, 0.2, 0.95, 0.5);
  legCent->SetBorderSize(0);
  legCent->SetFillStyle(0);

  for(int i = 0; i < NCLUS; i++){

    cMultCent[i] = new TCanvas(Form("cMultCent_%s", clusTrunc[i].data()), Form("cMultCent_%s", clusTrunc[i].data()), 500, 500);
    cMultCent[i]->SetLogy();

    for(int c = 0; c < NCENT; c++){

      hMult24[i][c]->SetLineColor(ccol[c]);
      hMult24[i][c]->SetMaximum(20.*hMult24[i][0]->GetMaximum());
      if(i==0) legCent->AddEntry(hMult24[i][c], Form("Cent %i-%i%s", cmin[c], cmax[c], "%"), "lp");

      cMultCent[i]->cd();
      if(c==0) hMult24[i][c]->Draw();
      else     hMult24[i][c]->Draw("same");

    }

    cMultCent[i]->cd();
    legCent->Draw("same");
    latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
    latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
    latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
    latex.DrawLatex(0.57, 0.70, Form("ClusterTrunc = %.1f", clusTruncVal[i]) );
    cMultCent[i]->SaveAs( Form("cMultCent_%s.png", clusTrunc[i].data()) );

  }

  //-- Group by clusTrunc
  TLegend * legClus = new TLegend(0.75, 0.2, 0.95, 0.5);
  legClus->SetBorderSize(0);
  legClus->SetFillStyle(0);

  for(int c = 0; c < NCENT; c++){

    cMultClus[c] = new TCanvas(Form("cMultClusc%i_%i", cmin[c], cmax[c]), Form("cMultClusc%i_%i", cmin[c], cmax[c]), 500, 500);
    cMultClus[c]->SetLogy();

    for(int i = 0; i < NCLUS; i++){

      hMult24[i][c]->SetLineColor(ccol[i]);
      if(c==0) legClus->AddEntry(hMult24[i][c], Form("ClusTrunc %.1f", clusTruncVal[i]), "lp");

      cMultClus[c]->cd();
      if(i==0) hMult24[i][c]->Draw();
      else     hMult24[i][c]->Draw("same");

    }

    cMultClus[c]->cd();
    legClus->Draw("same");
    latex.DrawLatex(0.65, 0.88, "Pixel reRECO");
    latex.DrawLatex(0.78, 0.82, "|#eta| < 2.4");
    latex.DrawLatex(0.53, 0.76, "0.3 < p_{T} < 3.0 GeV/c");
    latex.DrawLatex(0.7, 0.70, Form("Cent %i-%i%s", cmin[c], cmax[c], "%") );
    cMultClus[c]->SaveAs( Form("cMultClus_c%i_%i.png", cmin[c], cmax[c]) );

  }

  //-- Spit out tables

  //-- Total Events:
  std::cout << "Cut" << "\t" << "NEvents" << "\t" << "% Lost" << std::endl;
  for(int i = 0; i < NCLUS; i++){
    double pct = 100.*fabs( (double)NEvtTot[i] - (double)NEvtTot[0] ) / (double)NEvtTot[0];
    std::cout << clusTruncVal[i] << "\t" << NEvtTot[i] << "\t" << Form("%.1f", pct) << std::endl;
  }




}

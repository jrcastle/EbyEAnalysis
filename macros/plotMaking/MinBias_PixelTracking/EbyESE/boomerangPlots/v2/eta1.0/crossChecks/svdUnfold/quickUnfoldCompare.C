#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"

using namespace ebyese;

void quickUnfoldCompare(){


  TFile * fAna;
  TH1D * hObs[NCENT];

  TFile * fDAG;
  TH1D * hUnfoldDAG[NCENT][NITER];
  TH1D * hRefoldDAG[NCENT][NITER];

  int finalIter[NCENT];

  TFile * fSVD;
  TH1D * hKreg;
  TH1D * hUnfoldSVD[NCENT][NKREG];
  TH1D * hUnfoldSVD_RatioToDAG[NCENT][NKREG];
  TH1D * hRefoldSVD[NCENT][NKREG];

  double DAGchi2[NCENT];
  double SVDchi2[NCENT];

  TLatex latex;

  //
  // MAIN
  //
  setTDRStyle();
  latex.SetNDC();

  //-- Get files
  fAna = new TFile("../../AnalyzerResults/CastleEbyE.root");
  fDAG = new TFile("../../UnfoldResults/dataResp/data2.root");
  fSVD = new TFile("../../UnfoldResults/dataResp/data2_svd.root");

  hKreg = (TH1D*) fSVD->Get("hKreg");

  for(int icent = 0; icent < NCENT; icent++){

    //-- Observed
    hObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );

    //-- DAG
    for(int i = 0; i < NITER; i++){

      //--Un(Re)fold distns
      hUnfoldDAG[icent][i] = (TH1D*) fDAG->Get( Form("hreco%i_c%i", iter[i], icent) );
      hUnfoldDAG[icent][i]->GetXaxis()->SetNdivisions(508);

      hRefoldDAG[icent][i] = (TH1D*) fDAG->Get( Form("hrefold%i_c%i", iter[i], icent) );

      //-- Chi2
      double chi2 = hRefoldDAG[icent][i]->Chi2Test(hObs[icent], "CHI2/NDF");

      if(chi2 < 1.2){
	finalIter[icent] = i;
	DAGchi2[icent]   = chi2;
	break;
      }
      if(i == NITER - 1 ){
        finalIter[icent] = i;
	DAGchi2[icent]   = chi2;
        break;
      }

    } //-- End iter loop

    //-- SVD

    for(int ik = 0; ik < 9; ik++){

      hRefoldSVD[icent][ik] = (TH1D*) fSVD->Get( Form("hrefoldkreg%i_c%i", ik, icent) );
      double chi2 = hRefoldSVD[icent][ik]->Chi2Test(hObs[icent], "CHI2/NDF");
      std::cout<<chi2<<std::endl;
      if(ik == 0) SVDchi2[icent] = chi2;

      hUnfoldSVD[icent][ik] = (TH1D*) fSVD->Get( Form("hrecokreg%i_c%i", ik, icent) );
      hUnfoldSVD[icent][ik]->SetLineColor(2);
      hUnfoldSVD[icent][ik]->SetMarkerColor(2);
      hUnfoldSVD[icent][ik]->GetXaxis()->SetNdivisions(508);

      hUnfoldSVD_RatioToDAG[icent][ik] = (TH1D*) hUnfoldSVD[icent][ik]->Clone( Form("hUnfoldSVD_RatioToDAGkreg%i_c%i", ik, icent) );
      hUnfoldSVD_RatioToDAG[icent][ik]->Divide(hUnfoldDAG[icent][finalIter[icent]]);
      hUnfoldSVD_RatioToDAG[icent][ik]->GetYaxis()->SetTitle("Ratio: SVD/DAG");
      hUnfoldSVD_RatioToDAG[icent][ik]->SetMinimum(-10.);
      hUnfoldSVD_RatioToDAG[icent][ik]->SetMaximum(4.);
    }

  }

  for(int icent = 0; icent < NCENT; icent++) std::cout << "Cent " << icent <<"\t" << Form("%.6f", DAGchi2[icent]) << "\t" << Form("%.6f", SVDchi2[icent]) << "\t" << DAGchi2[icent] - SVDchi2[icent] << std::endl;

  //-- Draw

  /*
  //-- Draw distributions ---------------------------
  TCanvas * c = new TCanvas("c", "c", 2000, 1500);
  c->Divide(4,3);

  TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  legInit( leg );
  leg->AddEntry(hUnfoldDAG[0][finalIter[0]], "D'Agostini", "lp");
  leg->AddEntry(hUnfoldSVD[0][0],            "SVD",        "lp");

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)-1;

    c->cd(icent+1);
    //c->cd(icent+1)->SetLogy();
    hUnfoldSVD[icent][kreg]->Draw();
    hUnfoldDAG[icent][i]->Draw("same");
    //hUnfoldSVD[icent]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
  }
  c->SaveAs("QuickUnfoldCompare.pdf");

  //-- Draw log distns
  TLine * l4s_SVD[NCENT];
  TLine * l4s_DAG[NCENT];
  TCanvas * clog0 = new TCanvas("clog0", "clog0", 2000, 1500);
  clog0->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)-1;

    double m = hUnfoldSVD[icent][kreg]->GetMean() + 4.*hUnfoldSVD[icent][kreg]->GetRMS();
    l4s_SVD[icent] = new TLine(m, 0, m, hUnfoldSVD[icent][kreg]->GetMaximum());
    l4s_SVD[icent]->SetLineColor(2);

    m = hUnfoldDAG[icent][i]->GetMean() + 4.*hUnfoldDAG[icent][i]->GetRMS();
    l4s_DAG[icent] = new TLine(m, 0, m, hUnfoldDAG[icent][i]->GetMaximum());

    clog0->cd(icent+1);
    clog0->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    l4s_SVD[icent]->Draw("same");
    l4s_DAG[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog0->SaveAs("QuickUnfoldCompare_Logkreg.pdf");

  TCanvas * clog1 = new TCanvas("clog1", "clog1", 2000, 1500);
  clog1->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1);

    clog1->cd(icent+1);
    clog1->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog1->SaveAs("QuickUnfoldCompare_Logkregp1.pdf");


  TCanvas * clog2 = new TCanvas("clog2", "clog2", 2000, 1500);
  clog2->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+1;

    clog2->cd(icent+1);
    clog2->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog2->SaveAs("QuickUnfoldCompare_Logkregp2.pdf");

  TCanvas * clog3 = new TCanvas("clog3", "clog3", 2000, 1500);
  clog3->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+2;

    clog3->cd(icent+1);
    clog3->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog3->SaveAs("QuickUnfoldCompare_Logkregp3.pdf");

  TCanvas * clog4 = new TCanvas("clog4", "clog4", 2000, 1500);
  clog4->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+3;

    clog4->cd(icent+1);
    clog4->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog4->SaveAs("QuickUnfoldCompare_Logkregp4.pdf");

  TCanvas * clog5 = new TCanvas("clog5", "clog5", 2000, 1500);
  clog5->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+4;

    clog5->cd(icent+1);
    clog5->cd(icent+1)->SetLogy();
    hUnfoldDAG[icent][i]->Draw();
    hUnfoldSVD[icent][kreg]->Draw("same");
    if(icent == 0) leg->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  clog5->SaveAs("QuickUnfoldCompare_Logkregp5.pdf");


  //-- Draw kreg/DAG --------------------------------
  TLine * l0[NCENT];

  TCanvas * c2 = new TCanvas("c2", "c2", 2000, 1500);
  c2->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)-1;
    l0[icent] = new TLine(hUnfoldSVD_RatioToDAG[0][0]->GetBinLowEdge(1),0., hUnfoldSVD_RatioToDAG[0][0]->GetBinLowEdge(NBinsV2[icent])+binw, 0.);

    c2->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c2->SaveAs("RatioSVD_DAG_kreg.pdf");

  TCanvas * c3 = new TCanvas("c3", "c3", 2000, 1500);
  c3->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1);

    c3->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c3->SaveAs("RatioSVD_DAG_kregp1.pdf"); 

  TCanvas * c4 = new TCanvas("c4", "c4", 2000, 1500);
  c4->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+1;

    c4->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c4->SaveAs("RatioSVD_DAG_kregp2.pdf");

  TCanvas * c5 = new TCanvas("c5", "c5", 2000, 1500);
  c5->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+2;

    c5->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c5->SaveAs("RatioSVD_DAG_kregp3.pdf");

  TCanvas * c6 = new TCanvas("c6", "c6", 2000, 1500);
  c6->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+3;

    c6->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c6->SaveAs("RatioSVD_DAG_kregp4.pdf");

  TCanvas * c7 = new TCanvas("c7", "c7", 2000, 1500);
  c7->Divide(4,3);
  for(int icent = 0; icent < NCENT; icent++){
    int i = finalIter[icent];
    int kreg = hKreg->GetBinContent(icent+1)+4;

    c7->cd(icent+1);
    hUnfoldSVD_RatioToDAG[icent][kreg]->Draw();
    l0[icent]->Draw("same");
    latex.DrawLatex(0.2, 0.2, Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"));
    latex.DrawLatex(0.2, 0.25, Form("kreg = %i", kreg));
  }
  c7->SaveAs("RatioSVD_DAG_kregp5.pdf");
  */
}

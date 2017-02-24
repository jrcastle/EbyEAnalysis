#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void truncComp(){


  TFile * f;
  TH1D * hObs[NCENT];
  double truncVal[NCENT];

  TLine * lineTrunc[NCENT];
  TLine * line4Sig[NCENT];

  TCanvas * c;

  //
  // MAIN
  //
  setTDRStyle();

  f = new TFile("CastleEbyE.root");

  c = new TCanvas("c", "c", 2000, 1500);
  c->Divide(4,3);

  for(int icent = 0; icent < NCENT; icent++){

    hObs[icent] = (TH1D*) f->Get( Form("qwebye/hVnFull_c%i", icent) );
    double max = 2.*hObs[icent]->GetMaximum();
    hObs[icent]->SetMaximum(max);
    double sig4 = hObs[icent]->GetMean() + 4.*hObs[icent]->GetRMS();

    //-- Find the trunc val
    for(int i = 1; i < NBins; i++){
      double bc = hObs[icent]->GetBinContent(i);
      if(bc < 10){
	truncVal[icent] = hObs[icent]->GetBinCenter(i-1);
	break;
      }
    }

    std::cout << Form("======== Cent %i ========", icent) << std::endl;
    std::cout << "Trunc = " << truncVal[icent] << std::endl;
    std::cout << "4sig  = " << sig4 << std::endl;

    //-- Make lines
    lineTrunc[icent] = new TLine(truncVal[icent], 0., truncVal[icent], max);
    lineTrunc[icent]->SetLineColor(2);

    line4Sig[icent]  = new TLine(sig4, 0., sig4, max);
    line4Sig[icent]->SetLineColor(4);

    //-- Draw
    c->cd(icent+1);
    c->cd(icent+1)->SetLogy();
    hObs[icent]->Draw();
    lineTrunc[icent]->Draw("same");
    line4Sig[icent]->Draw("same");

  }



}

#include "TGraphErrors.h"
#include "TFile.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void makeTheoryGraphs(){

  string fname[12] = {
    "v2_005.txt",  "v2_510.txt",  "v2_1015.txt", "v2_1520.txt", "v2_2025.txt",
    "v2_2530.txt", "v2_3035.txt", "v2_3540.txt", "v2_4045.txt", "v2_4550.txt",
    "v2_5055.txt", "v2_5560.txt"
  };

  TFile * fOut;
  TGraphErrors * grAMPT[NCENT];
  TGraphErrors * grTRENTO[NCENT];

  //
  // MAIN
  //
  setTDRStyle();

  fOut = new TFile("scale_vn_distribution.root", "recreate");

  for(int icent = 0; icent < NCENT; icent++){

    grAMPT[icent] = new TGraphErrors( Form("ampt/%s", fname[icent].data()), "%lg %lg %lg" );
    grAMPT[icent]->SetName( Form("grAMPT_c%i_%i", cent_min[icent], cent_max[icent]) );
    grAMPT[icent]->GetXaxis()->SetTitle("v_{2}/#LTv_{2}#GT");
    grAMPT[icent]->GetYaxis()->SetTitle("p(v_{2}/#LTv_{2}#GT)");

    grTRENTO[icent] = new TGraphErrors( Form("trento/%s", fname[icent].data()), "%lg %lg %lg" );
    grTRENTO[icent]->SetName( Form("grTRENTO_c%i_%i", cent_min[icent], cent_max[icent]) );
    grTRENTO[icent]->GetXaxis()->SetTitle("v_{2}/#LTv_{2}#GT");
    grTRENTO[icent]->GetYaxis()->SetTitle("p(v_{2}/#LTv_{2}#GT)");

    fOut->cd();
    grAMPT[icent]->Write();
    grTRENTO[icent]->Write();

  }

  //fOut->Write();

}

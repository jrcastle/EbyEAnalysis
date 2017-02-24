#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void makeDDResp(){

  bool vnObsPrior   = 1;
  const int norder_ = 2;

  TFile * fAna;
  TH1D * hVnObs[NCENT];
  TH2D * hSmear2D[NCENT];
  TH1D * hSmearX[NCENT];
  TH1D * hSmearY[NCENT];

  TFile * fOut;
  TH1D * hPrior[NCENT];
  TH2D * hRespDD[NCENT];

  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  //
  // MAIN
  //

  setTDRStyle();

  //-- Get the Analyzer file
  fAna = new TFile("../CastleEbyE.root");

  //-- Create the output file
  fOut = new TFile("dataDrivenResponseAndPriors.root", "recreate");

  for(int icent = 0; icent < NCENT; icent++){

    //-- Get the VnObs histogram if that is the prior you wish to use
    hVnObs[icent] = 0;
    if( vnObsPrior ){
      hVnObs[icent] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_c%i", icent) );
      hPrior[icent] = (TH1D*) hVnObs[icent]->Clone( Form("hPrior_c%i", icent) );
    }

    //-- Grab the 2D smearing histogram and project onto the x and y axes
    hSmear2D[icent] = (TH2D*) fAna->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );
    hSmearX[icent]  = (TH1D*) hSmear2D[icent]->ProjectionX( Form("hSmearX_c%i", icent) );
    hSmearY[icent]  = (TH1D*) hSmear2D[icent]->ProjectionY( Form("hSmearY_c%i", icent) );

	//-- Initialize the DD response histogram
    hRespDD[icent] = new TH2D( Form("hresp_c%i", icent), Form("hresp_c%i", icent), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
    hRespDD[icent]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
    hRespDD[icent]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
    hRespDD[icent]->SetOption("colz");

    //-- Loop over the number of events in each cent, EP and qn bin and fill the 2D response fn
    int Nevt = hSmear2D[icent]->GetEntries();
    std::cout<< "Processing centBin = "<< icent << std::endl;
    for(int ievt = 0 ; ievt < Nevt; ievt++){

      double vnPrior = hPrior[icent]->GetRandom();
      double vnObsx = 0.;
      double vnObsy = 0.;
      double smearx = hSmearX[icent]->GetRandom();
      double smeary = hSmearY[icent]->GetRandom();

      vnObsx = vnPrior + smearx;
      vnObsy = smeary;

      double vnObs = TMath::Sqrt( pow(vnObsx,2) + pow(vnObsy,2) );
      hRespDD[icent]->Fill( vnObs, vnPrior );

    } //-- End event loop

    //-- Write the histograms to the output file
    fOut->cd();
    hPrior[icent]->Write();
    hRespDD[icent]->Write();

  } //-- End cent loop

}

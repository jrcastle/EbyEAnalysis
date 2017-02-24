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
  const int norder_ = 4;

  TFile * fAna[NSPLIT];
  TH1D * hVnObs[NCENT][NSPLIT];
  TH2D * hSmear2D[NCENT][NSPLIT];
  TH1D * hSmearX[NCENT][NSPLIT];
  TH1D * hSmearY[NCENT][NSPLIT];

  TFile * fOut[NSPLIT];
  TH1D * hPrior[NCENT][NSPLIT];
  TH2D * hRespDD[NCENT][NSPLIT];

  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  //
  // MAIN
  //

  setTDRStyle();

  //-- Get the Analyzer file
  for(int iS = 0; iS < NSPLIT; iS++){
    fAna[iS] = new TFile( Form("../CastleEbyE_Split%i.root", iS) );

    //-- Create the output file
    fOut[iS] = new TFile( Form("dataDrivenResponseAndPriors_Split%i.root", iS), "recreate");

    for(int icent = 0; icent < NCENT; icent++){

      //-- Get the VnObs histogram if that is the prior you wish to use
      hVnObs[icent][iS] = 0;
      if( vnObsPrior ){
	hVnObs[icent][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_c%i", icent) );
	hPrior[icent][iS] = (TH1D*) hVnObs[icent][iS]->Clone( Form("hPrior_c%i", icent) );
      }

      //-- Grab the 2D smearing histogram and project onto the x and y axes
      hSmear2D[icent][iS] = (TH2D*) fAna[iS]->Get( Form("qwebye/h2Vn2D0v1_c%i", icent) );
      hSmearX[icent][iS]  = (TH1D*) hSmear2D[icent][iS]->ProjectionX( Form("hSmearX_c%i", icent) );
      hSmearY[icent][iS]  = (TH1D*) hSmear2D[icent][iS]->ProjectionY( Form("hSmearY_c%i", icent) );

      //-- Initialize the DD response histogram
      hRespDD[icent][iS] = new TH2D( Form("hresp_c%i", icent), Form("hresp_c%i", icent), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
      hRespDD[icent][iS]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
      hRespDD[icent][iS]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
      hRespDD[icent][iS]->SetOption("colz");

      //-- Loop over the number of events in each cent, EP and qn bin and fill the 2D response fn
      int Nevt = hSmear2D[icent][iS]->GetEntries();
      std::cout<< "Processing centBin = "<< icent << std::endl;
      for(int ievt = 0 ; ievt < Nevt; ievt++){

	double vnPrior = hPrior[icent][iS]->GetRandom();
	double vnObsx  = 0.;
	double vnObsy  = 0.;
	double smearx  = hSmearX[icent][iS]->GetRandom();
	double smeary  = hSmearY[icent][iS]->GetRandom();

	vnObsx = vnPrior + smearx;
	vnObsy = smeary;

	double vnObs = TMath::Sqrt( pow(vnObsx,2) + pow(vnObsy,2) );
	hRespDD[icent][iS]->Fill( vnObs, vnPrior );

      } //-- End event loop

      //-- Write the histograms to the output file
      fOut[iS]->cd();
      hPrior[icent][iS]->Write();
      hRespDD[icent][iS]->Write();

    } //-- End cent loop
  } //-- End Split loop
}

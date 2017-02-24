#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void makeDDResp(){

  bool vnObsPrior   = 1;
  const int norder_ = 3;

  TFile * fAna;
  TH1D * hVnObs[NCENT][NEPSymm][NQN];
  TH2D * hSmear2D[NCENT][NEPSymm][NQN];
  TH1D * hSmearX[NCENT][NEPSymm][NQN];
  TH1D * hSmearY[NCENT][NEPSymm][NQN];

  TFile * fOut;
  TH1D * hPrior[NCENT][NEPSymm][NQN];
  TH2D * hRespDD[NCENT][NEPSymm][NQN];

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
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	//-- Get the VnObs histogram if that is the prior you wish to use
	hVnObs[icent][iEP][iqn] = 0;
	if( vnObsPrior ){
	  hVnObs[icent][iEP][iqn] = (TH1D*) fAna->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	  hPrior[icent][iEP][iqn] = (TH1D*) hVnObs[icent][iEP][iqn]->Clone( Form("hPrior_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	}

	//-- Grab the 2D smearing histogram and project onto the x and y axes
	hSmear2D[icent][iEP][iqn] = (TH2D*) fAna->Get( Form("qwebye/h2Vn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hSmearX[icent][iEP][iqn]  = (TH1D*) hSmear2D[icent][iEP][iqn]->ProjectionX( Form("hSmearX_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	hSmearY[icent][iEP][iqn]  = (TH1D*) hSmear2D[icent][iEP][iqn]->ProjectionY( Form("hSmearY_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	//-- Initialize the DD response histogram
	hRespDD[icent][iEP][iqn] = new TH2D(Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
	hRespDD[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
	hRespDD[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
	hRespDD[icent][iEP][iqn]->SetOption("colz");

	//-- Loop over the number of events in each cent, EP and qn bin and fill the 2D response fn
	int Nevt = hSmear2D[icent][iEP][iqn]->GetEntries();
	std::cout<< "Processing centBin = "<< icent << "\tEP = " << EPSymmNames[iEP] << "\tqBin = " << iqn << std::endl;
	for(int ievt = 0 ; ievt < Nevt; ievt++){

	  double vnPrior = hPrior[icent][iEP][iqn]->GetRandom();
	  double vnObsx = 0.;
	  double vnObsy = 0.;
	  double smearx = hSmearX[icent][iEP][iqn]->GetRandom();
	  double smeary = hSmearY[icent][iEP][iqn]->GetRandom();

	  vnObsx = vnPrior + smearx;
	  vnObsy = smeary;

	  double vnObs = TMath::Sqrt( pow(vnObsx,2) + pow(vnObsy,2) );
	  hRespDD[icent][iEP][iqn]->Fill( vnObs, vnPrior );

	} //-- End event loop

	//-- Write the histograms to the output file
	fOut->cd();
	hPrior[icent][iEP][iqn]->Write();
	hRespDD[icent][iEP][iqn]->Write();

      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop



}

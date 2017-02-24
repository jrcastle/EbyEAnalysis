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
  const int norder_ = 4;

  TFile * fAna[NSPLIT];
  TH1D * hVnObs[NCENT][NEPSymm][NQN][NSPLIT];
  TH2D * hSmear2D[NCENT][NEPSymm][NQN][NSPLIT];
  TH1D * hSmearX[NCENT][NEPSymm][NQN][NSPLIT];
  TH1D * hSmearY[NCENT][NEPSymm][NQN][NSPLIT];

  TFile * fOut[NSPLIT];
  TH1D * hPrior[NCENT][NEPSymm][NQN][NSPLIT];
  TH2D * hRespDD[NCENT][NEPSymm][NQN][NSPLIT];

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
      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;
	for(int iqn = 0; iqn < NQN; iqn++){

	  //-- Get the VnObs histogram if that is the prior you wish to use
	  hVnObs[icent][iEP][iqn][iS] = 0;
	  if( vnObsPrior ){
	    hVnObs[icent][iEP][iqn][iS] = (TH1D*) fAna[iS]->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	    hPrior[icent][iEP][iqn][iS] = (TH1D*) hVnObs[icent][iEP][iqn][iS]->Clone( Form("hPrior_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	  }

	  //-- Grab the 2D smearing histogram and project onto the x and y axes
	  hSmear2D[icent][iEP][iqn][iS] = (TH2D*) fAna[iS]->Get( Form("qwebye/h2Vn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	  hSmearX[icent][iEP][iqn][iS]  = (TH1D*) hSmear2D[icent][iEP][iqn][iS]->ProjectionX( Form("hSmearX_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );
	  hSmearY[icent][iEP][iqn][iS]  = (TH1D*) hSmear2D[icent][iEP][iqn][iS]->ProjectionY( Form("hSmearY_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	  //-- Initialize the DD response histogram
	  hRespDD[icent][iEP][iqn][iS] = new TH2D(Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
	  hRespDD[icent][iEP][iqn][iS]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
	  hRespDD[icent][iEP][iqn][iS]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
	  hRespDD[icent][iEP][iqn][iS]->SetOption("colz");

	  //-- Loop over the number of events in each cent, EP and qn bin and fill the 2D response fn
	  int Nevt = hSmear2D[icent][iEP][iqn][iS]->GetEntries();
	  std::cout<< "Processing Split = " << iS << "\tcentBin = " << icent << "\tEP = " << EPSymmNames[iEP] << "\tqBin = " << iqn << std::endl;
	  for(int ievt = 0 ; ievt < Nevt; ievt++){

	    double vnPrior = hPrior[icent][iEP][iqn][iS]->GetRandom();
	    double vnObsx = 0.;
	    double vnObsy = 0.;
	    double smearx = hSmearX[icent][iEP][iqn][iS]->GetRandom();
	    double smeary = hSmearY[icent][iEP][iqn][iS]->GetRandom();

	    vnObsx = vnPrior + smearx;
	    vnObsy = smeary;

	    double vnObs = TMath::Sqrt( pow(vnObsx,2) + pow(vnObsy,2) );
	    hRespDD[icent][iEP][iqn][iS]->Fill( vnObs, vnPrior );

	  } //-- End event loop

	  //-- Write the histograms to the output file
	  fOut[iS]->cd();
	  hPrior[icent][iEP][iqn][iS]->Write();
	  hRespDD[icent][iEP][iqn][iS]->Write();

	} //-- End qn loop
      } //-- End EP loop
    } //-- End cent loop
  }//-- End split loop

}

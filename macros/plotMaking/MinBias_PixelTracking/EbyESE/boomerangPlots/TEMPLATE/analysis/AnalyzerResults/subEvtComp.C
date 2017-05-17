#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

void subEvtComp(){

  TFile * fAna;
  TH1I * Mult0[NCENT];
  TH1I * Mult1[NCENT];
  TH1D * vnObs0[NCENT];
  TH1D * vnObs1[NCENT];
  TH2D * hVn2D0v1[NCENT];

  //
  // MAIN
  //

  fAna = new TFile("CastleEby_TEST.root");
  for(int icent = 0; icent < NCENT; icent++){

    Mult0[icent]    = (TH1I*) fAna->Get( Form("qwebye/Mult0_c%i", icent) );
    Mult1[icent]    = (TH1I*) fAna->Get( Form("qwebye/Mult1_c%i", icent) );
    vnObs0[icent]   = (TH1D*) fAna->Get( Form("qwebye/hVnSub0_c%i", icent) );
    vnObs1[icent]   = (TH1D*) fAna->Get( Form("qwebye/hVnSub1_c%i", icent) );
    hVn2D0v1[icent] = (TH2D*) fAna->Get( Form("qwebye/hVn2D0v1_c%i", icent) );

    double mMult0  = Mult0[icent]->GetMean();
    double mMult1  = Mult1[icent]->GetMean();
    double mVn0    = vnObs0[icent]->GetMean();
    double mVn1    = vnObs1[icent]->GetMean();
    double mSEdffx = hVn2D0v1[icent]->ProjectionX()->GetMean();
    double mSEdffy = hVn2D0v1[icent]->ProjectionY()->GetMean();

    std::cout << Form("\n========= Cent %i =========", icent) << std::endl;
    std::cout << "<Mult0> = " << mMult0  << "\t<Mult1> = " << mMult1  << "\tPctDiff = " << 100.*fabs(1. - mMult0/mMult1) << "%" << std::endl;
    std::cout << "<Vn0>   = " << mVn0    << "\t<Vn1>   = " << mVn1    << "\tPctDiff = " << 100.*fabs(1. - mVn0/mVn1)     << "%" << std::endl; 
    std::cout << "<SEDx>  = " << mSEdffx << "\t<SEDy>  = " << mSEdffy << std::endl;
  
  }

}

#ifndef HeavyIonsAnalysis_EbyEAnalysis_TrackEfficiency_h
#define HeavyIonsAnalysis_EbyEAnalysis_TrackEfficiency_h

#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"


static const int nc        = 5;
static const float cbins[] = {0, 5, 10, 30, 50, 100};
static const int cmn[]     = {0, 5, 10, 30, 50};
static const int cmx[]     = {5, 10, 30, 50, 100}; 

class TrackEfficiency{
  
 public:
  //-- Constructor/Destructor
  TrackEfficiency(TString effOption);
  ~TrackEfficiency();

  //-- Methods
  double GetEff(double centval, double pT, double eta);
  
 private:
  TFile * fEff;
  TH1F * hCentBinning;
  TH2F * hEff_ptEta[nc];

  
};

//-- ===================================|===================================
//--                               CONSTRUCTOR
//-- ===================================|===================================

inline TrackEfficiency::TrackEfficiency(TString effFile){    

  fEff = new TFile( effFile );
  for( int icent = 0; icent < nc; icent++) hEff_ptEta[icent] = (TH2F*) fEff->Get( Form("Eff_%i_%i", cmn[icent], cmx[icent]) );

  hCentBinning = new TH1F("hCentBinning", "hCentBinning", nc, cbins);

}

//-- ===================================|===================================
//--                                DESTRUCTOR
//-- ===================================|===================================

inline TrackEfficiency::~TrackEfficiency(){

  if( fEff ) delete fEff;
  if( hCentBinning ) delete hCentBinning;
  for( int icent = 0; icent < nc; icent++){
    if( hEff_ptEta[icent] ) delete hEff_ptEta[icent];
  }

}

//-- ===================================|===================================
//--                                  GetEff
//-- ===================================|===================================

inline double TrackEfficiency::GetEff(double centval, double pT, double eta){

  int icent  = hCentBinning->FindBin(centval) - 1;
  int etaBin = hEff_ptEta[icent]->GetXaxis()->FindBin(eta);
  int ptBin  = hEff_ptEta[icent]->GetYaxis()->FindBin(pT);

  double w = hEff_ptEta[icent]->GetBinContent(etaBin, ptBin);
  return w;

}

    
#endif

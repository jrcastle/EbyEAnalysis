#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLine.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void makeQnDet(){

  bool testrun     = 0;
  int    norder_   = 2;
  double vtxCut    = 15.;

  TFile * fIn;
  TTree * tree;
  double centval;
  double vtx;
  double qnHFx_EP[NumEPNames];
  double qnHFy_EP[NumEPNames];
  double sumET_EP[NumEPNames];

  TFile * fOut;
  TH1D * hqnHFDet_x[NumEPNames];
  TH1D * hqnHFDet_y[NumEPNames];

  double qnx_AllEvts[NCENT][NumEPNames];
  double qny_AllEvts[NCENT][NumEPNames];
  int NEvents[NCENT][NumEPNames];
  int NFails[NCENT][NumEPNames];

  //
  // Main
  //
  setTDRStyle();

  //-- Get the analyzer file
  fIn = new TFile(fAnaTreeName);

  //-- Get the analyzer tree and set the branches
  tree = (TTree*) fIn->Get("ebyeana/tree");
  tree->SetBranchAddress("Cent", &centval);
  tree->SetBranchAddress("Vtx",  &vtx);
  tree->SetBranchAddress("qnHFx_EP",  &qnHFx_EP);
  tree->SetBranchAddress("qnHFy_EP",  &qnHFy_EP);
  tree->SetBranchAddress("sumET_EP",  &sumET_EP);

  //-- Set up the output file
  fOut = new TFile( Form("Q%iDet.root", norder_), "recreate");

  //Initialize all the QnDet histograms
  for(int iEP = 0; iEP < NumEPNames; iEP++){

    fOut->cd();
    hqnHFDet_x[iEP] = new TH1D( Form("hqnHFDet_x_%s", EPNames[iEP].data()), Form("hqnHFDet_x_%s", EPNames[iEP].data()), NCENT, centbinsDefault );
    hqnHFDet_y[iEP] = new TH1D( Form("hqnHFDet_y_%s", EPNames[iEP].data()), Form("hqnHFDet_y_%s", EPNames[iEP].data()), NCENT, centbinsDefault );

    for(int icent = 0; icent < NCENT; icent++){

      qnx_AllEvts[icent][iEP] = 0.;
      qny_AllEvts[icent][iEP] = 0.;
      NEvents[icent][iEP]     = 0;
      NFails[icent][iEP]      = 0;

    }

  }

  //-- Begin Event loop
  int Nevt;
  if(!testrun) Nevt = tree->GetEntries();
  else         Nevt = 10000;
  std::cout<<"Begin QN DETECTOR loop, contains "<< Nevt << " events..." << std::endl;
  for(int ievt = 0; ievt < Nevt; ievt++){

    if((ievt+1)% 500000 == 0) cout<<"Processing Event "<<ievt+1<<"\t"<<(100.*(ievt+1)/Nevt)<<"% Completed"<<endl;
    tree->GetEntry(ievt);

    //-- Vertex Cut
    if(TMath::Abs(vtx) > vtxCut) continue;

    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    int icent = hCentBins.FindBin(centval)-1;

    //-- Begin EP loop
    for(int iEP = 0; iEP < NEP; iEP++){

      NEvents[icent][iEP]++;
      double qx = qnHFx_EP[iEP];
      double qy = qnHFy_EP[iEP];
      double sumET = sumET_EP[iEP];

      if(sumET != 0){
	qx /= sumET;
	qy /= sumET;
	qnx_AllEvts[icent][iEP] += qx;
	qny_AllEvts[icent][iEP] += qy;
      }
      else NFails[icent][iEP]++;

    } //-- End EP loop

  } //-- End event loop
  std::cout<<"End QN DETECTOR loop"<<std::endl;

  //-- Average over all events
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NumEPNames; iEP++){

      double qxDet = qnx_AllEvts[icent][iEP];
      double qyDet = qny_AllEvts[icent][iEP];
      double N = (double)NEvents[icent][iEP] - (double)NFails[icent][iEP];

      qxDet /= N;
      qyDet /= N;

      hqnHFDet_x[iEP]->SetBinContent(icent+1, qxDet);
      hqnHFDet_y[iEP]->SetBinContent(icent+1, qyDet);

    }
  }

  fOut->Write();

}

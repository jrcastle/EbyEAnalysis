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

  TFile * fSplit;
  TTree * treeSplit;
  int iSplit;

  TFile * fOut;
  TH1D * hqnHFDet_x[NumEPNames][NSPLIT];
  TH1D * hqnHFDet_y[NumEPNames][NSPLIT];

  double qnx_AllEvts[NCENT][NumEPNames][NSPLIT];
  double qny_AllEvts[NCENT][NumEPNames][NSPLIT];
  int NEvents[NCENT][NumEPNames][NSPLIT];
  int NFails[NCENT][NumEPNames][NSPLIT];

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

  //-- Set up the file splitter
  fSplit = new TFile(fileSplit);
  treeSplit = (TTree*) fSplit->Get("SplitTree");
  treeSplit->SetBranchAddress("iSplit", &iSplit);

  //-- Set up the output file
  fOut = new TFile( Form("Q%iDet.root", norder_), "recreate");

  //Initialize all the QnDet histograms
  for(int iEP = 0; iEP < NumEPNames; iEP++){
    for(int iS = 0; iS < NSPLIT; iS++){

      fOut->cd();
      hqnHFDet_x[iEP][iS] = new TH1D( Form("hqnHFDet_x_%s_Split%i", EPNames[iEP].data(), iS), Form("hqnHFDet_x_%s_Split%i", EPNames[iEP].data(), iS), NCENT, centbinsDefault );
      hqnHFDet_y[iEP][iS] = new TH1D( Form("hqnHFDet_y_%s_Split%i", EPNames[iEP].data(), iS), Form("hqnHFDet_y_%s_Split%i", EPNames[iEP].data(), iS), NCENT, centbinsDefault );

      for(int icent = 0; icent < NCENT; icent++){

	qnx_AllEvts[icent][iEP][iS] = 0.;
	qny_AllEvts[icent][iEP][iS] = 0.;
	NEvents[icent][iEP][iS]     = 0;
	NFails[icent][iEP][iS]      = 0;

      }
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
    treeSplit->GetEntry(ievt);

    //-- Vertex Cut
    if(TMath::Abs(vtx) > vtxCut) continue;

    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    int icent = hCentBins.FindBin(centval)-1;

    //-- Begin EP loop
    for(int iEP = 0; iEP < NEP; iEP++){

      NEvents[icent][iEP][iSplit]++;
      double qx = qnHFx_EP[iEP];
      double qy = qnHFy_EP[iEP];
      double sumET = sumET_EP[iEP];

      if(sumET != 0){
	qx /= sumET;
	qy /= sumET;
	qnx_AllEvts[icent][iEP][iSplit] += qx;
	qny_AllEvts[icent][iEP][iSplit] += qy;
      }
      else NFails[icent][iEP][iSplit]++;

    } //-- End EP loop

  } //-- End event loop
  std::cout<<"End QN DETECTOR loop"<<std::endl;

  //-- Average over all events
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NumEPNames; iEP++){
      for(int iS = 0; iS < NSPLIT; iS++){

	double qxDet = qnx_AllEvts[icent][iEP][iS];
	double qyDet = qny_AllEvts[icent][iEP][iS];
	double N = (double)NEvents[icent][iEP][iS] - (double)NFails[icent][iEP][iS];

	qxDet /= N;
	qyDet /= N;

	hqnHFDet_x[iEP][iS]->SetBinContent(icent+1, qxDet);
	hqnHFDet_y[iEP][iS]->SetBinContent(icent+1, qyDet);

      }
    }
  }

  fOut->Write();

}

#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void makeQnCuts(){

  bool testrun    = 0;
  const int NBins = 1000;

  double qn_min  = 0;
  double qn_max  = 0.61;
  int    qnselect_ = 2;
  double vtxCut  = 15.;

  TFile * fIn;
  TTree * tree;
  double centval;
  double vtx;
  double qnHFx_EP[NumEPNames];
  double qnHFy_EP[NumEPNames];
  double sumET_EP[NumEPNames];

  TFile * fQNDet;
  TH1D * hqnHFDet_x[NumEPNames];
  TH1D * hqnHFDet_y[NumEPNames];

  TFile * fOut;
  TH1D * hqnHF_EP[NCENT][NEPSymm];
  TCanvas * chqnHF_EP[NCENT][NEPSymm];

  double splitBinMin[NCENT][NEPSymm][NQN];
  double splitBinMax[NCENT][NEPSymm][NQN];

  TLine * splitLine[NCENT][NEPSymm][NQN-1];

  TH1D * hqbins[NCENT][NEPSymm];
  double qbins[NCENT][NEPSymm][NQN+1];

  TLatex latex;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();

  //-- Get the analyzer file
  fIn = new TFile(fAnaTreeName);

  //-- Get the analyzer tree and set the branches
  tree = (TTree*) fIn->Get("ebyeana/tree");
  tree->SetBranchAddress("Cent", &centval);
  tree->SetBranchAddress("Vtx",  &vtx);
  tree->SetBranchAddress("qnHFx_EP",  &qnHFx_EP);
  tree->SetBranchAddress("qnHFy_EP",  &qnHFy_EP);
  tree->SetBranchAddress("sumET_EP",  &sumET_EP);

  //-- Get the QN Detector histograms
  fQNDet = new TFile( Form("../../../../../../v%i/eta2.4/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/Q%iDet.root", qnselect_, qnselect_) );
  for(int iEP = 0; iEP < NumEPNames; iEP++){
    int EPbin  = EPSymmPartnerBin[iEP];
    if( EPbin != EPSymmBin ) continue;
    hqnHFDet_x[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_x_%s", EPNames[iEP].data()) );
    hqnHFDet_y[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_y_%s", EPNames[iEP].data()) );
  }

  //-- Make the output file
  fOut = new TFile( Form("q%iCuts.root", qnselect_), "recreate" );

  //-- Initialize the 1D histograms
  for(int iEP = 0; iEP < NEPSymm; iEP++){
    if( iEP != EPSymmBin ) continue;
    for(int icent = 0; icent < NCENT; icent++){
      fOut->cd();
      hqnHF_EP[icent][iEP] = new TH1D( Form("hqnHF_c%i_EP%i", icent, iEP), Form("hqnHF_c%i_EP%i", icent, iEP), NBins, qn_min, qn_max );
      hqnHF_EP[icent][iEP]->GetXaxis()->SetTitle( Form("q_{%i}", qnselect_) );
      hqnHF_EP[icent][iEP]->GetYaxis()->SetTitle( "Events" );
    }
  }


  //-- Begin Event loop
  int Nevt;
  if(!testrun) Nevt = tree->GetEntries();
  else         Nevt = 10000;
  std::cout<<"Begin QN SPLIT loop, contains "<< Nevt << " events..." << std::endl;
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

      int EPbin  = EPSymmPartnerBin[iEP];
      if( EPbin != EPSymmBin ) continue;

      double qx    = qnHFx_EP[iEP];
      double qy    = qnHFy_EP[iEP];
      double qxDet = hqnHFDet_x[iEP]->GetBinContent(icent+1);
      double qyDet = hqnHFDet_y[iEP]->GetBinContent(icent+1);
      double sumET = sumET_EP[iEP];
      if(sumET == 0) continue;

      qx -= qxDet;
      qy -= qyDet;
      qx /= sumET;
      qy /= sumET;
      double qn    = TMath::Sqrt( qx*qx + qy*qy );

      //-- Consolodate into symmetric HF bins
      hqnHF_EP[icent][EPbin]->Fill(qn);

    } //-- End EP loop

  } //-- End event loop
  std::cout<<"End QN loop"<<std::endl;

  //-- Divide each histogram into NQN equal areas
  bool splitSuccess = 1;

  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;

      double integral = hqnHF_EP[icent][iEP]->Integral(1, NBins);
      double binwidth = integral/(double)NQN;
      int startbin = 1;
      double sum = 0.;

      for(int iqn = 0; iqn < NQN-1; iqn++){
	splitBinMin[icent][iEP][iqn] = startbin;

	for(int ibin = startbin; ibin <= NBins; ibin++){
	  double A = hqnHF_EP[icent][iEP]->Integral(startbin, ibin);
	  if( A > binwidth ){
	    double B = hqnHF_EP[icent][iEP]->Integral(startbin, ibin-1);
	    if( fabs(A - binwidth) < fabs(B - binwidth) ){
	      startbin = ibin+1;
	      splitBinMax[icent][iEP][iqn] = ibin;
	      break;
	    }
	    else{
	      startbin = ibin;
	      splitBinMax[icent][iEP][iqn] = ibin-1;
	      break;
	    }
	  }

	}//-- End bin loop

	sum += hqnHF_EP[icent][iEP]->Integral(splitBinMin[icent][iEP][iqn], splitBinMax[icent][iEP][iqn]);

	double binc = hqnHF_EP[icent][iEP]->GetBinCenter(startbin);
	double max  = hqnHF_EP[icent][iEP]->GetBinContent(startbin);
	splitLine[icent][iEP][iqn] = new TLine(binc, 0, binc, max);
	splitLine[icent][iEP][iqn]->SetLineColor(2);
	splitLine[icent][iEP][iqn]->SetLineStyle(2);
	splitLine[icent][iEP][iqn]->SetLineWidth(2);

      } //-- End split loop

      splitBinMin[icent][iEP][NQN-1] = splitBinMax[icent][iEP][NQN-2] + 1;

      // Make the last bin edge to be the bin where the integral is full
      for(int ibin = 1; ibin <= NBins; ibin++){

	double integ = hqnHF_EP[icent][iEP]->Integral(1, ibin);
	if( integ == integral ){
	  splitBinMax[icent][iEP][NQN-1] = ibin;
	  break;
	}

      } //-- End bin loop

      sum += hqnHF_EP[icent][iEP]->Integral(splitBinMin[icent][iEP][NQN-1], splitBinMax[icent][iEP][NQN-1]);
      if(sum == integral) std::cout<<"Histogram division check..."<<"\tCent bin = "<<icent<<"\tEP bin = "<<iEP<<"\tPASS!"<<std::endl;
      else{
	std::cout<<"Histogram division check..."<<"\tCent bin = "<<icent<<"\tEP bin = "<<iEP<<"\tFAIL!"<<std::endl;
	splitSuccess = 0;
      }

    } //-- End EP loop
  } //-- End cent loop

  if( !splitSuccess ){
    std::cout << "Histogram area splitting failed in one or more bin!" << std::endl;
    std::cout << "Refusing to go futher, please fix and run again." <<std::endl;
    std::cout << "Exiting... Have a nice day." << std::endl;
    //exit(0);
  }

  //-- If all is good, proceed to fill the 2D arrays that contain the bin boundaries for qn
  for(int iqn = 0; iqn <= NQN; iqn++){
    for(int icent = 0; icent < NCENT; icent++){
      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;

	if( iqn < NQN ){
	  int bmin = splitBinMin[icent][iEP][iqn];
	  double qmin = hqnHF_EP[icent][iEP]->GetBinLowEdge(bmin);
	  qbins[icent][iEP][iqn] = qmin; 
	}
	else{
	  int bmax = splitBinMax[icent][iEP][NQN-1] + 1;
	  double qmax = hqnHF_EP[icent][iEP]->GetBinLowEdge(bmax);
	  qbins[icent][iEP][iqn] = qmax;
	}
      }
    }
  }

  //-- Set up histograms based on these arrays (so I can use the find bin option in the subsequent programs)
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;

      fOut->cd();
      hqbins[icent][iEP] = new TH1D( Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent), Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent), NQN, qbins[icent][iEP] );
      hqbins[icent][iEP]->GetXaxis()->SetTitle( Form("q_%i", qnselect_) );

      for(int iqn = 1; iqn <= NQN; iqn++) hqbins[icent][iEP]->SetBinContent( iqn, iqn );

   
    }
  }

  fOut->Write();

}

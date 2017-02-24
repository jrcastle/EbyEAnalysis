#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"

#include <iostream>

using namespace ebyese;

void makeVNDet(){

  bool testrun          = 0;
  const int norder_     = 2;
  const double vtxCut   = 15.;

  const double ptMin = 0.3;
  const double ptMax = 3.00;
  const double etaMax = 1.0;

  static const int ptBinMin  = 0;
  static const int ptBinMax  = nptbinsDefault-1;
  static const int etaBinMin = 0;  //0;
  static const int etaBinMax = netabinsDefault-1;

  TFile * fAna;
  TTree * tree;
  double centval;
  double vtx;
  TH2D * sumw;
  TH2D * sumwqx;
  TH2D * sumwqy;
  TH2I * hMult;

  TFile * fSplit;
  TTree * treeSplit;
  int iSplit;

  TFile * fOut;

  TH1D * hVNDetX_0[NSPLIT];
  TH1D * hVNDetY_0[NSPLIT];
  TH1D * hVNDetX_1[NSPLIT];
  TH1D * hVNDetY_1[NSPLIT];
  TH1D * hVNDetX_full[NSPLIT];
  TH1D * hVNDetY_full[NSPLIT];

  double VNDetX_0[NCENT][NSPLIT];
  double VNDetY_0[NCENT][NSPLIT];
  double VNDetX_1[NCENT][NSPLIT];
  double VNDetY_1[NCENT][NSPLIT];
  double VNDetX_full[NCENT][NSPLIT];
  double VNDetY_full[NCENT][NSPLIT];

  int Nevents[NCENT][NSPLIT];
  int NFails[NCENT][NSPLIT];

  //
  // MAIN
  //


  //-- Set up the analyzer objects
  fAna = new TFile( fAnaTreeNameVtx3_15 );

  tree = (TTree*) fAna->Get("ebyeana/tree");
  sumwqx = new TH2D(Form("sumwqx%i", norder_), Form("sumwqx%i", norder_), nptbinsDefault, ptbinsDefault, netabinsDefault, etabinsDefault);
  sumwqy = new TH2D(Form("sumwqy%i", norder_), Form("sumwqy%i", norder_), nptbinsDefault, ptbinsDefault, netabinsDefault, etabinsDefault);
  sumw   = new TH2D("sumw",                    "sumw",                    nptbinsDefault, ptbinsDefault, netabinsDefault, etabinsDefault);
  hMult  = new TH2I("hMult",                   "hMult",                   nptbinsDefault, ptbinsDefault, netabinsDefault, etabinsDefault);
    
  tree->SetBranchAddress("Cent",                    &centval);
  tree->SetBranchAddress("Vtx",                     &vtx);
  tree->SetBranchAddress("mult",                    &hMult);
  tree->SetBranchAddress(Form("sumwqx%i", norder_), &sumwqx);
  tree->SetBranchAddress(Form("sumwqy%i", norder_), &sumwqy);
  tree->SetBranchAddress("sumw",                    &sumw);

  //-- Set up the file splitter
  fSplit = new TFile(fileSplit);
  treeSplit = (TTree*) fSplit->Get("SplitTree");
  treeSplit->SetBranchAddress("iSplit", &iSplit);

  //-- Setup the output objects
  fOut     = new TFile( Form("V%iDet.root", norder_), "recreate" );
  for(int iS = 0; iS < NSPLIT; iS++){
    fOut->cd();
    hVNDetX_0[iS]    = new TH1D( Form("hVNDetX_0_Split%i", iS),    Form("hVNDetX_0_Split%i", iS),    NCENT, centbinsDefault);
    hVNDetX_0[iS]->GetXaxis()->SetTitle("Centrality %");
    
    hVNDetY_0[iS]    = new TH1D( Form("hVNDetY_0_Split%i", iS),    Form("hVNDetY_0_Split%i", iS),    NCENT, centbinsDefault);
    hVNDetY_0[iS]->GetXaxis()->SetTitle("Centrality %");
    
    hVNDetX_1[iS]    = new TH1D( Form("hVNDetX_1_Split%i", iS),    Form("hVNDetX_1_Split%i", iS),    NCENT, centbinsDefault);
    hVNDetX_1[iS]->GetXaxis()->SetTitle("Centrality %");
    
    hVNDetY_1[iS]    = new TH1D( Form("hVNDetY_1_Split%i", iS),    Form("hVNDetY_1_Split%i", iS),    NCENT, centbinsDefault);
    hVNDetY_1[iS]->GetXaxis()->SetTitle("Centrality %");
    
    hVNDetX_full[iS] = new TH1D( Form("hVNDetX_full_Split%i", iS), Form("hVNDetX_full_Split%i", iS), NCENT, centbinsDefault);
    hVNDetX_full[iS]->GetXaxis()->SetTitle("Centrality %");
    
    hVNDetY_full[iS] = new TH1D( Form("hVNDetY_full_Split%i", iS), Form("hVNDetY_full_Split%i", iS), NCENT, centbinsDefault);
    hVNDetY_full[iS]->GetXaxis()->SetTitle("Centrality %");
  }    

  //-- initialize all variables
  for(int icent = 0; icent<NCENT; icent++){
    for(int iS = 0; iS < NSPLIT; iS++){

      VNDetX_0[icent][iS]     = 0.;
      VNDetY_0[icent][iS]     = 0.;
      VNDetX_1[icent][iS]     = 0.;
      VNDetY_1[icent][iS]     = 0.;
      VNDetX_full[icent][iS]  = 0.;
      VNDetY_full[icent][iS]  = 0.;

      Nevents[icent][iS]      = 0;
      NFails[icent][iS]       = 0;

    }
  }
    
  //
  // Calculate Vn_det
  //
    
  cout<<"Begin DETECTOR loop, contains "<<tree->GetEntries()<<" Events"<<endl;
    
  int N;
  if(testrun) N = 10000; 
  else        N = tree->GetEntries();

  //-- Begin event loop
  for(int ievent = 0; ievent < N; ievent++) {
        
    if((ievent+1)% 500000 == 0) cout << "Processing Event " << ievent+1 << "\t" << (100.*(ievent+1)/N) << "% Completed" << endl;
        
    int itree = tree->GetEntry(ievent);
    if(itree < 0){
      std::cout<<"!!! BAD EVENT !!!"<<std::endl;
      continue;
    }

    treeSplit->GetEntry(ievent);

    //-- Vertex Cut
    if(TMath::Abs(vtx) < 3.0 || TMath::Abs(vtx) > 15.0) continue;
      
    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    int icent = (centval - cent_min[0]) / centBinWidth;

    //-- Reset Raw and sumw values
    double VNRawX_0     = 0.;
    double VNRawY_0     = 0.;
    double VNRawX_1     = 0.;
    double VNRawY_1     = 0.;
    double VNRawX_full  = 0.;
    double VNRawY_full  = 0.;            
	
    double sumw_0       = 0.;
    double sumw_1       = 0.;
    double sumw_full    = 0.;
      
    double evtMult_0    = 0;
    double evtMult_1    = 0;
    double evtMult_full = 0;
      
    Nevents[icent][iSplit]++;

    //-- Begin analyzer histogram loops
    for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++){
      for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++){
	  
	double pt  = sumw->GetXaxis()->GetBinCenter(ipt+1);
        double eta = sumw->GetYaxis()->GetBinCenter(ieta+1);
        if( pt < ptMin || pt > ptMax ) continue;
        if( fabs( eta ) > etaMax )     continue;

	if(sumw->GetBinContent(ipt+1,ieta+1) !=0){
	    
	  //-- Subevent 0 (eta >= 0)
	  if(etabinsDefault[ieta] >= 0){
	    VNRawX_0   += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VNRawY_0   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	    sumw_0     += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_0  += hMult->GetBinContent(ipt+1,ieta+1);
	  }
	  //-- Subevent 1 (eta < 0)
	  else{
	    VNRawX_1   += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VNRawY_1   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	    sumw_1     += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_1  += hMult->GetBinContent(ipt+1,ieta+1);
	  }
	  //-- Full Event
	  VNRawX_full  += sumwqx->GetBinContent(ipt+1,ieta+1);
	  VNRawY_full  += sumwqy->GetBinContent(ipt+1,ieta+1);
	    
	  sumw_full    += sumw->GetBinContent(ipt+1,ieta+1);
	  evtMult_full += hMult->GetBinContent(ipt+1,ieta+1);
	}
	  
      } //-- End eta loop
      
    } //-- End pt loop
    
    //-- Only use events that have tracks in all subevents AND subevents that have at least two tracks to determine VN
    if( sumw_0 == 0 || sumw_1 == 0 || evtMult_0 < 2 || evtMult_1 < 2 ){
	NFails[icent][iSplit]++;
    }
    else{      
      VNDetX_0[icent][iSplit]    += VNRawX_0    / sumw_0;
      VNDetY_0[icent][iSplit]    += VNRawY_0    / sumw_0;
      VNDetX_1[icent][iSplit]    += VNRawX_1    / sumw_1;
      VNDetY_1[icent][iSplit]    += VNRawY_1    / sumw_1;
      VNDetX_full[icent][iSplit] += VNRawX_full / sumw_full;
      VNDetY_full[icent][iSplit] += VNRawY_full / sumw_full;
    }

  } //-- End event loop
    
  std::cout<<"End DETECTOR loop"<<std::endl;

  //-- Average VNDet over all events for each centrality, EP and, qn bin
  for(int icent = 0; icent < NCENT; icent++){
    for(int iS = 0; iS < NSPLIT; iS++){
   
      double  Neffective = (double) Nevents[icent][iS] - (double) NFails[icent][iS];
      if( Neffective == 0 ) continue;
    
      VNDetX_0[icent][iS]    /= Neffective;
      VNDetY_0[icent][iS]    /= Neffective;
      VNDetX_1[icent][iS]    /= Neffective;
      VNDetY_1[icent][iS]    /= Neffective;
      VNDetX_full[icent][iS] /= Neffective;
      VNDetY_full[icent][iS] /= Neffective;

      //-- Populate histograms that will be used by ReadTree_normDet.C
      hVNDetX_0[iS]    -> SetBinContent(icent+1, VNDetX_0[icent][iS]);
      hVNDetX_1[iS]    -> SetBinContent(icent+1, VNDetX_1[icent][iS]);
      hVNDetX_full[iS] -> SetBinContent(icent+1, VNDetX_full[icent][iS]);
      
      hVNDetY_0[iS]    -> SetBinContent(icent+1, VNDetY_0[icent][iS]);
      hVNDetY_1[iS]    -> SetBinContent(icent+1, VNDetY_1[icent][iS]);
      hVNDetY_full[iS] -> SetBinContent(icent+1, VNDetY_full[icent][iS]);

    }
  }
  fOut->Write();
  cout<<"File written, process completed"<<endl;

}

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"

#include <iostream>

using namespace hi;
using namespace ebyese;

void makeVNDet(){

  bool testrun          = 0;
  const int norder_     = 2;
  const int QnBinOrder_ = 2;
  const double vtxCut   = 15.;

  static const int ptBinMin  = 0;
  static const int ptBinMax  = nptbinsDefault-1;
  static const int etaBinMin = 0; //0;
  static const int etaBinMax = netabinsDefault-1;

  TFile * fAna;
  TTree * tree;
  double centval;
  double vtx;
  TH2D * sumw;
  TH2D * sumwqx;
  TH2D * sumwqy;
  TH2I * hMult;
  double qnHFx_EP[NumEPNames];
  double qnHFy_EP[NumEPNames];
  double sumET_EP[NumEPNames];

  TFile * fQNDet;
  TH1D * hqnHFDet_x[NumEPNames];
  TH1D * hqnHFDet_y[NumEPNames];

  TFile * fQN;
  TH1D * hqbins[NCENT][NEPSymm];

  TFile * fOut;
  TDirectory * SubEvt_0;
  TDirectory * SubEvt_1;
  TDirectory * FullEvt;

  TH2D * hVNDetX_0[NQN];
  TH2D * hVNDetY_0[NQN];
  TH2D * hVNDetX_1[NQN];
  TH2D * hVNDetY_1[NQN];
  TH2D * hVNDetX_full[NQN];
  TH2D * hVNDetY_full[NQN];

  double VNRawX_0[NCENT][NEPSymm][NQN];
  double VNRawY_0[NCENT][NEPSymm][NQN];
  double VNRawX_1[NCENT][NEPSymm][NQN];
  double VNRawY_1[NCENT][NEPSymm][NQN];
  double VNRawX_full[NCENT][NEPSymm][NQN];
  double VNRawY_full[NCENT][NEPSymm][NQN];

  double sumw_0[NCENT][NEPSymm][NQN];
  double sumw_1[NCENT][NEPSymm][NQN];
  double sumw_full[NCENT][NEPSymm][NQN];

  int evtMult_0[NCENT][NEPSymm][NQN];
  int evtMult_1[NCENT][NEPSymm][NQN];
  int evtMult_full[NCENT][NEPSymm][NQN];
 
  double VNDetX_0[NCENT][NEPSymm][NQN];
  double VNDetY_0[NCENT][NEPSymm][NQN];
  double VNDetX_1[NCENT][NEPSymm][NQN];
  double VNDetY_1[NCENT][NEPSymm][NQN];
  double VNDetX_full[NCENT][NEPSymm][NQN];
  double VNDetY_full[NCENT][NEPSymm][NQN];

  int Nevents[NCENT][NEPSymm][NQN];
  int NFails[NCENT][NEPSymm][NQN];

  //
  // MAIN
  //


  //-- Set up the analyzer objects
  fAna = new TFile("/rfs/jcastle/PbPb2015/PixelTracking_MB2/EbyETree_pixel_tight.root");

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
  tree->SetBranchAddress("qnHFx_EP",                &qnHFx_EP);
  tree->SetBranchAddress("qnHFy_EP",                &qnHFy_EP);
  tree->SetBranchAddress("sumET_EP",                &sumET_EP);

  //-- Get the QN Detector histograms                                                                                                                                                                                             
  fQNDet = new TFile( Form("../../../../../../v%i/eta2.4/AnalyzerResults/Q%iDet.root", QnBinOrder_, QnBinOrder_) );
  for(int iEP = 0; iEP < NumEPNames; iEP++){
    int EPbin  = EPSymmPartnerBin[iEP];
    if( EPbin != EPSymmBin ) continue;
    hqnHFDet_x[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_x_%s", EPNames[iEP].data()) );
    hqnHFDet_y[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_y_%s", EPNames[iEP].data()) );
  }

  //-- Setup the QN binning objects
  fQN = new TFile( Form( "../../../../../../v%i/eta2.4/AnalyzerResults/q%iCuts.root", QnBinOrder_, QnBinOrder_) );
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      hqbins[icent][iEP] = (TH1D*) fQN->Get( Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent) );
    }
  }

  //-- Setup the output objects
  fOut     = new TFile( Form("V%iDet.root", norder_), "recreate" );
  SubEvt_0 = fOut->mkdir("SubEvt_0");
  SubEvt_1 = fOut->mkdir("SubEvt_1");
  FullEvt  = fOut->mkdir("FullEvt");
    
  for(int iqn = 0; iqn < NQN; iqn++){
    SubEvt_0->cd();
    hVNDetX_0[iqn]    = new TH2D( Form("hVNDetX_0_qbin%i", iqn),    Form("hVNDetX_0_qbin%i", iqn),    NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetX_0[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetX_0[iqn]->SetOption("colz");
    
    hVNDetY_0[iqn]    = new TH2D( Form("hVNDetY_0_qbin%i", iqn),    Form("hVNDetY_0_qbin%i", iqn),    NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetY_0[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetY_0[iqn]->SetOption("colz");
    
    SubEvt_1->cd();
    hVNDetX_1[iqn]    = new TH2D( Form("hVNDetX_1_qbin%i", iqn),    Form("hVNDetX_1_qbin%i", iqn),    NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetX_1[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetX_1[iqn]->SetOption("colz");
    
    hVNDetY_1[iqn]    = new TH2D( Form("hVNDetY_1_qbin%i", iqn),    Form("hVNDetY_1_qbin%i", iqn),    NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetY_1[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetY_1[iqn]->SetOption("colz");
  
    FullEvt->cd();
    hVNDetX_full[iqn] = new TH2D( Form("hVNDetX_full_qbin%i", iqn), Form("hVNDetX_full_qbin%i", iqn), NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetX_full[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetX_full[iqn]->SetOption("colz");
    
    hVNDetY_full[iqn] = new TH2D( Form("hVNDetY_full_qbin%i", iqn), Form("hVNDetY_full_qbin%i", iqn), NCENT, centbinsDefault, NEPSymm, epbinsDefault );
    hVNDetY_full[iqn]->GetXaxis()->SetTitle("Centrality %");
    hVNDetY_full[iqn]->SetOption("colz");
    
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      hVNDetX_0[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
      hVNDetY_0[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
      hVNDetX_1[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
      hVNDetY_1[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
      hVNDetX_full[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
      hVNDetY_full[iqn]->GetYaxis()->SetBinLabel(iEP+1, EPSymmNames[iEP].data());
    }

  }
    
  //-- initialize all variables
  for(int icent = 0; icent<NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){

	VNDetX_0[icent][iEP][iqn]     = 0.;
	VNDetY_0[icent][iEP][iqn]     = 0.;
	VNDetX_1[icent][iEP][iqn]     = 0.;
	VNDetY_1[icent][iEP][iqn]     = 0.;
	VNDetX_full[icent][iEP][iqn]  = 0.;
	VNDetY_full[icent][iEP][iqn]  = 0.;

	evtMult_0[icent][iEP][iqn]    = 0;
	evtMult_1[icent][iEP][iqn]    = 0;
	evtMult_full[icent][iEP][iqn] = 0;

	Nevents[icent][iEP][iqn]      = 0;
	NFails[icent][iEP][iqn]       = 0;

      }
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
        
    tree->GetEntry(ievent);
      
    //-- Vertex Cut
    if(TMath::Abs(vtx) > vtxCut) continue;
      
    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    int icent = hCentBins.FindBin(centval)-1;

    //-- begin EP loop
    for(int iEP = 0; iEP < NEP; iEP++){

      int EPbin  = EPSymmPartnerBin[iEP];
      if( EPbin != EPSymmBin ) continue;

      //-- Calculate qbin
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
      double qn = TMath::Sqrt( qx*qx + qy*qy );
      int   iqn = hqbins[icent][EPbin]->FindBin( qn ) - 1;
      if(iqn >= NQN) continue;

      //-- Reset Raw and sumw values
      VNRawX_0[icent][EPbin][iqn]     = 0.;
      VNRawY_0[icent][EPbin][iqn]     = 0.;
      VNRawX_1[icent][EPbin][iqn]     = 0.;
      VNRawY_1[icent][EPbin][iqn]     = 0.;
      VNRawX_full[icent][EPbin][iqn]  = 0.;
      VNRawY_full[icent][EPbin][iqn]  = 0.;            
	
      sumw_0[icent][EPbin][iqn]       = 0.;
      sumw_1[icent][EPbin][iqn]       = 0.;
      sumw_full[icent][EPbin][iqn]    = 0.;
      
      evtMult_0[icent][EPbin][iqn]    = 0;
      evtMult_1[icent][EPbin][iqn]    = 0;
      evtMult_full[icent][EPbin][iqn] = 0;
      
      Nevents[icent][EPbin][iqn]++;

      //-- Begin analyzer histogram loops
      for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++){
	for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++){
	  
	  if(sumw->GetBinContent(ipt+1,ieta+1) !=0){
	    
	    //-- Subevent 0 (eta >= 0)
	    if(etabinsDefault[ieta] >= 0){
	      VNRawX_0[icent][EPbin][iqn]   += sumwqx->GetBinContent(ipt+1,ieta+1);
	      VNRawY_0[icent][EPbin][iqn]   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	      sumw_0[icent][EPbin][iqn]     += sumw->GetBinContent(ipt+1,ieta+1);
	      evtMult_0[icent][EPbin][iqn]  += hMult->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Subevent 1 (eta < 0)
	    else{
	      VNRawX_1[icent][EPbin][iqn]   += sumwqx->GetBinContent(ipt+1,ieta+1);
	      VNRawY_1[icent][EPbin][iqn]   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	      sumw_1[icent][EPbin][iqn]     += sumw->GetBinContent(ipt+1,ieta+1);
	      evtMult_1[icent][EPbin][iqn]  += hMult->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Full Event
	    VNRawX_full[icent][EPbin][iqn]  += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VNRawY_full[icent][EPbin][iqn]  += sumwqy->GetBinContent(ipt+1,ieta+1);
	    
	    sumw_full[icent][EPbin][iqn]    += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_full[icent][EPbin][iqn] += hMult->GetBinContent(ipt+1,ieta+1);
	  }
	  
	} //-- End eta loop
	
      } //-- End pt loop
      
      //-- Only use events that have tracks in all subevents AND subevents that have at least two tracks to determine VN
      if( sumw_0[icent][EPbin][iqn] == 0 || sumw_1[icent][EPbin][iqn] == 0 || evtMult_0[icent][EPbin][iqn] < 2 || evtMult_1[icent][EPbin][iqn] < 2 ){
	NFails[icent][EPbin][iqn]++;
      }
      else{      
	VNDetX_0[icent][EPbin][iqn]    += VNRawX_0[icent][EPbin][iqn]    / sumw_0[icent][EPbin][iqn];
	VNDetY_0[icent][EPbin][iqn]    += VNRawY_0[icent][EPbin][iqn]    / sumw_0[icent][EPbin][iqn];
	VNDetX_1[icent][EPbin][iqn]    += VNRawX_1[icent][EPbin][iqn]    / sumw_1[icent][EPbin][iqn];
	VNDetY_1[icent][EPbin][iqn]    += VNRawY_1[icent][EPbin][iqn]    / sumw_1[icent][EPbin][iqn];
	VNDetX_full[icent][EPbin][iqn] += VNRawX_full[icent][EPbin][iqn] / sumw_full[icent][EPbin][iqn];
	VNDetY_full[icent][EPbin][iqn] += VNRawY_full[icent][EPbin][iqn] / sumw_full[icent][EPbin][iqn];
      }

    } //-- End EP loop

  } //-- End event loop
    
  std::cout<<"End DETECTOR loop"<<std::endl;

  //-- Average VNDet over all events for each centrality, EP and, qn bin
  for(int iqn = 0; iqn < NQN; iqn++){
    for(int icent = 0; icent < NCENT; icent++){
      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;

	double  Neffective = (double) Nevents[icent][iEP][iqn] - (double) NFails[icent][iEP][iqn];
	if( Neffective == 0) continue;

	VNDetX_0[icent][iEP][iqn]    /= Neffective;
	VNDetY_0[icent][iEP][iqn]    /= Neffective;
	VNDetX_1[icent][iEP][iqn]    /= Neffective;
	VNDetY_1[icent][iEP][iqn]    /= Neffective;
	VNDetX_full[icent][iEP][iqn] /= Neffective;
	VNDetY_full[icent][iEP][iqn] /= Neffective;

	//-- Populate histograms that will be used by ReadTree_normDet.C
	hVNDetX_0[iqn]    -> SetBinContent(icent+1, iEP+1, VNDetX_0[icent][iEP][iqn]);
	hVNDetX_1[iqn]    -> SetBinContent(icent+1, iEP+1, VNDetX_1[icent][iEP][iqn]);
	hVNDetX_full[iqn] -> SetBinContent(icent+1, iEP+1, VNDetX_full[icent][iEP][iqn]);
        
	hVNDetY_0[iqn]    -> SetBinContent(icent+1, iEP+1, VNDetY_0[icent][iEP][iqn]);
	hVNDetY_1[iqn]    -> SetBinContent(icent+1, iEP+1, VNDetY_1[icent][iEP][iqn]);
	hVNDetY_full[iqn] -> SetBinContent(icent+1, iEP+1, VNDetY_full[icent][iEP][iqn]);

      }
    }
  } //-- END TRIPLE LOOP

  fOut->Write();
  cout<<"File written, process completed"<<endl;

}

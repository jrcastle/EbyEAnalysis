#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void ReadTree_normDet(){

  bool testrun          = 0;
  const int norder_     = 4;
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

  TFile * fVNDet;
  TH2D * hVNDetX_0[NQN];
  TH2D * hVNDetY_0[NQN];
  TH2D * hVNDetX_1[NQN];
  TH2D * hVNDetY_1[NQN];
  TH2D * hVNDetX_full[NQN];
  TH2D * hVNDetY_full[NQN];

  TFile * fHists;
  TDirectory * qwebye;
  TH2D * hVn2Dfull[NCENT][NEPSymm][NQN];
  TH2D * hVn2Dsub0[NCENT][NEPSymm][NQN];
  TH2D * hVn2Dsub1[NCENT][NEPSymm][NQN];
  TH2D * hVn2D0v1[NCENT][NEPSymm][NQN];
  TH1D * hVnFull[NCENT][NEPSymm][NQN];
  TH1D * hVnSub0[NCENT][NEPSymm][NQN];
  TH1D * hVnSub1[NCENT][NEPSymm][NQN];
  TH1I * Mult[NCENT][NEPSymm][NQN];
  TH2D * h2Vn2D0v1[NCENT][NEPSymm][NQN];
  TH1D * h2Vn2D0v1Magnitude[NCENT][NEPSymm][NQN];

  double VnRaw_x_0;
  double VnRaw_y_0;
  double VnRaw_x_1;
  double VnRaw_y_1;
  double VnRaw_x_full;
  double VnRaw_y_full;

  double sumw_0;
  double sumw_1;
  double sumw_full;

  double VnCorrected_x_0;
  double VnCorrected_y_0;
  double VnCorrected_x_1;
  double VnCorrected_y_1;
  double VnCorrected_x_full;
  double VnCorrected_y_full;

  int evtMult_0;
  int evtMult_1;
  int evtMult_full;

  //
  // MAIN
  //
  setTDRStyle();
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  TH1I::SetDefaultSumw2();

  //-- Set up analyzer objects
  fAna = new TFile(fAnaTreeName);

  tree   = (TTree *) fAna->Get("ebyeana/tree");
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
  fQNDet = new TFile( Form("../../../../../../v%i/eta2.4/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/Q%iDet.root", QnBinOrder_, QnBinOrder_) );
  for(int iEP = 0; iEP < NumEPNames; iEP++){
    int EPbin  = EPSymmPartnerBin[iEP];
    if( EPbin != EPSymmBin ) continue;
    hqnHFDet_x[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_x_%s", EPNames[iEP].data()) );
    hqnHFDet_y[iEP] = (TH1D*) fQNDet->Get( Form("hqnHFDet_y_%s", EPNames[iEP].data()) );
  }

  //-- Setup the QN binning objects
  fQN = new TFile( Form( "../../../../../../v%i/eta2.4/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/q%iCuts.root", QnBinOrder_, QnBinOrder_) );
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      hqbins[icent][iEP] = (TH1D*) fQN->Get( Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent) );
    }
  }

  //-- Set up VN detector objects 
  fVNDet = new TFile( Form("V%iDet.root", norder_ ) );
  for(int iqn = 0; iqn < NQN; iqn++){
    hVNDetX_0[iqn]    = (TH2D*) fVNDet->Get( Form("SubEvt_0/hVNDetX_0_qbin%i",   iqn) );
    hVNDetY_0[iqn]    = (TH2D*) fVNDet->Get( Form("SubEvt_0/hVNDetY_0_qbin%i",   iqn) );
    hVNDetX_1[iqn]    = (TH2D*) fVNDet->Get( Form("SubEvt_1/hVNDetX_1_qbin%i",   iqn) );
    hVNDetY_1[iqn]    = (TH2D*) fVNDet->Get( Form("SubEvt_1/hVNDetY_1_qbin%i",   iqn) );
    hVNDetX_full[iqn] = (TH2D*) fVNDet->Get( Form("FullEvt/hVNDetX_full_qbin%i", iqn) );
    hVNDetY_full[iqn] = (TH2D*) fVNDet->Get( Form("FullEvt/hVNDetY_full_qbin%i", iqn) );
  }

  //-- Setup the output objects
  fHists = new TFile("CastleEbyE.root","recreate");
  qwebye = (TDirectory*) fHists->mkdir("qwebye");

  for(int iqn = 0; iqn < NQN; iqn++){
    for(int icent = 0; icent < NCENT; icent++){
      for(int iEP = 0; iEP < NEPSymm; iEP++){    
	if( iEP != EPSymmBin ) continue;

	qwebye->cd();
	hVn2Dfull[icent][iEP][iqn]          = new TH2D(Form("hVn2Dfull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVn2Dfull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
	hVn2Dfull[icent][iEP][iqn]->SetOption("colz");
	hVn2Dfull[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs}", norder_) );
	hVn2Dfull[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs}", norder_) );

	hVn2Dsub0[icent][iEP][iqn]          = new TH2D( Form("hVn2Dsub0_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVn2Dsub0_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
	hVn2Dsub0[icent][iEP][iqn]->SetOption("colz");
	hVn2Dsub0[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,a}", norder_) );
	hVn2Dsub0[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,a}", norder_) );

	hVn2Dsub1[icent][iEP][iqn]          = new TH2D( Form("hVn2Dsub1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVn2Dsub1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
	hVn2Dsub1[icent][iEP][iqn]->SetOption("colz");
	hVn2Dsub1[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,b}", norder_) );
	hVn2Dsub1[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,b}", norder_) );

	hVn2D0v1[icent][iEP][iqn]           = new TH2D( Form("hVn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
	hVn2D0v1[icent][iEP][iqn]->SetOption("colz");
	hVn2D0v1[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b}", norder_, norder_) );
	hVn2D0v1[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,a} - v_{%i,y}^{obs,b}", norder_, norder_) );

	hVnFull[icent][iEP][iqn]            = new TH1D( Form("hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_] );
	hVnFull[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
	hVnFull[icent][iEP][iqn]->GetYaxis()->SetTitle( "Events" );

	hVnSub0[icent][iEP][iqn]            = new TH1D( Form("hVnSub0_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVnSub0_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_] );
	hVnSub0[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}^{obs,a}", norder_) );

	hVnSub1[icent][iEP][iqn]            = new TH1D( Form("hVnSub1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("hVnSub1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_] );
	hVnSub1[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}^{obs,b}", norder_) );

	Mult[icent][iEP][iqn]               = new TH1I( Form("Mult_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("Mult_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 300, 1, 1500 );
	Mult[icent][iEP][iqn]->GetXaxis()->SetTitle("Multiplicity");

	h2Vn2D0v1[icent][iEP][iqn]          = new TH2D( Form("h2Vn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("h2Vn2D0v1_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
	h2Vn2D0v1[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("(v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b})/2", norder_, norder_) );
	h2Vn2D0v1[icent][iEP][iqn]->GetYaxis()->SetTitle( Form("(v_{%i,y}^{obs,a} - v_{%i,y}^{obs,b})/2", norder_, norder_) );
	h2Vn2D0v1[icent][iEP][iqn]->SetOption("colz");

	h2Vn2D0v1Magnitude[icent][iEP][iqn] = new TH1D( Form("h2Vn2D0v1Magnitude_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), Form("h2Vn2D0v1Magnitude_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn), NBins, 0., vnMax[norder_] );
	h2Vn2D0v1Magnitude[icent][iEP][iqn]->GetXaxis()->SetTitle( Form("|(v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b})/2|", norder_, norder_) );

      }
    } 
  }
    
  //
  // Tree Loop
  //
    
  cout<<"Begin FINAL loop, contains "<<tree->GetEntries()<<" Events"<<endl;
    
  int N;
  if( testrun ) N = 10000;
  else          N = tree->GetEntries();

  for(int ievent = 0; ievent < N; ievent++) {
        
    if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<100.*(ievent+1)/(double)N<<"% Completed"<<endl;
      
    tree->GetEntry(ievent);

    //-- Vertex Cut
    if(TMath::Abs(vtx) < 3. || TMath::Abs(vtx) > 15.) continue;

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

      //-- Reset raw and sumw values
      VnRaw_x_0    = 0;
      VnRaw_y_0    = 0;
      VnRaw_x_1    = 0;
      VnRaw_y_1    = 0;
      VnRaw_x_full = 0;
      VnRaw_y_full = 0;                

      sumw_0       = 0;
      sumw_1       = 0;
      sumw_full    = 0;

      evtMult_0    = 0;
      evtMult_1    = 0;
      evtMult_full = 0;

      //-- Begin analyzer histogram loops
      for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++){
	for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++){

	  if(sumw->GetBinContent(ipt+1,ieta+1) !=0){

	    //-- Subevent 0 (eta >= 0)
	    if(etabinsDefault[ieta] >= 0){
	      VnRaw_x_0     += sumwqx->GetBinContent(ipt+1,ieta+1);
	      VnRaw_y_0     += sumwqy->GetBinContent(ipt+1,ieta+1);
	      sumw_0        += sumw->GetBinContent(ipt+1,ieta+1);
	      evtMult_0     += hMult->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Subevent 1 (eta < 0)
	    else{
	      VnRaw_x_1     += sumwqx->GetBinContent(ipt+1,ieta+1);
	      VnRaw_y_1     += sumwqy->GetBinContent(ipt+1,ieta+1);
	      sumw_1        += sumw->GetBinContent(ipt+1,ieta+1);
	      evtMult_1     += hMult->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Full Event
	    VnRaw_x_full    += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VnRaw_y_full    += sumwqy->GetBinContent(ipt+1,ieta+1);
	    sumw_full       += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_full    += hMult->GetBinContent(ipt+1,ieta+1);

	  }

	} //-- End eta loop
      } //-- End pt loop

      //-- Fill Histograms, only use events that have at least two tracks in each SE
      if(sumw_0 == 0 || sumw_1 == 0 || evtMult_1 < 2 || evtMult_full < 2) continue;

      VnRaw_x_full /= sumw_full;
      VnRaw_y_full /= sumw_full;
	
      VnRaw_x_0 /= sumw_0;
      VnRaw_y_0 /= sumw_0;
	
      VnRaw_x_1 /= sumw_1;
      VnRaw_y_1 /= sumw_1;
                
      //-- Full Tracker
      VnCorrected_x_full = VnRaw_x_full - hVNDetX_full[iqn]->GetBinContent(icent+1, EPbin+1);
      VnCorrected_y_full = VnRaw_y_full - hVNDetY_full[iqn]->GetBinContent(icent+1, EPbin+1);
      double vn_full     = TMath::Sqrt( VnCorrected_x_full * VnCorrected_x_full +  VnCorrected_y_full * VnCorrected_y_full);

      hVnFull[icent][EPbin][iqn]   -> Fill( vn_full );
      hVn2Dfull[icent][EPbin][iqn] -> Fill(VnCorrected_x_full, VnCorrected_y_full);

      Mult[icent][EPbin][iqn]->Fill(evtMult_full);

      //-- SubEvt 0 (Eta > 0)
      VnCorrected_x_0 = VnRaw_x_0 - hVNDetX_0[iqn]->GetBinContent(icent+1, EPbin+1);
      VnCorrected_y_0 = VnRaw_y_0 - hVNDetY_0[iqn]->GetBinContent(icent+1, EPbin+1);
      double vn_0     = TMath::Sqrt( VnCorrected_x_0 * VnCorrected_x_0 + VnCorrected_y_0 * VnCorrected_y_0 );

      hVnSub0[icent][EPbin][iqn]   -> Fill( vn_0 );
      hVn2Dsub0[icent][EPbin][iqn] -> Fill(VnCorrected_x_0, VnCorrected_y_0);
                
      //-- SubEvt 1 (Eta < 0)
      VnCorrected_x_1 = VnRaw_x_1 - hVNDetX_1[iqn]->GetBinContent(icent+1, EPbin+1);
      VnCorrected_y_1 = VnRaw_y_1 - hVNDetY_1[iqn]->GetBinContent(icent+1, EPbin+1);
      double vn_1     = TMath::Sqrt( VnCorrected_x_1 * VnCorrected_x_1 + VnCorrected_y_1 *VnCorrected_y_1 );

      hVnSub1[icent][EPbin][iqn]   -> Fill( vn_1 );
      hVn2Dsub1[icent][EPbin][iqn] -> Fill(VnCorrected_x_1, VnCorrected_y_1);

      //-- SubEvt Difference
      double vn0m1_x = VnCorrected_x_0 - VnCorrected_x_1;
      double vn0m1_y = VnCorrected_y_0 - VnCorrected_y_1;
      hVn2D0v1[icent][EPbin][iqn]->Fill(vn0m1_x, vn0m1_y);	

      //-- SubEvt Difference for DD response
      double vn0m1_x2 = (VnCorrected_x_0 - VnCorrected_x_1) / 2.;
      double vn0m1_y2 = (VnCorrected_y_0 - VnCorrected_y_1) / 2.;
      h2Vn2D0v1[icent][EPbin][iqn]->Fill(vn0m1_x2, vn0m1_y2);

      double vn0m12 = TMath::Sqrt( pow(vn0m1_x2, 2) + pow(vn0m1_y2, 2) );
      h2Vn2D0v1Magnitude[icent][EPbin][iqn]->Fill(vn0m12);

    } //-- End EP loop        

  } //-- End Event loop

  cout<<"End FINAL loop!"<<endl;
        
  fHists->Write();
  cout<<"File written, process completed"<<endl;
    
}
    
    
    
    
    
    
    
    
    
    
    

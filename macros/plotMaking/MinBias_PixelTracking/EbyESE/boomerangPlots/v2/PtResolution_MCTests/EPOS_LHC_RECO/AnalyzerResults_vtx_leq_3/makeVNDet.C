#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
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

  TFile * fOut;
  TDirectory * SubEvt_0;
  TDirectory * SubEvt_1;
  TDirectory * FullEvt;

  TH1D * hVNDetX_0;
  TH1D * hVNDetY_0;
  TH1D * hVNDetX_1;
  TH1D * hVNDetY_1;
  TH1D * hVNDetX_full;
  TH1D * hVNDetY_full;

  double VNRawX_0[NCENT];
  double VNRawY_0[NCENT];
  double VNRawX_1[NCENT];
  double VNRawY_1[NCENT];
  double VNRawX_full[NCENT];
  double VNRawY_full[NCENT];

  double sumw_0[NCENT];
  double sumw_1[NCENT];
  double sumw_full[NCENT];

  int evtMult_0[NCENT];
  int evtMult_1[NCENT];
  int evtMult_full[NCENT];
 
  double VNDetX_0[NCENT];
  double VNDetY_0[NCENT];
  double VNDetX_1[NCENT];
  double VNDetY_1[NCENT];
  double VNDetX_full[NCENT];
  double VNDetY_full[NCENT];

  int Nevents[NCENT];
  int NFails[NCENT];

  //
  // MAIN
  //


  //-- Set up the analyzer objects
  fAna = new TFile("/rfs/jcastle/PbPb2015/PixelTracking_MB2/EbyETree_EPOS_LHC_RECO.root");

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

  //-- Setup the output objects
  fOut     = new TFile( Form("V%iDet.root", norder_), "recreate" );
  SubEvt_0 = fOut->mkdir("SubEvt_0");
  SubEvt_1 = fOut->mkdir("SubEvt_1");
  FullEvt  = fOut->mkdir("FullEvt");
    
  SubEvt_0->cd();
  hVNDetX_0    = new TH1D( "hVNDetX_0",    "hVNDetX_0",    NCENT, centbinsDefault);
  hVNDetX_0->GetXaxis()->SetTitle("Centrality %");
    
  hVNDetY_0    = new TH1D( "hVNDetY_0",    "hVNDetY_0",    NCENT, centbinsDefault);
  hVNDetY_0->GetXaxis()->SetTitle("Centrality %");
    
  SubEvt_1->cd();
  hVNDetX_1    = new TH1D( "hVNDetX_1",    "hVNDetX_1",    NCENT, centbinsDefault);
  hVNDetX_1->GetXaxis()->SetTitle("Centrality %");
    
  hVNDetY_1    = new TH1D( "hVNDetY_1",    "hVNDetY_1",    NCENT, centbinsDefault);
  hVNDetY_1->GetXaxis()->SetTitle("Centrality %");
  
  FullEvt->cd();
  hVNDetX_full = new TH1D( "hVNDetX_full", "hVNDetX_full", NCENT, centbinsDefault);
  hVNDetX_full->GetXaxis()->SetTitle("Centrality %");
    
  hVNDetY_full = new TH1D( "hVNDetY_full", "hVNDetY_full", NCENT, centbinsDefault);
  hVNDetY_full->GetXaxis()->SetTitle("Centrality %");
    
  //-- initialize all variables
  for(int icent = 0; icent<NCENT; icent++){

    VNDetX_0[icent]     = 0.;
    VNDetY_0[icent]     = 0.;
    VNDetX_1[icent]     = 0.;
    VNDetY_1[icent]     = 0.;
    VNDetX_full[icent]  = 0.;
    VNDetY_full[icent]  = 0.;

    evtMult_0[icent]    = 0;
    evtMult_1[icent]    = 0;
    evtMult_full[icent] = 0;

    Nevents[icent]      = 0;
    NFails[icent]       = 0;

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
    if(TMath::Abs(vtx) > 3) continue;
      
    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    //int icent = hCentBins.FindBin(centval)-1;
    int icent = 0;
    //-- Reset Raw and sumw values
    VNRawX_0[icent]     = 0.;
    VNRawY_0[icent]     = 0.;
    VNRawX_1[icent]     = 0.;
    VNRawY_1[icent]     = 0.;
    VNRawX_full[icent]  = 0.;
    VNRawY_full[icent]  = 0.;            
	
    sumw_0[icent]       = 0.;
    sumw_1[icent]       = 0.;
    sumw_full[icent]    = 0.;
      
    evtMult_0[icent]    = 0;
    evtMult_1[icent]    = 0;
    evtMult_full[icent] = 0;
      
    Nevents[icent]++;

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
	    VNRawX_0[icent]   += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VNRawY_0[icent]   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	    sumw_0[icent]     += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_0[icent]  += hMult->GetBinContent(ipt+1,ieta+1);
	  }
	  //-- Subevent 1 (eta < 0)
	  else{
	    VNRawX_1[icent]   += sumwqx->GetBinContent(ipt+1,ieta+1);
	    VNRawY_1[icent]   += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	    sumw_1[icent]     += sumw->GetBinContent(ipt+1,ieta+1);
	    evtMult_1[icent]  += hMult->GetBinContent(ipt+1,ieta+1);
	  }
	  //-- Full Event
	  VNRawX_full[icent]  += sumwqx->GetBinContent(ipt+1,ieta+1);
	  VNRawY_full[icent]  += sumwqy->GetBinContent(ipt+1,ieta+1);
	    
	  sumw_full[icent]    += sumw->GetBinContent(ipt+1,ieta+1);
	  evtMult_full[icent] += hMult->GetBinContent(ipt+1,ieta+1);
	}
	  
      } //-- End eta loop
      
    } //-- End pt loop
    
    //-- Only use events that have tracks in all subevents AND subevents that have at least two tracks to determine VN
    if( sumw_0[icent] == 0 || sumw_1[icent] == 0 || evtMult_0[icent] < 2 || evtMult_1[icent] < 2 ){
	NFails[icent]++;
    }
    else{      
      VNDetX_0[icent]    += VNRawX_0[icent]    / sumw_0[icent];
      VNDetY_0[icent]    += VNRawY_0[icent]    / sumw_0[icent];
      VNDetX_1[icent]    += VNRawX_1[icent]    / sumw_1[icent];
      VNDetY_1[icent]    += VNRawY_1[icent]    / sumw_1[icent];
      VNDetX_full[icent] += VNRawX_full[icent] / sumw_full[icent];
      VNDetY_full[icent] += VNRawY_full[icent] / sumw_full[icent];
    }

  } //-- End event loop
    
  std::cout<<"End DETECTOR loop"<<std::endl;

  //-- Average VNDet over all events for each centrality, EP and, qn bin
  for(int icent = 0; icent < NCENT; icent++){
  
    double  Neffective = (double) Nevents[icent] - (double) NFails[icent];
    if( Neffective == 0 ) continue;
    
    VNDetX_0[icent]    /= Neffective;
    VNDetY_0[icent]    /= Neffective;
    VNDetX_1[icent]    /= Neffective;
    VNDetY_1[icent]    /= Neffective;
    VNDetX_full[icent] /= Neffective;
    VNDetY_full[icent] /= Neffective;

    //-- Populate histograms that will be used by ReadTree_normDet.C
    hVNDetX_0    -> SetBinContent(icent+1, VNDetX_0[icent]);
    hVNDetX_1    -> SetBinContent(icent+1, VNDetX_1[icent]);
    hVNDetX_full -> SetBinContent(icent+1, VNDetX_full[icent]);
        
    hVNDetY_0    -> SetBinContent(icent+1, VNDetY_0[icent]);
    hVNDetY_1    -> SetBinContent(icent+1, VNDetY_1[icent]);
    hVNDetY_full -> SetBinContent(icent+1, VNDetY_full[icent]);

  }

  fOut->Write();
  cout<<"File written, process completed"<<endl;

}

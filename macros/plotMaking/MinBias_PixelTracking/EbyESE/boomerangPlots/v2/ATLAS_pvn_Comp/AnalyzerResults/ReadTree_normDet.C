#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void ReadTree_normDet(){

  bool testrun          = 0;
  const int norder_     = 2;
  const double vtxCut   = 15.;

  const double ptMin = 0.5;
  const double ptMax = 20.0;
  const double etaMax = 2.4;

  const int ptBinMin  = 0;
  const int ptBinMax  = nptbinsDefault-1;
  const int etaBinMin = 0; //0;
  const int etaBinMax = netabinsDefault-1;

  const double PTBINSDEFAULT[] = {
    0.50,    1.00,    20.0
  };

  TFile * fAna;
  TTree * tree;
  double centval;
  double vtx;
  TH2D * sumw;
  TH2D * sumwqx;
  TH2D * sumwqy;
  TH2I * hMult;

  TFile * fVNDet;
  TH1D * hVNDetX_0;
  TH1D * hVNDetY_0;
  TH1D * hVNDetX_1;
  TH1D * hVNDetY_1;
  TH1D * hVNDetX_full;
  TH1D * hVNDetY_full;

  TFile * fHists;
  TDirectory * qwebye;
  TH2D * hVn2Dfull[NCENT];
  TH2D * hVn2Dsub0[NCENT];
  TH2D * hVn2Dsub1[NCENT];
  TH2D * hVn2D0v1[NCENT];
  TH1D * hVnFull[NCENT];
  TH1D * hVnSub0[NCENT];
  TH1D * hVnSub1[NCENT];
  TH1I * Mult[NCENT];
  TH2D * h2Vn2D0v1[NCENT];
  TH1D * h2Vn2D0v1Magnitude[NCENT];

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
  fAna = new TFile( "/rfs/jcastle/PbPb2015/PixelTracking_MB2/EbyETree_ATLASComp.root" );

  tree   = (TTree *) fAna->Get("ebyeana/tree");
  sumwqx = new TH2D(Form("sumwqx%i", norder_), Form("sumwqx%i", norder_), nptbinsDefault, PTBINSDEFAULT, netabinsDefault, etabinsDefault);
  sumwqy = new TH2D(Form("sumwqy%i", norder_), Form("sumwqy%i", norder_), nptbinsDefault, PTBINSDEFAULT, netabinsDefault, etabinsDefault);
  sumw   = new TH2D("sumw",                    "sumw",                    nptbinsDefault, PTBINSDEFAULT, netabinsDefault, etabinsDefault);
  hMult  = new TH2I("hMult",                   "hMult",                   nptbinsDefault, PTBINSDEFAULT, netabinsDefault, etabinsDefault);

  tree->SetBranchAddress("Cent",                    &centval);
  tree->SetBranchAddress("Vtx",                     &vtx);
  tree->SetBranchAddress("mult",                    &hMult);
  tree->SetBranchAddress(Form("sumwqx%i", norder_), &sumwqx);
  tree->SetBranchAddress(Form("sumwqy%i", norder_), &sumwqy);
  tree->SetBranchAddress("sumw",                    &sumw);

  //-- Set up VN detector objects 
  fVNDet = new TFile( Form("V%iDet.root", norder_ ) );
  hVNDetX_0    = (TH1D*) fVNDet->Get("SubEvt_0/hVNDetX_0");
  hVNDetY_0    = (TH1D*) fVNDet->Get("SubEvt_0/hVNDetY_0");
  hVNDetX_1    = (TH1D*) fVNDet->Get("SubEvt_1/hVNDetX_1");
  hVNDetY_1    = (TH1D*) fVNDet->Get("SubEvt_1/hVNDetY_1");
  hVNDetX_full = (TH1D*) fVNDet->Get("FullEvt/hVNDetX_full");
  hVNDetY_full = (TH1D*) fVNDet->Get("FullEvt/hVNDetY_full");

  //-- Setup the output objects
  fHists = new TFile("CastleEbyE.root","recreate");
  qwebye = (TDirectory*) fHists->mkdir("qwebye");

  for(int icent = 0; icent < NCENT; icent++){

    qwebye->cd();
    hVn2Dfull[icent]          = new TH2D(Form("hVn2Dfull_c%i", icent), Form("hVn2Dfull_c%i", icent), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
    hVn2Dfull[icent]->SetOption("colz");
    hVn2Dfull[icent]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs}", norder_) );
    hVn2Dfull[icent]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs}", norder_) );

    hVn2Dsub0[icent]          = new TH2D( Form("hVn2Dsub0_c%i", icent), Form("hVn2Dsub0_c%i", icent), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
    hVn2Dsub0[icent]->SetOption("colz");
    hVn2Dsub0[icent]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,a}", norder_) );
    hVn2Dsub0[icent]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,a}", norder_) );

    hVn2Dsub1[icent]          = new TH2D( Form("hVn2Dsub1_c%i", icent), Form("hVn2Dsub1_c%i", icent), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
    hVn2Dsub1[icent]->SetOption("colz");
    hVn2Dsub1[icent]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,b}", norder_) );
    hVn2Dsub1[icent]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,b}", norder_) );

    hVn2D0v1[icent]           = new TH2D( Form("hVn2D0v1_c%i", icent), Form("hVn2D0v1_c%i", icent), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
    hVn2D0v1[icent]->SetOption("colz");
    hVn2D0v1[icent]->GetXaxis()->SetTitle( Form("v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b}", norder_, norder_) );
    hVn2D0v1[icent]->GetYaxis()->SetTitle( Form("v_{%i,y}^{obs,a} - v_{%i,y}^{obs,b}", norder_, norder_) );

    hVnFull[icent]            = new TH1D( Form("hVnFull_c%i", icent), Form("hVnFull_c%i", icent), NBins, 0., vnMax[norder_] );
    hVnFull[icent]->GetXaxis()->SetTitle( Form("v_{%i}", norder_) );
    hVnFull[icent]->GetYaxis()->SetTitle( "Events" );

    hVnSub0[icent]            = new TH1D( Form("hVnSub0_c%i", icent), Form("hVnSub0_c%i", icent), NBins, 0., vnMax[norder_] );
    hVnSub0[icent]->GetXaxis()->SetTitle( Form("v_{%i}^{obs,a}", norder_) );

    hVnSub1[icent]            = new TH1D( Form("hVnSub1_c%i", icent), Form("hVnSub1_c%i", icent), NBins, 0., vnMax[norder_] );
    hVnSub1[icent]->GetXaxis()->SetTitle( Form("v_{%i}^{obs,b}", norder_) );
    
    Mult[icent]               = new TH1I( Form("Mult_c%i", icent), Form("Mult_c%i", icent), 250, 1, 5000 );
    Mult[icent]->GetXaxis()->SetTitle("Multiplicity");
	
    h2Vn2D0v1[icent]          = new TH2D( Form("h2Vn2D0v1_c%i", icent), Form("h2Vn2D0v1_c%i", icent), 2*NBins, -vnMax[norder_], vnMax[norder_], 2*NBins, -vnMax[norder_], vnMax[norder_] );
    h2Vn2D0v1[icent]->GetXaxis()->SetTitle( Form("(v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b})/2", norder_, norder_) );
    h2Vn2D0v1[icent]->GetYaxis()->SetTitle( Form("(v_{%i,y}^{obs,a} - v_{%i,y}^{obs,b})/2", norder_, norder_) );
    h2Vn2D0v1[icent]->SetOption("colz");

    h2Vn2D0v1Magnitude[icent] = new TH1D( Form("h2Vn2D0v1Magnitude_c%i", icent), Form("h2Vn2D0v1Magnitude_c%i", icent), NBins, 0., vnMax[norder_] );
    h2Vn2D0v1Magnitude[icent]->GetXaxis()->SetTitle( Form("|(v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b})/2|", norder_, norder_) );

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
    if(TMath::Abs(vtx) > vtxCut) continue;

    //-- Calculate centbin
    if( centval > cent_max[NCENT-1]) continue;
    int icent = hCentBins.FindBin(centval)-1;

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

	double pt  = sumw->GetXaxis()->GetBinCenter(ipt+1);
        double eta = sumw->GetYaxis()->GetBinCenter(ieta+1);
        if( pt < ptMin || pt > ptMax ) continue;
        if( fabs( eta ) > etaMax )     continue;

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
    VnCorrected_x_full = VnRaw_x_full - hVNDetX_full->GetBinContent(icent+1);
    VnCorrected_y_full = VnRaw_y_full - hVNDetY_full->GetBinContent(icent+1);
    double vn_full     = TMath::Sqrt( VnCorrected_x_full * VnCorrected_x_full +  VnCorrected_y_full * VnCorrected_y_full);

    hVnFull[icent]   -> Fill( vn_full );
    hVn2Dfull[icent] -> Fill(VnCorrected_x_full, VnCorrected_y_full);

    Mult[icent]->Fill(evtMult_full);

    //-- SubEvt 0 (Eta > 0)
    VnCorrected_x_0 = VnRaw_x_0 - hVNDetX_0->GetBinContent(icent+1);
    VnCorrected_y_0 = VnRaw_y_0 - hVNDetY_0->GetBinContent(icent+1);
    double vn_0     = TMath::Sqrt( VnCorrected_x_0 * VnCorrected_x_0 + VnCorrected_y_0 * VnCorrected_y_0 );

    hVnSub0[icent]   -> Fill( vn_0 );
    hVn2Dsub0[icent] -> Fill(VnCorrected_x_0, VnCorrected_y_0);
                
    //-- SubEvt 1 (Eta < 0)
    VnCorrected_x_1 = VnRaw_x_1 - hVNDetX_1->GetBinContent(icent+1);
    VnCorrected_y_1 = VnRaw_y_1 - hVNDetY_1->GetBinContent(icent+1);
    double vn_1     = TMath::Sqrt( VnCorrected_x_1 * VnCorrected_x_1 + VnCorrected_y_1 *VnCorrected_y_1 );

    hVnSub1[icent]   -> Fill( vn_1 );
    hVn2Dsub1[icent] -> Fill(VnCorrected_x_1, VnCorrected_y_1);

    //-- SubEvt Difference
    double vn0m1_x = VnCorrected_x_0 - VnCorrected_x_1;
    double vn0m1_y = VnCorrected_y_0 - VnCorrected_y_1;
    hVn2D0v1[icent]->Fill(vn0m1_x, vn0m1_y);	

    //-- SubEvt Difference for DD response
    double vn0m1_x2 = (VnCorrected_x_0 - VnCorrected_x_1) / 2.;
    double vn0m1_y2 = (VnCorrected_y_0 - VnCorrected_y_1) / 2.;
    h2Vn2D0v1[icent]->Fill(vn0m1_x2, vn0m1_y2);

    double vn0m12 = TMath::Sqrt( pow(vn0m1_x2, 2) + pow(vn0m1_y2, 2) );
    h2Vn2D0v1Magnitude[icent]->Fill(vn0m12);

  } //-- End Event loop

  cout<<"End FINAL loop!"<<endl;
        
  fHists->Write();
  cout<<"File written, process completed"<<endl;
   
}
    
    
    
    
    
    
    
    
    
    
    

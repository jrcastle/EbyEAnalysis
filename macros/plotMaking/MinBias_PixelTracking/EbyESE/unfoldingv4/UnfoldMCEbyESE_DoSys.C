#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/differentialBinning.h"
#include "/home/j550c590/tdrstyle.C"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfold.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>

//int centbin = 0;  // 0 - 10 %
//int centbin = 1;  // 10 - 20 %
//int centbin = 2;  // 20 - 30 %
//int centbin = 3;  // 30 - 40 %
//int centbin = 4;  // 40 - 50 %
//int centbin = 5;  // 50 - 60 %
//int centbin = 6;  // 60 - 70 %
//int centbin = 7;  // 70 - 80 %
//int centbin = 8;  // 80 - 90 %
//int centbin = 9;  // 90 - 100 %

//
//-- MAIN
//

void UnfoldMCEbyE_DoSys( int ipt ){

  int ptBin = ipt;

  bool Resp_Data     = 1;
  bool StudTresp     = 0;
  bool GaussResp     = 0;

  bool overshoot20   = 0;
  bool undershoot20  = 0;
  bool weighRespObs  = 0;

  bool smoothIter    = 0;
  bool smoothResp    = 1;
  bool isMC          = 1;
  bool sw            = 1;

  double sigmaFit = 2.;
  int VN          = 2;

  static const int NITER  = 8;
  static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};

  RooUnfoldResponse * response[NCENT];

  TLatex latex;

  TFile * fsave;
  TFile * fData;
  TFile * fPrior;

  //-- Unfolding Objects
  TH2D * h2D[NCENT];
  TH2D * h2Dsub0[NCENT];
  TH2D * h2Dsub1[NCENT];
  TH2D * h2Dx[NCENT];

  TH1D * h1D[NCENT];
  TH1D * h1Dsub0[NCENT];
  TH1D * h1Dsub1[NCENT];
  TH1D * h1Dx[NCENT];
  TH1D * h1Dy[NCENT];
  TH1I * hMult[NCENT];
  TH2D * h22D01[NCENT];
  TH1D * h22D01Magnitude[NCENT];
  TH1D * hw[NCENT];
  TF1 * f1Dx[NCENT];
  TF1 * f1Dy[NCENT];
  TF1 * f1Dx_StudT[NCENT];
  TF1 * f1Dy_StudT[NCENT];
  TF1 * fRespStudT[NCENT];

  //-- Data objects for fitting analytic response functions
  double sigma_x[NCENT];
  double sigma_xe[NCENT];
  double mean_x[NCENT];
  double chi2_x[NCENT];
  double ndf_x[NCENT];
  double sigma_y[NCENT];
  double sigma_ye[NCENT];
  double mean_y[NCENT];
  double chi2_y[NCENT];
  double ndf_y[NCENT];
  double sigma_2SE[NCENT];
  double Vn_mean[NCENT];
  double Vn_rms[NCENT];

  TH2D * hresp[NCENT];
  TH1D * hreco[NCENT][NITER];
  TH1D * hrefold[NCENT][NITER];

  //-- Data Histos
  static const int NVn = 7;
  TH2D * hVn2Dfull[NVn][NCENT_PBPB];
  TH2D * hVn2Dsub0[NVn][NCENT_PBPB];
  TH2D * hVn2Dsub1[NVn][NCENT_PBPB];
  TH2D * hVn2D0v1[NVn][NCENT_PBPB];
  TH1D * hVnFull[NVn][NCENT_PBPB];
  TH1D * hVnSub0[NVn][NCENT_PBPB];
  TH1D * hVnSub1[NVn][NCENT_PBPB];
  TH1I * hMultRaw[NVn][NCENT_PBPB];

  TH2D * h2Vn2D0v1[NVn][NCENT_PBPB];
  TH1D * h2Vn2D0v1Magnitude[NVn][NCENT_PBPB];
  TH2D * hResp[NVn][NCENT_PBPB];


  setTDRStyle();
  latex.SetNDC();

  fsave = 0;
  if( isMC ) fsave = new TFile(Form("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/txt/PbPb_2015/MC/data%i.root", VN), "recreate");
  else       fsave = new TFile(Form("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/txt/PbPb_2015/data/data%i.root", VN), "recreate");

  bool R1 = Resp_Data && StudTresp && GaussResp;
  bool R2 = Resp_Data && StudTresp;
  bool R3 = Resp_Data && GaussResp;
  bool R4 = StudTresp && GaussResp;

  if( R1 || R2 || R3 || R4){
    std::cout<<"WARNING! More than one response function defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  //-- Get histos from datafile
  fData = 0;
  if( isMC ) fData = new TFile("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/data/PbPb_2015/MC/CastleEbyE.root");
  else       fData = new TFile("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/data/PbPb_2015/data/CastleEbyE.root");

  //-- If using the data-driven response function, grab the file that has them
  fPrior = 0;
  if( Resp_Data ) fPrior = new TFile("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/dataDrivenRespAndPriors.root");


  for(int ivn = 1; ivn <=4; ivn++){
    if(ivn != VN) continue;

    for(int icent = 0; icent < NCENT_PBPB; icent++){

      if(icent == NCENT40_cutoff) break;

      hVn2Dfull[ivn][icent]          = (TH2D*) fData->Get(Form("qwebye/hVn2Dfull_%i_%i",ivn,icent));
      hVn2Dsub0[ivn][icent]          = (TH2D*) fData->Get(Form("qwebye/hVn2Dsub0_%i_%i",ivn,icent));
      hVn2Dsub1[ivn][icent]          = (TH2D*) fData->Get(Form("qwebye/hVn2Dsub1_%i_%i",ivn,icent));
      hVn2D0v1[ivn][icent]           = (TH2D*) fData->Get(Form("qwebye/hVn2D0v1_%i_%i",ivn,icent));
      hVnFull[ivn][icent]            = (TH1D*) fData->Get(Form("qwebye/hVnFull_%i_%i",ivn,icent));
      hVnSub0[ivn][icent]            = (TH1D*) fData->Get(Form("qwebye/hVnSub0_%i_%i",ivn,icent));
      hVnSub1[ivn][icent]            = (TH1D*) fData->Get(Form("qwebye/hVnSub1_%i_%i",ivn,icent));
      hMultRaw[ivn][icent]           = (TH1I*) fData->Get(Form("qwebye/Mult_%i_%i",ivn,icent));
      h2Vn2D0v1[ivn][icent]          = (TH2D*) fData->Get(Form("qwebye/h2Vn2D0v1_%i_%i",ivn,icent));
      h2Vn2D0v1Magnitude[ivn][icent] = (TH1D*) fData->Get(Form("qwebye/h2Vn2D0v1Magnitude_%i_%i",ivn,icent));
      hResp[ivn][icent]              = (TH2D*) fData->Get(Form("qwebye/hResp_%i_%i",ivn,icent));;

    }

  }

  //-- Begin Centrality loop
  for(int c = 0; c < NCENT; c++){
    
    if(c > 5) break;
    std::cout<<"!! Processing Cent = "<<c<<std::endl;

    //-- Initiate unfolding histos and condense from 40 to 10 cent bins
    h2D[c] = (TH2D*)hVn2Dfull[VN][4*c]->Clone(Form("h2D_%i_%i", VN, c));
    h2D[c]->Add(hVn2Dfull[VN][4*c+1]);
    h2D[c]->Add(hVn2Dfull[VN][4*c+2]);
    h2D[c]->Add(hVn2Dfull[VN][4*c+3]);
    if ( h2D[c]->GetEntries() < 1000 ) continue;

    h2Dsub0[c] = (TH2D*)hVn2Dsub0[VN][4*c]->Clone(Form("h2Dsub0_%i_%i", VN, c));
    h2Dsub0[c]->Add(hVn2Dsub0[VN][4*c+1]);
    h2Dsub0[c]->Add(hVn2Dsub0[VN][4*c+2]);
    h2Dsub0[c]->Add(hVn2Dsub0[VN][4*c+3]);

    h2Dsub1[c] = (TH2D*)hVn2Dsub1[VN][4*c]->Clone(Form("h2Dsub1_%i_%i", VN, c));
    h2Dsub1[c]->Add(hVn2Dsub1[VN][4*c+1]);
    h2Dsub1[c]->Add(hVn2Dsub1[VN][4*c+2]);
    h2Dsub1[c]->Add(hVn2Dsub1[VN][4*c+3]);

    h2Dx[c] = (TH2D*)hVn2D0v1[VN][4*c]->Clone(Form("h2Dx_%i_%i", VN, c));
    h2Dx[c]->Add(hVn2D0v1[VN][4*c+1]);
    h2Dx[c]->Add(hVn2D0v1[VN][4*c+2]);
    h2Dx[c]->Add(hVn2D0v1[VN][4*c+3]);

    h1D[c] = (TH1D*)hVnFull[VN][4*c]->Clone(Form("h1D_%i_%i", VN, c));
    h1D[c]->Add(hVnFull[VN][4*c+1]);
    h1D[c]->Add(hVnFull[VN][4*c+2]);
    h1D[c]->Add(hVnFull[VN][4*c+3]);

    h1Dsub0[c] = (TH1D*)hVnSub0[VN][4*c]->Clone(Form("h1Dsub0_%i_%i", VN, c));
    h1Dsub0[c]->Add(hVnSub0[VN][4*c+1]);
    h1Dsub0[c]->Add(hVnSub0[VN][4*c+2]);
    h1Dsub0[c]->Add(hVnSub0[VN][4*c+3]);

    h1Dsub1[c] = (TH1D*)hVnSub1[VN][4*c]->Clone(Form("h1Dsub1_%i_%i", VN, c));
    h1Dsub1[c]->Add(hVnSub1[VN][4*c+1]);
    h1Dsub1[c]->Add(hVnSub1[VN][4*c+2]);
    h1Dsub1[c]->Add(hVnSub1[VN][4*c+3]);

    hMult[c] = (TH1I*) hMultRaw[VN][4*c]->Clone(Form("hMult_%i_%i", VN, c));
    hMult[c]->Add(hMultRaw[VN][4*c+1]);
    hMult[c]->Add(hMultRaw[VN][4*c+2]);
    hMult[c]->Add(hMultRaw[VN][4*c+3]);

    h22D01[c] = (TH2D*) h2Vn2D0v1[VN][4*c]->Clone(Form("h22D01_%i_%i", VN, c));
    h22D01[c]->Add(h2Vn2D0v1[VN][4*c+1]);
    h22D01[c]->Add(h2Vn2D0v1[VN][4*c+2]);
    h22D01[c]->Add(h2Vn2D0v1[VN][4*c+3]);

    h22D01Magnitude[c] = (TH1D*) h2Vn2D0v1Magnitude[VN][4*c]->Clone(Form("h22D01Magnitude_%i_%i", VN, c));
    h22D01Magnitude[c]->Add(h2Vn2D0v1Magnitude[VN][4*c+1]);
    h22D01Magnitude[c]->Add(h2Vn2D0v1Magnitude[VN][4*c+2]);
    h22D01Magnitude[c]->Add(h2Vn2D0v1Magnitude[VN][4*c+3]);

    if(Resp_Data){
      hresp[c] = 0;
      if( isMC ) hresp[c] = (TH2D*) fPrior->Get( Form("MC/hrespDD_MC_p%i_c%i", ptBin, c) );
      else       hresp[c] = (TH2D*) fPrior->Get( Form("DATA/hrespDD_DATA_p%i_c%i", ptBin, c) );
      hresp[c]->SetName( Form("hresp_%i_%i", VN, c) );

      if( weighRespObs ){
	TH1D * px = hresp[c]->ProjectionX();
	for(int x = 1; x <= NBins[ptBin][c]; x++){
	  for(int y = 1; y <= NBins[ptBin][c]; y++){

	    double projx = px->GetBinContent(x);
	    double truth = h1D[c]->GetBinContent(x);
	    double resp  = hresp[c]->GetBinContent(x,y);

	    double projxe = px->GetBinError(x);
            double truthe = h1D[c]->GetBinError(x);
            double respe  = hresp[c]->GetBinError(x,y);

	    double w;
	    double we;
	    if(projx == 0){
	      w  = 0;
	      we = 0;
	    }
	    else{
	      w = truth/projx;
	      we = TMath::Sqrt( pow(truthe*resp/projx, 2) + pow( projxe*truth*resp/projx/projx,2) + pow(respe*truth/projx, 2) );
	    }

	    hresp[c]->SetBinContent(x, y, w*resp);
	    hresp[c]->SetBinError(x, y, we);

	  }
	}
	delete px;
      }

    }

    h1Dx[c] = h2Dx[c]->ProjectionX(Form("h1Dx_%i_%i", VN, c));
    h1Dy[c] = h2Dx[c]->ProjectionY(Form("h1Dy_%i_%i", VN, c));

    f1Dx[c] = new TF1(Form("f1Dx_%i_%i", VN, c), "gaus", -v2Max[ptBin][c], v2Max[ptBin][c]);
    f1Dy[c] = new TF1(Form("f1Dy_%i_%i", VN, c), "gaus", -v2Max[ptBin][c], v2Max[ptBin][c]);

    int ndf = hMult[c]->GetMean() - 1;

    f1Dx_StudT[c] = new TF1(Form("f1Dx_StudT_%i_%i", VN, c), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -v2Max[ptBin][c],v2Max[ptBin][c]);
    f1Dx_StudT[c]->SetParName(0,"Norm");
    f1Dx_StudT[c]->SetParName(1,"Mean");
    f1Dx_StudT[c]->SetParName(2,"Sigma");
    f1Dx_StudT[c]->SetParName(3,"nu");
    f1Dx_StudT[c]->SetParameters(h1Dx[c]->GetMaximum(), 0, ((ndf-2.)/ndf)*h1Dx[c]->GetRMS(),ndf);
    f1Dx_StudT[c]->FixParameter(3,ndf);
    f1Dx_StudT[c]->SetLineColor(4);
    f1Dx_StudT[c]->SetLineWidth(2);

    f1Dy_StudT[c] = new TF1(Form("f1Dy_StudT_%i_%i", VN, c), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -v2Max[ptBin][c],v2Max[ptBin][c]);
    f1Dy_StudT[c]->SetParName(0,"Norm");
    f1Dy_StudT[c]->SetParName(1,"Mean");
    f1Dy_StudT[c]->SetParName(2,"Sigma");
    f1Dy_StudT[c]->SetParName(3,"nu");
    f1Dy_StudT[c]->SetParameters(h1Dy[c]->GetMaximum(), 0, ((ndf-2.)/ndf)*h1Dy[c]->GetRMS(),ndf);
    f1Dy_StudT[c]->FixParameter(3,ndf);
    f1Dy_StudT[c]->SetLineColor(4);
    f1Dy_StudT[c]->SetLineWidth(2);

    //-- [0] = vn_meas
    //-- [1] = vn_true
    //-- [2] = delta_vn (obtained from 2SE fits)
    //-- [3] = nu
    fRespStudT[c] = new TF1(Form("fRespStudT_%i_%i", VN, c),"[0] * pow(1 + ( pow([0],2) + pow([1],2) - 2 * [0] * [1] * TMath::Cos(x) ) / ( [3] * pow([2],2) ), -0.5*([3]+1) )", 0, 2*TMath::Pi());

    //-- Gaussian 2SE Fit
    if(GaussResp){

      //-- Set up response function
      hresp[c] = new TH2D(Form("hresp_%i_%i", VN, c), "hresp", NBins[ptBin][c], 0., v2Max[ptBin][c], NBins[ptBin][c], 0., v2Max[ptBin][c]);
      hresp[c]->SetOption("colz");

      h1Dx[c]->Fit(f1Dx[c], "NLM", "", h1Dx[c]->GetMean() - sigmaFit * h1Dx[c]->GetRMS(), h1Dx[c]->GetMean() + sigmaFit * h1Dx[c]->GetRMS());
      h1Dy[c]->Fit(f1Dy[c], "NLM", "", h1Dy[c]->GetMean() - sigmaFit * h1Dy[c]->GetRMS(), h1Dy[c]->GetMean() + sigmaFit * h1Dy[c]->GetRMS());

      sigma_x[c]   = f1Dx[c]->GetParameter("Sigma");
      sigma_xe[c]  = f1Dx[c]->GetParError(f1Dx[c]->GetParNumber("Sigma"));
      mean_x[c]    = f1Dx[c]->GetParameter("Mean");
      chi2_x[c]    = f1Dx[c]->GetChisquare();

      sigma_y[c]   = f1Dy[c]->GetParameter("Sigma");
      sigma_ye[c]  = f1Dy[c]->GetParError(f1Dy[c]->GetParNumber("Sigma"));
      mean_y[c]    = f1Dy[c]->GetParameter("Mean");
      chi2_y[c]    = f1Dy[c]->GetChisquare();

      sigma_2SE[c] = 0.5*(sigma_x[c]+sigma_y[c]);
      ndf_x[c]     = f1Dx[c]->GetNDF();
      ndf_y[c]     = f1Dy[c]->GetNDF();

      double sigma = sigma_2SE[c]/2.;

      //-- Gaussian Response Function
      if ( sw ) {
	hw[c] = (TH1D*) h1D[c]->Clone(Form("hw_%i", c));
      }
      for ( int i = 1; i <= NBins[ptBin][c]; i++ ) {
	for ( int j = 1; j <= NBins[ptBin][c]; j++ ) {
	  double w = 1.;
	  double v_mess = hresp[c]->GetXaxis()->GetBinCenter(i); 
	  double v_true = hresp[c]->GetYaxis()->GetBinCenter(j);
	  if ( sw ) {
            w = hw[c]->GetBinContent(j);
          }
	  double resp = v_mess * TMath::Gaus(sqrt(v_mess*v_mess + v_true*v_true), 0, sigma) * TMath::BesselI0( v_mess*v_true/sigma/sigma );
	  //if ( i == 1 ) cout << "!!! i = " << i << "\t j = " << j << "\t resp = " << resp << endl;                                                                                                                     
	  if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
	  hresp[c]->SetBinContent(i, j, resp*w);
	}
      }
    }

    //-- Student's T response function (numerically integrated)
    if(StudTresp){

      //-- Set up response function
      hresp[c] = new TH2D(Form("hresp_%i_%i", VN, c), "hresp", NBins[ptBin][c], 0., v2Max[ptBin][c], NBins[ptBin][c], 0., v2Max[ptBin][c]);
      hresp[c]->SetOption("colz");

      h1Dx[c]->Fit(f1Dx_StudT[c], "NLM", "", h1Dx[c]->GetMean() - sigmaFit * h1Dx[c]->GetRMS(), h1Dx[c]->GetMean() + sigmaFit * h1Dx[c]->GetRMS());
      h1Dy[c]->Fit(f1Dy_StudT[c], "NLM", "", h1Dy[c]->GetMean() - sigmaFit * h1Dy[c]->GetRMS(), h1Dy[c]->GetMean() + sigmaFit * h1Dy[c]->GetRMS());

      sigma_x[c]   = f1Dx_StudT[c]->GetParameter(2);
      sigma_xe[c]  = f1Dx_StudT[c]->GetParError(f1Dx_StudT[c]->GetParNumber("Sigma"));
      mean_x[c]    = f1Dx_StudT[c]->GetParameter(1);
      chi2_x[c]    = f1Dx_StudT[c]->GetChisquare();

      sigma_y[c]   = f1Dy_StudT[c]->GetParameter(2);
      sigma_ye[c]  = f1Dy_StudT[c]->GetParError(f1Dy_StudT[c]->GetParNumber("Sigma"));
      mean_y[c]    = f1Dy_StudT[c]->GetParameter(1);
      chi2_y[c]    = f1Dy_StudT[c]->GetChisquare();

      sigma_2SE[c] = 0.5*(sigma_x[c]+sigma_y[c]);
      ndf_x[c]     = f1Dx_StudT[c]->GetNDF();
      ndf_y[c]     = f1Dy_StudT[c]->GetNDF();

      double sigma = sigma_2SE[c]/2.;

      if ( sw ) {
        hw[c] = (TH1D*) h1D[c]->Clone(Form("hw_%i", c));
      }
      for ( int i = 1; i <= NBins[ptBin][c]; i++ ) {
        for ( int j = 1; j <= NBins[ptBin][c]; j++ ) {
	  double w = 1.;
          double v_mess = hresp[c]->GetXaxis()->GetBinCenter(i);
	  double v_true = hresp[c]->GetYaxis()->GetBinCenter(j);
	  if ( sw ) {
            w = hw[c]->GetBinContent(j);
          }
	  //-- [0] = vn_meas
	  //-- [1] = vn_true
	  //-- [2] = delta_vn (obtained from 2SE fits)
	  //-- [3] = nu 
	  fRespStudT[c]->SetParameters(v_mess, v_true, sigma, ndf);
          double resp = fRespStudT[c]->Integral(0,2*TMath::Pi());
	  if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
          hresp[c]->SetBinContent(i, j, resp*w);
        }
      } 
    }

    if(!hresp[c]){
      std::cout<<"Response function has not built, please check the code..."<<std::endl;
      break;
    }

    //-- Unfold!
    response[c] = new RooUnfoldResponse(0, 0, hresp[c], Form("response_%i_%i", VN, c));
    for(int i = 0; i < NITER; i++){

      RooUnfoldBayes unfolder(response[c], h1D[c], iter[i] );
      hreco[c][i]   = (TH1D*) unfolder.Hreco();
      hreco[c][i]->SetName(Form("hreco%i_%i_%i", iter[i], VN, c) );
      hrefold[c][i] = (TH1D*) response[c]->ApplyToTruth( hreco[c][i] );
      hrefold[c][i]->SetName( Form("hrefold%i_%i_%i", iter[i], VN, c) );
    }  

    //Save Files
    fsave->cd();
    h2D[c]->Write();
    h2Dsub0[c]->Write();
    h2Dsub1[c]->Write();
    h2Dx[c]->Write();

    h1D[c]->Write();
    h1Dsub0[c]->Write();
    h1Dsub1[c]->Write();
    h1Dx[c]->Write();
    h1Dy[c]->Write();

    if(GaussResp){
      f1Dx[c]->Write();
      f1Dy[c]->Write();
    }
    if(StudTresp){
      f1Dx_StudT[c]->Write();
      f1Dy_StudT[c]->Write();
    }

    if(Resp_Data){
      h22D01[c]->Write();
      h22D01Magnitude[c]->Write();
    }

    hresp[c]->Write();

    for(int i = 0; i < NITER; i++){
      hreco[c][i]->Write();
      hrefold[c][i]->Write();
    }

  }

}

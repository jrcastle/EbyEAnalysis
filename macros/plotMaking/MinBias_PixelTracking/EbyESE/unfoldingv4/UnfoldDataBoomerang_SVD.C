#include "TSVDUnfold.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldSvd.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldBayes.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace ebyese;

//int centbin = 0;  //  0 -  5 %
//int centbin = 1;  //  5 - 10 %
//int centbin = 2;  // 10 - 15 %
//int centbin = 3;  // 15 - 20 %
//int centbin = 4;  // 20 - 25 %
//int centbin = 5;  // 25 - 30 %
//int centbin = 6;  // 30 - 35 %
//int centbin = 7;  // 35 - 40 %
//int centbin = 8;  // 40 - 45 %
//int centbin = 9;  // 45 - 50 %
//int centbin = 10; // 50 - 55 %
//int centbin = 11; // 55 - 60 %

//
//-- MAIN
//

void UnfoldDataBoomerang_SVD( int nord ){

  int norder_ = nord;

  bool Resp_Data     = 1;
  bool StudTresp     = 0;
  bool GaussResp     = 0;

  bool smoothIter    = 0;
  bool smoothResp    = 0;
  bool isMC          = 0;
  bool sw            = 1;

  double sigmaFit = 2.;

  TLatex latex;

  TFile * fsave;
  TFile * fData;
  TFile * fPrior;

  //-- Unfolding Objects
  TH2D * h2D[NCENT];
  TH2D * h2Dx[NCENT];

  TH1D * h1D[NCENT];
  TH1D * h1Dsub0[NCENT];
  TH1D * h1Dsub1[NCENT];
  TH1D * h1Dx[NCENT];
  TH1D * h1Dy[NCENT];
  TH1I * hMult[NCENT];
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

  //-- SVD Stuff
  TH1D * hdi[NCENT];
  RooUnfoldResponse * response[NCENT];
  RooUnfoldSvd * unfold0[NCENT];
  RooUnfoldSvd * unfoldkreg[9][NCENT];
  TSVDUnfold * svdUnf[NCENT];
  TH2D * hCovMat[NCENT][9];

  TH1D * hreco0[NCENT];
  TH1D * hrecokreg[9][NCENT];

  TH1D * hrefold[NCENT][9];
  TH1D * hKreg;

  setTDRStyle();
  latex.SetNDC();

  fsave = 0;
  if( isMC ) fsave = new TFile(Form("txt/PbPb_2015/MC/data%i.root", norder_), "recreate");
  else       fsave = new TFile(Form("txt/PbPb_2015/data/data%i.root", norder_), "recreate");

  bool R1 = Resp_Data && StudTresp && GaussResp;
  bool R2 = Resp_Data && StudTresp;
  bool R3 = Resp_Data && GaussResp;
  bool R4 = StudTresp && GaussResp;

  if( R1 || R2 || R3 || R4){
    std::cout<<"WARNING! More than one response function defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  hKreg = new TH1D("hKreg", "hKreg", NCENT, centbinsDefault);
  hKreg->GetXaxis()->SetTitle("Centrality %");
  hKreg->GetYaxis()->SetTitle("kreg");

  //-- Get histos from datafile
  fData = 0;
  if( isMC ) fData = new TFile("data/PbPb_2015/MC/CastleEbyE.root");
  else       fData = new TFile("data/PbPb_2015/data/CastleEbyE.root");

  //-- If using the data-driven response function, grab the file that has them
  fPrior = 0;
  if( Resp_Data ) fPrior = new TFile("DDResp/dataDrivenResponseAndPriors.root");

  //-- Begin loop over centrality, EP, and qn
  for(int c = 0; c < NCENT; c++){

    std::cout<<"!! Processing Cent = "<<c<<std::endl;

    h2D[c]             = (TH2D*) fData->Get( Form("qwebye/hVn2Dfull_c%i", c) );
    h2Dx[c]            = (TH2D*) fData->Get( Form("qwebye/hVn2D0v1_c%i",  c) );
    h1D[c]             = (TH1D*) fData->Get( Form("qwebye/hVnFull_c%i",   c) );
    hMult[c]           = (TH1I*) fData->Get( Form("qwebye/Mult_c%i",      c) );

    if ( h2D[c]->GetEntries() < 1000 ) continue;

    //-- Data-driven response
    hresp[c] = 0;
    if(Resp_Data) hresp[c] = (TH2D*) fPrior->Get( Form("hresp_c%i", c) );

    h1Dx[c] = h2Dx[c]->ProjectionX( Form("h1Dx_c%i", c) );
    h1Dy[c] = h2Dx[c]->ProjectionY( Form("h1Dy_c%i", c) );

    //-- Gaussian 2SE Fit
    if(GaussResp){

      //-- Set up response function
      hresp[c] = new TH2D( Form("hresp_c%i", c), Form("hresp_c%i", c), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
      hresp[c]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
      hresp[c]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
      hresp[c]->SetOption("colz");

      f1Dx[c] = new TF1( Form("f1Dx_c%i", c), "gaus", -vnMax[norder_], vnMax[norder_] );
      f1Dy[c] = new TF1( Form("f1Dy_c%i", c), "gaus", -vnMax[norder_], vnMax[norder_] );

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
      for ( int i = 1; i <= NBins; i++ ) {
	for ( int j = 1; j <= NBins; j++ ) {
	  double w = 1.;
	  double v_mess = hresp[c]->GetXaxis()->GetBinCenter(i); 
	  double v_true = hresp[c]->GetYaxis()->GetBinCenter(j);
	  if ( sw ) {
	    w = h1D[c]->GetBinContent(j);
	  }
	  double resp = v_mess * TMath::Gaus(sqrt(v_mess*v_mess + v_true*v_true), 0, sigma) * TMath::BesselI0( v_mess*v_true/sigma/sigma );
	  //if ( i == 1 ) cout << "!!! i = " << i << "\t j = " << j << "\t resp = " << resp << endl;                                                                                                                     
	  if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
	  hresp[c]->SetBinContent(i, j, resp*w);
	}
      }
    } //-- End if(GaussResp)

    //-- Student's T response function (numerically integrated)
    if(StudTresp){

      //-- Set up response function
      hresp[c] = new TH2D( Form("hresp_c%i", c), Form("hresp_c%i", c), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
      hresp[c]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
      hresp[c]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
      hresp[c]->SetOption("colz");

      int ndf = hMult[c]->GetMean() - 1;

      f1Dx_StudT[c] = new TF1(Form("f1Dx_StudT_c%i", c), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
      f1Dx_StudT[c]->SetParName(0,"Norm");
      f1Dx_StudT[c]->SetParName(1,"Mean");
      f1Dx_StudT[c]->SetParName(2,"Sigma");
      f1Dx_StudT[c]->SetParName(3,"nu");
      f1Dx_StudT[c]->SetParameters(h1Dx[c]->GetMaximum(), 0, ((ndf-2.)/ndf)*h1Dx[c]->GetRMS(),ndf);
      f1Dx_StudT[c]->FixParameter(3,ndf);
      f1Dx_StudT[c]->SetLineColor(4);
      f1Dx_StudT[c]->SetLineWidth(2);

      f1Dy_StudT[c] = new TF1(Form("f1Dy_StudT_c%i", c), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
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
      fRespStudT[c] = new TF1(Form("fRespStudT_c%i", c),"[0] * pow(1 + ( pow([0],2) + pow([1],2) - 2 * [0] * [1] * TMath::Cos(x) ) / ( [3] * pow([2],2) ), -0.5*([3]+1) )", 0, 2*TMath::Pi());
      
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

      for ( int i = 1; i <= NBins; i++ ) {
	for ( int j = 1; j <= NBins; j++ ) {
	  double w = 1.;
	  double v_mess = hresp[c]->GetXaxis()->GetBinCenter(i);
	  double v_true = hresp[c]->GetYaxis()->GetBinCenter(j);
	  if ( sw ) {
	    w = h1D[c]->GetBinContent(j);
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
    } // End if(StudTresp)

    if(!hresp[c]){
      std::cout<<"Response function has not been built, please check the code..."<<std::endl;
      break;
    }

    //-- Unfold!
    response[c] = new RooUnfoldResponse( 0, 0, hresp[c], Form("response_c%i", c) );

    /////////////////////////////////////////////////////////////////////////////
    unfold0[c] = new RooUnfoldSvd( response[c], h1D[c], 0 );
    hreco0[c] = (TH1D*) unfold0[c]->Hreco();
    hreco0[c]->SetName(Form("hreco0_c%i", c));

    svdUnf[c]  = (TSVDUnfold*) unfold0[c]->Impl();
    hdi[c]     = (TH1D*) svdUnf[c]->GetD();
    hdi[c]->SetName(Form("hdi_c%i", c));

    int kreg = 0;
    for(int i = 1; i<= hdi[c]->GetNbinsX(); i++){
      double diMag = hdi[c]->GetBinContent(i);
      if(diMag < 1.){
        kreg = i-1;
        break;
      }
    }

    hKreg->SetBinContent(c+1, kreg);

    int dummy = 0;
    for(int ik = kreg-4; ik <= kreg+4; ik++){
      unfoldkreg[dummy][c] = new RooUnfoldSvd( response[c], h1D[c], ik );
      hrecokreg[dummy][c] = (TH1D*) unfoldkreg[dummy][c]->Hreco();
      hrecokreg[dummy][c]->SetName(Form("hrecokreg%i_c%i", dummy, c));
      hrefold[dummy][c] = (TH1D*) response[c]->ApplyToTruth( hrecokreg[dummy][c], Form("hrefoldkreg%i_c%i", dummy, c) );

      TMatrixD covMat = unfoldkreg[dummy][c]->Ereco();
      hCovMat[c][dummy] = new TH2D(Form("hCovMatkreg%i_c%i", dummy, c), Form("hCovMatkreg%i_c%i", dummy, c), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
      hCovMat[c][dummy]->GetXaxis()->SetTitle( Form("v_{%i}^{Unfold}", norder_) );
      hCovMat[c][dummy]->GetYaxis()->SetTitle( Form("v_{%i}^{Unfold}", norder_) );
      hCovMat[c][dummy]->SetOption("colz");
      M2H(covMat, hCovMat[c][dummy]);

      dummy++;
    }
    /////////////////////////////////////////////////////////////////////////////

    //-- Save Files
    fsave->cd();
    if(GaussResp){
      f1Dx[c]->Write();
      f1Dy[c]->Write();
    }
    if(StudTresp){
      f1Dx_StudT[c]->Write();
      f1Dy_StudT[c]->Write();
    }

    h1D[c]->Write();
    
    h1Dx[c]->Write();
    h1Dy[c]->Write();

    hresp[c]->Write();
    
    hdi[c]->Write();
    hreco0[c]->Write();
    for(int i=0; i<9; i++){
      hrecokreg[i][c]->Write();
      hrefold[i][c]->Write();
      hCovMat[c][i]->Write();
    }
	
  } //-- End cent loop

  hKreg->Write();

}

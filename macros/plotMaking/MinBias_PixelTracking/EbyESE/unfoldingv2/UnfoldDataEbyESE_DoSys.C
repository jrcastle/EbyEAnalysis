#include "TLatex.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldBayes.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

//
//-- MAIN
//

void UnfoldDataEbyESE_DoSys( int nord ){

  int norder_ = nord;

  bool Resp_Data     = 0;
  bool StudTresp     = 0;
  bool GaussResp     = 1;

  bool smoothIter    = 0;
  bool smoothResp    = 0;
  bool isMC          = 0;
  bool sw            = 1;

  double sigmaFit = 2.;

  static const int NITER  = 8;
  static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};

  RooUnfoldResponse * response[NCENT][NEPSymm][NQN];

  TLatex latex;

  TFile * fsave;
  TFile * fData;
  TFile * fPrior;

  //-- Unfolding Objects
  TH2D * h2D[NCENT][NEPSymm][NQN];
  TH2D * h2Dx[NCENT][NEPSymm][NQN];

  TH1D * h1D[NCENT][NEPSymm][NQN];
  TH1D * h1Dsub0[NCENT][NEPSymm][NQN];
  TH1D * h1Dsub1[NCENT][NEPSymm][NQN];
  TH1D * h1Dx[NCENT][NEPSymm][NQN];
  TH1D * h1Dy[NCENT][NEPSymm][NQN];
  TH1I * hMult[NCENT][NEPSymm][NQN];
  TF1 * f1Dx[NCENT][NEPSymm][NQN];
  TF1 * f1Dy[NCENT][NEPSymm][NQN];
  TF1 * f1Dx_StudT[NCENT][NEPSymm][NQN];
  TF1 * f1Dy_StudT[NCENT][NEPSymm][NQN];
  TF1 * fRespStudT[NCENT][NEPSymm][NQN];

  //-- Data objects for fitting analytic response functions
  double sigma_x[NCENT][NEPSymm][NQN];
  double sigma_xe[NCENT][NEPSymm][NQN];
  double mean_x[NCENT][NEPSymm][NQN];
  double chi2_x[NCENT][NEPSymm][NQN];
  double ndf_x[NCENT][NEPSymm][NQN];
  double sigma_y[NCENT][NEPSymm][NQN];
  double sigma_ye[NCENT][NEPSymm][NQN];
  double mean_y[NCENT][NEPSymm][NQN];
  double chi2_y[NCENT][NEPSymm][NQN];
  double ndf_y[NCENT][NEPSymm][NQN];
  double sigma_2SE[NCENT][NEPSymm][NQN];
  double Vn_mean[NCENT][NEPSymm][NQN];
  double Vn_rms[NCENT][NEPSymm][NQN];

  TH2D * hresp[NCENT][NEPSymm][NQN];
  TH1D * hreco[NCENT][NEPSymm][NQN][NITER];
  TH1D * hrefold[NCENT][NEPSymm][NQN][NITER];

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

  //-- Get histos from datafile
  fData = 0;
  if( isMC ) fData = new TFile("data/PbPb_2015/MC/CastleEbyE.root");
  else       fData = new TFile("data/PbPb_2015/data/CastleEbyE.root");

  //-- If using the data-driven response function, grab the file that has them
  fPrior = 0;
  if( Resp_Data ) fPrior = new TFile("DDResp/dataDrivenResponseAndPriors.root");

  //-- Begin loop over centrality, EP, and qn
  for(int c = 0; c < NCENT; c++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      for(int iqn = 0; iqn < NQN; iqn++){    

	if(iEP != EPSymmBin) continue;
	std::cout<<"!! Processing Cent = "<<c<<std::endl;

	h2D[c][iEP][iqn]             = (TH2D*) fData->Get( Form("qwebye/hVn2Dfull_%s_c%i_qbin%i",          EPSymmNames[iEP].data(), c, iqn) );
	h2Dx[c][iEP][iqn]            = (TH2D*) fData->Get( Form("qwebye/hVn2D0v1_%s_c%i_qbin%i",           EPSymmNames[iEP].data(), c, iqn) );
	h1D[c][iEP][iqn]             = (TH1D*) fData->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i",            EPSymmNames[iEP].data(), c, iqn) );
	//h1Dsub0[c][iEP][iqn]         = (TH1D*) fData->Get( Form("qwebye/hVnSub0_%s_c%i_qbin%i",            EPSymmNames[iEP].data(), c, iqn) );
	//h1Dsub1[c][iEP][iqn]         = (TH1D*) fData->Get( Form("qwebye/hVnSub1_%s_c%i_qbin%i",            EPSymmNames[iEP].data(), c, iqn) );
	hMult[c][iEP][iqn]           = (TH1I*) fData->Get( Form("qwebye/Mult_%s_c%i_qbin%i",               EPSymmNames[iEP].data(), c, iqn) );

	if ( h2D[c][iEP][iqn]->GetEntries() < 1000 ) continue;

	hresp[c][iEP][iqn] = 0;
        if(Resp_Data) hresp[c][iEP][iqn] = (TH2D*) fPrior->Get( Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn) );

	h1Dx[c][iEP][iqn] = h2Dx[c][iEP][iqn]->ProjectionX( Form("h1Dx_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn) );
	h1Dy[c][iEP][iqn] = h2Dx[c][iEP][iqn]->ProjectionY( Form("h1Dy_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn) );

	//-- Gaussian 2SE Fit
	if(GaussResp){

	  //-- Set up response function
	  hresp[c][iEP][iqn] = new TH2D( Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
	  hresp[c][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
          hresp[c][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
	  hresp[c][iEP][iqn]->SetOption("colz");

	  f1Dx[c][iEP][iqn] = new TF1( Form("f1Dx_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), "gaus", -vnMax[norder_], vnMax[norder_] );
	  f1Dy[c][iEP][iqn] = new TF1( Form("f1Dy_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), "gaus", -vnMax[norder_], vnMax[norder_] );

	  h1Dx[c][iEP][iqn]->Fit(f1Dx[c][iEP][iqn], "NLM", "", h1Dx[c][iEP][iqn]->GetMean() - sigmaFit * h1Dx[c][iEP][iqn]->GetRMS(), h1Dx[c][iEP][iqn]->GetMean() + sigmaFit * h1Dx[c][iEP][iqn]->GetRMS());
	  h1Dy[c][iEP][iqn]->Fit(f1Dy[c][iEP][iqn], "NLM", "", h1Dy[c][iEP][iqn]->GetMean() - sigmaFit * h1Dy[c][iEP][iqn]->GetRMS(), h1Dy[c][iEP][iqn]->GetMean() + sigmaFit * h1Dy[c][iEP][iqn]->GetRMS());

	  sigma_x[c][iEP][iqn]   = f1Dx[c][iEP][iqn]->GetParameter("Sigma");
	  sigma_xe[c][iEP][iqn]  = f1Dx[c][iEP][iqn]->GetParError(f1Dx[c][iEP][iqn]->GetParNumber("Sigma"));
	  mean_x[c][iEP][iqn]    = f1Dx[c][iEP][iqn]->GetParameter("Mean");
	  chi2_x[c][iEP][iqn]    = f1Dx[c][iEP][iqn]->GetChisquare();
	  
	  sigma_y[c][iEP][iqn]   = f1Dy[c][iEP][iqn]->GetParameter("Sigma");
	  sigma_ye[c][iEP][iqn]  = f1Dy[c][iEP][iqn]->GetParError(f1Dy[c][iEP][iqn]->GetParNumber("Sigma"));
	  mean_y[c][iEP][iqn]    = f1Dy[c][iEP][iqn]->GetParameter("Mean");
	  chi2_y[c][iEP][iqn]    = f1Dy[c][iEP][iqn]->GetChisquare();

	  sigma_2SE[c][iEP][iqn] = 0.5*(sigma_x[c][iEP][iqn]+sigma_y[c][iEP][iqn]);
	  ndf_x[c][iEP][iqn]     = f1Dx[c][iEP][iqn]->GetNDF();
	  ndf_y[c][iEP][iqn]     = f1Dy[c][iEP][iqn]->GetNDF();

	  double sigma = sigma_2SE[c][iEP][iqn]/2.;

	  //-- Gaussian Response Function
	  for ( int i = 1; i <= NBins; i++ ) {
	    for ( int j = 1; j <= NBins; j++ ) {
	      double w = 1.;
	      double v_mess = hresp[c][iEP][iqn]->GetXaxis()->GetBinCenter(i); 
	      double v_true = hresp[c][iEP][iqn]->GetYaxis()->GetBinCenter(j);
	      if ( sw ) {
		w = h1D[c][iEP][iqn]->GetBinContent(j);
	      }
	      double resp = v_mess * TMath::Gaus(sqrt(v_mess*v_mess + v_true*v_true), 0, sigma) * TMath::BesselI0( v_mess*v_true/sigma/sigma );
	      //if ( i == 1 ) cout << "!!! i = " << i << "\t j = " << j << "\t resp = " << resp << endl;                                                                                                                     
	      if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
	      hresp[c][iEP][iqn]->SetBinContent(i, j, resp*w);
	    }
	  }
	} //-- End if(GaussResp)

	//-- Student's T response function (numerically integrated)
	if(StudTresp){

	  //-- Set up response function
	  hresp[c][iEP][iqn] = new TH2D( Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), Form("hresp_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), NBins, 0., vnMax[norder_], NBins, 0., vnMax[norder_]);
	  hresp[c][iEP][iqn]->GetXaxis()->SetTitle( Form("v_{%i}^{obs}", norder_) );
          hresp[c][iEP][iqn]->GetYaxis()->SetTitle( Form("v_{%i}^{true}", norder_) );
	  hresp[c][iEP][iqn]->SetOption("colz");

	  int ndf = hMult[c][iEP][iqn]->GetMean() - 1;

	  f1Dx_StudT[c][iEP][iqn] = new TF1(Form("f1Dx_StudT_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
	  f1Dx_StudT[c][iEP][iqn]->SetParName(0,"Norm");
	  f1Dx_StudT[c][iEP][iqn]->SetParName(1,"Mean");
	  f1Dx_StudT[c][iEP][iqn]->SetParName(2,"Sigma");
	  f1Dx_StudT[c][iEP][iqn]->SetParName(3,"nu");
	  f1Dx_StudT[c][iEP][iqn]->SetParameters(h1Dx[c][iEP][iqn]->GetMaximum(), 0, ((ndf-2.)/ndf)*h1Dx[c][iEP][iqn]->GetRMS(),ndf);
	  f1Dx_StudT[c][iEP][iqn]->FixParameter(3,ndf);
	  f1Dx_StudT[c][iEP][iqn]->SetLineColor(4);
	  f1Dx_StudT[c][iEP][iqn]->SetLineWidth(2);

	  f1Dy_StudT[c][iEP][iqn] = new TF1(Form("f1Dy_StudT_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn), "[0] * pow([3]/( [3]+pow( (x-[1])/[2], 2) ), (1.+[3])/2.)  / ( TMath::Sqrt([3]) * [2] * TMath::Beta([3]/2.,0.5) )", -vnMax[norder_], vnMax[norder_]);
	  f1Dy_StudT[c][iEP][iqn]->SetParName(0,"Norm");
	  f1Dy_StudT[c][iEP][iqn]->SetParName(1,"Mean");
	  f1Dy_StudT[c][iEP][iqn]->SetParName(2,"Sigma");
	  f1Dy_StudT[c][iEP][iqn]->SetParName(3,"nu");
	  f1Dy_StudT[c][iEP][iqn]->SetParameters(h1Dy[c][iEP][iqn]->GetMaximum(), 0, ((ndf-2.)/ndf)*h1Dy[c][iEP][iqn]->GetRMS(),ndf);
	  f1Dy_StudT[c][iEP][iqn]->FixParameter(3,ndf);
	  f1Dy_StudT[c][iEP][iqn]->SetLineColor(4);
	  f1Dy_StudT[c][iEP][iqn]->SetLineWidth(2);

	  //-- [0] = vn_meas
	  //-- [1] = vn_true
	  //-- [2] = delta_vn (obtained from 2SE fits)
	  //-- [3] = nu
	  fRespStudT[c][iEP][iqn] = new TF1(Form("fRespStudT_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn),"[0] * pow(1 + ( pow([0],2) + pow([1],2) - 2 * [0] * [1] * TMath::Cos(x) ) / ( [3] * pow([2],2) ), -0.5*([3]+1) )", 0, 2*TMath::Pi());

	  h1Dx[c][iEP][iqn]->Fit(f1Dx_StudT[c][iEP][iqn], "NLM", "", h1Dx[c][iEP][iqn]->GetMean() - sigmaFit * h1Dx[c][iEP][iqn]->GetRMS(), h1Dx[c][iEP][iqn]->GetMean() + sigmaFit * h1Dx[c][iEP][iqn]->GetRMS());
	  h1Dy[c][iEP][iqn]->Fit(f1Dy_StudT[c][iEP][iqn], "NLM", "", h1Dy[c][iEP][iqn]->GetMean() - sigmaFit * h1Dy[c][iEP][iqn]->GetRMS(), h1Dy[c][iEP][iqn]->GetMean() + sigmaFit * h1Dy[c][iEP][iqn]->GetRMS());

	  sigma_x[c][iEP][iqn]   = f1Dx_StudT[c][iEP][iqn]->GetParameter(2);
	  sigma_xe[c][iEP][iqn]  = f1Dx_StudT[c][iEP][iqn]->GetParError(f1Dx_StudT[c][iEP][iqn]->GetParNumber("Sigma"));
	  mean_x[c][iEP][iqn]    = f1Dx_StudT[c][iEP][iqn]->GetParameter(1);
	  chi2_x[c][iEP][iqn]    = f1Dx_StudT[c][iEP][iqn]->GetChisquare();
	  
	  sigma_y[c][iEP][iqn]   = f1Dy_StudT[c][iEP][iqn]->GetParameter(2);
	  sigma_ye[c][iEP][iqn]  = f1Dy_StudT[c][iEP][iqn]->GetParError(f1Dy_StudT[c][iEP][iqn]->GetParNumber("Sigma"));
	  mean_y[c][iEP][iqn]    = f1Dy_StudT[c][iEP][iqn]->GetParameter(1);
	  chi2_y[c][iEP][iqn]    = f1Dy_StudT[c][iEP][iqn]->GetChisquare();

	  sigma_2SE[c][iEP][iqn] = 0.5*(sigma_x[c][iEP][iqn]+sigma_y[c][iEP][iqn]);
	  ndf_x[c][iEP][iqn]     = f1Dx_StudT[c][iEP][iqn]->GetNDF();
	  ndf_y[c][iEP][iqn]     = f1Dy_StudT[c][iEP][iqn]->GetNDF();

	  double sigma = sigma_2SE[c][iEP][iqn]/2.;

	  for ( int i = 1; i <= NBins; i++ ) {
	    for ( int j = 1; j <= NBins; j++ ) {
	      double w = 1.;
	      double v_mess = hresp[c][iEP][iqn]->GetXaxis()->GetBinCenter(i);
	      double v_true = hresp[c][iEP][iqn]->GetYaxis()->GetBinCenter(j);
	      if ( sw ) {
		w = h1D[c][iEP][iqn]->GetBinContent(j);
	      }
	      //-- [0] = vn_meas
	      //-- [1] = vn_true
	      //-- [2] = delta_vn (obtained from 2SE fits)
	      //-- [3] = nu 
	      fRespStudT[c][iEP][iqn]->SetParameters(v_mess, v_true, sigma, ndf);
	      double resp = fRespStudT[c][iEP][iqn]->Integral(0,2*TMath::Pi());
	      if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
	      hresp[c][iEP][iqn]->SetBinContent(i, j, resp*w);
	    }
	  } 
	} // End if(StudTresp)

	if(!hresp[c][iEP][iqn]){
	  std::cout<<"Response function has not built, please check the code..."<<std::endl;
	  break;
	}

	//-- Unfold!
	response[c][iEP][iqn] = new RooUnfoldResponse( 0, 0, hresp[c][iEP][iqn], Form("response_%s_c%i_qbin%i", EPSymmNames[iEP].data(), c, iqn) );
	for(int i = 0; i < NITER; i++){

	  RooUnfoldBayes unfolder(response[c][iEP][iqn], h1D[c][iEP][iqn], iter[i] );
	  hreco[c][iEP][iqn][i]   = (TH1D*) unfolder.Hreco();
	  hreco[c][iEP][iqn][i]->SetName(Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), c, iqn) );
	  hrefold[c][iEP][iqn][i] = (TH1D*) response[c][iEP][iqn]->ApplyToTruth( hreco[c][iEP][iqn][i] );
	  hrefold[c][iEP][iqn][i]->SetName( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), c, iqn) );
	}  

	//-- Save Files
	fsave->cd();
	if(GaussResp){
	  f1Dx[c][iEP][iqn]->Write();
	  f1Dy[c][iEP][iqn]->Write();
	}
	if(StudTresp){
	  f1Dx_StudT[c][iEP][iqn]->Write();
	  f1Dy_StudT[c][iEP][iqn]->Write();
	}

	h1Dx[c][iEP][iqn]->Write();
	h1Dy[c][iEP][iqn]->Write();
	
	hresp[c][iEP][iqn]->Write();
	
	for(int i = 0; i < NITER; i++){
	  hreco[c][iEP][iqn][i]->Write();
	  hrefold[c][iEP][iqn][i]->Write();
	}
	
      } //-- End qn loop
    } //-- End EP loop
  } //-- End cent loop

}

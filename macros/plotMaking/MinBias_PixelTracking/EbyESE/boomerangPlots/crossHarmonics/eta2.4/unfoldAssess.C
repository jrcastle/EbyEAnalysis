#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;

void unfoldAssess(){

  double tkEta = 2.4;
  bool dosys   = 0;

  bool looseChi2IterCut   = 0;
  bool nominalChi2IterCut = 1;
  bool tightChi2IterCut   = 0;

  bool gaussResp = 0;
  bool studTResp = 0;
  bool dataResp  = 1;

  double RMSv2Min   = 0.0;
  double RMSv2Max   = 0.26;
  double MEANv2Min  = 0.0;
  double MEANv2Max  = 0.26;
  double STDEVv2Min = 0.0;
  double STDEVv2Max = 0.085;

  double RMSv3Min   = 0.0;
  double RMSv3Max   = 0.11;
  double MEANv3Min  = 0.0;
  double MEANv3Max  = 0.11;
  double STDEVv3Min = 0.0;
  double STDEVv3Max = 0.048;

  double RMSv4Min   = 0.0;
  double RMSv4Max   = 0.11;
  double MEANv4Min  = 0.0;
  double MEANv4Max  = 0.11;
  double STDEVv4Min = 0.0;
  double STDEVv4Max = 0.048;

  double momentTitleSize = 0.05;
  double momentTitleOffS = 1.5;
  double momentAxisLabS  = 0.04;

  TLatex latex;
  TLatex latex2;

  //-- Analyzer Output
  TFile * fAna[NVN];
  TH1D * hObs[NVN][NCENT];

  //-- Unfolding output
  TFile * fUnfold[NVN];
  TH2D * hResp[NVN][NCENT];
  TH1D * hUnfold[NVN][NCENT][NITER];
  TH1D * hRefold[NVN][NCENT][NITER];

  //-- v3 vs v2 with q2 selection
  TGraphErrors * grRMSV3V2;
  TGraphErrors * grMeanV3V2;
  TGraphErrors * grStDevV3V2;

  //-- v4 vs v2 with q2 selection
  TGraphErrors * grRMSV4V2;
  TGraphErrors * grMeanV4V2;
  TGraphErrors * grStDevV4V2;

  double unfold_meanVN[NVN][NCENT];
  double unfold_StDevVN[NVN][NCENT];
  double unfold_RMSVN[NVN][NCENT];

  double unfold_meanVNe[NVN][NCENT];
  double unfold_StDevVNe[NVN][NCENT];
  double unfold_RMSVNe[NVN][NCENT];

  //-- Canvases
  TCanvas * cMoment2x2_V3V2;
  TCanvas * cMoment2x2_V4V2;

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();
  latex2.SetNDC();
  latex2.SetTextFont(43);
  latex2.SetTextSize(15);

  //-- Set the chi2 cutoff
  bool c2c1 = looseChi2IterCut   && nominalChi2IterCut;
  bool c2c2 = looseChi2IterCut   && tightChi2IterCut;
  bool c2c3 = nominalChi2IterCut && tightChi2IterCut;
  bool c2c4 = looseChi2IterCut   && nominalChi2IterCut && tightChi2IterCut;

  if( c2c1 || c2c2 || c2c3 || c2c4){
    std::cout<<"WARNING! More than one chi2 cutoff scenario defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }
  double chi2Cut;
  if( looseChi2IterCut )   chi2Cut = 1.5;
  if( nominalChi2IterCut ) chi2Cut = 1.2;
  if( tightChi2IterCut )   chi2Cut = 1.0;

  bool R1 = dataResp && studTResp && gaussResp;
  bool R2 = dataResp && studTResp;
  bool R3 = dataResp && gaussResp;
  bool R4 = studTResp && gaussResp;

  if( R1 || R2 || R3 || R4){
    std::cout<<"WARNING! More than one response function defined for unfolding.  Check the flags at the beginning of this macro and fix your mistake."<<std::endl;
    std::cout<<"Exiting macro now...  Have a nice day!"<<std::endl;
    exit(0);
  }

  //-- Start looping over the data...
  for(int ivn = 0; ivn < NVN; ivn++){

    //-- Get the analyzer files
    fAna[ivn] = new TFile( Form("../../v%i/eta%.1f/AnalyzerResults/CastleEbyE.root", vn_[ivn], tkEta) );

    //-- Get the unfold data files
    fUnfold[ivn] = 0;
    if( gaussResp ){
      if( !dosys ) fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/gaussResp/data%i.root", vn_[ivn], tkEta, vn_[ivn]) );
      else         fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/gaussResp/data%i_dosys.root", vn_[ivn], tkEta, vn_[ivn]) );
  }
    if( studTResp ){
      if( !dosys ) fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/studTResp/data%i.root", vn_[ivn], tkEta, vn_[ivn]) );
      else         fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/studTResp/data%i_dosys.root", vn_[ivn], tkEta, vn_[ivn]) );
  }
    if( dataResp ){
      if( !dosys ) fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/dataResp/data%i.root", vn_[ivn], tkEta, vn_[ivn]) );
      else         fUnfold[ivn] = new TFile( Form("../../v%i/eta%.1f/UnfoldResults/dataResp/data%i_dosys.root", vn_[ivn], tkEta, vn_[ivn]) );
    }

    for(int icent = 0; icent < NCENT; icent++){

      //-- Get the VN observed histogram
      hObs[ivn][icent] = (TH1D*) fAna[ivn]->Get( Form("qwebye/hVnFull_c%i", icent) );
      hObs[ivn][icent]->SetName( Form("qwebye/hV%iFull_c%i", vn_[ivn], icent) );
      hObs[ivn][icent]->SetMaximum( 10.*hObs[ivn][icent]->GetMaximum() );

      for(int i = 0; i < NITER; i++){

	//-- Get the unfolded histograms
	hUnfold[ivn][icent][i] = (TH1D*) fUnfold[ivn]->Get( Form("hreco%i_c%i", iter[i], icent) );
	hUnfold[ivn][icent][i]->SetName( Form("hreco%i_v%i_c%i", iter[i], vn_[ivn], icent) );
	hUnfold[ivn][icent][i]->SetLineColor(col[i]);
	hUnfold[ivn][icent][i]->SetMarkerColor(col[i]);

	//-- Get the refolded histograms
	hRefold[ivn][icent][i] = (TH1D*) fUnfold[ivn]->Get( Form("hrefold%i_c%i", iter[i], icent) );
        hRefold[ivn][icent][i]->SetName( Form("hreco%i_v%i_c%i", iter[i], vn_[ivn], icent) );
        hRefold[ivn][icent][i]->SetLineColor(col[i]);
        hRefold[ivn][icent][i]->SetMarkerColor(col[i]); 

	//-- Calculate the moments of each unfolded histogram
	double uMu   = hUnfold[ivn][icent][i]->GetMean();
	double uMue  = hUnfold[ivn][icent][i]->GetMeanError();
	double uSig  = hUnfold[ivn][icent][i]->GetRMS();
	double uSige = hUnfold[ivn][icent][i]->GetRMSError();
	double uRMS  = TMath::Sqrt( pow(uMu, 2) + pow(uSig, 2) );
	double uRMSe = TMath::Sqrt( ( pow(uMu * uMue, 2) + pow(uSig * uSige, 2) ) / ( pow(uMu,2) + pow(uSig,2) ) );

	//-- Chi squares
	double chi2NDF_Iteration;

	chi2NDF_Iteration = hRefold[ivn][icent][i]->Chi2Test(hObs[ivn][icent], "CHI2/NDF");

	if( chi2NDF_Iteration <= 1. ){
	  unfold_meanVN[ivn][icent]   = uMu;
	  unfold_StDevVN[ivn][icent]  = uSig;
	  unfold_RMSVN[ivn][icent]    = uRMS;
		
	  unfold_meanVNe[ivn][icent]  = uMue;
	  unfold_StDevVNe[ivn][icent] = uSige;
	  unfold_RMSVNe[ivn][icent]   = uRMSe;
	  std::cout<<"v" << vn_[ivn] << "\tcentBin = "<< icent <<"\titer = " << iter[i] <<std::endl;
	  break;
	}
	if( i == NITER-1 ){
	  unfold_meanVN[ivn][icent]   = uMu;
	  unfold_StDevVN[ivn][icent]  = uSig;
	  unfold_RMSVN[ivn][icent]    = uRMS;

	  unfold_meanVNe[ivn][icent]  = uMue;
	  unfold_StDevVNe[ivn][icent] = uSige;
	  unfold_RMSVNe[ivn][icent]   = uRMSe;
	  std::cout<<"v" << vn_[ivn] << "\tcentBin = "<< icent <<"\titer = " << iter[i] <<std::endl;
	}
	    
      } //-- End unfold iteration loop

    } //-- End cent loop

    //-- v3 VS v2
    if( vn_[ivn] == 3){

      grMeanV3V2  = new TGraphErrors(NCENT, unfold_meanVN[0], unfold_meanVN[ivn], unfold_meanVNe[0], unfold_meanVNe[ivn]);
      grMeanV3V2->GetXaxis()->SetLimits(MEANv2Min, MEANv2Max);
      grMeanV3V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grMeanV3V2->GetXaxis()->SetTitle("#LTv_{2}#GT");
      grMeanV3V2->GetYaxis()->SetRangeUser(MEANv3Min, MEANv3Max);
      grMeanV3V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grMeanV3V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grMeanV3V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grMeanV3V2->GetYaxis()->SetTitle(Form("#LTv_{%i}#GT",vn_[ivn]));
      grMeanV3V2->SetLineColor(8);
      grMeanV3V2->SetMarkerColor(8);
      grMeanV3V2->SetMarkerStyle(20);
      grMeanV3V2->GetXaxis()->SetNdivisions(509);
      grMeanV3V2->GetYaxis()->SetNdivisions(509);

      grStDevV3V2  = new TGraphErrors(NCENT, unfold_StDevVN[0], unfold_StDevVN[ivn], unfold_StDevVNe[0], unfold_StDevVNe[ivn]);
      grStDevV3V2->GetXaxis()->SetLimits(STDEVv2Min, STDEVv2Max);
      grStDevV3V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grStDevV3V2->GetXaxis()->SetTitle("#sigma_{v_{2}}");
      grStDevV3V2->GetYaxis()->SetRangeUser(STDEVv3Min, STDEVv3Max);
      grStDevV3V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grStDevV3V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grStDevV3V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grStDevV3V2->GetYaxis()->SetTitle(Form("#sigma_{v_{%i}}", vn_[ivn]));
      grStDevV3V2->SetLineColor(46);
      grStDevV3V2->SetMarkerColor(46);
      grStDevV3V2->SetMarkerStyle(20);
      grStDevV3V2->GetXaxis()->SetNdivisions(509);
      grStDevV3V2->GetYaxis()->SetNdivisions(509);

      grRMSV3V2  = new TGraphErrors(NCENT, unfold_RMSVN[0], unfold_RMSVN[ivn], unfold_RMSVNe[0], unfold_RMSVNe[ivn]);
      grRMSV3V2->GetXaxis()->SetLimits(RMSv2Min, RMSv2Max);
      grRMSV3V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grRMSV3V2->GetXaxis()->SetTitle("#sqrt{#LTv_{2}^{2}#GT}");
      grRMSV3V2->GetYaxis()->SetRangeUser(RMSv3Min, RMSv3Max);
      grRMSV3V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grRMSV3V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grRMSV3V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grRMSV3V2->GetYaxis()->SetTitle(Form("#sqrt{#LTv_{%i}^{2}#GT}", vn_[ivn]));
      grRMSV3V2->SetLineColor(9);
      grRMSV3V2->SetMarkerColor(9);
      grRMSV3V2->SetMarkerStyle(20);
      grRMSV3V2->GetXaxis()->SetNdivisions(509);
      grRMSV3V2->GetYaxis()->SetNdivisions(509);

    }

    //-- v4 vs v2
    if( vn_[ivn] == 4){

      grMeanV4V2  = new TGraphErrors(NCENT, unfold_meanVN[0], unfold_meanVN[ivn], unfold_meanVNe[0], unfold_meanVNe[ivn]);
      grMeanV4V2->GetXaxis()->SetLimits(MEANv2Min, MEANv2Max);
      grMeanV4V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grMeanV4V2->GetXaxis()->SetTitle("#LTv_{2}#GT");
      grMeanV4V2->GetYaxis()->SetRangeUser(MEANv4Min, MEANv4Max);
      grMeanV4V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grMeanV4V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grMeanV4V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grMeanV4V2->GetYaxis()->SetTitle(Form("#LTv_{%i}#GT",vn_[ivn]));
      grMeanV4V2->SetLineColor(8);
      grMeanV4V2->SetMarkerColor(8);
      grMeanV4V2->SetMarkerStyle(20);
      grMeanV4V2->GetXaxis()->SetNdivisions(509);
      grMeanV4V2->GetYaxis()->SetNdivisions(509);
      
      grStDevV4V2  = new TGraphErrors(NCENT, unfold_StDevVN[0], unfold_StDevVN[ivn], unfold_StDevVNe[0], unfold_StDevVNe[ivn]);
      grStDevV4V2->GetXaxis()->SetLimits(STDEVv2Min, STDEVv2Max);
      grStDevV4V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grStDevV4V2->GetXaxis()->SetTitle("#sigma_{v_{2}}");
      grStDevV4V2->GetYaxis()->SetRangeUser(STDEVv4Min, STDEVv4Max);
      grStDevV4V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grStDevV4V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grStDevV4V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grStDevV4V2->GetYaxis()->SetTitle(Form("#sigma_{v_{%i}}", vn_[ivn]));
      grStDevV4V2->SetLineColor(46);
      grStDevV4V2->SetMarkerColor(46);
      grStDevV4V2->SetMarkerStyle(20);
      grStDevV4V2->GetXaxis()->SetNdivisions(509);
      grStDevV4V2->GetYaxis()->SetNdivisions(509);

      grRMSV4V2  = new TGraphErrors(NCENT, unfold_RMSVN[0], unfold_RMSVN[ivn], unfold_RMSVNe[0], unfold_RMSVNe[ivn]);
      grRMSV4V2->GetXaxis()->SetLimits(RMSv2Min, RMSv2Max);
      grRMSV4V2->GetXaxis()->SetLabelSize(momentAxisLabS);
      grRMSV4V2->GetXaxis()->SetTitle("#sqrt{#LTv_{2}^{2}#GT}");
      grRMSV4V2->GetYaxis()->SetRangeUser(RMSv4Min, RMSv4Max);
      grRMSV4V2->GetYaxis()->SetTitleSize( momentTitleSize );
      grRMSV4V2->GetYaxis()->SetTitleOffset( momentTitleOffS );
      grRMSV4V2->GetYaxis()->SetLabelSize( momentAxisLabS );
      grRMSV4V2->GetYaxis()->SetTitle(Form("#sqrt{#LTv_{%i}^{2}#GT}", vn_[ivn]));
      grRMSV4V2->SetLineColor(9);
      grRMSV4V2->SetMarkerColor(9);
      grRMSV4V2->SetMarkerStyle(20);
      grRMSV4V2->GetXaxis()->SetNdivisions(509);
      grRMSV4V2->GetXaxis()->SetNdivisions(509);

    }

  } //-- End vn loop

  //-- 2x2 v3v2 correlations
  cMoment2x2_V3V2 = new TCanvas( "cMoment2x2_V3V2", "cMoment2x2_V3V2", 1000, 1000);
  cMoment2x2_V3V2->Divide(2,2);
  cMoment2x2_V3V2->cd(1);
  grRMSV3V2->Draw("alp");
  cMoment2x2_V3V2->cd(2);
  grMeanV3V2->Draw("alp");
  cMoment2x2_V3V2->cd(3);
  grStDevV3V2->Draw("alp");
  cMoment2x2_V3V2->Update();
  cMoment2x2_V3V2->SaveAs( "plots/cMoment2x2_V3V2.pdf" );

  //-- 2x2 v4v2 correlations
  cMoment2x2_V4V2 = new TCanvas( "cMoment2x2_V4V2", "cMoment2x2_V4V2", 1000, 1000);
  cMoment2x2_V4V2->Divide(2,2);
  cMoment2x2_V4V2->cd(1);
  grRMSV4V2->Draw("alp");
  cMoment2x2_V4V2->cd(2);
  grMeanV4V2->Draw("alp");
  cMoment2x2_V4V2->cd(3);
  grStDevV4V2->Draw("alp");
  grStDevV4V2->Draw("lpsame");
  cMoment2x2_V4V2->Update();
  cMoment2x2_V4V2->SaveAs( "plots/cMoment2x2_V4V2.pdf" );


  for(int icent = 0; icent < NCENT; icent++) std::cout << "cent = " << icent << "\tv2 = "<< unfold_meanVN[0][icent] << "\tv4 = "<< unfold_meanVN[2][icent] << std::endl;


}

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyECumu.h"
#include "/home/j550c590/tdrstyle.C"

#include <iostream>

using namespace hi;
using namespace ebyese;


void v4RMS_v2RMSsq(){

  int norder_     = 2;
  double tkEta    = 2.4;
  int QnBinOrder_ = 2;

  bool dosys_      = 0;
  double sysWidth_ = 0.1;

  double mar   = 0.2;
  double offsx = 1.2;
  double offsy = 1.5;

  bool looseChi2IterCut   = 0;
  bool nominalChi2IterCut = 1;
  bool tightChi2IterCut   = 0;

  double qnmin = 0;
  double qnmax = 0.26;
  double v4Rv2RsMin = 0.;
  double v4Rv2RsMax = 25.;

  TLatex latex;

  TFile * fQn;
  TH1D * hqbins[NCENT][NEPSymm];
  TH1D * hqnHF_EP[NCENT][NEPSymm];
  double qnBinCenter[NCENT][NEPSymm][NQN];
  double qnBinCentere[NCENT][NEPSymm][NQN];

  //-- Analyzer Output
  TFile * fAna[NVN];
  TH1D * hObs[NVN][NCENT][NEPSymm][NQN];

  //-- Unfolding output
  TFile * fUnfold[NVN];
  TH1D * hUnfold[NVN][NCENT][NEPSymm][NQN][NITER];
  TH1D * hRefold[NVN][NCENT][NEPSymm][NQN][NITER];

  //-- Statististical Errors
  TFile * fStatErr[NVN];
  TH1D * hVarianceOfMean_RMSVn[NVN][NEPSymm][NQN];
  TH1D * hVarianceOfMean_RelFluctVn[NVN][NEPSymm][NQN];

  //-- Final Unfolded iteration
  int iterCut[NVN][NCENT][NEPSymm][NQN];

  //-- Moment vs vm 
  double rmsVn_vs_vm[NVN][NCENT][NEPSymm][NQN];
  double rmsVn_vs_vm_statErr[NVN][NCENT][NEPSymm][NQN];
  double relFluctVn_vs_vm[NVN][NCENT][NEPSymm][NQN];
  double relFluctVn_vs_vm_statErr[NVN][NCENT][NEPSymm][NQN];

  double v4RMS_v2RMSsq[NCENT][NEPSymm][NQN];
  double v4RMS_v2RMSsqe[NCENT][NEPSymm][NQN];

  double v4RMS_v2RMSsq_IdealRat[NCENT][NEPSymm][NQN];
  double v4RMS_v2RMSsq_IdealRate[NCENT][NEPSymm][NQN];

  TGraphErrors * grV4RMSV2RMSsq[NCENT][NEPSymm];
  TGraphErrors * grV4RMSV2RMSsq_IdealRat[NCENT][NEPSymm];
  TH1D * h[NCENT][NEP];

  //
  // MAIN
  //

  setTDRStyle();
  latex.SetNDC();

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

  //double sysWidth[NQN];
  //-- Widths for systematic error bars
  //for(intiqn = 0; iqn < NQN; iqn++) sysWidth[iqn] = sysWidth_;

  fQn = new TFile( Form( "../../v%i/eta2.4/AnalyzerResults/q%iCuts.root", QnBinOrder_, QnBinOrder_) );
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      hqbins[icent][iEP]   = (TH1D*) fQn->Get( Form("hqbins_%s_c%i", EPSymmNames[iEP].data(), icent) );
      hqnHF_EP[icent][iEP] = (TH1D*) fQn->Get( Form("hqnHF_c%i_EP%i", icent, iEP) );
      for(int iqn = 0; iqn < NQN; iqn++){
	//-- Determine the SW qn values
	double qMin = hqbins[icent][iEP]->GetBinLowEdge(iqn+1);
	double qMax = hqbins[icent][iEP]->GetBinLowEdge(iqn+2);
	int binMin  = hqnHF_EP[icent][iEP]->FindBin(qMin);
	int binMax  = hqnHF_EP[icent][iEP]->FindBin(qMax);
	double qnSW = 0.;
	double sumw = 0;
	for(int ib = binMin; ib <= binMax; ib++){
	  double w  = hqnHF_EP[icent][iEP]->GetBinContent(ib);
	  double qn = hqnHF_EP[icent][iEP]->GetBinCenter(ib);
	  qnSW += w * qn;
	  sumw += w;
	}
	qnSW /= sumw;
	qnBinCenter[icent][iEP][iqn]  = qnSW;
	qnBinCentere[icent][iEP][iqn] = 0;
      }
    }
  }



  for(int ivn = 0; ivn < NVN; ivn++){

    //-- Get the Analyzer output file
    fAna[ivn] = new TFile( Form("../../v%i/eta2.4/AnalyzerResults/CastleEbyE.root", vn_[ivn]) );

    //-- Get Unfolding output
    fUnfold[ivn] = new TFile( Form("../../v%i/eta2.4/UnfoldResults/dataResp/data%i.root", vn_[ivn], vn_[ivn]) );

    //-- Get Statististical Errors
    fStatErr[ivn] = new TFile( Form("../../statErrorHandle/v%i/eta2.4/StatisticalUncertainties_v%i.root", vn_[ivn], vn_[ivn]) );
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      for(int iqn = 0; iqn < NQN; iqn++){
	hVarianceOfMean_RMSVn[ivn][iEP][iqn]      = (TH1D*) fStatErr[ivn]->Get( Form("hVarianceOfMean_RMSVn_%s_qbin%i",      EPSymmNames[iEP].data(), iqn) );
	hVarianceOfMean_RelFluctVn[ivn][iEP][iqn] = (TH1D*) fStatErr[ivn]->Get( Form("hVarianceOfMean_RelFluctVn_%s_qbin%i", EPSymmNames[iEP].data(), iqn) );
      }
    }

    //-- Start looping over the data...
    for(int icent = 0; icent < NCENT; icent++){

      for(int iEP = 0; iEP < NEPSymm; iEP++){
	if( iEP != EPSymmBin ) continue;

	for(int iqn = 0; iqn < NQN; iqn++){

	  //-- Get the VN observed histogram
	  hObs[ivn][icent][iEP][iqn] = (TH1D*) fAna[ivn]->Get( Form("qwebye/hVnFull_%s_c%i_qbin%i", EPSymmNames[iEP].data(), icent, iqn) );

	  iterCut[ivn][icent][iEP][iqn] = 0;
	  for(int i = 0; i < NITER; i++){

	    //-- Get the unfolded histograms
	    hUnfold[ivn][icent][iEP][iqn][i] = (TH1D*) fUnfold[ivn]->Get( Form("hreco%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	    //-- Get the Refolded histograms
	    hRefold[ivn][icent][iEP][iqn][i] = (TH1D*) fUnfold[ivn]->Get( Form("hrefold%i_%s_c%i_qbin%i", iter[i], EPSymmNames[iEP].data(), icent, iqn) );

	    //-- Chi2 Cut
	    double chi2NDF_Refold = hRefold[ivn][icent][iEP][iqn][i]->Chi2Test(hObs[ivn][icent][iEP][iqn], "CHI2/NDF");

	    if( chi2NDF_Refold <= chi2Cut ){
	      iterCut[ivn][icent][iEP][iqn] = i;
	      break;
	    }
	    if( i == NITER-1 ){
	      iterCut[ivn][icent][iEP][iqn] = i;
	      break;
	    }

	  } //-- End unfold iteration loop


	  //-- Fill Arrays for plots
	  int it = iterCut[ivn][icent][iEP][iqn];
	  FixUnfold( hUnfold[ivn][icent][iEP][iqn][it] );

	  //-- Values
	  double mean     = hUnfold[ivn][icent][iEP][iqn][it]->GetMean();
	  double stdev    = hUnfold[ivn][icent][iEP][iqn][it]->GetRMS();
	  double rms      = sqrt( mean*mean + stdev*stdev );
	  double relfluct = stdev / mean;

	  //-- Stat Errors
	  double rmsStatErr      = sqrt( hVarianceOfMean_RMSVn[ivn][iEP][iqn]->GetBinContent(icent+1) );
	  double relFluctStatErr = sqrt( hVarianceOfMean_RelFluctVn[ivn][iEP][iqn]->GetBinContent(icent+1) );

	  rmsVn_vs_vm[ivn][icent][iEP][iqn]         = rms;
	  rmsVn_vs_vm_statErr[ivn][icent][iEP][iqn] = rmsStatErr;

	  relFluctVn_vs_vm[ivn][icent][iEP][iqn]         = relfluct;
	  relFluctVn_vs_vm_statErr[ivn][icent][iEP][iqn] = relFluctStatErr;

	} //-- End QN loop

	//-- Initialize TGraphErrors
	if( ivn == NVN-1 ){

	  h[icent][iEP] = new TH1D(Form("h%i%i", icent, iEP), Form("h%i%i", icent, iEP), 100, 0, 25);

	  for(int iqn = 0; iqn < NQN; iqn++){


	    // v4RMS / v2RMS^2
	    double v4RMS  = rmsVn_vs_vm[2][icent][iEP][iqn];
	    double v4RMSe = rmsVn_vs_vm_statErr[2][icent][iEP][iqn];

	    double v2RMS  = rmsVn_vs_vm[0][icent][iEP][iqn];
            double v2RMSe = rmsVn_vs_vm_statErr[0][icent][iEP][iqn];

	    double v4Rv2R2 = v4RMS / pow(v2RMS,2);
	    double v4Rv2R2e = sqrt( pow(v4RMSe,2)/pow(v2RMS,2) + pow(2.*v4RMS*v2RMSe,2)/pow(v2RMS,6) );

	    v4RMS_v2RMSsq[icent][iEP][iqn]  = v4Rv2R2;
	    v4RMS_v2RMSsqe[icent][iEP][iqn] = v4Rv2R2e;

	    h[icent][iEP]->Fill(v4Rv2R2);


	    //-- (v4RMS / v2RMS^2) / (1 + 4*(relFluctV2)^2)
	    double relFluctV2 = relFluctVn_vs_vm[0][icent][iEP][iqn];
	    double relFluctV2e = relFluctVn_vs_vm_statErr[ivn][icent][iEP][iqn];

	    double RHS = 1. + 4. * pow(relFluctV2,2);
	    double RHSe = 8.*relFluctV2*relFluctV2e; 

	    double IdealRat = v4Rv2R2 / RHS;
	    double IdealRate = sqrt( pow(v4Rv2R2e/RHS, 2) + pow(v4Rv2R2*RHSe/RHS/RHS, 2) );


	    v4RMS_v2RMSsq_IdealRat[icent][iEP][iqn]  = IdealRat;
	    v4RMS_v2RMSsq_IdealRate[icent][iEP][iqn] = IdealRate;

	  }

	  grV4RMSV2RMSsq[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], v4RMS_v2RMSsq[icent][iEP], qnBinCentere[icent][iEP], v4RMS_v2RMSsqe[icent][iEP]);
	  grV4RMSV2RMSsq_IdealRat[icent][iEP] = new TGraphErrors(NQN, qnBinCenter[icent][iEP], v4RMS_v2RMSsq_IdealRat[icent][iEP], qnBinCentere[icent][iEP], v4RMS_v2RMSsq_IdealRate[icent][iEP]);

	  //-- Format Graphs

	  //-- V3 vs V2
	  grV4RMSV2RMSsq[icent][iEP]->GetXaxis()->SetTitle( "HF q_{2}" );
	  grV4RMSV2RMSsq[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
	  grV4RMSV2RMSsq[icent][iEP]->GetXaxis()->SetNdivisions(507);
	  grV4RMSV2RMSsq[icent][iEP]->GetYaxis()->SetTitle( "#sqrt{#LTv_{4}^{2}#GT} / #LTv_{2}^{2}#GT" );
	  grV4RMSV2RMSsq[icent][iEP]->GetYaxis()->SetRangeUser(v4Rv2RsMin, v4Rv2RsMax);
	  grV4RMSV2RMSsq[icent][iEP]->GetYaxis()->SetNdivisions(509);
	  grV4RMSV2RMSsq[icent][iEP]->SetLineColor( centCol[icent] );
	  grV4RMSV2RMSsq[icent][iEP]->SetMarkerColor( centCol[icent] );
	  grV4RMSV2RMSsq[icent][iEP]->SetMarkerStyle( centMark[icent] );

	  grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetXaxis()->SetTitle( "HF q_{2}" );
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetXaxis()->SetLimits(qnmin, qnmax);
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetXaxis()->SetNdivisions(507);
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetYaxis()->SetTitle( "#frac{#sqrt{#LTv_{4}^{2}#GT}}{#LTv_{2}^{2}#GT} / #left[1 + 4#left(#frac{#sigma_{v_{2}}}{#LT v_{2} #GT}#right)^{2}#right]" );
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetYaxis()->SetRangeUser(0, 10);
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetYaxis()->SetNdivisions(509);
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->SetLineColor( centCol[icent] );
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->SetMarkerColor( centCol[icent] );
          grV4RMSV2RMSsq_IdealRat[icent][iEP]->SetMarkerStyle( centMark[icent] );

	} //-- End if( ivn == NVN-1 )

      } //-- End EP loop
    } //-- End cent loop
  } //-- End vn loop

  double aveQ2[NCENT];
  double aveOverQ2[NCENT];
  for(int icent = 0; icent < NCENT; icent++){
    aveQ2[icent] = hqnHF_EP[icent][EPSymmBin]->GetMean();
    aveOverQ2[icent] = h[icent][EPSymmBin]->GetMean();
  }
  TGraph * g = new TGraph(NCENT, aveQ2, aveOverQ2);
  g->SetLineWidth(3);

  //-- DRAW!!!!!
  TLegend * legCent = new TLegend(0.41, 0.40, 0.97, 0.74);
  legCent->SetFillStyle(0);
  legCent->SetBorderSize(0);
  legCent->SetNColumns(2);
  for(int icent = 0; icent < NCENT; icent++) legCent->AddEntry(grV4RMSV2RMSsq[icent][EPSymmBin], Form("Cent %i - %i%s", cent_min[icent], cent_max[icent], "%"), "lp");

  TLegend * l2 = new TLegend(0.53, 0.28, 0.95, 0.46);
  legInit(l2);
  l2->AddEntry(g, "Average centrality behavior", "l");

  TLine * lp5 = new TLine(qnmin, 0.5, qnmax, 0.5);
  lp5->SetLineStyle(2);

  TCanvas * c = new TCanvas("c", "c", 1000, 500);
  c->Divide(2,1);
  c->cd(1)->SetLeftMargin(0.25);
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      grV4RMSV2RMSsq[icent][iEP]->GetXaxis()->SetTitleOffset(offsx);
      grV4RMSV2RMSsq[icent][iEP]->GetYaxis()->SetTitleOffset(offsy);
      if(icent == 0) grV4RMSV2RMSsq[icent][iEP]->Draw("alp");
      else           grV4RMSV2RMSsq[icent][iEP]->Draw("lpsame");
    }
  }
  legCent->Draw("same");
  //g->Draw("lsame");
  //l2->Draw("same");
  latex.DrawLatex(0.53, 0.88, "CMS #it{Preliminary}");
  latex.DrawLatex(0.73, 0.82, Form("|#eta| < %.1f", tkEta) );
  latex.DrawLatex(0.48, 0.76, Form("%.1f < p_{T} < %.1f GeV/c", pt_min[0], pt_max[NPT-1]) );

  c->cd(2)->SetLeftMargin(0.25);
  for(int icent = 0; icent < NCENT; icent++){
    for(int iEP = 0; iEP < NEPSymm; iEP++){
      if( iEP != EPSymmBin ) continue;
      grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetXaxis()->SetTitleOffset(offsx);
      grV4RMSV2RMSsq_IdealRat[icent][iEP]->GetYaxis()->SetTitleOffset(offsy);
      if(icent == 0) grV4RMSV2RMSsq_IdealRat[icent][iEP]->Draw("alp");
      else           grV4RMSV2RMSsq_IdealRat[icent][iEP]->Draw("lpsame");
    }
  }
  lp5->Draw("same");
  c->SaveAs("plots/grV4RMSV2RMSsq.pdf");

}

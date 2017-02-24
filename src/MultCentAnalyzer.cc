// -*- C++ -*-
//
// Package:    MultCentAnalyzer
// Class:      MultCentAnalyzer
//
/**\class MultCentAnalyzer MultCentAnalyzer.cc RecoHI/MultCentAnalyzer/src/MultCentAnalyzer.cc
 
   Description: <one line class summary>
 
   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Sergey Petrushanko
//         Created:  Fri Jul 11 10:05:00 2008
// $Id: MultCentAnalyzer.cc,v 1.18 2011/10/07 09:41:29 yilmaz Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
//#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
//#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
//#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/TrackEfficiency.h"
#include <DataFormats/HeavyIonEvent/interface/ClusterCompatibility.h>
#include "/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"

#include "TROOT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>

#include <vector>
#include <iostream>
using std::vector;
using std::rand;
using namespace std;
using namespace hi;
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

/*
  static const int nptbinsDefault = 16;
  static const double ptbinsDefault[]={
  0.2,  0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0,
  2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0};
  static const int netabinsDefault = 12;
  static const double etabinsDefault[]={-2.4, -2.0, -1.6, -1.2,
  -0.8, -0.4, 0.0,  0.4,  0.8,
  1.2,  1.6,  2.0,  2.4};
*/

//-- NO PIXEL TRACKING
//static const int nptbinsDefault = 5;
//static const double ptbinsDefault[]={1.00, 1.25, 1.50, 2.00, 2.50, 3.00};

//-- ADD PIXEL TRACKING 
static const int nptbinsDefault = 2;
static const double ptbinsDefault[]={0.30, 1.00, 3.00};

static const int netabinsDefault = 4;
static const double etabinsDefault[]= {-2.4, -1.0, 0.0, 1.0, 2.4};

static const int nptspectrumbins  = 10;
static const double ptspectrumbins[] = {0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.50, 2.0, 2.5, 3.0};

static const int NCENT = 20;
static const double centbinsdefault[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
static const int cmin[NCENT] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95};
static const int cmax[NCENT] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};

//
// class decleration
//

class MultCentAnalyzer : public edm::EDAnalyzer {
public:
  explicit MultCentAnalyzer(const edm::ParameterSet&);
  ~MultCentAnalyzer();
    
private:
  //edm::Service<TFileService> fs;
    
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  static double determineQuality(const reco::ClusterCompatibility & cc,
		   double minZ, double maxZ) ;    
  // ----------member data ---------------------------
  int nOrder_;

  unsigned int runno_;
    
  edm::Service<TFileService> fs;
    
  edm::InputTag vertexTag_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag caloTag_;
  edm::EDGetTokenT<CaloTowerCollection> caloToken_;
  edm::Handle<CaloTowerCollection> caloCollection_;

  edm::EDGetTokenT<reco::Centrality> CentralityTag_;
  edm::EDGetTokenT<int> CentralityBinTag_;

  edm::InputTag clusterCompatibilityTag_;
  edm::EDGetTokenT<reco::ClusterCompatibility> clusterCompatibilityToken_;
  edm::Handle<reco::ClusterCompatibility> clusterCompatibilityCollection_;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  double vtx;
    
  double minpt_;
  double maxpt_;
  double minet_;
  double maxet_;
  double etaMax_;
  double minvz_;
  double maxvz_;
  int nvtx_;

  int trackQualityCuts_;
  double dzCut_;
  double d0Cut_;
  double ptResCut_;

  bool usePixelTeff_;
  TrackEfficiency * pxTeff;
    
  bool FirstEvent_;
  string effTable_;

  TH2D * hMultCent_eta24_pt03_3_hiPixelAndGeneral;
  TH2D * hMultCent_eta24_pt03_1_hiPixelAndGeneral;
  TH2D * hMultCent_eta24_pt1_3_hiPixelAndGeneral;
  TH2D * hMultCent_eta10_pt03_3_hiPixelAndGeneral;
  TH2D * hMultCent_eta10_pt03_1_hiPixelAndGeneral;
  TH2D * hMultCent_eta10_pt1_3_hiPixelAndGeneral;

  TH2D * hMultCent_eta24_pt03_3_hiPixel;
  TH2D * hMultCent_eta24_pt03_1_hiPixel;
  TH2D * hMultCent_eta24_pt1_3_hiPixel;
  TH2D * hMultCent_eta10_pt03_3_hiPixel;
  TH2D * hMultCent_eta10_pt03_1_hiPixel;
  TH2D * hMultCent_eta10_pt1_3_hiPixel;

  TH2D * hMultCent_eta24_pt03_3_hiGeneral;
  TH2D * hMultCent_eta24_pt03_1_hiGeneral;
  TH2D * hMultCent_eta24_pt1_3_hiGeneral;
  TH2D * hMultCent_eta10_pt03_3_hiGeneral;
  TH2D * hMultCent_eta10_pt03_1_hiGeneral;
  TH2D * hMultCent_eta10_pt1_3_hiGeneral;

  TH1D * hPtRes;

  TH1D * hPtSpectrum_24;
  TH1D * hPtSpectrum_10;

  TH1D * hCent;
  TH2D * hClusVtxCompInteg;
  TH2D * hClusVtxComp[NCENT];

  TTree * tree;
  double centval;
  int Mult_pt03_3_eta24;
  int Mult_pt03_3_eta10;
  double clusVtxQual;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MultCentAnalyzer::MultCentAnalyzer(const edm::ParameterSet& iConfig):
  runno_(0),
  vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag_")) ),
  trackToken_( consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackTag_")) ),
  CentralityTag_(consumes<reco::Centrality>(iConfig.getUntrackedParameter<edm::InputTag>("CentralityTag_"))),
  CentralityBinTag_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("CentralityBinTag_"))),
  clusterCompatibilityToken_( consumes<reco::ClusterCompatibility>(iConfig.getUntrackedParameter<edm::InputTag>("ClusterCompatibilityTag_")) )
{
  runno_      = 0;
  FirstEvent_ = kTRUE;
    
  vertexTag_               = iConfig.getUntrackedParameter<edm::InputTag>("vertexTag_");
  trackTag_                = iConfig.getUntrackedParameter<edm::InputTag>("trackTag_");
  clusterCompatibilityTag_ = iConfig.getUntrackedParameter<edm::InputTag>("ClusterCompatibilityTag_");
  usePixelTeff_            = iConfig.getUntrackedParameter<bool>("usePixelTeff_",false);
  effTable_                = iConfig.getUntrackedParameter<std::string>("effTable_");
  trackQualityCuts_        = iConfig.getUntrackedParameter<int>("trackQualityCuts_",1);
  pxTeff                   = 0;

  if( trackQualityCuts_ == 1){
    //-- Nominal Cuts
    dzCut_    = 3.0;
    d0Cut_    = 3.0;
    ptResCut_ = 0.1;
  }
  if( trackQualityCuts_ == 2){
    //-- Loose Cuts
    dzCut_    = 5.0;
    d0Cut_    = 5.0;
    ptResCut_ = 0.1;
  }
  if( trackQualityCuts_ == 3){
    //-- Tight Cuts
    dzCut_    = 2.0;
    d0Cut_    = 2.0;
    ptResCut_ = 0.05;
  }

  if( usePixelTeff_ ) pxTeff = new TrackEfficiency(effTable_);

  tree = fs->make<TTree>("tree","tree");
  tree->Branch("Cent",              &centval,           "cent/D");
  tree->Branch("Mult_pt03_3_eta24", &Mult_pt03_3_eta24, "Mult_pt03_3_eta24/i");
  tree->Branch("Mult_pt03_3_eta10", &Mult_pt03_3_eta10, "Mult_pt03_3_eta10/i");
  tree->Branch("clusVtxQual",       &clusVtxQual,       "clusVtxQual/D");

  /*
  hPtSpectrum_24 = fs->make<TH1D>("hPtSpectrum_24", "hPtSpectrum_24", nptspectrumbins, ptspectrumbins);
  hPtSpectrum_10 = fs->make<TH1D>("hPtSpectrum_10", "hPtSpectrum_10", nptspectrumbins, ptspectrumbins);

  // |eta| < 2.4, 0.3 < pT < 3.0
  hMultCent_eta24_pt03_3_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta24_pt03_3_hiPixelAndGeneral", "hMultCent_eta24_pt03_3_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_3_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_3_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_3_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_3_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt03_3_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_3_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta24_pt03_3_hiPixel = fs->make<TH2D>("hMultCent_eta24_pt03_3_hiPixel", "hMultCent_eta24_pt03_3_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_3_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_3_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_3_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_3_hiPixel->GetNbinsY()-1);
  hMultCent_eta24_pt03_3_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_3_hiPixel->SetOption("colz");

  hMultCent_eta24_pt03_3_hiGeneral = fs->make<TH2D>("hMultCent_eta24_pt03_3_hiGeneral", "hMultCent_eta24_pt03_3_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_3_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_3_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_3_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_3_hiGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt03_3_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_3_hiGeneral->SetOption("colz");

  // |eta| < 2.4, 0.3 < pT < 1.0
  hMultCent_eta24_pt03_1_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta24_pt03_1_hiPixelAndGeneral", "hMultCent_eta24_pt03_1_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_1_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_1_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_1_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_1_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt03_1_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_1_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta24_pt03_1_hiPixel = fs->make<TH2D>("hMultCent_eta24_pt03_1_hiPixel", "hMultCent_eta24_pt03_1_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_1_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_1_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_1_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_1_hiPixel->GetNbinsY()-1);
  hMultCent_eta24_pt03_1_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_1_hiPixel->SetOption("colz");

  hMultCent_eta24_pt03_1_hiGeneral = fs->make<TH2D>("hMultCent_eta24_pt03_1_hiGeneral", "hMultCent_eta24_pt03_1_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt03_1_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt03_1_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt03_1_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt03_1_hiGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt03_1_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt03_1_hiGeneral->SetOption("colz");

  // |eta| < 2.4, 1.0 < pT < 3.0
  hMultCent_eta24_pt1_3_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta24_pt1_3_hiPixelAndGeneral", "hMultCent_eta24_pt1_3_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt1_3_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt1_3_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt1_3_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt1_3_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt1_3_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt1_3_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta24_pt1_3_hiPixel = fs->make<TH2D>("hMultCent_eta24_pt1_3_hiPixel", "hMultCent_eta24_pt1_3_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt1_3_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt1_3_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt1_3_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta24_pt1_3_hiPixel->GetNbinsY()-1);
  hMultCent_eta24_pt1_3_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt1_3_hiPixel->SetOption("colz");

  hMultCent_eta24_pt1_3_hiGeneral = fs->make<TH2D>("hMultCent_eta24_pt1_3_hiGeneral", "hMultCent_eta24_pt1_3_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta24_pt1_3_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta24_pt1_3_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta24_pt1_3_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta24_pt1_3_hiGeneral->GetNbinsY()-1);
  hMultCent_eta24_pt1_3_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta24_pt1_3_hiGeneral->SetOption("colz");

  // |eta| < 1.0, 0.3 < pT < 3.0
  hMultCent_eta10_pt03_3_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta10_pt03_3_hiPixelAndGeneral", "hMultCent_eta10_pt03_3_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_3_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_3_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_3_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_3_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt03_3_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_3_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta10_pt03_3_hiPixel = fs->make<TH2D>("hMultCent_eta10_pt03_3_hiPixel", "hMultCent_eta10_pt03_3_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_3_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_3_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_3_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_3_hiPixel->GetNbinsY()-1);
  hMultCent_eta10_pt03_3_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_3_hiPixel->SetOption("colz");

  hMultCent_eta10_pt03_3_hiGeneral = fs->make<TH2D>("hMultCent_eta10_pt03_3_hiGeneral", "hMultCent_eta10_pt03_3_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_3_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_3_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_3_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_3_hiGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt03_3_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_3_hiGeneral->SetOption("colz");

  // |eta| < 1.0, 0.3 < pT < 1.0
  hMultCent_eta10_pt03_1_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta10_pt03_1_hiPixelAndGeneral", "hMultCent_eta10_pt03_1_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_1_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_1_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_1_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_1_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt03_1_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_1_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta10_pt03_1_hiPixel = fs->make<TH2D>("hMultCent_eta10_pt03_1_hiPixel", "hMultCent_eta10_pt03_1_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_1_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_1_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_1_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_1_hiPixel->GetNbinsY()-1);
  hMultCent_eta10_pt03_1_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_1_hiPixel->SetOption("colz");

  hMultCent_eta10_pt03_1_hiGeneral = fs->make<TH2D>("hMultCent_eta10_pt03_1_hiGeneral", "hMultCent_eta10_pt03_1_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt03_1_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt03_1_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt03_1_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt03_1_hiGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt03_1_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt03_1_hiGeneral->SetOption("colz");

  // |eta| < 1.0, 1.0 < pT < 3.0
  hMultCent_eta10_pt1_3_hiPixelAndGeneral = fs->make<TH2D>("hMultCent_eta10_pt1_3_hiPixelAndGeneral", "hMultCent_eta10_pt1_3_hiPixelAndGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt1_3_hiPixelAndGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt1_3_hiPixelAndGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt1_3_hiPixelAndGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt1_3_hiPixelAndGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt1_3_hiPixelAndGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt1_3_hiPixelAndGeneral->SetOption("colz");

  hMultCent_eta10_pt1_3_hiPixel = fs->make<TH2D>("hMultCent_eta10_pt1_3_hiPixel", "hMultCent_eta10_pt1_3_hiPixel", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt1_3_hiPixel->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt1_3_hiPixel->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt1_3_hiPixel->GetYaxis()->SetRange(1, hMultCent_eta10_pt1_3_hiPixel->GetNbinsY()-1);
  hMultCent_eta10_pt1_3_hiPixel->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt1_3_hiPixel->SetOption("colz");

  hMultCent_eta10_pt1_3_hiGeneral = fs->make<TH2D>("hMultCent_eta10_pt1_3_hiGeneral", "hMultCent_eta10_pt1_3_hiGeneral", 200, 0, 100, 250, 0, 4000);
  hMultCent_eta10_pt1_3_hiGeneral->GetYaxis()->SetTitle("Event Multiplicity");
  hMultCent_eta10_pt1_3_hiGeneral->GetYaxis()->SetNdivisions(505);
  hMultCent_eta10_pt1_3_hiGeneral->GetYaxis()->SetRange(1, hMultCent_eta10_pt1_3_hiGeneral->GetNbinsY()-1);
  hMultCent_eta10_pt1_3_hiGeneral->GetXaxis()->SetTitle("Event Centrality [%]");
  hMultCent_eta10_pt1_3_hiGeneral->SetOption("colz");


  hCent = fs->make<TH1D>("hCent", "hCent", NCENT, centbinsdefault);
  hCent->GetXaxis()->SetTitle("Centrality [%]");
  */
  hPtRes= fs->make<TH1D>("hPtRes", "hPtRes", 200, 0., 1. ); 

  // hClusVtxComp
  hClusVtxCompInteg = fs->make<TH2D>("hClusVtxCompInteg", "hClusVtxCompInteg", 4000, 0, 100000, 200, 0, 20);
  hClusVtxCompInteg->GetXaxis()->SetTitle("N_{pixels}");
  hClusVtxCompInteg->GetYaxis()->SetTitle("clus-vtx compatibility");
  hClusVtxCompInteg->SetOption("colz");
  /*
  for(int c = 0; c < NCENT; c++){
    hClusVtxComp[c] = fs->make<TH2D>(Form("hClusVtxComp_c%i_%i", cmin[c], cmax[c]), Form("hClusVtxComp_c%i_%i", cmin[c], cmax[c]), 300, 0, 60000, 10000, 0, 20);
    hClusVtxComp[c]->GetXaxis()->SetTitle("N_{pixels}");
    hClusVtxComp[c]->GetYaxis()->SetTitle("clus-vtx compatibility");
    hClusVtxComp[c]->SetOption("colz");
  }
  */
  minpt_  = iConfig.getUntrackedParameter<double>("minpt_",  0.0);
  maxpt_  = iConfig.getUntrackedParameter<double>("maxpt_",  8.0);
  minet_  = iConfig.getUntrackedParameter<double>("minet_",  -1.);
  maxet_  = iConfig.getUntrackedParameter<double>("maxet_",  -1.);
  etaMax_ = iConfig.getUntrackedParameter<double>("etaMax_", 2.4);
  minvz_  = iConfig.getUntrackedParameter<double>("minvz_",  -15.);
  maxvz_  = iConfig.getUntrackedParameter<double>("maxvz_",  15.);
  nvtx_   = iConfig.getUntrackedParameter<int>("nvtx_", 100);

}


MultCentAnalyzer::~MultCentAnalyzer()
{
    
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called to produce and analyze the data  ------------
void
MultCentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
    
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
  }//First event
    
  //
  //Get Centrality
  //
  edm::Handle<int> cbin_;
  iEvent.getByToken(CentralityBinTag_,cbin_);
  int hiBin = *cbin_; //HF tower centrality
  int hiBinHF = hiBin;
  centval = 0.5 * (hiBinHF + 0.5);
  //hCent->Fill(centval);

  //
  //ClusterCompatibility
  //

  iEvent.getByToken(clusterCompatibilityToken_, clusterCompatibilityCollection_);

  double minZ_    = -20.0;
  double maxZ_    = 20.05;
  clusVtxQual     = determineQuality(*clusterCompatibilityCollection_, minZ_, maxZ_);
  double nPxlHits = clusterCompatibilityCollection_->nValidPixelHits();

  //int c = hCent->FindBin(centval)-1;
  //hClusVtxComp[c]->Fill(nPxlHits, clusVtxQual);
  hClusVtxCompInteg->Fill(nPxlHits, clusVtxQual);
  
  //
  //Get Vertex
  //
  iEvent.getByToken(vertexToken_, vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else vzr_sell = -999.9;
        
  vtx = vzr_sell;

  VertexCollection recoV = *vertex_;

  int primaryvtx = 0;
  math::XYZPoint vv1( recoV[primaryvtx].position().x(), recoV[primaryvtx].position().y(), recoV[primaryvtx].position().z() );
  double vz      = recoV[primaryvtx].z();

  bool vSize    = true;
  bool vZaccept = true;
  bool vTrkSize = true;

  if( (int) recoV.size() > nvtx_ )         vSize    = false;
  if( vz < minvz_ || vz > maxvz_ )         vZaccept = false;
  if( recoV[primaryvtx].tracksSize() < 1 ) vTrkSize = false;
    
  if( vSize && vZaccept && vTrkSize ){

    //
    //Tracking part
    //
    iEvent.getByToken(trackToken_, trackCollection_);

    int evtMult_eta24_pt03_3_hiPixelAndGeneral = 0;
    int evtMult_eta24_pt03_1_hiPixelAndGeneral = 0;
    int evtMult_eta24_pt1_3_hiPixelAndGeneral  = 0;
    int evtMult_eta10_pt03_3_hiPixelAndGeneral = 0;
    int evtMult_eta10_pt03_1_hiPixelAndGeneral = 0;
    int evtMult_eta10_pt1_3_hiPixelAndGeneral  = 0;

    int evtMult_eta24_pt03_3_hiPixel = 0;
    int evtMult_eta24_pt03_1_hiPixel = 0;
    int evtMult_eta24_pt1_3_hiPixel  = 0;
    int evtMult_eta10_pt03_3_hiPixel = 0;
    int evtMult_eta10_pt03_1_hiPixel = 0;
    int evtMult_eta10_pt1_3_hiPixel  = 0;

    int evtMult_eta24_pt03_3_hiGeneral = 0;
    int evtMult_eta24_pt03_1_hiGeneral = 0;
    int evtMult_eta24_pt1_3_hiGeneral  = 0;
    int evtMult_eta10_pt03_3_hiGeneral = 0;
    int evtMult_eta10_pt03_1_hiGeneral = 0;
    int evtMult_eta10_pt1_3_hiGeneral  = 0;

    Mult_pt03_3_eta24 = 0;
    Mult_pt03_3_eta10 = 0;

    double vxError = recoV[primaryvtx].xError();
    double vyError = recoV[primaryvtx].yError();
    double vzError = recoV[primaryvtx].zError();

    for(TrackCollection::const_iterator itTrack = trackCollection_->begin();
	itTrack != trackCollection_->end();
	++itTrack) {

      double d0          = -1.* itTrack->dxy(vv1);
      double derror      = sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
      double dz          = itTrack->dz(vv1);
      double dzerror     = sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
      double algoOffline = itTrack->originalAlgo();
      bool pixTrax       = 0;
      int nHits          = itTrack->numberOfValidHits();

      if( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
      if( itTrack->charge() == 0 )                         continue;
      if( itTrack->pt() < minpt_ )                         continue;
      if( itTrack->pt() > maxpt_ )                         continue;
      if( fabs(itTrack->eta()) > etaMax_ )                 continue;

      //-- Pixel Tracks
      if(itTrack->pt() < 2.4 && (nHits==3 || nHits==4 || nHits==5 || nHits==6) ) pixTrax = 1;

      //-- Not Pixel Tracks
      if( !pixTrax ){
	if( nHits < 11 )                                                                                      continue;
	if( fabs( dz/dzerror ) > dzCut_ )                                                                     continue;
	if( fabs( d0/derror ) > d0Cut_ )                                                                      continue;
	if( itTrack->ptError()/itTrack->pt() > ptResCut_ )                                                    continue;
	if( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 )         continue;
	if( itTrack->pt() > 2.4 && !(algoOffline==4 || algoOffline==5 || algoOffline==6 || algoOffline==7 ) ) continue;
      }

      hPtRes->Fill( itTrack->ptError()/itTrack->pt() );
            
      double w = 1.;
      if(usePixelTeff_) w = pxTeff->GetEff(centval, itTrack->pt(), itTrack->eta()); 
      if(w > 0.0) {

	//int ptbin       = hPtSpectrum_24->GetXaxis()->FindBin( itTrack->pt() );
	//double binwidth = hPtSpectrum_24->GetBinWidth( ptbin );
	//if( abs(itTrack->eta()) < 2.4) hPtSpectrum_24->Fill( itTrack->pt(), 1./w/binwidth);
	//if( abs(itTrack->eta()) < 1.0) hPtSpectrum_10->Fill( itTrack->pt(), 1./w/binwidth);

	if( abs(itTrack->eta()) <= 1.0){
	  if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta10_pt03_3_hiPixelAndGeneral++;
	  if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) Mult_pt03_3_eta10++;
	  if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta10_pt03_1_hiPixelAndGeneral++;
	  if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta10_pt1_3_hiPixelAndGeneral++;
	}
	if( abs(itTrack->eta()) <= 2.4){
	  if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta24_pt03_3_hiPixelAndGeneral++;
	  if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) Mult_pt03_3_eta24++;
          if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta24_pt03_1_hiPixelAndGeneral++;
          if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta24_pt1_3_hiPixelAndGeneral++;
	}

	if(pixTrax){
	  if( abs(itTrack->eta()) <= 1.0){
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta10_pt03_3_hiPixel++;
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta10_pt03_1_hiPixel++;
	    if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta10_pt1_3_hiPixel++;
	  }
	  if( abs(itTrack->eta()) <= 2.4){
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta24_pt03_3_hiPixel++;
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta24_pt03_1_hiPixel++;
	    if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta24_pt1_3_hiPixel++;
	  }
	}
	else{
	  if( abs(itTrack->eta()) <= 1.0){
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta10_pt03_3_hiGeneral++;
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta10_pt03_1_hiGeneral++;
	    if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta10_pt1_3_hiGeneral++;
	  }
	  if( abs(itTrack->eta()) <= 2.4){
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 3.0) evtMult_eta24_pt03_3_hiGeneral++;
	    if( itTrack->pt()>= 0.3 && itTrack->pt()< 1.0) evtMult_eta24_pt03_1_hiGeneral++;
	    if( itTrack->pt()>= 1.0 && itTrack->pt()< 3.0) evtMult_eta24_pt1_3_hiGeneral++;
	  }
	}

      }
                
    } //-- End track loop
    /*
    hMultCent_eta24_pt03_3_hiPixelAndGeneral ->Fill(centval, evtMult_eta24_pt03_3_hiPixelAndGeneral);
    hMultCent_eta24_pt03_1_hiPixelAndGeneral ->Fill(centval, evtMult_eta24_pt03_1_hiPixelAndGeneral);
    hMultCent_eta24_pt1_3_hiPixelAndGeneral  ->Fill(centval, evtMult_eta24_pt1_3_hiPixelAndGeneral);
    hMultCent_eta10_pt03_3_hiPixelAndGeneral ->Fill(centval, evtMult_eta10_pt03_3_hiPixelAndGeneral);
    hMultCent_eta10_pt03_1_hiPixelAndGeneral ->Fill(centval, evtMult_eta10_pt03_1_hiPixelAndGeneral);
    hMultCent_eta10_pt1_3_hiPixelAndGeneral  ->Fill(centval, evtMult_eta10_pt1_3_hiPixelAndGeneral);

    hMultCent_eta24_pt03_3_hiPixel ->Fill(centval, evtMult_eta24_pt03_3_hiPixel);
    hMultCent_eta24_pt03_1_hiPixel ->Fill(centval, evtMult_eta24_pt03_1_hiPixel);
    hMultCent_eta24_pt1_3_hiPixel  ->Fill(centval, evtMult_eta24_pt1_3_hiPixel);
    hMultCent_eta10_pt03_3_hiPixel ->Fill(centval, evtMult_eta10_pt03_3_hiPixel);
    hMultCent_eta10_pt03_1_hiPixel ->Fill(centval, evtMult_eta10_pt03_1_hiPixel);
    hMultCent_eta10_pt1_3_hiPixel  ->Fill(centval, evtMult_eta10_pt1_3_hiPixel);

    hMultCent_eta24_pt03_3_hiGeneral ->Fill(centval, evtMult_eta24_pt03_3_hiGeneral);
    hMultCent_eta24_pt03_1_hiGeneral ->Fill(centval, evtMult_eta24_pt03_1_hiGeneral);
    hMultCent_eta24_pt1_3_hiGeneral  ->Fill(centval, evtMult_eta24_pt1_3_hiGeneral);
    hMultCent_eta10_pt03_3_hiGeneral ->Fill(centval, evtMult_eta10_pt03_3_hiGeneral);
    hMultCent_eta10_pt03_1_hiGeneral ->Fill(centval, evtMult_eta10_pt03_1_hiGeneral);
    hMultCent_eta10_pt1_3_hiGeneral  ->Fill(centval, evtMult_eta10_pt1_3_hiGeneral);
    */
    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void
MultCentAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MultCentAnalyzer::endJob() {
}

double
MultCentAnalyzer::determineQuality(const reco::ClusterCompatibility & cc,
                                               double minZ, double maxZ) 
{
  // will compare cluster compatibility at a determined best 
  // z position to + and - 10 cm from the best position
  float best_z = 0.;
  int best_n= 0.,low_n = 0.,high_n = 0.;


  // look for best vertex z position within zMin to zMax range
  // best position is determined by maximum nHit with 
  // chi used for breaking a tie
  int nhits_max = 0;
  double chi_max = 1e+9;
  for( int i=0; i<cc.size(); i++ )
    {
      if( cc.z0(i) > maxZ || cc.z0(i) < minZ ) continue;
      if(cc.nHit(i) == 0) continue;
      if(cc.nHit(i) > nhits_max) {
	chi_max = 1e+9;
	nhits_max = cc.nHit(i);
      }
      if(cc.nHit(i) >= nhits_max && cc.chi(i) < chi_max) {
	chi_max = cc.chi(i);
	best_z = cc.z0(i); best_n = cc.nHit(i);
      }
    }

  // find compatible clusters at + or - 10 cm of the best, 
  // or get as close as possible in terms of z position.
  double low_target = best_z - 10.0;
  double high_target = best_z + 10.0;
  double low_match = 1000., high_match = 1000.;
  for( int i=0; i<cc.size(); i++ )
    {  
      if( fabs(cc.z0(i)-low_target) < low_match )
	{
	  low_n = cc.nHit(i); 
	  low_match = fabs(cc.z0(i)-low_target);
	}
      if( fabs(cc.z0(i)-high_target) < high_match )
	{
	  high_n = cc.nHit(i); 
	  high_match = fabs(cc.z0(i)-high_target);
	}
    }

  // determine vertex compatibility quality score
  double clusVtxQual=0.0;
  if ((low_n+high_n)> 0)
    clusVtxQual = (2.0*best_n)/(low_n+high_n);  // A/B
  else if (best_n > 0)
    clusVtxQual = 1000.0;                      // A/0 (set to arbitrarily large number)
  else
    clusVtxQual = 0;   

  return clusVtxQual;

}



//define this as a plug-in
DEFINE_FWK_MODULE(MultCentAnalyzer);

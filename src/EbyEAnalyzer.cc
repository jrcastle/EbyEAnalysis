// -*- C++ -*-
//
// Package:    EbyEAnalyzer
// Class:      EbyEAnalyzer
//
/**\class EbyEAnalyzer EbyEAnalyzer.cc HeavyIonsAnalysis/EbyEAnalysis/src/EbyEAnalyzer.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  James Castle
//         Created:  Feb. 2017

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/EbyEAnalysis/interface/TrackEfficiency.h"
#include "HeavyIonsAnalysis/TrackingCode/HIRun2015Ana/macros/TrackCorrector3D.h"
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

//
// class decleration
//

class EbyEAnalyzer : public edm::EDAnalyzer {
public:
  explicit EbyEAnalyzer(const edm::ParameterSet&);
  ~EbyEAnalyzer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual double SmearPt(double pt, double smearpct);
  virtual void HFQVectors(const edm::Event& iEvent);
  virtual void GENTrackerQVectors(const edm::Event& iEvent);
  virtual void RECOTrackerQVectors(const edm::Event& iEvent, math::XYZPoint vv1);
  virtual void endJob() ;

  // ----------member data ---------------------------
  unsigned int runno_;

  edm::Service<TFileService> fs;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag genTrackTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genTrackToken_;
  edm::Handle<reco::GenParticleCollection> genTrackCollection_;

  edm::InputTag caloTag_;
  edm::EDGetTokenT<CaloTowerCollection> caloToken_;
  edm::Handle<CaloTowerCollection> caloCollection_;

  edm::EDGetTokenT<reco::Centrality> CentralityTag_;
  edm::EDGetTokenT<int> CentralityBinTag_;

  edm::InputTag HepMCEvtTag_;
  edm::EDGetTokenT<edm::HepMCProduct> HepMCEvtToken_;
  edm::Handle<edm::HepMCProduct> HepMCEvt_ ;

  edm::InputTag V2TrueTag_;
  edm::EDGetTokenT<double> V2TrueToken_;
  edm::Handle<double> V2True_ ;

  bool genAna_;
  bool smearPt_;
  double smearPtPct_;

  double V2True;

  int hiBinHF;
  double centval;
  double centmax_;
  double b;

  int vs_sell;
  float vzr_sell;
  double vtx;

  double vxError;
  double vyError;
  double vzError;

  int nptbinsDefault_;
  std::vector<double> ptbinsDefault_;
  int netabinsDefault_;
  std::vector<double> etabinsDefault_;

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
  double pixChi2NLayerCut_;
  double pixd0Cut_;

  bool useTeff_;
  TrackCorrector3D * teff;

  bool usePixelTeff_;
  TrackEfficiency * pxTeff;

  TTree * tree;

  TRandom3 * ran;

  bool FirstEvent_;
  string effTable_;

  bool Branch_V2True;
  bool Branch_Vtx;
  bool Branch_sumw;
  bool Branch_Run;
  bool Branch_mult;
  bool Branch_HFqn;

  bool Subevent_Standard;

  TH2D * sumw;
  TH2D * sumwqx2;
  TH2D * sumwqy2;
  TH2D * sumwqx3;
  TH2D * sumwqy3;
  TH2D * sumwqx4;
  TH2D * sumwqy4;

  TH2I * mult;

  double qnHFx_EP[NumEPNames];
  double qnHFy_EP[NumEPNames];
  double sumET_EP[NumEPNames];

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
EbyEAnalyzer::EbyEAnalyzer(const edm::ParameterSet& iConfig):
  runno_(0),
  vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag_")) ),
  trackToken_( consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackTag_")) ),
  genTrackToken_( consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genTrackTag_")) ),
  caloToken_( consumes<CaloTowerCollection>(iConfig.getUntrackedParameter<edm::InputTag>("caloTag_")) ),
  CentralityTag_(consumes<reco::Centrality>(iConfig.getUntrackedParameter<edm::InputTag>("CentralityTag_"))),
  CentralityBinTag_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("CentralityBinTag_"))),
  HepMCEvtToken_( consumes<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("HepMCEvtTag_")) ),
  V2TrueToken_( consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("V2TrueTag_")) )
{
  runno_      = 0;
  FirstEvent_ = kTRUE;

  ran = new TRandom3(0);

  vertexTag_          = iConfig.getUntrackedParameter<edm::InputTag>("vertexTag_");
  genAna_             = iConfig.getUntrackedParameter<bool>("GenAna_", false);
  caloTag_            = iConfig.getUntrackedParameter<edm::InputTag>("caloTag_");
  centmax_            = iConfig.getUntrackedParameter<double>("centmax_", 70.);
  useTeff_            = iConfig.getUntrackedParameter<bool>("useTeff_", true);
  usePixelTeff_       = iConfig.getUntrackedParameter<bool>("usePixelTeff_", false);
  effTable_           = iConfig.getUntrackedParameter<std::string>("effTable_");
  trackQualityCuts_   = iConfig.getUntrackedParameter<int>("trackQualityCuts_", 1);
  smearPt_            = iConfig.getUntrackedParameter<bool>("smearPt_", false);
  smearPtPct_         = iConfig.getUntrackedParameter<double>("smearPtPct_", 0.1);
  teff                = 0;
  pxTeff              = 0;

  if( smearPt_ ) std::cout << Form("Smearing pt with a %.2f resolution", smearPtPct_) << std::endl;

  //-- Set genTrack tags
  if(genAna_) genTrackTag_ = iConfig.getUntrackedParameter<edm::InputTag>("genTrackTag_");
  else        trackTag_    = iConfig.getUntrackedParameter<edm::InputTag>("trackTag_");

  //-- Set Track quality cuts
  if( trackQualityCuts_ == 1){
    //-- Nominal Cuts
    dzCut_            = 3.0;
    d0Cut_            = 3.0;
    ptResCut_         = 0.1;
    pixChi2NLayerCut_ = 12;
    pixd0Cut_         = 8;
  }
  if( trackQualityCuts_ == 2){
    //-- Loose Cuts
    dzCut_            = 5.0;
    d0Cut_            = 5.0;
    ptResCut_         = 0.1;
    pixChi2NLayerCut_ = 18;
    pixd0Cut_         = 10;
  }
  if( trackQualityCuts_ == 3){
    //-- Tight Cuts
    dzCut_            = 2.0;
    d0Cut_            = 2.0;
    ptResCut_         = 0.05;
    pixChi2NLayerCut_ = 9;
    pixd0Cut_         = 6;
  }

  //-- Set efficiency table
  if(useTeff_){
    teff = new TrackCorrector3D(effTable_);
    teff->load("HITrackCorrections");
  }
  if( usePixelTeff_ ) pxTeff = new TrackEfficiency(effTable_);

  //-- Choose which branches will be stored in the tree
  Branch_V2True = iConfig.getUntrackedParameter<bool>("Branch_V2True", false);
  Branch_Vtx    = iConfig.getUntrackedParameter<bool>("Branch_Vtx",    true);
  Branch_sumw   = iConfig.getUntrackedParameter<bool>("Branch_sumw",   true);
  Branch_Run    = iConfig.getUntrackedParameter<bool>("Branch_Run",    true);
  Branch_mult   = iConfig.getUntrackedParameter<bool>("Branch_mult",   true);
  Branch_HFqn   = iConfig.getUntrackedParameter<bool>("Branch_HFqn",   true);

  //-- Choose how subevents will be divided
  Subevent_Standard = iConfig.getUntrackedParameter<bool>("Subevent_Standard",true);
  if(Subevent_Standard) std::cout<<"Standard subevent selection will be used"<<std::endl;

  //-- Initilize q-vector and weight sum histograms
  nptbinsDefault_  = iConfig.getUntrackedParameter<int>("nptbinsDefault_");
  ptbinsDefault_   = iConfig.getUntrackedParameter< std::vector<double> >("ptbinsDefault_");
  netabinsDefault_ = iConfig.getUntrackedParameter<int>("netabinsDefault_");
  etabinsDefault_  = iConfig.getUntrackedParameter< std::vector<double> >("etabinsDefault_");

  tree    = fs->make<TTree>("tree", "EP tree");
  sumwqx2 = fs->make<TH2D>( "sumwqx2", "sumwqx2", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumwqy2 = fs->make<TH2D>( "sumwqy2", "sumwqy2", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumwqx3 = fs->make<TH2D>( "sumwqx3", "sumwqx3", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumwqy3 = fs->make<TH2D>( "sumwqy3", "sumwqy3", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumwqx4 = fs->make<TH2D>( "sumwqx4", "sumwqx4", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumwqy4 = fs->make<TH2D>( "sumwqy4", "sumwqy4", nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  sumw    = fs->make<TH2D>( "sumw",    "sumw",    nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );
  mult    = fs->make<TH2I>( "mult",    "mult",    nptbinsDefault_, &(ptbinsDefault_[0]), netabinsDefault_, &(etabinsDefault_[0]) );

  sumwqx2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqx2->GetYaxis()->SetTitle("#eta");
  sumwqx2->SetOption("colz");

  sumwqx3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqx3->GetYaxis()->SetTitle("#eta");
  sumwqx3->SetOption("colz");

  sumwqx4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqx4->GetYaxis()->SetTitle("#eta");
  sumwqx4->SetOption("colz");

  sumwqy2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqy2->GetYaxis()->SetTitle("#eta");
  sumwqy2->SetOption("colz");

  sumwqy3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqy3->GetYaxis()->SetTitle("#eta");
  sumwqy3->SetOption("colz");

  sumwqy4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqy4->GetYaxis()->SetTitle("#eta");
  sumwqy4->SetOption("colz");

  sumw->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumw->GetYaxis()->SetTitle("#eta");
  sumw->SetOption("colz");

  mult->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  mult->GetYaxis()->SetTitle("#eta");
  mult->SetOption("colz");

  //-- Set Tree branches
  if(!genAna_)      tree->Branch("Cent",    &centval, "cent/D");
  else              tree->Branch("b",       &b,       "b/D");
  if(Branch_V2True) tree->Branch("V2True",  &V2True,  "V2True/D");
  if(Branch_Vtx)    tree->Branch("Vtx",     &vtx,     "vtx/D");
  if(Branch_Run)    tree->Branch("Run",     &runno_,  "run/i");
  if(Branch_mult)   tree->Branch("mult",    "TH2I",   &mult,    128000, 0);
  if(Branch_sumw){
    tree->Branch("sumw",    "TH2D",   &sumw,    128000, 0);
    tree->Branch("sumwqx2", "TH2D",   &sumwqx2, 128000, 0);
    tree->Branch("sumwqy2", "TH2D",   &sumwqy2, 128000, 0);
    tree->Branch("sumwqx3", "TH2D",   &sumwqx3, 128000, 0);
    tree->Branch("sumwqy3", "TH2D",   &sumwqy3, 128000, 0);
    tree->Branch("sumwqx4", "TH2D",   &sumwqx4, 128000, 0);
    tree->Branch("sumwqy4", "TH2D",   &sumwqy4, 128000, 0);
  }
  if(Branch_HFqn){
    TString epnames = EPNames[0].data();
    epnames = epnames+"/D";
    for(int i = 0; i<NumEPNames; i++) if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    tree->Branch("qnHFx_EP",  &qnHFx_EP,  epnames.Data() );
    tree->Branch("qnHFy_EP",  &qnHFy_EP,  epnames.Data() );
    tree->Branch("sumET_EP",  &sumET_EP,  epnames.Data() );
  }

  //-- min/max kinematic varibales
  minpt_  = ptbinsDefault_[0];
  maxpt_  = ptbinsDefault_[nptbinsDefault_];
  etaMax_ = etabinsDefault_[netabinsDefault_];

  minet_  = iConfig.getUntrackedParameter<double>("minet_",  -1.);
  maxet_  = iConfig.getUntrackedParameter<double>("maxet_",  -1.);
  minvz_  = iConfig.getUntrackedParameter<double>("minvz_",  -15.);
  maxvz_  = iConfig.getUntrackedParameter<double>("maxvz_",  15.);
  nvtx_   = iConfig.getUntrackedParameter<int>("nvtx_", 100);

}


EbyEAnalyzer::~EbyEAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce and analyze the data  ------------
double
EbyEAnalyzer::SmearPt( double pt, double smearpct )
{
  double ptmin = (1. - smearpct) * pt;
  double ptmax = (1. + smearpct) * pt;
  double smearpt = ran->Uniform(ptmin, ptmax);
  return smearpt;
}

void
EbyEAnalyzer::HFQVectors(const edm::Event& iEvent)
{
  //-- Reset event values
    for(int i = 0; i<NumEPNames; i++) {
    qnHFx_EP[i]  = 0.;
    qnHFy_EP[i]  = 0.;
    sumET_EP[i]  = 0.;
  }

  double tower_eta, tower_phi;
  double tower_energyet, tower_energyet_e, tower_energyet_h;
  iEvent.getByToken(caloToken_,caloCollection_);

  if(caloCollection_.isValid()){
    for (CaloTowerCollection::const_iterator j = caloCollection_->begin();
        j !=caloCollection_->end();
        j++) {

      tower_eta          = j->eta();
      tower_phi          = j->phi();
      tower_energyet_e   = j->emEt();
      tower_energyet_h   = j->hadEt();
      tower_energyet     = tower_energyet_e + tower_energyet_h;
      double minet       = minet_;
      double maxet       = maxet_;

      for(int i = 0; i<NumEPNames; i++) {
        if(minet_<0) minet = minTransverse[i];
        if(maxet_<0) maxet = maxTransverse[i];
        if(tower_energyet<minet) continue;
        if(tower_energyet>maxet) continue;
        double w = tower_energyet;

        if( tower_eta >= EPEtaMin1[i] && tower_eta < EPEtaMax1[i] ){
          qnHFx_EP[i]  += w * TMath::Cos( 2. * tower_phi ); 
          qnHFy_EP[i]  += w * TMath::Sin( 2. * tower_phi ); 
          sumET_EP[i]  += w;
        }

      } //-- End EP loop

    } //-- End calo tower loop

    //int evt = iEvent.id().event();
    //int lumi = iEvent.id().luminosityBlock();
    //double q2x = (qnHFx_EP[0] + qnHFx_EP[1]) / (sumET_EP[0] + sumET_EP[1]);
    //double q2y = (qnHFy_EP[0] + qnHFy_EP[1]) / (sumET_EP[0] + sumET_EP[1]);
    //double q2 = sqrt(q2x*q2x + q2y*q2y);
    //if(q2 > 0.2) std::cout << Form("Run = %i\t", runno_) << Form("Lumi = %i\t", lumi) << Form("Event = %i\t", evt) << Form("q2 = %.2f", q2) << std::endl;

  } //-- End if(caloCollection_.isValid())

}

void
EbyEAnalyzer::GENTrackerQVectors(const edm::Event& iEvent)
{
  using namespace reco;
  iEvent.getByToken(genTrackToken_, genTrackCollection_);
  for (GenParticleCollection::const_iterator itTrack = genTrackCollection_->begin();
      itTrack != genTrackCollection_->end();
      ++itTrack
      ){

    if( itTrack->status() != 1 )         continue;
    if( itTrack->charge() == 0 )         continue;
    if( itTrack->pt() < minpt_ )         continue;
    if( itTrack->pt() > maxpt_ )         continue;
    if( fabs(itTrack->eta()) > etaMax_ ) continue;

    double w = 1.;
    sumwqx2 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(2.*itTrack->phi())/w);
    sumwqy2 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(2.*itTrack->phi())/w);
    sumwqx3 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(3.*itTrack->phi())/w);
    sumwqy3 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(3.*itTrack->phi())/w);
    sumwqx4 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(4.*itTrack->phi())/w);
    sumwqy4 -> Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(4.*itTrack->phi())/w);

    mult->Fill(itTrack->pt(), itTrack->eta(), 1 );
    sumw->Fill(itTrack->pt(), itTrack->eta(), 1./w);

  }

}

void
EbyEAnalyzer::RECOTrackerQVectors(const edm::Event& iEvent, math::XYZPoint vv1)
{
  using namespace reco;
  iEvent.getByToken(trackToken_, trackCollection_);
  for(TrackCollection::const_iterator itTrack = trackCollection_->begin();
      itTrack != trackCollection_->end();
      ++itTrack){

    double d0          = -1.* itTrack->dxy(vv1);
    double derror      = sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
    double dz          = itTrack->dz(vv1);
    double dzerror     = sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
    double algoOffline = itTrack->originalAlgo();
    int nHits          = itTrack->numberOfValidHits();
    int nLayers        = itTrack->hitPattern().trackerLayersWithMeasurement();
    double chi2NDF     = itTrack->normalizedChi2();

    bool pixTrax       = 0;
    bool highPurity    = itTrack->quality(reco::TrackBase::highPurity);
    double tracketa    = itTrack->eta();
    double trackpt     = itTrack->pt();
    double ptErr       = itTrack->ptError();
    int trackCharge    = itTrack->charge();

    if( smearPt_ ) trackpt = SmearPt(trackpt, smearPtPct_);

    if( !highPurity )              continue;
    if( trackCharge == 0 )         continue;
    if( trackpt < minpt_ )         continue;
    if( trackpt > maxpt_ )         continue;
    if( fabs(tracketa) > etaMax_ ) continue;

    if(trackpt < 2.4 && (nHits==3 || nHits==4 || nHits==5 || nHits==6) ) pixTrax = true;

    if( pixTrax ){
      //-- Pixel track cuts
      if( (chi2NDF/(double)nLayers) > pixChi2NLayerCut_) continue;
      if( fabs( dz/dzerror ) > pixd0Cut_)                continue;
    }
    else{
      //-- General track cuts
      if( nHits < 11 )                                                                                continue;
      if( fabs( dz/dzerror ) > dzCut_ )                                                               continue;
      if( fabs( d0/derror ) > d0Cut_ )                                                                continue;
      if( ptErr/trackpt > ptResCut_ )                                                                 continue;
      if( (chi2NDF/(double)nLayers) > 0.15 )                                                          continue;
      if( trackpt > 2.4 && !(algoOffline==4 || algoOffline==5 || algoOffline==6 || algoOffline==7 ) ) continue;
    }

    double w = 1.;
    if(useTeff_) w = teff->getWeight(trackpt, tracketa, hiBinHF);
    if(usePixelTeff_) w = pxTeff->GetEff(centval, trackpt, tracketa);
    if(w == 0.0) continue;
    if(Subevent_Standard){
      sumwqx2 -> Fill(trackpt, tracketa, TMath::Cos(2.*itTrack->phi())/w);
      sumwqy2 -> Fill(trackpt, tracketa, TMath::Sin(2.*itTrack->phi())/w);
      sumwqx3 -> Fill(trackpt, tracketa, TMath::Cos(3.*itTrack->phi())/w);
      sumwqy3 -> Fill(trackpt, tracketa, TMath::Sin(3.*itTrack->phi())/w);
      sumwqx4 -> Fill(trackpt, tracketa, TMath::Cos(4.*itTrack->phi())/w);
      sumwqy4 -> Fill(trackpt, tracketa, TMath::Sin(4.*itTrack->phi())/w);

      mult->Fill(trackpt, tracketa, 1 );
      sumw->Fill(trackpt, tracketa, 1./w);
    }

  } //-- end track loop

}

void
EbyEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  if( !genAna_ ){
    edm::Handle<int> cbin_;
    iEvent.getByToken(CentralityBinTag_,cbin_);
    int hiBin = *cbin_; //HF tower centrality
    hiBinHF = hiBin;
    centval = 0.5 * (hiBinHF + 0.5);
    if(centval > centmax_) return;
  }
  else{
    iEvent.getByToken( HepMCEvtToken_ , HepMCEvt_ ) ;
    HepMC::GenEvent * genevt = (HepMC::GenEvent *)HepMCEvt_->GetEvent();
    HepMC::HeavyIon * hi = genevt->heavy_ion();
    b = hi->impact_parameter();
  }


  if(genAna_){

    if(Branch_V2True){
      iEvent.getByToken(V2TrueToken_, V2True_);
      V2True = *V2True_;
    }

    sumw->Reset();
    sumwqx2->Reset();
    sumwqy2->Reset();
    sumwqx3->Reset();
    sumwqy3->Reset();
    sumwqx4->Reset();
    sumwqy4->Reset();
    mult->Reset();
    GENTrackerQVectors(iEvent);
  }
  else{
    //
    // Get Vertex
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
    vxError = recoV[primaryvtx].xError();
    vyError = recoV[primaryvtx].yError();
    vzError = recoV[primaryvtx].zError();
    double vz = recoV[primaryvtx].z();

    bool vSize    = true;
    bool vZaccept = true;
    bool vTrkSize = true;

    if( (int) recoV.size() > nvtx_ )             vSize    = false;
    if( fabs(vz) < minvz_ || fabs(vz) > maxvz_ ) vZaccept = false;
    if( recoV[primaryvtx].tracksSize() < 1 )     vTrkSize = false;

    if( !vSize || !vZaccept || !vTrkSize ) return;

    //
    // Tracking part
    //
    sumw->Reset();
    sumwqx2->Reset();
    sumwqy2->Reset();
    sumwqx3->Reset();
    sumwqy3->Reset();
    sumwqx4->Reset();
    sumwqy4->Reset();
    mult->Reset();
    RECOTrackerQVectors(iEvent, vv1);
  }

  //
  // Calo part
  //
  if(Branch_HFqn) HFQVectors(iEvent);

  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void
EbyEAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EbyEAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EbyEAnalyzer);

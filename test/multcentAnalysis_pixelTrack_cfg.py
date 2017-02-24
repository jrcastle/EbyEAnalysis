import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("MultCentAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('HeavyIonsAnalysis.EbyEAnalysis.multcentana_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v12', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )  
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/0039105D-4775-E611-898A-F01FAFE0F396.root' ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/26C53197-D777-E611-B7FD-782BCB4FBF56.root', ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/3067D9EB-D977-E611-9626-F01FAFD9CF48.root', ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/A23EB3F2-E177-E611-9638-003048F2E66E.root', ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/F01B4F9A-D777-E611-9A58-F01FAFD9C778.root', ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/F69A28F0-FA7C-E611-90DC-003048F1B960.root'  ##Pixel Reco
#        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias3-PromptReco-AOD/02748005-39A7-E511-8A33-02163E0126D8.root' ##Prompt Reco
  )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("MultCent_NoCCTune.root")
)

# Track Quality Cuts
# ==1 => Nominal
# ==2 => Loose
# ==3 => Tight
process.multcentana.trackQualityCuts_ = cms.untracked.int32(1) #####PAY ATTENTION!!!!!!! 

process.multcentana.usePixelTeff_ = cms.untracked.bool(True) #####PAY ATTENTION!!!!!!!
process.multcentana.effTable_     = cms.untracked.string("EffCorrectionsPixel_TT_pt_0_10_v2.root") ##-- CRAB Jobs

process.multcentana.etaMax_ = cms.untracked.double(1.0)
process.multcentana.minpt_  = cms.untracked.double(0.3)
process.multcentana.maxpt_  = cms.untracked.double(3.0)

process.multcentana.trackTag_ = cms.untracked.InputTag('hiGeneralAndPixelTracks')

# MinBias trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.minbiasHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.minbiasHLT.HLTPaths = [
    "HLT_HIL1MinimumBiasHF1AND_*",
    "HLT_HIL1MinimumBiasHF2AND_*"
]

# Event Selection
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.HIClusterCompatibilityFilter_cfi')
process.clusterCompatibilityFilter.clusterPars  = cms.vdouble(0.0,0.006)
process.clusterCompatibilityFilter.clusterTrunc = cms.double(2.0)

# Pixel ReRECO CC Tune
#process.clusterCompatibilityFilter.clusterTrunc             = cms.double(3.76)
#process.clusterCompatibilityFilter.pixelTune                = cms.untracked.bool(True)
#process.clusterCompatibilityFilter.nhitsLineTrunc           = cms.untracked.int32(1000)
#process.clusterCompatibilityFilter.pixelTuneLineClusterPars = cms.untracked.vdouble(2.14116, 0.000928176)
#process.clusterCompatibilityFilter.pixelTunePolyClusterPars = cms.untracked.vdouble(3.01672, 5.63152e-05, -3.82712e-09, 1.31546e-13, -2.2803e-18, 1.57971e-23)


process.eventSelection = cms.Sequence(
    process.hfCoincFilter3
    + process.primaryVertexFilter
    + process.clusterCompatibilityFilter
)

process.p = cms.Path(process.minbiasHLT*
                     process.eventSelection*
                     process.centralityBin*
                     process.multcentana)

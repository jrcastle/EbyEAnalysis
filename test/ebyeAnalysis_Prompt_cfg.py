import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("EbyEAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('HeavyIonsAnalysis.EbyEAnalysis.ebyeana_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v12', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )  
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/j/jcastle/HIMinimumBias2-PixelTracking/0039105D-4775-E611-898A-F01FAFE0F396.root'
#        'file:/afs/cern.ch/user/q/qwang/work/public/PbPb2015_tracking/PbPb_HIMB3_AOD.root'
#        'file:/afs/cern.ch/user/q/qwang/work/public/PbPb2015_tracking/pixeltracking_1.root'
#        'file:/afs/cern.ch/user/q/qwang/work/public/PbPb2015_tracking/hydjet_pixel.root'
  )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("EbyETree_NewCCTune_3pct.root")
)

# Track Quality Cuts
# ==1 => Nominal
# ==2 => Loose
# ==3 => Tight
process.ebyeana.trackQualityCuts_ = cms.untracked.int32(1) #####PAY ATTENTION!!!!!!! 

process.ebyeana.Subevent_Standard = cms.untracked.bool(True)
process.ebyeana.Subevent_Charge   = cms.untracked.bool(False)

process.ebyeana.useTeff_      = cms.untracked.bool(False)
process.ebyeana.usePixelTeff_ = cms.untracked.bool(True) #####PAY ATTENTION!!!!!!!
process.ebyeana.effTable_     = cms.untracked.string("EffCorrectionsPixel_TT_pt_0_10_v2.root") ##-- CRAB Jobs

process.ebyeana.Branch_Run = cms.untracked.bool(False)

process.ebyeana.etaMax_ = cms.untracked.double(2.4)
process.ebyeana.minpt_  = cms.untracked.double(0.3)
process.ebyeana.maxpt_  = cms.untracked.double(3.0)

process.ebyeana.trackTag_ = cms.untracked.InputTag('hiGeneralAndPixelTracks')

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
#process.clusterCompatibilityFilter.clusterPars = cms.vdouble(0.0,0.006)
#process.clusterCompatibilityFilter.clusterTrunc = cms.double(2.0)

# Pixel ReRECO CC Tune
process.clusterCompatibilityFilter.pixelTune                = cms.untracked.bool(True)

# 0.5%
#process.clusterCompatibilityFilter.pixelTunePolyClusterPars = cms.untracked.vdouble(2.80464, 1.31997e-05, -6.65202e-10, 3.39093e-14, -9.93328e-19, 1.01562e-23)
#process.clusterCompatibilityFilter.pixelTuneLineClusterPars = cms.untracked.vdouble(1.82726, 0.000989945)
#process.clusterCompatibilityFilter.clusterTrunc             = cms.double(3.64)
#process.clusterCompatibilityFilter.nhitsLineTrunc           = cms.untracked.int32(1000)
# 1%
#process.clusterCompatibilityFilter.pixelTunePolyClusterPars = cms.untracked.vdouble(3.00911, 2.10645e-05, -2.63219e-09, 1.25211e-13, -2.53694e-18, 1.88619e-23)
#process.clusterCompatibilityFilter.pixelTuneLineClusterPars = cms.untracked.vdouble(2.04286, 0.000476187)
#process.clusterCompatibilityFilter.clusterTrunc             = cms.double(3.71)
#process.clusterCompatibilityFilter.nhitsLineTrunc           = cms.untracked.int32(2100)
# 2%
#process.clusterCompatibilityFilter.pixelTunePolyClusterPars = cms.untracked.vdouble(3.01672, 5.63152e-05, -3.82712e-09, 1.31546e-13, -2.2803e-18, 1.57971e-23)
#process.clusterCompatibilityFilter.pixelTuneLineClusterPars = cms.untracked.vdouble(2.14116, 0.000928176)
#process.clusterCompatibilityFilter.clusterTrunc             = cms.double(3.76)
#process.clusterCompatibilityFilter.nhitsLineTrunc           = cms.untracked.int32(1000)
# 3%
process.clusterCompatibilityFilter.pixelTunePolyClusterPars = cms.untracked.vdouble(3.12211, 5.04532e-05, -2.46064e-09, 7.08956e-14, -1.22693e-18, 9.17866e-24)
process.clusterCompatibilityFilter.pixelTuneLineClusterPars = cms.untracked.vdouble(2.24092, 0.000929259)
process.clusterCompatibilityFilter.clusterTrunc             = cms.double(3.79)
process.clusterCompatibilityFilter.nhitsLineTrunc           = cms.untracked.int32(1000)


process.eventSelection = cms.Sequence(
    process.hfCoincFilter3
    + process.primaryVertexFilter
    + process.clusterCompatibilityFilter
)

process.p = cms.Path(process.minbiasHLT*
                     process.eventSelection*
                     process.centralityBin*
                     process.ebyeana)

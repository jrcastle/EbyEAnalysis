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
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v14', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )  
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/j/jcastle/HYDJETMBPixel/02E1D6EB-CC90-E611-878D-44A842B46A98.root', ## HYDJET MC 1
        'file:/afs/cern.ch/work/j/jcastle/HYDJETMBPixel/06508579-3691-E611-99DB-549F35AC7DFB.root', ## HYDJET MC 2
        'file:/afs/cern.ch/work/j/jcastle/HYDJETMBPixel/0AABAEBB-3491-E611-A5EC-1418774117C7.root'  ## HYDJET MC 3
  )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("MultCent.root")
)

# Track Quality Cuts
# ==1 => Nominal
# ==2 => Loose
# ==3 => Tight
process.multcentana.trackQualityCuts_ = cms.untracked.int32(1) #####PAY ATTENTION!!!!!!! 

process.multcentana.usePixelTeff_ = cms.untracked.bool(True) #####PAY ATTENTION!!!!!!!
process.multcentana.effTable_     = cms.untracked.string("EffCorrectionsPixel_TT_pt_0_10_v2.root") ##-- CRAB Jobs

process.multcentana.etaMax_ = cms.untracked.double(2.4)
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
process.clusterCompatibilityFilter.clusterTrunc = cms.double(3.8)

process.eventSelection = cms.Sequence(
    process.hfCoincFilter3
    + process.primaryVertexFilter
    + process.clusterCompatibilityFilter
)

process.p = cms.Path(#process.minbiasHLT*
                     #process.eventSelection*
                     process.centralityBin*
                     process.multcentana)

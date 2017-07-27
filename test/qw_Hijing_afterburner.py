import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.HiGenCommon.AfterBurnerGenerator_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("EmptySource")

process.generator = cms.EDFilter("HijingGeneratorFilter",
                                 rotateEventPlane = cms.bool(True),
                                 frame = cms.string('CMS     '),
                                 targ = cms.string('A       '),
                                 izp = cms.int32(82),
                                 bMin = cms.double(0),
                                 izt = cms.int32(82),
                                 proj = cms.string('A       '),
                                 comEnergy = cms.double(5023.0),
                                 iat = cms.int32(208),
                                 bMax = cms.double(15),
                                 iap = cms.int32(208)
                                 )

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('HIJING generator'),
    name = cms.untracked.string('$Source: QW $')
    )

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('Hijing.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

process.AftBurner.modv1 = cms.InputTag("0.0")
#process.AftBurner.modv2 = cms.InputTag("fEllP_8pct_v2")
#process.AftBurner.modv2 = cms.InputTag("pBG_4pct_v2")
process.AftBurner.modv2 = cms.InputTag("0.0")
process.AftBurner.modv3 = cms.InputTag("0.0")
process.AftBurner.fluct_v3 = cms.double(0.0)
#process.AftBurner.modv3 = cms.InputTag("pBG_2pct_v3")
process.AftBurner.modmethod = cms.int32(2) 
#process.AftBurner.modmethod = cms.int32(1) # Turn off non-flow

################ comment this out to run without QWNtrkOfflineProducer
#process.QWGenEvent = cms.EDProducer("QWGenEventProducer",
#                                    trackSrc  = cms.untracked.InputTag("genParticles")
#                                    )
#
#process.vectPhi = cms.EDAnalyzer('QWVectorAnalyzer',
#                                 src = cms.untracked.InputTag("QWGenEvent", "phi"),
#                                 hNbins = cms.untracked.int32(5000),
#                                 hstart = cms.untracked.double(0),
#                                 hend = cms.untracked.double(5000),
#                                 cNbins = cms.untracked.int32(1000),
#                                 cstart = cms.untracked.double(-3.14159265358979323846),
#                                 cend = cms.untracked.double(3.14159265358979323846),
#                                 )
#
#process.vectEta = process.vectPhi.clone(
#		src = cms.untracked.InputTag("QWGenEvent", "eta")
#		)
#
#process.vectPt = process.vectPhi.clone(
#		src = cms.untracked.InputTag("QWGenEvent", "pt")
#		)
#
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('cumu.root')
#                                   )
#
#process.pgen = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.AfterBurner+process.GeneInfo + process.QWGenEvent + process.vectPhi + process.vectEta)
############### comment END
#process.pgen = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.AfterBurner+process.GeneInfo)
process.pgen = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.GeneInfo)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

#################### EbyE ####################
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

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("EbyETree_HIJING_NoFlow.root")
)

# Track Quality Cuts
# ==1 => Nominal
# ==2 => Loose
# ==3 => Tight
process.ebyeana.trackQualityCuts_ = cms.untracked.int32(1) #####PAY ATTENTION!!!!!!! 

process.ebyeana.Subevent_Standard = cms.untracked.bool(True)
process.ebyeana.Subevent_Charge   = cms.untracked.bool(False)

process.ebyeana.useTeff_      = cms.untracked.bool(False)
process.ebyeana.usePixelTeff_ = cms.untracked.bool(False) #####PAY ATTENTION!!!!!!!
process.ebyeana.effTable_     = cms.untracked.string("EffCorrectionsPixel_TT_pt_0_10_v2.root")

process.ebyeana.minvz_ = cms.untracked.double(-15.0)
process.ebyeana.maxvz_ = cms.untracked.double(15.0)

process.ebyeana.Branch_V2True = cms.untracked.bool(False)
process.ebyeana.Branch_Run    = cms.untracked.bool(False)

process.ebyeana.etaMax_     = cms.untracked.double(2.4)  #####PAY ATTENTION!!!!!!! 
process.ebyeana.minpt_      = cms.untracked.double(0.3)  #####PAY ATTENTION!!!!!!! 
process.ebyeana.maxpt_      = cms.untracked.double(3.0 ) #####PAY ATTENTION!!!!!!! 
process.ebyeana.smearPt_    = cms.untracked.bool(False)  #####PAY ATTENTION!!!!!!!
process.ebyeana.smearPtPct_ = cms.untracked.double(0.05)

process.ebyeana.trackTag_    = cms.untracked.InputTag('hiGeneralAndPixelTracks')

process.ebyeana.GenAna_      = cms.untracked.bool(True)  #####PAY ATTENTION!!!!!!!
process.ebyeana.genTrackTag_ = cms.untracked.InputTag('genParticles')

process.ebyeana.overrideCent_ = cms.untracked.int32(1)

##############################################


process.generation_step   = cms.Path(process.pgen)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.ebye_step         = cms.EndPath(process.ebyeana)

process.schedule = cms.Schedule(
    process.generation_step,
    process.RAWSIMoutput_step,
    process.ebye_step
)


for path in process.paths:
                getattr(process,path)._seq = process.generator * getattr(process,path)._seq


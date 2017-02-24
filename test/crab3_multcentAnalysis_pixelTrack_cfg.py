from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB

config.section_('General')
config.General.requestName            = 'EbyEPbPb2015_PixelMB2_MultCentAna_NoCCTune'
config.General.transferOutputs        = True
config.General.transferLogs           = False

config.section_('JobType')
config.JobType.outputFiles            = ['MultCent_NoCCTune.root']
config.JobType.pyCfgParams            = ['noprint']
config.JobType.pluginName             = 'Analysis'
config.JobType.psetName               = 'multcentAnalysis_pixelTrack_cfg.py'
config.JobType.maxJobRuntimeMin       = 1315
config.JobType.inputFiles             = ['EffCorrectionsPixel_TT_pt_0_10_v2.root']

config.section_('Data')
#config.Data.inputDBS                  = 'phys03'
config.Data.inputDataset              = '/HIMinimumBias2/HIRun2015-25Aug2016-v1/AOD'
#config.Data.inputDataset              = '/HIMinimumBias2/qwang-crab_HIMB2_PixelTracking_v3-f359806bdfa543922b20b1cc8503759a/USER'
config.Data.allowNonValidInputDataset = True
config.Data.runRange                  = '262548-263757'
config.Data.unitsPerJob               = 50
config.Data.publication               = False
config.Data.splitting                 = 'LumiBased'
config.Data.outLFNDirBase             = '/store/user/jcastle/EbyEPbPb2015_PixelMB2_MultCentAna_NoCCTune'
#config.Data.lumiMask                  = 'Cert_262548-263757_PromptReco_HICollisions15_JSON_v2.txt'
config.Data.lumiMask                  = 'smallSampleGoldenJSON.txt'
#config.Data.lumiMask                  = 'run263614_JSON.txt'

config.section_('User')

config.section_('Site')
config.Site.whitelist                 = ['T2_US_Vanderbilt', 'T2_US_Purdue']
config.Site.storageSite               = 'T2_US_Vanderbilt'

if __name__ == '__main__':

   from CRABAPI.RawCommand import crabCommand
   from httplib import HTTPException

   # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
   # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
   config.General.workArea = 'crab_projects'

   def submit(config):
       try:
           crabCommand('submit', config = config)
       except HTTPException, hte:
           print hte.headers

   #############################################################################################
   ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
   #############################################################################################
submit(config)

#config.General.requestName = 'EbyEPbPb2015_HIMinimumBias3_Check'
#config.Data.inputDataset = '/HIMinimumBias3/HIRun2015-PromptReco-v1/AOD'
#config.Data.outLFNDirBase = '/store/user/jcastle/EbyEPbPb2015_HIMinimumBias3_Check'
#submit(config)

#config.General.requestName = 'EbyEPbPb2015_HIMinimumBias4_Check'
#config.Data.inputDataset = '/HIMinimumBias4/HIRun2015-PromptReco-v1/AOD'
#config.Data.outLFNDirBase = '/store/user/jcastle/EbyEPbPb2015_HIMinimumBias4_Check'
#submit(config)

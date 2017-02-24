from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB

config.section_('General')
config.General.requestName            = 'EbyEPbPb2015_MC'
config.General.transferOutputs        = True
config.General.transferLogs           = False

config.section_('JobType')
config.JobType.outputFiles            = ['MultCent.root']
config.JobType.pyCfgParams            = ['noprint']
config.JobType.pluginName             = 'Analysis'
config.JobType.psetName               = 'multcentAnalysis_mc_cfg.py'
config.JobType.maxJobRuntimeMin       = 1315
config.JobType.inputFiles             = ['EffCorrectionsPixel_TT_pt_0_10_v2.root']

config.section_('Data')
#config.Data.inputDBS                  = 'phys03'
config.Data.inputDataset              = '/Hydjet_Quenched_MinBias_5020GeV_750/HINPbPbWinter16DR-NoPU_75X_mcRun2_HeavyIon_v2_75X_mcRun2_HeavyIon_v14-v1/AODSIM'
config.Data.allowNonValidInputDataset = True
config.Data.unitsPerJob               = 25
config.Data.publication               = False
config.Data.splitting                 = 'FileBased'
config.Data.outLFNDirBase             = '/store/user/jcastle/EbyEPbPb2015_MC'

config.section_('User')

config.section_('Site')
config.Site.whitelist                 = ['T2_US_Vanderbilt', 'T2_US_Purdue', 'T2_US_MIT']
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

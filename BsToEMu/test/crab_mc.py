from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'mc_RunIISummerUL2018'
config.General.workArea = 'crab_Apr01'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/confFile_mc_cfg.py'

#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 2000
#config.JobType.numCores = 8
config.Data.inputDataset =  '/BsToEMu-pythia8-evtgen/sakumar-crab_bs2EMu_MiniAOD_UL2018-07bb2832fd9cf08ee8da01c42829422a/USER'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10

config.Data.outLFNDirBase = '/store/user/sakumar/bs2emu/RunIISummerUL2018/'
config.Data.publication = False 
config.Site.whitelist = ['T2_IN_TIFR']
config.Site.storageSite = 'T2_IN_TIFR'

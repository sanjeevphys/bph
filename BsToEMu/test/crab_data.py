from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ParkingBPH3_Run2018D_UL2018'
config.General.workArea = 'crab_Apr01'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/confFile_data_cfg.py'

config.Data.inputDBS = 'global'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 2000
#config.JobType.numCores = 8
config.Data.inputDataset =  '/ParkingBPH3/Run2018D-UL2018_MiniAODv2-v1/MINIAOD'

config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

config.Data.outLFNDirBase = '/store/user/sakumar/bs2emu/data/'
config.Data.publication = False 
config.Site.storageSite = 'T2_IN_TIFR'

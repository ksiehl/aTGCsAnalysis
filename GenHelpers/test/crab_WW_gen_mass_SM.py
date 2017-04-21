from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'WW_gen_SM'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/aTGCRun2/CMSSW_7_4_14/src/aTGCsAnalysis/GenHelpers/test/makeMgenWW.py'

config.section_("Data")
config.Data.inputDataset = '/WWToLNuQQ_13TeV-powheg/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

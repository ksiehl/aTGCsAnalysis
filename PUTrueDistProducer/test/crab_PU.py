from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PU_MC'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/RunII_preparation/CMSSW_7_4_7/src/WW-analysis/PUTrueDistProducer/test/PU.py'

config.section_("Data")
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

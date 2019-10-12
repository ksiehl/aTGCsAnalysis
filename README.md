aTGC Analysis
========


This is the analysis code for anomalous triple gauge couplings at 13 TeV using CMSSW framework.

Setup Instructions
------------------

```
# ONLY GAURANTEED TO WORK IN SLC6 AT THE MOMENT
# make sure the architectrure is correct:
export SCRAM_ARCH=slc6_amd64_gcc530

# Setup CMSSW
cmsrel CMSSW_8_0_28
cd CMSSW_8_0_28/src
cmsenv

# Necessary for doing MET corrections, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
git cms-merge-topic -u cms-met:METRecipe_8020_for80Xintegration

# HEEP electron ID, see https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2
# brings in HEEP V70 into VID:
git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout
# for other E/gamma IDs in VID if you wish to have them:
git cms-merge-topic ikrav:egm_id_80X_v3
# we need this for the mva weights which runs in VID regardless if you need it or not:
mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/ 
# we need this for the mva weights which runs in VID regardless if you need it or not:
git clone git@github.com:cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data

# Only if we want to include L1ECAL pre-firing, include L1ECALPrefiring from https://github.com/lathomas/cmssw/tree/L1Prefiring_8_0_32/L1Prefiring/EventWeightProducer and uncomment the lines in the analysis code afterwards. Copy the input histogram in files directory to aTGCsAnalysis/TreeMaker/data.

# Checkout aTGC analysis code
git clone -b 80X git@github.com:ksiehl/aTGCsAnalysis.git

# Compile
scram b -j 10

# PDF variation map may need to be updated
vi aTGCsAnalysis/TreeMaker/plugins/PDFVariationMap.h

# if running from lpc, the eos fuse mount shouldn't be directly accesed, 
# so the pathname must be changed in the analysis python templates
cd aTGCsAnalysis/Common/
vi sedscript.sh
bash sedscript.sh

# Jobs can be submit and retrieved using python scripts
cd test/

# First, it may be necesarry to change the config.Site.storeageSite
# ...as well as possibly uncommenting and setting config.Data.outLFNDirBase
# Also, if running again, it may be prudent to update config.Data.outputDatasetTag
vi templates/template.txt

# setup crab
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -valid 192:00

# execute scripts; <name_of_crabjob> should be changed every time it submits to crab;
# there must always be a second argument; if the second argument is anything other than 'for-real',
# it will go through a dryrun and not actually submit to crab
python submit_jobs.py '<name_of_crabjob>' 'for-real'

# check the status of the jobs
unbuffer sh status_jobs.sh '<name_of_crabjob>' | tee outputfile.log

# I like to choose on a case by case basis what I resubmit, and do it manually,
# so I don't run retreive_jobs.py or addExtendedStates.sh
bash get_jobs.sh 'name_of_crabjob'
bash rename.sh # you'll need to rename the storage location here
bash hadd_files.sh # again, rename storage locaton

# or you can use the retrieve jobs script
python retrieve_jobs.py
bash addExtendedStats.sh # once again, you'll have to change the source directory

# Compile the plotting code
cd aTGCsAnalysis/Common/test/Plotting/
make

# Add weight samples after putting in right luminosity
vi CMSLumi.cpp
cd ..
root -l addWeightSamplesWithPU.cpp+ # again, here, you will have to modify the input directory to match the output destination of the get_jobs.sh script

# If this step is in some script, I can't find it, and it wasn't in the orignial instructions
hadd WW-aTGC_mu.root WW-aTGC_MWW-*_mu.root
hadd WZ-aTGC_mu.root WW-aTGC_MWZ-*_mu.root
hadd WW-aTGC_ele.root WW-aTGC_MWW-*_ele.root
hadd WZ-aTGC_ele.root WW-aTGC_MWZ-*_ele.root

# Draw
cd Plotting
./draw <arguments>

# Allowed options:
  --help                produce help message
  --channel arg         channel to be used: ele or mu
  --CR arg              control region to make plots with
  --output arg          tag for the output directory
  --input arg           name of input directory
  --withMC              use Monte Carlo or not
  --withData            use data or not
  --withSignal          draw signal or not
  --withSystematics     calculate systematics or not. If not statistical uncertainties are calculated and drawn.
  --wantToWriteHists    to write histograms to the local file

control region values: ttbar, WJets, TTBarEnrichedInclusive, TTBarEnrichedBTagVeto, signal, signalWW, signalWZ

An example below makes plots in the ttbar control region in the electron channel with data, Monte-Carlo, signal and no systematics :
./draw --CR ttbar --channel {ele,mu} --output ttbar_CR --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData


```

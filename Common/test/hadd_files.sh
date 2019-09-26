#! /bin/bash

cd /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage

# hadd all of the files that need hadding
hadd WJets_Ht100To200_ele.root WJets_Ht100To200-ext*_ele.root
hadd WJets_Ht200To400_ele.root WJets_Ht200To400-ext*_ele.root
hadd WJets_Ht400To600_ele.root WJets_Ht400To600-ext*_ele.root
hadd WJets_Ht600To800_ele.root WJets_Ht600To800-ext*_ele.root
hadd WJets_Ht800To1200_ele.root WJets_Ht800To1200-ext*_ele.root
hadd WJets_Ht1200To2500_ele.root WJets_Ht1200To2500-ext*_ele.root
hadd WJets_Ht2500ToInf_ele.root WJets_Ht2500ToInf-ext*_ele.root

hadd WW_ele.root WW-ext*_ele.root
hadd data_ele.root data-Run*_ele.root

#################################################################
hadd WJets_Ht100To200_mu.root WJets_Ht100To200-ext*_mu.root
hadd WJets_Ht200To400_mu.root WJets_Ht200To400-ext*_mu.root
hadd WJets_Ht400To600_mu.root WJets_Ht400To600-ext*_mu.root
hadd WJets_Ht600To800_mu.root WJets_Ht600To800-ext*_mu.root
hadd WJets_Ht800To1200_mu.root WJets_Ht800To1200-ext*_mu.root
hadd WJets_Ht1200To2500_mu.root WJets_Ht1200To2500-ext*_mu.root
hadd WJets_Ht2500ToInf_mu.root WJets_Ht2500ToInf-ext*_mu.root

hadd WW_mu.root WW-ext*_mu.root
hadd data_mu.root data-Run*_mu.root

# now remove files no longer needed
# this step can be commented out if you're worried about hadding errors

rm *-ext*_ele.root
rm *-ext*_mu.root

rm data-Run*_ele.root
rm data-Run*_mu.root

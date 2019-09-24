#! /bin/bash

set -e
for file in `cat crablist`
do
	echo "acting on $file today."
	crab getoutput crab_projects/crab_${file}_$1
	cd crab_projects/crab_${file}_${1}/results
	hadd ${file}.root tree_*.root
	rm tree_*.root
	mv ${file}.root /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/
	cd -
done

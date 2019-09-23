#! /bin/bash

set -e
for file in `cat crablist`
do
	echo "acting on $file today."
	crab getoutput crab_projects/crab_$file
	cd crab_projects/crab_$file/results
	hadd ${file}.root tree_*.root
	rm tree_*.root
	mv ${file}.root /afs/cern.ch/work/k/ksiehl/public/ntuple_output_storage/
	cd -
done

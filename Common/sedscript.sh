#/bin/bash 

#XXX:if the source command is used, cp will still be aliased to cp -i, and noclobber will still be set
#replace file:///eos/
#with root://cms-xrd-global.cern.ch//eos/

for x in `ls analysis_*.py`
do
	echo $x
	sed 's_file:///eos/_root://cms-xrd-global.cern.ch//eos/_' < $x > file
	#sed 's_#input = cms.untracked.int32(1000) #test_#input = cms.untracked.int32(100) #test_' < $x > file
	cp file $x
done

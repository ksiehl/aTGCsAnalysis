#!/bin/sh

for file in `ls -d crab_projects/*$1`
do
	echo "=========================STATUS CHECK FOR CRAB JOB: $file======="
	echo "================================================================"
	crab status "$file"
	echo -e "=========================================================\n\n"
done


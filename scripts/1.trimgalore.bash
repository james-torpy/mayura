#!/bin/bash

#make modules to load here

#make directory hierachy
projectname="mayura"
samplename="subset"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
rawDir="$projectDir/raw_files"
resultsDir="$projectDir/results"

#input/output directories
inDir="$rawDir/$samplename"
outType="trimgalore"
outPath="$resultsDir/$samplename.$outType"

#scripts/log directories
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#extension of files to be used
inExt="fastq.gz"


# fetch inFiles:

files=( $(ls $inDir/*.$inExt) )
for inFile in ${files[@]} ;do
	echo -e
	echo The file used is: $file

	uniqueID=`basename $inFile | sed s/_R.*//g`
	outDir=$outPath/$uniqueID/

	mkdir -p $outDir
	
	echo -e
	echo "This is inFile:"
	echo $inFile
	echo -e
	echo "This is the outDir:"
	echo $outDir

	trimgalore_line="trim_galore $inFile --fastqc --length 16 -o $outDir"
	echo -e
	echo The trimgalore_line is:
	echo $trimgalore_line
	echo -e

#submit the job to the cluster:
	qsub -N TRIMGALORE_$uniqueID -wd $logDir -b y -j y -R y -pe smp 1 -V $trimgalore_line

	j=$(($j+2))

done;



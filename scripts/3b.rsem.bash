#!/bin/bash

module load gi/boost/1.53.0

numcores=6

#make directory hierachy
projectname="mayura"
samplename="subset"
genomeName="GRCm38_p4"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"
genomeDir="$homeDir/genomes/$genomeName"

#genome reference file
genome_refFile="$genomeDir/$genomeName"

#input/output directories
inType="star"
outType="rsem"

inPath="$resultsDir/$samplename.$inType"
outPath="$resultsDir/$samplename.$outType"

#scripts/log directory
scriptsDir="$projectDir/scripts"
logDir="$scriptsDir/logs"
mkdir -p $logDir

#fetch directory names for files:
i=0
directory_names=`ls $inPath | grep -v "_trimmed" | grep -v "filters" | grep -v "test"`
for directory in ${directory_names[@]};do
	echo The directory used is: $directory;
	echo -e
	directoriesTotal[i]=$directory
	let i++
done;

#set up conditions to perform analysis on bam files from STAR alignment output
j=0
echo The total number of directories is: ${#directoriesTotal[@]}
echo -e

#set up directories specific to each file being analysed
while [ $j -lt ${#directoriesTotal[@]} ]; do

	unique_inDir=$inPath/${directoriesTotal[$j]}
	inFile=$unique_inDir/Aligned.toTranscriptome.out.bam
	outDir=$outPath/${directoriesTotal[$j]}
	out_filePrefix=$outDir/${directoriesTotal[$j]}
	uniqueID=`basename $unique_inDir`

	 mkdir -p $outDir

	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is the inFile:
	echo $inFile
	echo -e
	echo This is the genome_refFile:
	echo $genome_refFile
	echo -e
	echo This is the out_filePrefix:
	echo $out_filePrefix
	echo -e
	
	#perform analysis on bam files using the reference genome
	rsem_line="rsem-calculate-expression \
	-p $numcores \
	--bam $inFile \
	$genome_refFile \
	$out_filePrefix"
	
	echo This is the rsem_line:
	echo $rsem_line
	echo -e
	
	#submit job with name 'RSEM_count_$samplename' to 10 cluster cores:
	qsub -N RSEM_$uniqueID -hold_jid STAR_$uniqueID -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_line

        j=$(($j+1))

done;



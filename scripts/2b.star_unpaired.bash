#!/bin/bash

module load gi/samtools/1.2

numcores=6

#make directory hierachy
projectname="mayura"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome directory
genome="GRCm38.p4.gencode.pa"

genomeDir="$homeDir/genomes/$genome"

echo -e
echo This is the genomeDir:
echo $genomeDir
echo -e

#input/output types
samplenames=( "subset" )
inType="trimgalore"
outType="star"

#extension of files to be used:
inExt=".fq.gz"

#scripts/logs directory
scriptsPath="$projectDir/scripts"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename.$inType"
	outPath="$resultsDir/$samplename.$outType"
	
	echo This is the inPath:
	echo $inPath
	echo -e

#fetch inFiles:

	files=( $(ls $inPath/**/*$inExt | grep -v unpaired) )
    for inFile in ${files[@]}; do
	
		echo The file used is: $inFile
	
		uniqueID=`basename $inFile | sed s/_trimmed$inExt//`
		outDir=$outPath/$uniqueID/
			
		mkdir -p $outDir

		echo -e
		echo This is the uniqueID:
		echo $uniqueID
		echo -e
		echo This is the outDir:
		echo $outDir
		echo -e

#align reads of input files with STAR, output into .bam files:
		starJobName="star."$uniqueID
        bamJobName="bam."$uniqueID
        sortJobName="sort."$uniqueID

        indexJobName="index."$uniqueID
        indexStatsJobName="indexstats."$uniqueID
        outSam=$outDir"Aligned.out.sam"
        outBam=$outDir"$uniqueID.bam"
        outSortedBam=$outDir"$uniqueID.sorted.bam"

      	star_line="STAR --runMode alignReads \
     --genomeDir $genomeDir \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD \
     --outFilterMultimapNmax 50 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.06 \
     --alignIntronMin 20 \
     --alignIntronMax 1500000 \
     --alignSJoverhangMin 6 \
     --alignSJDBoverhangMin 1 \
     --readFilesIn $inFile \
     --outFileNamePrefix $outDir \
     --runThreadN $numcores \
     --quantMode TranscriptomeSAM \
     --outFilterMatchNmin 45 \
     --outSAMtype BAM SortedByCoordinate \
     --outWigType wiggle \
     --outWigNorm None \
     --limitBAMsortRAM 80000000000"

        echo $star_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
        qsub -N STAR_$uniqueID -hold_jid TRIMGALORE_$UniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $star_line

	done;
done;

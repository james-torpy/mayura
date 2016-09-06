#!/bin/bash

numcores=6

module load gi/boost/1.53.0

#directory hierachy
projectname="mayura"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"

#genome/annotation file directories
genomeName="GRCm38_p4"
annotationName="gencode.vM4.annotation"

genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
annotationFile="$genomeDir/$annotationName.gtf"
outDir="$genomeDir"

mkdir -p $outDir

outFile="$outDir/$genomeName"

#scripts/log directory
scriptsPath="$projectDir/scripts"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir:
echo $outDir
echo -e

#Set up conditions
rsem_ref_line="rsem-prepare-reference --gtf $annotationFile $genomeFile $outFile"

echo This is the rsem_ref_line:
echo $rsem_ref_line

qsub -N RSEM_ref_$genome -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_ref_line

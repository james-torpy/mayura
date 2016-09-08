#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Download the reads
#
dx-download-all-inputs --except genomeindex

#
# Fetch and uncompress the genome
#
mkdir ./genomeIndex
dx cat "$genomeindex" | tar xvf - -C ./genomeIndex

STAR --runMode alignReads \
     --genomeDir ./genomeIndex/ \
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
     --readFilesIn ./in/reads/* \
     --runThreadN `nproc`  \
     --quantMode TranscriptomeSAM \
     --outFilterMatchNmin 45 \
     --outSAMtype BAM SortedByCoordinate \
     --outWigType wiggle \
     --outWigNorm None \
     --limitBAMsortRAM 80000000000

#
# Upload results
#
name="$reads_prefix"
name="${name%_1}"
if [ "$sample_name" != "" ]; then
  name="$sample_name"
fi

mkdir -p ~/out/transcriptome/star
mkdir -p ~/out/genome/star
mkdir -p ~/out/logs/logs

mv Aligned.toTranscriptome.out.bam ~/out/transcriptome/star/"$name".transcriptome.bam
mv Aligned.sortedByCoord.out.bam ~/out/genome/star/"$name".genome.sorted.bam
mv Log.final.out ~/out/logs/"$name".log

dx-upload-all-outputs

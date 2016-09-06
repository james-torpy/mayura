#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# Download the inputs (./in/genome_fastagz/* and optionally ./in/genes_gtf/*):
dx-download-all-inputs --parallel
gunzip ./in/genome_fastagz/*

# Set up options:
if [ -d ./in/genes_gtf ]
then
  opts="--sjdbGTFfile ./in/genes_gtf/* --sjdbOverhang $overhang"
fi

mkdir ./genome
STAR --runMode genomeGenerate --genomeDir ./genome/ --genomeFastaFiles ./in/genome_fastagz/* --runThreadN `nproc` $opts $extra_options

# Upload results:
cd ./genome
mkdir -p ~/out/index_tar
tar cf ~/out/index_tar/"$genome_fastagz_prefix".star-index.tar *

dx-upload-all-outputs

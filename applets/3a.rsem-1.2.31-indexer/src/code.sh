#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Download the inputs (./in/genome_fastagz/* and ./in/genes_gtf/*)
#
dx-download-all-inputs --parallel
gunzip ./in/genome_fastagz/*

#
# Run RSEM
#
mkdir ./genome
rsem-prepare-reference $extra_options --gtf ./in/genes_gtf/* ./in/genome_fastagz/* ./genome/genome

#
# Upload results
#
cd ./genome
mkdir -p ~/out/index_tar
tar cf ~/out/index_tar/"$genome_fastagz_prefix"-"$genes_gtf_prefix".rsem-index.tar *
dx-upload-all-outputs


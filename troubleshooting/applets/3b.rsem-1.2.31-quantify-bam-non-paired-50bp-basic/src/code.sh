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
mkdir ./genome
dx cat "$genomeindex" | tar xvf - -C ./genome
genome=$(ls ./genome/*.grp)
genome="${genome%.grp}"

#
# Run RSEM
#
mkdir -p ~/out/results/rsem/basic
rsem-calculate-expression --bam -p `nproc` --no-bam-output ./in/name_sorted_bam/* "$genome" out/results/rsem/basic/"$name_sorted_bam_prefix"


# Upload results
#
dx-upload-all-outputs

rm -r ~/out/results/rsem/basic/*.stat
#

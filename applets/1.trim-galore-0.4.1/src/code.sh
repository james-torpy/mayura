#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Download the reads
#
dx-download-all-inputs --parallel

#
# Run trim_galore
#

# Requires fastqc and cutadapt in the path
PATH="/FastQC/:$HOME/.local/bin:$PATH"

mkdir results
trim_galore ./in/reads/* $extra_options -o results

#
# Upload results
#
name="$reads_prefix"
name="${name%_1}"
name="${name%_R1}"
if [ "$sample_name" != "" ]; then
  name="$sample_name"
fi

mkdir -p ~/out/reads/trimgalore ~/out/others/trimgalore/other_files
mv results/*.fq.gz ~/out/reads/trimgalore/"$name".fq.gz
mv results/* ~/out/others/trimgalore/other_files

dx-upload-all-outputs --parallel

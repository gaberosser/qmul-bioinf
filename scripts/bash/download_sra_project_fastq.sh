#!/usr/bin/env bash

# Call with a single variable, the project ID (e.g. PRJNAxxxx or SRPxxxx)
# This downloads the data in parallel and extracts the fastq files to the current directory

PROJ_ID="$1"

parallel -j 4 'fastq-dump --split-files {}' ::: $(esearch -db sra -query $1 | efetch --format runinfo | cut -d ',' -f 1 | grep 'SRR')
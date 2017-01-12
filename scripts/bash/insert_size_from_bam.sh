#!/bin/bash
# Compute the mean and std of the aligned insert from a bam file
# Input: EITHER path to a BAM or SAM file OR path to a directory containing BAM files

# SAM flag 66 = 64 + 2: "read mapped in proper pair" AND "first in pair"
# Extract the 9th column = TLEN (length of aligned section)

# Define number of reads to use. The estimate is unstable if this number is too small
n=1000000

function readone {
    # get filename
    filename=$(basename "$1")

    res=`
    samtools view -f67 $1 |
    cut -f 9 |
    head -n $n |
    sed 's/^-//' |
    awk '{sum+=$1;sum2+=($1*$1)} END {print sum/NR,sqrt((sum2 - (sum*sum) / NR) / (NR - 1))}' |
    tr ' ' ','
    `

    echo "${filename},${res}"
}



if [[ -d $1 ]]; then
    # iterate over all files in directory
    echo "file,mean,stdev"
    filelist=`ls $1/*.{bam,BAM,sam,SAM} 2>/dev/null`
    for f in $filelist; do
        readone $f
    done
elif [[ -f $1 ]]; then
    echo "file,mean,stdev"
    readone $f
else
    echo "$1 is not a recognised file or directory"
fi
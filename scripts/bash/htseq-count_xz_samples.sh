#!/bin/bash

indir="/media/gabriel/raid1_4tb/data/Xinyu_Zhang_NEB_mRNASeq_GC-XZ-2499_270115-19825815/alignment_and_counting/"

declare -a arr=(
"XZ-3-21077646/alignments/XZ-3.alignments.sorted.bam"
"XZ-4-21099422/alignments/XZ-4.alignments.sorted.bam"
"XZ-5-21101414/alignments/XZ-5.alignments.sorted.bam"
"XZ-6-21077644/alignments/XZ-6.alignments.sorted.bam"
"XZ-7-21099419/alignments/XZ-7.alignments.sorted.bam"
"XZ-8-21108311/alignments/XZ-8.alignments.sorted.bam"
)

for i in "${arr[@]}"
do
    echo "$i"
    outfn=$(expr substr $i 1 4).count
    echo "$outfn"
    python ~/bioinf/rnaseq/htseq_count.py -d -f bam -s reverse "$indir$i" ~/bioinf/data/ensemble_human_genome/gr37/htseq.features_counts.chr.dill > $outfn &
done
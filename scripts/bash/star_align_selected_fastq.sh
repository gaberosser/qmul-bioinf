#!/bin/bash
OUTDIR="star_alignment_mouse"
for VAR in WTCHG_339858_219127 WTCHG_339858_220139 WTCHG_339858_202113
do
CMD="STAR --runThreadN 12 --genomeDir /home/gabriel/local_data/reference_genomes/mouse/ensembl/GRCm38.p5/star/ --readFilesIn ./${VAR}_1.fastq.gz ./${VAR}_2.fastq.gz --outFileNamePrefix ${OUTDIR}/${VAR} --outSAMtype ""BAM Unsorted SortedByCoordinate"" --outSAMstrandField intronMotif --quantMode GeneCounts --readFilesCommand zcat 2> ${OUTDIR}/${VAR}.stderr"
echo $CMD
eval $CMD
done


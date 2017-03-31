#!/bin/bash
ALIGN_PY="/home/gabriel/bioinf/scripts/star_alignment.py"
OUTDIR="star_alignment"
SRADIR="sra"
REF="/home/gabriel/data/reference_genomes/ensembl/GRCh38/star"
VENV="/home/gabriel/bioinf/venv"

source "${VENV}/bin/activate"

for i in $SRADIR/*.sra; do
	CMD="fastq-dump --split-files --outdir ./ $i"
	echo $CMD
	eval $CMD
done

CMD="${ALIGN_PY} ./ $REF $OUTDIR --runThreadN 8"

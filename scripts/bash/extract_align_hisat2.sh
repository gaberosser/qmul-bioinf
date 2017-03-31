#!/bin/bash
ALIGN_PY="/home/gabriel/bioinf/scripts/hisat2_alignment.py"
REF="/home/gabriel/data/reference_genomes/ensembl/GRCh38/ht2/Homo_sapiens.GRCh38.dna.primary_assembly"
VENV="/home/gabriel/bioinf/venv"

source "${VENV}/bin/activate"

for i in *.sra; do
	bam=`echo $i | sed 's/sra$/bam/'`
	if [ -s $bam ]; then
		echo "File named ${bam} already exists. Skipping"
	else
		# extract fastq files: will create two files, one per read
		CMD="fastq-dump --split-files --outdir ./ $i"
		echo $CMD
		eval $CMD
		# align them, convert to bam, delete sam
		CMD="python ${ALIGN_PY} ./ $REF ./ -p 12 --no-discordant"
		echo $CMD
		# eval $CMD && rm *.fastq  # rm only run if the previous command executes successfully
		eval $CMD
	fi
done

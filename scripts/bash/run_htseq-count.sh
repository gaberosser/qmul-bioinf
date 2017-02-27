#!/usr/bin/env bash

# $1: input directory. All bam files found will be run. The BAM files must be SORTED by coordinate.
# $2: gtf/gff reference.
# $3: output directory. Count files, logs and stderr will be written here. Created if doesn't exist.
# $4 onwards: sent to htseq-count

called="$@"

indir="$1"
shift

if [ ! -f $1 ]; then
    echo "Error: cannot find GTF file $2."
    exit 1
fi
gtf="$1"
shift

if [ ! -d $1 ]; then
	echo "Creating output path $1"
	mkdir "$1"
fi
outdir="$1"
shift

logfile="${outdir}/htseq-count.log"
echo "All logging goes to $logfile."

echo "$called" > $logfile

for f in $indir/*.{b,s}am; do
    if [ -e $f ]; then
        bn=`basename $f | sed 's/\.[bs]am$//'`
        ext=`echo $f | sed 's/.*\.\([bs]am\)$/\1/'`
        stderr="$outdir/${bn}.stderr"
        stdout="$outdir/${bn}.count"
        CMD="htseq-count -f $ext $@ $f $gtf 1> $stdout 2> $stderr"
        echo $CMD > $logfile
        echo $CMD
        eval $CMD
    fi
done

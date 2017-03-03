#!/bin/bash

# the first input argument is the directory to operate upon if supplied
# otherwise, use the current directory
echo "$@"

# first input arg is the dir containing the input bams

d="$1"
echo "Called `basename $0` with input directory $d"
shift

if [ ! -d $d ]; then
	echo "Error: $d is not a directory"
	exit 1
fi

outd="$1"
shift

if [ ! -d $outd ]; then
	echo "Creating output dir $outd"
	mkdir $outd
fi

for i in $d/*.bam; do
	bn=`basename $i | sed 's/\.bam$//'`
	outf="$outd/$bn"
	# piping stderr results in a huge file (app emits incremental updates)
	# stderr="$outd/${bn}.stderr"
	CMD="cufflinks $@ -o $outf $i"
	echo $CMD
	eval $CMD
done

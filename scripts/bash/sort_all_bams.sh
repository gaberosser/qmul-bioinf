#!/bin/bash

# the first input argument is the directory to operate upon if supplied
# otherwise, use the current directory
echo "$@"

if [ ! -n "$1" ]; then
	echo "No input dir supplied, using current dir."
	d="./"
else
	d="$1"
	echo "Called `basename $0` with input directory $d"
fi

if [ ! -d $d ]; then
	echo "Error: $d is not a directory"
	exit 1
fi

sortd="${d}/sorted/"
if [ ! -d $sortd ]; then
	echo "Creating output directory ${sortd}."
	mkdir $sortd
fi

# pop the first input arg
shift

for i in $d/*.bam; do
	outf="$sortd/`basename $i`"
	echo "Sorting input file $i -> ${outf}."
	CMD="samtools sort -o ${outf} $i $@"
	echo $CMD
	eval $CMD
done

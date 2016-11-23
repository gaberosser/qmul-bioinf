#!/bin/bash

# get filename and extension
filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"
dirname=$(dirname "$1")

# output filename
outfn="$dirname/$filename.sorted.$extension"

echo "Sorting file $1 by name and outputting to $outfn..."
samtools sort -n --threads 4 -o "$outfn" "$1"

#!/usr/bin/env bash

# This kind of script is useful for selectively symlinking fastq files when samples from different species have been
# sequenced in the same run.

# set cutoff based onthe dividing point - may need two?

CUTOFF="267171"

# assume we are inside the destination directory (e.g. 'human') and fastq files reside one level up

for i in ../*.fastq.gz;
    # this stem needs to be changed based on the filenames
    do t=$(echo $i | sed 's/..\/WTCHG_430074_\([0-9]*\)_.*/\1/');
    if (("$t" > "$CUTOFF"));
        then CMD="ln -s $i ./";
        echo $CMD;
        eval $CMD;
    fi;
done
#!/usr/bin/env bash
# $1: the 'completed' file
# $2: the 'submitted' file
# $3: the params file
echo $1
echo $2

lines=$(comm -13 <(sort $1) <(sort $2) | cut -d , -f 1)

if [ -z "$lines" ]; then
    echo "All submitted IDs were completed. Good job!"
else
    echo "Found some runs that were not completed. Parameters follow."
    for i in $lines; do
        sed -n "$i p" $3
    done
fi

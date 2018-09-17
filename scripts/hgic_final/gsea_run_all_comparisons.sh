#!/usr/bin/env bash

# Update this to point to the relevant GMX file
GMX_FILE="$HOME/Dropbox/research/qmul/hGIC_project/current/core_pipeline/rnaseq/gsea/input_data/msigdb.v6.1.symbols.gmt"

# Update this to point to the GSEA java file
GSEA_PATH="/opt/gsea/gsea-3.0.jar"


runGsea() {
    b=$(basename $1 .params)

    if [[ $b = *"nsc"* ]]; then
        c="#GBM_versus_NSC";
    else
        c="#GBM_versus_iNSC";
    fi;

    java -Xmx16192m -cp $GSEA_PATH xtools.gsea.Gsea \
    -res ${b}.gct -cls ${b}.cls${c} -gmx GMX_FILE \
    -out ./${b} \
    -param_file ${b}.params

    STATUS=$?

    if [ $STATUS == 0 ]; then
        echo "SUCCESS: ${b}"
        mv ${b}.gct ${b}.cls ${b}.params ${b}/
    else
        echo "FAILED: ${b}"
    fi
}

# export so that the parallel env has the function in scope
export -f runGsea

parallel -j 3 runGsea ::: *.params
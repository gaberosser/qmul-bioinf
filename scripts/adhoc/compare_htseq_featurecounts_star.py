from matplotlib import pyplot as plt
import pandas as pd
import os
from settings import DATA_DIR_NON_GIT
from load_data import rnaseq_data


if __name__ == "__main__":
    INDIR = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE73721')

    # featureCounts with hisat2
    fcount_files = [os.path.join(INDIR, 'hisat2_alignment', 'featureCounts_unstr', 'counts.txt')]
    metafile = [os.path.join(INDIR, 'sources.csv')]
    samples = ['SRR2557089']
    dat_h2_fc, met_h2_fc = rnaseq_data.featurecounts(
        fcount_files,
        [metafile],
        samples=samples,
        annotate_by='Approved Symbol'
    )

    # htseq-count with hisat2
    htseq_files = [os.path.join(INDIR, 'hisat2_alignment', 'htseq-count', 'SRR2557089.counts')]
    dat_ht2_ht, met_h2_ht = rnaseq_data.htseqcounts(
        count_files=htseq_files,
        metafile=metafile,
        annotate_by='Approved Symbol',
    )
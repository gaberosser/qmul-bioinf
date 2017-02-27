from matplotlib import pyplot as plt
import pandas as pd
import os
from settings import DATA_DIR_NON_GIT
from load_data import rnaseq_data


if __name__ == "__main__":
    INDIR = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE73721')

    # featureCounts with hisat2
    fcount_files = [os.path.join(INDIR, 'hisat2_alignment', 'featureCounts_unstr', 'counts.txt')]
    metafile = os.path.join(INDIR, 'sources.csv')
    samples = ['SRR2557090']
    dat_h2_fc, met_h2_fc = rnaseq_data.featurecounts(
        fcount_files,
        [metafile],
        samples=samples,
        annotate_by='Approved Symbol'
    )

    # htseq-count with hisat2
    htseq_files = [os.path.join(INDIR, 'hisat2_alignment', 'htseq-count', 'SRR2557090.count')]
    dat_h2_ht, met_h2_ht = rnaseq_data.htseqcounts(
        count_files=htseq_files,
        metafile=metafile,
        annotate_by='Approved Symbol',
    )

    # count by star
    starcount_files = [os.path.join(INDIR, 'star_alignment', 'SRR2557090ReadsPerGene.out.tab')]
    dat_star = pd.read_csv(starcount_files[0], header=None, index_col=0, sep='\t').iloc[:, 0]  # unstranded
    dat_star = rnaseq_data.annotate(dat_star, annotate_by='Approved Symbol')

    all_dat = pd.concat((dat_h2_fc, dat_h2_ht, dat_star), axis=1)
    all_dat.columns = ['fc', 'htseq', 'STAR']

    all_dat.corr()


    ## WTCHG run 1
    INDIR = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704/161219_K00198_0151_BHGYHTBBXX')
    fcount_files = [os.path.join(INDIR, 'featureCounts', 'counts.txt')]
    metafile = os.path.join(INDIR, 'sources.csv')

    dat_h2_fc, met_h2_fc = rnaseq_data.featurecounts(
        fcount_files,
        [metafile],
        samples=['GBM018'],
        annotate_by='Approved Symbol'
    )


    # htseq-count: INCOMPLETE

    # count by star
    starcount_files = [os.path.join(INDIR, 'star_alignment', 'WTCHG_338493_201101ReadsPerGene.out.tab')]
    dat_star = pd.read_csv(starcount_files[0], header=None, index_col=0, sep='\t').iloc[:, 2]  # rev stranded
    dat_star = rnaseq_data.annotate(dat_star, annotate_by='Approved Symbol')

    all_dat = pd.concat((dat_h2_fc, dat_star), axis=1)
    all_dat.columns = ['fc', 'STAR']
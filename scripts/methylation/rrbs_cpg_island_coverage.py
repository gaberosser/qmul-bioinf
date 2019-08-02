import sys
import os
import pandas as pd
import numpy as np
from glob import glob
from settings import DATA_DIR, GIT_LFS_DATA_DIR
from utils import genomics, output, log
from matplotlib import pyplot as plt
import seaborn as sns

logger = log.get_console_logger(__name__)


def plot_one_hist(dat, ax, *args, **kwargs):
    mval = dat[:, 0] / dat.sum(axis=1).astype(float) * 100.
    ax.hist(mval, *args, **kwargs)


if __name__ == "__main__":

    min_coverage = 10
    outdir = output.unique_output_dir("rrbs_methylation_cpg_islands", reuse_empty=True)

    cpg_island_tsv = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.tsv')
    cpg_regions = pd.read_csv(cpg_island_tsv, sep='\t', header=0)

    indir = os.path.join(DATA_DIR, 'rrbseq', 'GC-CV-7163', 'trim_galore', 'mouse', 'bismark')
    subdir = "GC-CV-7163-{i}_S{i}"
    flist = glob(os.path.join(indir, "*.bismark.cov.gz"))

    chroms = [str(t) for t in range(1, 20)]

    chrom_lengths = genomics.reference_genome_chrom_lengths(tax_id=10090)
    # discard unplaced scaffolds, MT, X, Y
    chrom_lengths = chrom_lengths.loc[chroms]

    res_cpg = {}
    res_outside = {}

    for f in flist:
        fstem = os.path.split(f)[-1].replace(".bismark.cov.gz", "")

        logger.info("Starting analysis of %s.", f)
        dat = pd.read_csv(f, sep='\t', header=None, index_col=None)
        dat.columns = ['chr', 'coord', 'end', 'meth_pct', 'n_m', 'n_u']
        # strange behaviour: some chrom names are int, others string
        # fix now by casting all to string
        dat.loc[:, 'chr'] = dat.chr.astype(str)

        # pick only the data and regions we are interested in
        dat = dat.loc[dat.chr.isin(chroms)]

        # require a minimum coverage
        dat = dat.loc[(dat.n_u + dat.n_m) >= min_coverage]

        in_cpg = []
        not_in_cpg = []

        for c, the_dat in dat.groupby('chr'):
            logger.info("Chromosome %s" % c)
            the_dat.insert(the_dat.shape[1], 'in_cpg_island', False)
            coord = the_dat.coord.values
            idx = np.zeros(len(coord)).astype(bool)

            the_cpg = cpg_regions.loc[cpg_regions.loc[:, 'chrom'] == c]

            for i, row in the_cpg.iterrows():
                st = row.chromStart
                en = row.chromEnd
                idx[(coord >= st) & (coord <= en)] = True

            in_cpg.append(the_dat.loc[idx, ['n_m', 'n_u']].values)
            not_in_cpg.append(the_dat.loc[~idx, ['n_m', 'n_u']].values)

        res_cpg[fstem] = np.concatenate(in_cpg, axis=0)
        res_outside[fstem] = np.concatenate(not_in_cpg, axis=0)

    for f in flist:
        fstem = os.path.split(f)[-1].replace(".bismark.cov.gz", "")
        outfn = os.path.join(outdir, "%s_m_values.png" % fstem)

        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
        plot_one_hist(res_cpg[fstem], axs[0], 50)
        # axs[0].hist(res_cpg[fstem], 50)
        axs[0].set_title("CpG Islands")

        plot_one_hist(res_outside[fstem], axs[1], 50)
        # axs[1].hist(res_outside[fstem], 50)
        axs[1].set_title("Non-CpG Islands")
        axs[0].set_xlabel('Percent methylated')
        axs[0].set_ylabel('Frequency')
        axs[1].set_ylabel('Frequency')
        fig.tight_layout()
        fig.savefig(outfn)
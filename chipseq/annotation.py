import pandas as pd
import numpy as np
from . import feature_enrichment
from . import loader
from settings import LOCAL_DATA_DIR
import os


if __name__ == '__main__':
    # recreate the essential operation of Homer suite's annotatePeaks.pl
    peaks = loader.load_macs2_by_patient('050', run_type='default')
    dat = peaks.data.loc[:, 'GBM050H3K27ac'].dropna()
    fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.release87',
        'gtf',
        'Homo_sapiens.GRCh38.87.gtf.gz'
    )

    reg, names = feature_enrichment.get_gene_tss_from_gtf(fn, distance=0)
    chroms = [str(t) for t in range(23)]
    closest_tss = {}
    for c in chroms:
        this_tss = [t for t in names if t[0] == c]
        this_tss_loc = np.array([t[1] for t in this_tss])
        closest_tss[c] = {}
        this_dat = dat.loc[dat.chrom == c]
        search_ix = np.searchsorted(this_tss_loc, this_dat.start + this_dat.rel_peak_pos, side='left')
        ix = search_ix
        ix[ix == this_tss_loc.shape[0]] -= 1
        tss_coords = this_tss_loc[ix]

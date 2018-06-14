"""
The HipSci open access data are available in zip files, each containing a collection of .txt files.
These give intensity values, beta values and detection P values for each probe, as exported from the
GenomeStudio software.
We also have a sources.csv metadata file.
This script loads the files, filters out probes that have a high detection P value, then assembles them into a
single pandas DataFrame.
Finally, the dataframe is exported to a gzipped CSV file.
"""

import os
import pandas as pd
from settings import DATA_DIR_NON_GIT
from glob import glob
import re
import numpy as np
from utils import output


if __name__ == "__main__":
    basedir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'hipsci_ipsc')
    indir = os.path.join(basedir, 'raw')
    meta_fn = os.path.join(basedir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    meta.index = meta.index.str.replace('.', '-')

    flist = glob(os.path.join(indir, "*.txt"))
    raw = {}
    snames = []
    array_types = {}
    for fn in flist:
        sn = re.sub(r"\.mtarray.*", "", os.path.split(fn)[1])
        if sn not in meta.index:
            print "Sample %s is not in the meta file. Skipping" % sn
            continue

        if 'HumanMethylation450' in fn:
            at = '450K'
        elif 'EPIC' in fn:
            at = 'EPIC'
        else:
            print "Unknown array type for file %s" % fn
            at = None
        this_dat = pd.read_csv(fn, sep='\t', header=0, index_col=1)
        pval_ix = this_dat.columns.str.lower().str.contains('pval')
        beta_ix = this_dat.columns.str.lower().str.contains('beta')
        if pval_ix.sum() == 0:
            print "Sample %s has no pval col" % sn
        elif pval_ix.sum() == 1:
            pvals = this_dat.loc[:, pval_ix].squeeze()
            this_dat = this_dat.loc[pvals <= 0.01]
        else:
            print "Sample %s has >1 pval col" % sn

        if beta_ix.sum() != 1:
            print "Sample %s has no beta col!"
        else:
            this_dat = this_dat.loc[:, beta_ix].squeeze()
            this_dat.name = sn
            this_dat.index.name = 'probe_id'
            snames.append(sn)
            raw[sn] = this_dat
            array_types[sn] = at

    # amass all results and leave NaN in place for now
    # null values seem to be indicated by a single space (great...)
    res = pd.DataFrame(raw).replace(' ', np.nan).astype(float)

    # reorder array_types accordingly
    array_types = pd.Series(array_types).loc[res.columns]

    # separate by array type
    epic_res = res.loc[:, array_types == 'EPIC'].dropna()
    four50k_res = res.loc[:, array_types == '450K'].dropna()
    print "EPIC arrays: we retain %d probes" % epic_res.shape[0]
    print "450K arrays: we retain %d probes" % four50k_res.shape[0]

    # output all results, dropping NaN
    all_res = res.dropna()

    outdir = output.unique_output_dir("hipsci_methylation")
    all_res.to_csv(os.path.join(outdir, 'beta_raw.csv.gz'), compression='gzip')

    # output split results
    subdir = os.path.join(outdir, 'EPIC')
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    epic_res.to_csv(os.path.join(subdir, 'beta_raw.csv.gz'), compression='gzip')

    subdir = os.path.join(outdir, '450K')
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    four50k_res.to_csv(os.path.join(subdir, 'beta_raw.csv.gz'), compression='gzip')
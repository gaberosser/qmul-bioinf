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
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'hipsci_ipsc', 'raw')
    flist = glob(os.path.join(indir, "*.txt"))
    raw = {}
    snames = []
    for fn in flist:
        sn = re.sub(r"\.mtarray.*", "", os.path.split(fn)[1])
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

    # null values seem to be indicated by a single space (great...)
    res = pd.DataFrame(raw).replace(' ', np.nan).astype(float)
    outdir = output.unique_output_dir("hipsci_methylation")
    res.to_csv(os.path.join(outdir, 'beta.csv.gz'), compression='gzip')

import pandas as pd
import os
import re
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR, RNASEQ_DIR
from matplotlib import pyplot as plt
import seaborn as sns
from utils import log, output
logger = log.get_console_logger()


if __name__ == '__main__':
    xcell_tcga_fn = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'xCell_TCGA_RSEM.txt')
    xcell_tcga = pd.read_csv(xcell_tcga_fn, sep='\t', header=0, index_col=0)

    # convert sample name format
    xcell_tcga.columns = xcell_tcga.columns.str.replace(r'\.[0-9]{2}$', '')

    # load IPA signatures
    ## TODO: extract these!

    # rnaseq GBM data
    tcga_gbm_meta_fn = os.path.join(RNASEQ_DIR, 'tcga_gbm', 'gliovis_tcga_gbmlgg_meta.csv')
    tcga_gbm_meta = pd.read_csv(tcga_gbm_meta_fn, header=0, index_col=0)
    # filter to include only GBM IDH1wt
    tcga_gbm_meta = tcga_gbm_meta.loc[
        (tcga_gbm_meta.Histology == 'GBM') & (tcga_gbm_meta['IDH.status'] == 'WT')
    ]

    # reduce to intersection, reporting drop outs
    ix = tcga_gbm_meta.index.intersection(xcell_tcga.columns)
    dropped = tcga_gbm_meta.index.difference(xcell_tcga.columns)
    if len(dropped):
        logger.warning("Dropping %d RNA-Seq samples that are not in the xCell data. %d remain.", len(dropped), len(ix))
    xcell_tcga_rnaseq = xcell_tcga.loc[:, ix]


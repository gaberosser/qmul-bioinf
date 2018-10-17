import pandas as pd
import os
import re
import csv
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR, RNASEQ_DIR
from matplotlib import pyplot as plt
import seaborn as sns
from utils import log, output
logger = log.get_console_logger()


def load_ipa_signatures(fn):
    res = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f)
        for row in c:
            res[row[0]] = row[1:]
    return res


if __name__ == '__main__':
    xcell_tcga_fn = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'xCell_TCGA_RSEM.txt')
    xcell_tcga = pd.read_csv(xcell_tcga_fn, sep='\t', header=0, index_col=0)

    # convert sample name format
    xcell_tcga.columns = xcell_tcga.columns.str.replace(r'\.[0-9]{2}$', '')

    # load IPA signatures
    ipa_sign_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'ipa_pathways', 'ipa_exported_pathways.csv')
    ipa_signatures = load_ipa_signatures(ipa_sign_fn)

    # rnaseq GBM data and metadata
    tcga_gbm_meta_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'tcga', 'GlioVis_TCGA_GBMLGG.meta.xlsx')
    tcga_gbm_meta = pd.read_excel(tcga_gbm_meta_fn, header=0, index_col=0)
    # filter to include only GBM IDH1wt
    tcga_gbm_meta = tcga_gbm_meta.loc[
        (tcga_gbm_meta.Histology == 'GBM') & (tcga_gbm_meta['IDH.status'] == 'WT')
    ]

    # reduce to intersection, reporting drop outs
    ix = tcga_gbm_meta.index.intersection(xcell_tcga.columns)
    dropped = tcga_gbm_meta.index.difference(xcell_tcga.columns)
    if len(dropped):
        logger.warning("Dropping %d GlioVis samples that are not in the xCell data. %d remain.", len(dropped), len(ix))
    xcell_tcga_rnaseq = xcell_tcga.loc[:, ix]
    tcga_gbm_rnaseq_meta = tcga_gbm_meta.loc[ix]

    tcga_rnaseq_data_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'tcga', 'rnaseq.xlsx')
    tcga_rnaseq_data = pd.read_excel(tcga_rnaseq_data_fn, header=0, index_col=0, sheet_name='fpkm')


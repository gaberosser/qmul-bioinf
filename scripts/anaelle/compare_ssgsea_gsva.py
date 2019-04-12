from rnaseq.gsea import ssgsea, run_one_ssgsea
import pandas as pd
from settings import DATA_DIR_NON_GIT, GIT_LFS_DATA_DIR, HGIC_LOCAL_DIR
import os
import csv
import references
import datetime
from matplotlib import pyplot as plt
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import seaborn as sns
import numpy as np
import collections
from scipy import stats

from utils.output import unique_output_dir


if __name__ == '__main__':
    # Step 1: compare RNA-Seq count data and some known gene signatures
    # We just want to demonstrate that ssGSEA and GSVA are similar
    ipa_pathway_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/ipa_pathways',
        'ipa_exported_pathways_ensembl_ids.csv'
    )
    ipa_pathways = {}
    with open(ipa_pathway_fn, 'rb') as f:
        c = csv.reader(f)
        for row in c:
            ipa_pathways[row[0]] = row[2:]
    pass
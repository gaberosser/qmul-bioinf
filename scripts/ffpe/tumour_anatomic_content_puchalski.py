from settings import GIT_LFS_DATA_DIR, HGIC_LOCAL_DIR
from stats import nht
from utils import output
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP

import os
from scipy import stats
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def nh_id_to_patient_id(arr):
    nh_id = pd.Index(arr).str.replace(r'(_?)(DEF|SP).*', '')
    return [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    rna_ff_samples = [
        'NH15_1661DEF2C',
        'NH15_1877_SP1C',
        'NH15_2101_DEF1A',
        'NH16_270_DEF1Ereplacement',
        'NH16_616DEF1B',
        'NH16_677_SP1A',
        'NH16_2063_DEF1Areplacement',
        'NH16_2214DEF1A',
        'NH16_2255DEF1B2',
        'NH16_2806DEF3A1'
    ]

    xcell_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/characterisation/ffpe_cell_type_determination/xcell/xcell_results_salmon_tpm.xlsx'
    )
    xcell_res = pd.read_excel(xcell_fn, sheet_name=None)
    xcell_prop = xcell_res['Proportions'][rna_ff_samples]
    xcell_prop.columns = nh_id_to_patient_id(xcell_prop.columns)

    xcell_pval = xcell_res['P values'][rna_ff_samples]
    xcell_pval.columns = nh_id_to_patient_id(xcell_pval.columns)


    puch_fn = os.path.join(
        GIT_LFS_DATA_DIR,
        'puchalski_2018_tumour_anatomy',
        'ffpe_salmon_tpm_cibersort.csv'
    )
    puch_res = pd.read_csv(puch_fn, header=0, index_col=0).loc[rna_ff_samples]
    puch_res.index = nh_id_to_patient_id(puch_res.index)

    # correlation: Tregs : tumour core proportion
    a = xcell_prop.loc['Tregs']
    b = puch_res.loc[a.index, 'CT']

    lr = stats.linregress(a.values.astype(float), b.values.astype(float))

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.scatter(a.values, b.values, label=None)
    x = np.array([a.min(), a.max()])
    ax.plot(x, lr.intercept + lr.slope * x, 'k--', label='p = %.2f' % lr.pvalue)
    ax.set_xlabel('Proportion Tregs')
    ax.set_ylabel('Proportion core tumour')
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "all_samples_ct_vs_tregs.png"), dpi=200)

    print "Full dataset"
    print "Pearson r: %.3f (p=%.3f)" % stats.pearsonr(a, b)
    print "Spearman r: %.3f (p=%.3f)" % stats.spearmanr(a, b)

    # remove low quality FFPE samples
    a = a.drop(['026', '052'])
    b = b.loc[a.index]

    lr = stats.linregress(a.values.astype(float), b.values.astype(float))

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.scatter(a.values, b.values, label=None)
    x = np.array([a.min(), a.max()])
    ax.plot(x, lr.intercept + lr.slope * x, 'k--', label='p = %.2f' % lr.pvalue)
    ax.set_xlabel('Proportion Tregs')
    ax.set_ylabel('Proportion core tumour')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "good_samples_ct_vs_tregs.png"), dpi=200)

    print "Good quality sequencing only"
    print "Pearson r: %.3f (p=%.3f)" % stats.pearsonr(a, b)
    print "Spearman r: %.3f (p=%.3f)" % stats.spearmanr(a, b)
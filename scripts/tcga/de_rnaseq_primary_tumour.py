from rnaseq import loader, differential_expression, general, filter
import os
import pandas as pd
from utils import output
from settings import RNASEQ_DIR
import numpy as np


"""
Aim:

Carry out DE analysis on the TCGA RNA-Seq data, all primary tumour samples vs all solid healthy tissue
"""

if __name__ == "__main__":
    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'QLGLM'
    }

    outdir = output.unique_output_dir("tcga_de")
    indir_pt = os.path.join(RNASEQ_DIR, 'tcga_gbm', 'primary_tumour')
    indir_hn = os.path.join(RNASEQ_DIR, 'tcga_gbm', 'solid_tissue_normal')


    dat_fn = os.path.join(indir_pt, 'rnaseq.htseq.csv.gz')
    dat_hn_fn = os.path.join(indir_hn, 'rnaseq_normal.htseq.csv.gz')
    meta_fn = os.path.join(indir_pt, 'brennan_s7.csv')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    dat_pt = pd.read_csv(dat_fn, header=0, index_col=0)
    dat_hn = pd.read_csv(dat_hn_fn, header=0, index_col=0)

    dat_pt.columns = [t[:12] for t in dat_pt.columns]
    dat_pt = dat_pt.loc[:, ~dat_pt.columns.duplicated()]
    meta = meta.loc[dat_pt.columns]
    meta = meta.loc[~meta.index.duplicated()]
    dat_hn.columns = ["%s_normal" % t[:12] for t in dat_hn.columns]

    # filter samples
    keep = (meta['gcimp_methylation'] == 'non-G-CIMP') & (meta['idh1_status'] == 'WT')
    dat_pt = dat_pt.loc[:, keep]

    all_dat = pd.concat((dat_pt, dat_hn), axis=1)
    all_dat = filter.filter_by_cpm(all_dat, min_n_samples=3, min_cpm=1)

    groups = ['GBM'] * dat_pt.shape[1] + ['control'] * dat_hn.shape[1]

    de_res = differential_expression.run_one_de(
        all_dat,
        groups,
        ('GBM', 'control'),
        **de_params
    )

    # de_res.to_excel(os.path.join(outdir, 'tcga_primary_tumour_vs_normal_solid.xlsx'))

    de_res_full = differential_expression.run_one_de(
        all_dat,
        groups,
        ('GBM', 'control'),
        return_full=True,
        **de_params
    )
    de_res_full.to_excel(os.path.join(outdir, 'tcga_primary_tumour_vs_normal_solid.xlsx'))

    from matplotlib import pyplot as plt
    import seaborn as sns

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axvline(0., c='k', ls='--')
    # ax.axhline(0., c='k', ls='-')
    ax.scatter(de_res_full.logFC, -np.log10(de_res_full.FDR), c='gray', alpha=0.2, label=None)
    ix_up = de_res.logFC > 0
    ix_down = de_res.logFC < 0
    lup = ax.scatter(
        de_res.logFC.loc[ix_up],
        -np.log10(de_res.FDR).loc[ix_up],
        c='r',
        alpha=0.8,
        label="%d genes upregulated" % ix_up.sum()
    )
    ldown = ax.scatter(
        de_res.logFC.loc[ix_down],
        -np.log10(de_res.FDR).loc[ix_down],
        c='g',
        alpha=0.8,
        label="%d genes downregulated" % ix_down.sum()
    )
    ax.set_xlim([-12.5, 12.5])
    ax.set_xlabel('logFC')
    ax.set_ylabel('-log10(FDR)')
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "volcano_plot_de.png"), dpi=200)





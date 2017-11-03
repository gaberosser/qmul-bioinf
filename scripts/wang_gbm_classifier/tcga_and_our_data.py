import pandas as pd
import numpy as np
from rnaseq import gsea
from load_data import rnaseq_data
from utils.output import unique_output_dir
from settings import OUTPUT_DIR
import os
import references
from matplotlib import pyplot as plt
import seaborn as sns


def prepare_gct_files(outdir=None):
    """
    Prepare the GCT files required to perform classification:
    - Our GBM FFPE and cell culture samples
    - TCGA RNA-Seq cohort
    - Both combined
    In all cases, use FPKM units and gene symbols, as these are used by Wang
    """
    if outdir is None:
        outdir = unique_output_dir("gct_files_for_wang")

    infiles = []

    # 1) Our data
    obj_ffpe = rnaseq_data.load_by_patient('all', type='ffpe')
    dat_ffpe = obj_ffpe.get_fpkm()
    dat_ffpe.columns = ['%s_FFPE' % t for t in obj_ffpe.meta.reference_id]
    obj_cc = rnaseq_data.load_by_patient(patient_ids='all')
    dat_cc = obj_cc.get_fpkm()
    dat_cc = dat_cc.loc[:, obj_cc.meta.type == 'GBM']
    dat_all = pd.concat((dat_cc, dat_ffpe), axis=1)
    idx = references.ensembl_to_gene_symbol(dat_all.index).dropna()
    dat_all = dat_all.loc[idx.index]
    dat_all.index = idx
    fn = os.path.join(outdir, "gbm_ffpe_cc_fpkm.gct")
    gsea.data_to_gct(dat_all, fn)
    infiles.append(fn)

    # 2) TCGA (IDH1 WT only)
    tcga_dat, tcga_meta = rnaseq_data.tcga_primary_gbm(units='fpkm')
    tcga_dat = tcga_dat.loc[:, tcga_meta.idh1_status == 'WT']
    idx = references.ensembl_to_gene_symbol(tcga_dat.index).dropna()
    idx = idx.loc[~idx.index.duplicated()]
    tcga_dat = tcga_dat.loc[idx.index]
    tcga_dat.index = idx
    fn = os.path.join(outdir, "tcga_idh1_wt_fpkm.gct")
    gsea.data_to_gct(tcga_dat, fn)
    infiles.append(fn)

    # 3) Combined
    dat = gsea.combine_gct_files(*infiles)
    fn = os.path.join(outdir, "tcga_idh1_wt_and_gbm_ffpe_cc_fpkm.gct")
    gsea.data_to_gct(dat, fn)


def load_pvalue_results(fn):
    dat = pd.read_csv(fn, header=0, index_col=0, delimiter='\t')
    # only keep the p values
    ncol = dat.columns.size
    dat = dat.iloc[:, (ncol / 2):]
    dat.columns = dat.columns.str.replace('_pval', '')
    return dat


def simplicity_score(pvals):
    """
    For each sample (row), compute the simplicity score defined in Wang et al.
    :param pvals:
    :return:
    """
    # Rank the pvalues. This method chooses the first column it encounters in the event of a tie. This is fine as it
    # doesn't affect the outcome.
    n_cls = pvals.columns.size
    if n_cls < 2:
        raise AttributeError("Cannot compute a simplicity score with fewer than 2 classes")
    rnk = pvals.rank(axis=1, method='first')
    adds = pd.Series(index=pvals.index)
    adns = pd.Series(index=pvals.index)
    rng = pd.Series(index=pvals.index)
    for ix in pvals.index:
        p = pvals.loc[ix].values
        r = rnk.loc[ix].values
        p0 = p[r == 1]
        adds.loc[ix] = (p[r > 1] - p0).sum()
        this_adns = 0.
        for i in range(2, n_cls + 1):
            for j in range(i, n_cls + 1):
                this_adns += (p[r == j] - p[r == i])
        adns.loc[ix] = this_adns
        rng.loc[ix] = p[r == n_cls] - p0
    return (adds - adns) * rng / float(n_cls - 1)


if __name__ == '__main__':
    alpha = 0.04
    outdir = unique_output_dir("wang_classification")

    # our data

    indir = os.path.join(OUTPUT_DIR, 'wang_classification')
    fn = os.path.join(indir, 'p_result_gbm_cc_and_ffpe_fpkm.gct.txt')
    pvals_ours = load_pvalue_results(fn)
    ss_ours = simplicity_score(pvals_ours)
    nm_ours = (pvals_ours < alpha).sum(axis=1)

    # compare with methylation subtype
    heidelberg_subtype_ffpe = pd.Series({
        '017': 'RTK II',
        '018': 'RTK I',
        '019': 'RTK I',
        '026': 'MES',
        '030': 'RTK I',
        '031': 'RTK I',
        '044': 'MES',
        '049': 'None',
        '050': 'RTK II',
        '052': 'MES',
        '054': 'RTK II',
        '061': 'RTK II',
    })

    # obj_ours_cc = rnaseq_data.load_by_patient('all', include_control=False)
    # obj_ours_ffpe = rnaseq_data.load_by_patient('all', type='ffpe', include_control=False)


    # TCGA data (primary tumours)

    fn = os.path.join(indir, 'p_result_tcga_idh1_wt_fpkm.gct.txt')
    pvals_tcga = load_pvalue_results(fn)
    ss_tcga = simplicity_score(pvals_tcga)
    nm_tcga = (pvals_tcga < alpha).sum(axis=1)
    cls_tcga = pd.Series(index=pvals_tcga.index)
    min_idx = np.argmin(pvals_tcga.values, axis=1)
    cls_tcga.loc[nm_tcga == 1] = pvals_tcga.columns[min_idx[nm_tcga == 1]]
    cls_tcga.loc[nm_tcga > 1] = 'Multi'

    # load TCGA meta so we have the previous classification
    _, tcga_meta = rnaseq_data.tcga_primary_gbm(units='fpkm')
    # change index to match
    tcga_meta.index = tcga_meta.index.str.replace('-', '.')
    # extract only the same samples we have here
    tcga_meta = tcga_meta.loc[pvals_tcga.index].dropna(how='all')

    # contingency table and simplicity of each element
    previous = tcga_meta.loc[:, 'expression_subclass']
    new = cls_tcga
    ctg_tcga = pd.DataFrame(
        index=['Proneural', 'Mesenchymal', 'Classical', 'Neural', 'G-CIMP', 'None'],
        columns=['Proneural', 'Mesenchymal', 'Classical', 'Multi', 'None']
    )
    ctg_ss_tcga = pd.DataFrame(ctg_tcga, copy=True)

    for ix in ctg_tcga.index:
        if ix == "None":
            the_ids = tcga_meta.loc[previous.isnull()].index
        else:
            the_ids = tcga_meta.loc[previous == ix].index
        for col in ctg_tcga.columns:
            the_match = new.loc[the_ids]
            if col == "None":
                this_ix = the_match.isnull()
            else:
                this_ix = (the_match == col)
            ctg_tcga.loc[ix, col] = this_ix.sum()
            ctg_ss_tcga.loc[ix, col] = ss_tcga.loc[this_ix.index[this_ix]].mean()

    ax = sns.heatmap(ctg_ss_tcga.fillna(0), vmin=0, vmax=1, cmap='RdBu_r', annot=True, cbar=False)
    # turn the axis labels
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "simplicity_scores_tcga.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "simplicity_scores_tcga.pdf"))

    ax = sns.heatmap(ctg_tcga.astype(float), cmap='RdBu_r', annot=True, cbar=False)
    # turn the axis labels
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "contingency_tcga.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "contingency_tcga.pdf"))

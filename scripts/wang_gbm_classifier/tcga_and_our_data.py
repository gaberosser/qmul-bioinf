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


def contingency_table(new, previous, vals=None, val_func=np.mean):
    """
    Previous values go on the INDEX, new values on the COLUMNS
    :param new:
    :param previous:
    :param vals:
    :return:
    """
    _, new_cls = new.factorize()
    new_cls = set(new_cls).difference({'None', 'Multi'})

    _, prev_cls = previous.factorize()
    prev_cls = set(prev_cls).difference({'None', 'Multi'})

    new_idx = list(new_cls) + ['Multi', 'None']
    prev_idx = list(prev_cls) + ['Multi', 'None']

    ctg = pd.DataFrame(
        index=prev_idx,
        columns=new_idx
    )

    for ix in ctg.index:

        if ix == "None":
            the_ids = previous.loc[previous.isnull()].index
        else:
            the_ids = previous.loc[previous == ix].index
        if len(the_ids) == 0:
            continue

        for col in ctg.columns:
            the_match = new.loc[the_ids]
            if col == "None":
                this_ix = the_match.isnull()
            else:
                this_ix = (the_match == col)
            if vals is None:
                # just count
                ctg.loc[ix, col] = this_ix.sum()
            else:
                # store values
                if val_func is None:
                    ctg.loc[ix, col] = vals.loc[this_ix.index[this_ix]].tolist()
                else:
                    ctg.loc[ix, col] = val_func(vals.loc[this_ix.index[this_ix]])
                
    return ctg


if __name__ == '__main__':
    alpha = 0.05
    outdir = unique_output_dir()

    # our data

    indir = os.path.join(OUTPUT_DIR, 'wang_classification')
    fn = os.path.join(indir, 'p_result_gbm_cc_and_ffpe_fpkm.gct.txt')
    pvals_ours = load_pvalue_results(fn)
    ss_ours = simplicity_score(pvals_ours)
    nm_ours = (pvals_ours < alpha).sum(axis=1)
    cls_ours = pd.Series(index=pvals_ours.index)
    min_idx = np.argmin(pvals_ours.values, axis=1)
    cls_ours.loc[nm_ours == 1] = pvals_ours.columns[min_idx[nm_ours == 1]]
    cls_ours.loc[nm_ours > 1] = 'Multi'

    # easy to read table
    our_summary = pd.concat((cls_ours, ss_ours, pvals_ours), axis=1).sort_index()

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
    cls_tcga_export = pd.DataFrame(index=pvals_tcga.index, columns=['Class'])

    min_idx = np.argmin(pvals_tcga.values, axis=1)
    cls_tcga.loc[nm_tcga == 1] = pvals_tcga.columns[min_idx[nm_tcga == 1]]
    cls_tcga_export.loc[nm_tcga == 1, 'Class'] = pvals_tcga.columns[min_idx[nm_tcga == 1]]

    cls_tcga.loc[nm_tcga > 1] = 'Multi'
    for row in nm_tcga.index[nm_tcga > 1]:
        p = pvals_tcga.loc[row]
        cls_tcga_export.loc[row, 'Class'] = ','.join(p.index[p < alpha])
    cls_tcga_export.insert(1, 'Simplicity score', ss_tcga)

    cls_tcga_export.to_csv(os.path.join(outdir, "tcga_wang_classification.csv"))

    # load TCGA meta so we have the previous classification
    _, tcga_meta = rnaseq_data.tcga_primary_gbm(units='fpkm')
    # change index to match
    tcga_meta.index = tcga_meta.index.str.replace('-', '.')
    # extract only the same samples we have here
    tcga_meta = tcga_meta.loc[pvals_tcga.index].dropna(how='all')

    # contingency table and simplicity of each element
    idx = ['Proneural', 'Mesenchymal', 'Classical', 'Neural', 'G-CIMP', 'None']
    cols = ['Proneural', 'Mesenchymal', 'Classical', 'Multi', 'None']

    ctg_tcga = contingency_table(cls_tcga, tcga_meta.loc[:, 'expression_subclass'])
    ctg_ss_tcga = contingency_table(cls_tcga, tcga_meta.loc[:, 'expression_subclass'], ss_tcga)

    ctg_tcga = ctg_tcga.loc[idx, cols]
    ctg_ss_tcga = ctg_ss_tcga.loc[idx, cols]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(ctg_ss_tcga.fillna(0), vmin=0, vmax=1, cmap='RdBu_r', annot=True, cbar=False, ax=ax)
    # turn the axis labels
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "simplicity_scores_tcga.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "simplicity_scores_tcga.pdf"))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(ctg_tcga.astype(float), cmap='RdBu_r', annot=True, cbar=False, ax=ax)
    # turn the axis labels
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "contingency_tcga.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "contingency_tcga.pdf"))

    # how is the classification affected if we combine our data and TCGA data?
    fn = os.path.join(indir, 'p_result_tcga_idh1_wt_and_gbm_cc_ffpe_fpkm.gct.txt')
    pvals_both = load_pvalue_results(fn)
    ss_both = simplicity_score(pvals_both)
    nm_both = (pvals_both < alpha).sum(axis=1)
    cls_both = pd.Series(index=pvals_both.index)
    min_idx = np.argmin(pvals_both.values, axis=1)
    cls_both.loc[nm_both== 1] = pvals_both.columns[min_idx[nm_both == 1]]
    cls_both.loc[nm_both > 1] = 'Multi'

    # our data contingency
    previous = cls_ours
    new = cls_both.loc[cls_ours.index]
    cols = ['Proneural', 'Mesenchymal', 'Classical', 'Multi', 'None']
    ctg_ours_both = contingency_table(new, previous).loc[cols, cols]

    previous = cls_tcga
    new = cls_both.loc[cls_tcga.index]
    ctg_tcga_both = contingency_table(new, previous).loc[cols, cols]

    # display as a pie/bar-chart
    pids = heidelberg_subtype_ffpe.sort_index().index
    classes =  ['Proneural', 'Mesenchymal', 'Classical']
    ind = np.arange(len(pids))

    p_inv = 1. - pvals_ours

    ffpe_vals = []
    culture_p1_vals = []
    culture_p2_vals = []

    for p in pids:
        this_ffpe = p_inv.loc["GBM%s_FFPE" % p, classes]
        ffpe_vals.append(this_ffpe.values)
        this_cc = np.where(p_inv.index.str.contains("GBM%s_P" % p))[0]
        culture_p1_vals.append(p_inv.iloc[this_cc[0]].loc[classes].values)
        if len(this_cc) > 1:
            culture_p2_vals.append(p_inv.iloc[this_cc[1]].loc[classes].values)
        else:
            culture_p2_vals.append(np.zeros(len(classes)))

    ffpe_vals = np.concatenate((np.zeros((len(pids), 1)), ffpe_vals), axis=1)
    ccp1_vals = np.concatenate((np.zeros((len(pids), 1)), culture_p1_vals), axis=1)
    ccp2_vals = np.concatenate((np.zeros((len(pids), 1)), culture_p2_vals), axis=1)

    ffpe_cs = ffpe_vals.cumsum(axis=1)
    ccp1_cs = ccp1_vals.cumsum(axis=1)
    ccp2_cs = ccp2_vals.cumsum(axis=1)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = .2
    offset = 0.05
    colours = {
        'Proneural': '#7030A0',
        'Mesenchymal': '#22D05E',
        'Classical': '#F79646',
    }

    for i in range(len(classes)):
        c = colours[classes[i]]
        x_ffpe = ind
        x_cc1 = ind + width + offset
        x_cc2 = ind + 2 * (width + offset)

        the_ffpe = ffpe_vals[:, i + 1]
        the_ffpe_bottom = ffpe_cs[:, i]
        ax.bar(x_ffpe, the_ffpe, width, bottom=the_ffpe_bottom, color=c, edgecolor=c, lw=1.5)

        the_cc1 = ccp1_vals[:, i + 1]
        the_cc1_bottom = ccp1_cs[:, i]
        ax.bar(x_cc1, the_cc1, width, bottom=the_cc1_bottom, color=c, edgecolor=c, lw=1.5)

        the_cc2 = ccp2_vals[:, i + 1]
        the_cc2_bottom = ccp2_cs[:, i]
        ax.bar(x_cc2, the_cc2, width, bottom=the_cc2_bottom, color=c, edgecolor=c, lw=1.5)

        idx = the_ffpe > (1 - alpha)
        ax.bar(x_ffpe[idx], the_ffpe[idx], 0.25, bottom=the_ffpe_bottom[idx],
               color='none', edgecolor='k', lw=1.5, zorder=999)

        idx = the_cc1 > (1 - alpha)
        ax.bar(x_cc1[idx], the_cc1[idx], 0.25, bottom=the_cc1_bottom[idx],
               color='none', edgecolor='k', lw=1.5, zorder=999)

        idx = the_cc2 > (1 - alpha)
        ax.bar(x_cc2[idx], the_cc2[idx], 0.25, bottom=the_cc2_bottom[idx],
               color='none', edgecolor='k', lw=1.5, zorder=999)

    ax.set_xticks(ind + 1.5 * (width + offset))
    ax.set_xticklabels(pids, horizontalalignment='center')

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'wang_pvalues_chart.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'wang_pvalues_chart.pdf'))

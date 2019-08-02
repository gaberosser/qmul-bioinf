from rnaseq import gsva, loader
import pandas as pd
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR, DATA_DIR
from plotting import venn, common

import os
import references
import datetime
from matplotlib import pyplot as plt, colors
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import seaborn as sns
import numpy as np
import collections
from scipy import stats
from utils import setops, output, log

logger = log.get_console_logger()

def line_plot_pvalues_slope(
        pvals,
        slopes,
        cmap=plt.get_cmap('Reds'),
        alpha=None,
        log_scale=True,
        vmin=None,
        vmax=None,
):
    """
    Plot summarising pvalue and correlation slope simultaneously using position (slope) and colour (pval).
    Hard coded into vertical orientation (TODO: make this a parameter??)
    :param pvals: DataFrame. Index and columns will be used for ticklabels. Suggest using -log10
    :param slopes: DataFrame, must match pvals
    :param pct_to_size_func:
    :param cmap: cmap used to represent
    :param alpha: If supplied, highlight all results with p < alpha.
    :param log_scale: If True (default), convert p values to -log10 scale.
    :param vmin: For pvalue shading. If not supplied, min of data will be used
    :param vmax: For pvalue shading. If not supplied, max of data will be used
    :return:
    """
    if sorted(pvals.index) != sorted(slopes.index):
        raise AttributeError("Index of pvals and concords must match")
    if sorted(pvals.columns) != sorted(slopes.columns):
        raise AttributeError("Columns of pvals and concords must match")

    if alpha is None:
        alpha = -1.

    signif = pvals < alpha

    if log_scale:
        pvals = -np.log10(pvals.astype(float))

    slopes = slopes.loc[pvals.index, pvals.columns]

    if vmin is None:
        vmin = pvals.values.min()
    if vmax is None:
        vmax = pvals.values.max()

    ny, nx = pvals.shape
    markers = common.get_best_marker_map(nx)

    gs = plt.GridSpec(nrows=2, ncols=1, height_ratios=[1, 19])
    fig = plt.figure(figsize=(1.5 * nx, .5 * ny))
    ax = fig.add_subplot(gs[1])
    ax.invert_yaxis()

    cax = fig.add_subplot(gs[0])

    plogp_norm = colors.Normalize(vmin=vmin, vmax=vmax)
    plogp_sm = plt.cm.ScalarMappable(cmap=cmap, norm=plogp_norm)

    for i, col in enumerate(slopes.columns):
        ew = [1.5 if t else 0.7 for t in signif[col].values]
        ax.scatter(
            slopes[col],
            range(ny),
            c=[plogp_sm.to_rgba(t) for t in pvals[col].values],
            s=60,
            edgecolor='k',
            linewidths=ew,
            marker=markers[i]
        )

    ax.set_xlabel('Slope', fontsize=12)
    ax.set_yticks(range(ny))
    ax.set_yticklabels(pvals.index, fontsize=12)
    # ax.set_xlim([50, 100])
    plogp_sm.set_array(pvals)
    fig.colorbar(plogp_sm, cax=cax, orientation='horizontal')
    cax.xaxis.set_label_position('top')
    cax.set_xlabel(r'$-\log_{10}(p)$')

    type_attrs = {
        'class': 'line',
        'linestyle': 'none',
        'markeredgecolor': 'k',
        'markeredgewidth': 1.,
        'markerfacecolor': 'none',
        'markersize': 8
    }

    leg_dict = {}
    for i, col in enumerate(pvals.columns):
        leg_dict[col] = dict(type_attrs)
        leg_dict[col].update({'marker': markers[i]})

    common.add_custom_legend(ax, leg_dict, loc_outside=True, loc_outside_horiz='right')
    gs.update(left=0.2, bottom=0.1, right=0.72, top=0.95, hspace=0.12)

    return {
        'fig': fig,
        'ax': ax,
        'gs': gs,
        'cax': cax
    }


kegg_mtor_from_msigdb = [
    "AKT3", "EIF4B", "EIF4E", "EIF4EBP1", "AKT1", "AKT2", "FIGF", "PIK3R5", "MTOR", "RICTOR", "EIF4E1B", "ULK3",
    "RPS6KA6", "HIF1A", "IGF1", "INS", "PDPK1", "CAB39", "PGF", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1",
    "PIK3R2", "DDIT4", "PRKAA1", "PRKAA2", "MAPK1", "MAPK3", "RPTOR", "RHEB", "RPS6", "RPS6KA1", "RPS6KA2",
    "RPS6KA3", "RPS6KB1", "RPS6KB2", "MLST8", "BRAF", "STK11", "TSC1", "TSC2", "VEGFA", "VEGFB", "VEGFC", "CAB39L",
    "ULK1", "PIK3R3", "STRADA", "EIF4E2", "ULK2"
]

# downloaded directly from KEGG website (hsa04150)
kegg_mtor_from_kegg = [
    "SLC7A5", "SLC3A2", "SLC38A9", "ATP6V1A", "ATP6V1B1", "ATP6V1B2", "ATP6V1C2", "ATP6V1C1", "ATP6V1D", "ATP6V1E2",
    "ATP6V1E1", "ATP6V1F", "ATP6V1G1", "ATP6V1G3", "ATP6V1G2", "ATP6V1H", "LAMTOR1", "LAMTOR2", "LAMTOR3",
    "LAMTOR4", "LAMTOR5", "FLCN", "FNIP1", "FNIP2", "RRAGA", "RRAGB", "RRAGC", "RRAGD", "SESN2", "CASTOR1",
    "CASTOR2", "MIOS", "SEH1L", "WDR24", "WDR59", "SEC13", "DEPDC5", "NPRL2", "NPRL3", "SKP2", "RNF152", "RPTOR",
    "AKT1S1", "MTOR", "DEPTOR", "MLST8", "TELO2", "TTI1", "CLIP1", "GRB10", "ULK1", "ULK2", "EIF4EBP1", "EIF4E",
    "EIF4E2", "EIF4E1B", "RPS6KB1", "RPS6KB2", "EIF4B", "RPS6", "STRADA", "STRADB", "STK11", "CAB39", "CAB39L",
    "PRKAA1", "PRKAA2", "TSC1", "TSC2", "TBC1D7", "TBC1D7", "RHEB", "DDIT4", "WNT1", "WNT2", "WNT2B", "WNT3",
    "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10B",
    "WNT10A", "WNT11", "WNT16", "FZD1", "FZD7", "FZD2", "FZD3", "FZD4", "FZD5", "FZD8", "FZD6", "FZD10", "FZD9",
    "LRP5", "LRP6", "DVL3", "DVL2", "DVL1", "GSK3B", "TNF", "TNFRSF1A", "IKBKB", "INS", "IGF1", "INSR", "IGF1R",
    "GRB2", "SOS1", "SOS2", "HRAS", "KRAS", "NRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA3",
    "RPS6KA1", "RPS6KA2", "RPS6KA6", "IRS1", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3CA", "PIK3CD", "PIK3CB", "PTEN",
    "PDPK1", "AKT1", "AKT2", "AKT3", "CHUK", "MAPKAP1", "RICTOR", "PRR5", "RHOA", "PRKCA", "PRKCB", "PRKCG", "SGK1",
    "LPIN1"
]

pid_mtor = [
    "SSPO", "SGK1", "EEF2K", "IKBKB", "PLD2", "PDPK1", "ATG13", "ULK1", "NRAS", "HRAS", "KRAS", "RAF1", "EIF4E",
    "EEF2", "BRAF", "PRKCA", "RPS6KB1", "EIF4B", "CCNE1", "CDK2", "YY1", "YWHAQ", "MAPK3", "MAPK1", "PML", "CLIP1",
    "AKT1", "YWHAB", "SFN", "IRS1", "MAP2K2", "SREBF1", "MTOR", "PXN", "TSC2", "EIF4A1", "RHOA", "YWHAG", "YWHAE",
    "RAC1", "YWHAZ", "PRR5", "CYCS", "MAP2K1", "YWHAH", "BNIP3", "PLD1", "EIF4EBP1", "RHEB", "RPS6KA1", "PDCD4",
    "RRAGB", "RICTOR", "RRAGA", "ULK2", "RPTOR", "DEPTOR", "RB1CC1", "TSC1", "AKT1S1", "MAPKAP1", "MLST8",
    "POLDIP3", "RRAGC", "RRAGD", "DDIT4", "RRN3", "PPARGC1A", "FBXW11"
]

biocarta_mtor = [
    "EIF4A1", "EIF4A2", "EIF4B", "EIF4E", "EIF4EBP1", "EIF4G1", "EIF4G2", "AKT1", "FKBP1A", "MTOR", "PDK2", "PDPK1",
    "PIK3CA", "PIK3R1", "PPP2CA", "PTEN", "RPS6", "RPS6KB1", "TSC1", "TSC2", "MKNK1", "EIF3A", "EIF4G3"
]

manual_gene_name_correction = {
    'ATG13': 'KIAA0652',
    'CASTOR1': 'GATSL3',
    'CASTOR2': 'GATSL2',  # might also be GATSL1?
    'DEPTOR': 'DEPDC6',
    'LAMTOR1': 'C11orf59',
    'LAMTOR2': 'ROBLD3',
    'LAMTOR3': 'MAPKSP1',
    'LAMTOR4': 'C7orf59',
    'LAMTOR5': 'HBXIP',
    'TTI1': 'KIAA0406',
    'VEGFD': 'FIGF',
    'HLA-DMB': 'HLA.DMB',
    'HLA-DQA1': 'HLA.DQA1',
    'HLA-DRB5': 'HLA.DRB5'
}


def z_transform(df, axis=None):
    if axis is None:
        return (df - df.values.flatten().mean()) / df.values.flatten().std()
    elif axis == 0:
        i = 0
        j = 1
    elif axis == 1:
        i = 1
        j = 0
    else:
        raise NotImplementedError("axis argument must be None, 0 or 1.")
    return df.subtract(df.mean(axis=i), axis=j).divide(df.std(axis=i), axis=j)


def ols_plot(y, x, add_intercept=True, alpha=0.05, xlim=None, ax=None):
    """
    Generate a scatter plot with OLS prediction plus confidence intervals
    :param y:
    :param x:
    :param add_intercept:
    :param alpha:
    :param ax:
    :return:
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    try:
        x = x.astype(float)
    except Exception:
        pass

    if add_intercept:
        X = sm.add_constant(x)
    else:
        X = x

    model = sm.OLS(y, X)
    res = model.fit()

    # plot data
    ax.scatter(x, y, marker='o')
    if xlim is None:
        xlim = np.array(ax.get_xlim())

    xx = np.linspace(xlim[0], xlim[1], 100)

    # compute prediction and confidence intervals
    if add_intercept:
        b0, b1 = res.params
        sdev, lower, upper = wls_prediction_std(res, sm.add_constant(xx), alpha=alpha)
        # b0_min, b0_max = res.conf_int(alpha=alpha)[0]
        # b1_min, b1_max = res.conf_int(alpha=alpha)[1]

    else:
        b1 = res.params[0]
        b0 = 0.
        sdev, lower, upper = wls_prediction_std(res, xx, alpha=alpha)
        # b0 = b0_min = b0_max = 0.
        # b1_min, b1_max = res.conf_int(alpha=alpha)[0]

    ax.plot(xx, b0 + b1 * xx, 'k-', lw=1.5)
    ax.fill_between(xx, lower, upper, edgecolor='b', facecolor='b', alpha=0.4)

    # lower = b0_min + b1_min * xlim
    # upper = b0_max + b1_max * xlim
    # ax.fill_between(xlim, lower, upper, edgecolor='b', facecolor='b', alpha=0.4)

    ax.set_xlim(xlim)
    return res, ax


def scatter_plot_with_linregress(x, y, group_list=None, groups=None):
    nrow = 2
    reduce_before_return = False
    if groups is None:
        # no subgroups: just run with 'all'
        nrow = 1
        ncol = 1
        group_list = ['foo']
        groups = pd.Series('foo', index=x.index)
        reduce_before_return = True
    else:
        if group_list is None:
            group_list = groups.unique()
        # add 'all' to the group list
        group_list = list(group_list) + [None]
        ncol = int(np.ceil(len(group_list) * 0.5))
    res = pd.DataFrame(index=group_list, columns=['slope', 'intercept', 'rvalue', 'pvalue', 'stderr'])
    sm_res = {}

    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True)
    if reduce_before_return:
        axs = np.array([axs])

    axs_seen = set(axs.flat)

    for i, sg in enumerate(group_list):
        if sg is None:
            sg_idx = pd.Series(True, index=groups.index)
            ttl = 'All'
        else:
            sg_idx = (groups == sg)
            ttl = sg
        this_x = x.loc[sg_idx].values.astype(float)
        this_y = y.loc[sg_idx].values.astype(float)
        lr = stats.linregress(this_x, this_y)
        res.loc[sg] = lr

        ax = axs.flat[i]
        axs_seen.remove(ax)

        sm_res[ttl], _ = ols_plot(
            this_y,
            this_x,
            xlim=[-3.5, 3.5],
            ax=ax
        )

        rsq = lr.rvalue ** 2
        sl = lr.slope
        pval = lr.pvalue

        if np.abs(sl - sm_res[ttl].params[-1]) > 1e-3:
            logger.warn("Subgroup %s. stats.linregress slope doesn't agree with statsmodels OLS.", sg)

        if pval < 0.05:
            lbl = "$R^2 = %.2f$\n$\mathrm{slope}=%.2f$\n$p=\mathbf{%.3e}$" % (rsq, sl, pval)
        else:
            lbl = "$R^2 = %.2f$\n$\mathrm{slope}=%.2f$\n$p=%.3e$" % (rsq, sl, pval)
        ax.text(
            1.,
            0.,
            lbl,
            bbox={'facecolor': 'w', 'alpha': 0.3},
            verticalalignment='bottom',
            horizontalalignment='right',
            transform=ax.transAxes
        )
        ax.set_ylim([-4, 4])
        ax.set_xlabel(x.name)
        ax.set_ylabel(y.name)
        if not reduce_before_return:
            ax.set_title(ttl)

    if reduce_before_return:
        res = res.loc['foo']

    for ax in axs_seen:
        ax.set_visible(False)

    return {
        'fig': fig,
        'axs': axs,
        'linregress': res,
        'statsmodels': sm_res
    }


def get_slope_and_pval(plot_dict, col_order=None):
    """
    Extract the slope and pvalue of the linear regression results
    :param col_order: if supplied, this gives the order of the index in the returned DataFrames.
    :param plot_dict: Dictionary of input data.
    Keys will be used in the results. Values are the dict output from scatter_plot_with_linregress
    :return: Two pd.DataFrame objects: s, p
    """
    if col_order is None:
        col_order = plot_dict.values()[0]['statsmodels'].keys()
    s = pd.DataFrame(index=plot_dict.keys(), columns=col_order)
    p = s.copy()

    for k, v in plot_dict.items():
        s.loc[k] = [v['statsmodels'][t].params[-1] for t in col_order]
        p.loc[k] = [v['statsmodels'][t].f_pvalue for t in col_order]

    return s, p


def plot_signature_vs_gene(dat, es, the_gene, geneset_name, ax=None):
    the_expr = dat.loc[the_gene]
    # Z transform the signature scores for this gene set
    the_signature = es.loc[geneset_name]
    the_signature = (the_signature - the_signature.mean()) / the_signature.std()
    # ensure the ordering is the same
    the_signature = the_signature.loc[the_expr.index]
    lr = stats.linregress(the_signature.astype(float), np.log2(the_expr + 1))
    x_lr = np.array([the_signature.min(), the_signature.max()])
    y_lr = lr.intercept + lr.slope * x_lr

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = None

    ax.scatter(the_signature, np.log2(the_expr + 1))
    ax.plot(x_lr, y_lr, 'k--')
    ax.set_xlabel('Normalised ssGSEA score')
    ax.set_ylabel('log2(%s)' % the_gene)
    if fig is not None:
        fig.tight_layout()
    return ax


def mtor_signature_dict():
    # kegg mtor from msigdb (hsa04150)
    for arr in [kegg_mtor_from_msigdb, kegg_mtor_from_kegg, pid_mtor, biocarta_mtor]:
        for i, t in enumerate(arr):
            if t in manual_gene_name_correction:
                arr[i] = manual_gene_name_correction[t]

    return {
        'kegg_msigdb': kegg_mtor_from_msigdb,
        'kegg': kegg_mtor_from_kegg,
        'biocarta': biocarta_mtor,
        'pid': pid_mtor
    }


def tam_signature_dict():
    # Human-specific by Muller et al.
    muller_tam_signature_fn = os.path.join(GIT_LFS_DATA_DIR, 'muller_2017_tam', '13059_2017_1362_MOESM5_ESM.xlsx')
    muller_tam_signatures = pd.read_excel(muller_tam_signature_fn, header=0, index_col=None)
    muller_tam_signatures = {
        'MG': muller_tam_signatures['MG_Markers'].dropna().values,
        'BMDM': muller_tam_signatures['Mac_Markers'].dropna().values,
    }

    # Mouse version (needs translating) by Bowman et al.
    bowman_tam_signature_fn = os.path.join(DATA_DIR, 'rnaseq', 'GSE86573', 'table_S2.csv')
    bowman_tam_signatures = pd.read_csv(bowman_tam_signature_fn, header=0, index_col=None)

    # generate series of orthologs of the relevant gene signatures
    orth = references.homologs_table(references.mouse_tid, references.human_tid)
    orth = orth.set_index('gene_symbol_10090').squeeze()

    bowman_tam_signatures = {
        'MG': orth.reindex(bowman_tam_signatures['MG'].dropna().values).dropna().values,
        'BMDM': orth.reindex(bowman_tam_signatures['BMDM'].dropna().values).dropna().values,
    }

    return {
        'bowman': bowman_tam_signatures,
        'muller': muller_tam_signatures,
    }


if __name__ == "__main__":
    """
    Use the TCGA cohort to validate a link between the mTOR pathway and the proportion of microglial and macrophage
    immune infiltrate in the bulk samples.
    
    mTOR is assessed using a known set of genes.
    
    Tumour-associated bone marrow-derived macrophages (TAM-BMDM) and microglia (TAM-MG) are distinguished using a 
    human signature from Muller et al. (Genome Biol 2017) or a converted mouse signature from Bowman et al. (???).
    """
    rnaseq_type = 'gliovis'
    remove_idh1 = True
    tam_signature_source = 'bowman'
    # tam_signature_source = 'muller'
    mtor_source = 'kegg_msigdb'  # ('kegg', 'pid', 'biocarta')
    # class_method = 'wang'
    class_method = 'verhaak'
    # toggle allowing more than one class to be used
    allow_multiple_classes = False

    # load mTOR signatures
    mtor_gs_dict = mtor_signature_dict()

    # load MG/BMDM signatures
    tam_gs_dict = tam_signature_dict()

    mtor_geneset = mtor_gs_dict[mtor_source]
    tam_genesets = tam_gs_dict[tam_signature_source]

    genesets = dict(tam_genesets)
    genesets['mTOR'] = mtor_geneset

    outdir = output.unique_output_dir()

    # export all signatures to a file
    from utils import dictionary

    all_gs_dict = dictionary.nested_dict_to_flat(tam_gs_dict)
    all_gs_dict[('mTOR',)] = mtor_geneset

    for_export = pd.DataFrame(index=range(max((len(t) for t in all_gs_dict.values()))), columns=[])
    for k, v in all_gs_dict.items():
        the_key = '_'.join(k)
        for_export.loc[range(len(v)), the_key] = sorted(v)
    for_export.fillna('', inplace=True)
    for_export = for_export.sort_index(axis=1)
    for_export.to_excel(os.path.join(outdir, "gene_sets.xlsx"), index=False)

    # Venn diagram showing various mTOR signature options and overlap between them
    fig, ax = plt.subplots()
    venn.venn_diagram(*mtor_gs_dict.values(), set_labels=mtor_gs_dict.keys(), ax=ax)
    fig.tight_layout()
    ax.set_facecolor('w')
    fig.savefig(os.path.join(outdir, "venn_mtor_genesets.png"), dpi=200)

    basedir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/tcga'
    )

    brennan_s7_fn = os.path.join(basedir, "brennan_s7.csv")
    brennan_s7 = pd.read_csv(brennan_s7_fn, header=0, index_col=0)

    if rnaseq_type == 'counts':
        rnaseq_dat_fn = os.path.join(basedir, 'rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(basedir, 'rnaseq.meta.xlsx')
        sheet_name = 'htseq'
        wang_fn = os.path.join(basedir, 'wang_classification', 'tcga_counts_wang_classification.csv')
    elif rnaseq_type == 'fpkm':
        rnaseq_dat_fn = os.path.join(basedir, 'rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(basedir, 'rnaseq.meta.xlsx')
        sheet_name = 'fpkm'
        wang_fn = os.path.join(basedir, 'wang_classification', 'tcga_fpkm_wang_classification.csv')
    elif rnaseq_type == 'gliovis':
        rnaseq_dat_fn = os.path.join(basedir, 'gliovis', 'gliovis_tcga_gbm_rnaseq.xlsx')
        wang_fn = os.path.join(basedir, 'gliovis', 'wang_classification', 'tcga_gliovis_wang_classification.csv')
        sheet_name = 0
    else:
        raise NotImplementedError("Unrecognised rnaseq data type")


    rnaseq_dat_raw = pd.read_excel(rnaseq_dat_fn, header=0, index_col=0, sheet_name=sheet_name)
    wang_classes = pd.read_csv(wang_fn, header=0, index_col=0)

    if rnaseq_type == 'gliovis':
        rnaseq_meta_fn = os.path.join(basedir, 'gliovis', 'GlioVis_TCGA_GBMLGG.meta.xlsx')
        rnaseq_meta = pd.read_excel(rnaseq_meta_fn, header=0, index_col=0)
        # filter only GBM
        rnaseq_meta = rnaseq_meta.loc[rnaseq_meta.Histology == 'GBM']
        rnaseq_dat_raw = rnaseq_dat_raw.loc[:, rnaseq_meta.index]

        rnaseq_meta.rename(
            columns={'IDH.status': 'idh1_status', 'Subtype.original': 'expression_subclass'},
            inplace=True
        )

    else:
        # simplify sample naming
        new_cols = rnaseq_dat_raw.columns.str.replace(r'(?P<x>TCGA-[0-9]{2}-[0-9]{4})-.*', r'\g<x>')

        # rnaseq_meta = rnaseq_meta.loc[~new_cols.duplicated()]
        rnaseq_dat_raw = rnaseq_dat_raw.loc[:, ~new_cols.duplicated()]
        # rnaseq_meta.index = new_cols[~new_cols.duplicated()]
        rnaseq_dat_raw.columns = new_cols[~new_cols.duplicated()]
        rnaseq_meta = brennan_s7.reindex(rnaseq_dat_raw.columns)


    if remove_idh1:
        # filter IDH1 mutants
        idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')

        rnaseq_meta = rnaseq_meta.loc[idh1_wt]
        rnaseq_dat = rnaseq_dat_raw.loc[:, rnaseq_meta.index]
    else:
        rnaseq_dat = rnaseq_dat_raw.loc[:, rnaseq_dat_raw.columns.str.contains('TCGA')]

    if rnaseq_type != 'gliovis':
        # add gene symbols for gene signature scoring?
        gs = references.ensembl_to_gene_symbol(rnaseq_dat.index).dropna()
        rnaseq_dat = rnaseq_dat.loc[gs.index]
        rnaseq_dat.index = gs.values

    if rnaseq_type == 'counts':
        # convert to CPM
        rnaseq_dat = rnaseq_dat.divide(rnaseq_dat.sum(axis=0), axis=1) * 1e6

    rnaseq_meta.insert(0, 'wang_classification_simplicity', wang_classes.loc[rnaseq_meta.index, 'Simplicity score'])
    rnaseq_meta.insert(0, 'wang_classification_num_matches', wang_classes.loc[rnaseq_meta.index, 'Number of matches'])
    rnaseq_meta.insert(0, 'wang_classification', wang_classes.loc[rnaseq_meta.index, 'Wang subclass'])

    # check that signature genes are all found in the data
    for k, v in genesets.items():
        for i, t in enumerate(v):
            if t in manual_gene_name_correction:
                v[i] = manual_gene_name_correction[t]
        g_in = rnaseq_dat.index.intersection(v)
        if set(g_in) != set(v):
            missing = set(v).difference(rnaseq_dat.index)
            logger.warn(
                "%d genes in the %s signature do not match with the data index and will be dropped: %s.",
                len(missing),
                k,
                ', '.join(missing)
            )
            genesets[k] = g_in

    # check here whether there is any overlap
    vs, vc = setops.venn_from_arrays(*genesets.values())
    n_overlap = sum([vc[t] for t in setops.binary_combinations_sum_gte(len(genesets), 2)])
    if n_overlap > 0:
        logger.warn(
            "The %d gene signatures used here have %d overlapping genes - please check this is OK.",
            len(genesets),
            n_overlap
        )

    # run ssGSEA then Z transform the results
    es = gsva.ssgsea(rnaseq_dat, genesets)
    es_z = z_transform(es, axis=1)

    # export
    for_export = es_z.transpose()
    for_export.insert(for_export.shape[1], 'Verhaak classification', rnaseq_meta.loc[for_export.index, 'expression_subclass'])
    for_export.insert(for_export.shape[1], 'Wang classification', rnaseq_meta.loc[for_export.index, 'wang_classification'])
    for_export.to_excel(os.path.join(outdir, "tcga_signature_scores_and_subgroups.xlsx"))

    # boxplot by subgroup
    if class_method == 'verhaak':
        groups = rnaseq_meta.expression_subclass
        # remove any small groups (e.g. a single G-CIMP instance)
        group_list = groups.value_counts()
        group_list = group_list.index[group_list > 2]
    elif class_method == 'wang':
        groups = rnaseq_meta.wang_classification
        if not allow_multiple_classes:
            groups[rnaseq_meta.wang_classification_num_matches != 1] = None
        group_list = groups.dropna().unique()
    else:
        raise NotImplementedError("Subgroup type not recognised.")

    groups = groups.fillna('NONE')
    group_list = sorted(group_list)

    bplot = {}
    for k in genesets:
        the_data = es_z.loc[k]
        bplot[k] = collections.OrderedDict()
        for sg in group_list:
            bplot[k][sg] = the_data.loc[groups.fillna('').str.contains(sg)].values

        lbl, tmp = zip(*bplot[k].items())
        tmp = [list(t) for t in tmp]
        fig = plt.figure(num=k, figsize=(5, 4))
        ax = fig.add_subplot(111)
        sns.boxplot(data=tmp, orient='v', ax=ax, color='0.5')
        ax.set_xticklabels(lbl, rotation=45)
        ax.set_ylabel("Normalised ssGSEA score")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.png' % k.lower()), dpi=200)
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.pdf' % k.lower()))

    # is the correlation between MG / BMDM and mTOR higher in a given subgroup?
    gs = plt.GridSpec(6, 3)
    fig = plt.figure(figsize=(9, 6))
    # left panel is 2 x 2, comprising all 4 subgroups

    dict_mg = scatter_plot_with_linregress(es_z.loc['mTOR'], es_z.loc['MG'], group_list, groups)
    dict_bmdm = scatter_plot_with_linregress(es_z.loc['mTOR'], es_z.loc['BMDM'], group_list, groups)

    dict_mg['fig'].savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_mg['fig'].savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.pdf"))

    dict_bmdm['fig'].savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_bmdm['fig'].savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.pdf"))

    # check for MG / BMDM correlation
    dict_both = scatter_plot_with_linregress(es_z.loc['MG'],  es_z.loc['BMDM'], group_list, groups)
    dict_both['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_both['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_subgroup.pdf"))

    # again but across groups
    dict_both_uniform = scatter_plot_with_linregress(es_z.loc['MG'],  es_z.loc['BMDM'])
    dict_both_uniform['fig'].set_size_inches([6, 4])
    dict_both_uniform['fig'].tight_layout()
    dict_both_uniform['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation.png"), dpi=200)
    dict_both_uniform['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation.pdf"))

    # summary plot with all information
    alpha = 0.01
    slope_cmap = plt.get_cmap('RdBu_r')
    slope_norm = common.MidpointNormalize(vmin=-.5, vmax=1.2, midpoint=0.)
    slope_sm = plt.cm.ScalarMappable(cmap=slope_cmap, norm=slope_norm)

    group_list_extended = list(group_list) + ['All']

    for_plot = collections.OrderedDict([
        # ('MG-BMDM', dict_both),
        ('mTOR-MF', dict_bmdm),
        ('mTOR-MG', dict_mg)
    ])

    s, p = get_slope_and_pval(
        for_plot,
        col_order=group_list_extended,
    )

    # p_to_size = lambda t: min(150., 45 - 15 * np.log10(t))
    p_to_size = lambda t: min(500., 100 + 10 * np.log10(t) ** 2)

    x = range(len(group_list_extended))
    y_fun = lambda t: [t] * len(group_list_extended)

    fig, ax = plt.subplots(figsize=(6, 2.4))
    for i, k in enumerate(s.index):
        ax.scatter(
            x,
            y_fun(i),
            color=[slope_sm.to_rgba(t) for t in s.loc[k]],
            s=[p_to_size(t) for t in p.loc[k]],
            edgecolor='k',
            linewidth=[.5 if t > alpha else 1.5 for t in p.loc[k]]
        )
    ax.grid('off')
    ax.set_facecolor('w')
    ax.set_xticks(x)
    ax.set_xticklabels(group_list_extended)
    ax.set_yticks(range(p.shape[0]))
    ax.set_yticklabels(p.index)
    ax.set_ylim([-.5, p.shape[0] - 0.5])

    slope_sm.set_array(s.values)
    cbar = fig.colorbar(slope_sm)
    cbar.set_label('Slope')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "tcga_correlation_summary.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "tcga_correlation_summary.pdf"))

    plt_dict = line_plot_pvalues_slope(p, s, alpha=.05)
    plt_dict['fig'].set_size_inches([6., 2.5])
    plt_dict['gs'].update(bottom=0.23, top=0.9, hspace=0.4, right=0.75)
    plt_dict['fig'].savefig(os.path.join(outdir, "tcga_correlation_summary_line.png"), dpi=200)
    plt_dict['fig'].savefig(os.path.join(outdir, "tcga_correlation_summary_line.pdf"))
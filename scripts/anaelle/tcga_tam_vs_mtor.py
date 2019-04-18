from rnaseq import gsva, loader
import pandas as pd
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR
from plotting import venn

import os
import references
import datetime
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import seaborn as sns
import numpy as np
import collections
from scipy import stats
from utils import setops, output, log
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP

logger = log.get_console_logger()


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


if __name__ == "__main__":
    """
    Use the TCGA cohort to validate a link between the mTOR pathway and the proportion of microglial and macrophage
    immune infiltrate in the bulk samples.
    
    mTOR is assessed using a known set of genes.
    
    Tumour-associated bone marrow-derived macrophages (TAM-BMDM) and microglia (TAM-MG) are distinguished using a 
    human signature from Muller et al. (Genome Biol 2017). 
    """
    rnaseq_type = 'gliovis'
    remove_idh1 = True
    mtor_source = 'kegg_msigdb'  # ('kegg', 'pid', 'biocarta')
    class_method = 'wang'  # ('verhaak')
    # toggle allowing more than one class to be used
    allow_multiple_classes = True

    # kegg mtor from msigdb (hsa04150)
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

    for arr in [kegg_mtor_from_msigdb, kegg_mtor_from_kegg, pid_mtor, biocarta_mtor]:
        for i, t in enumerate(arr):
            if t in manual_gene_name_correction:
                arr[i] = manual_gene_name_correction[t]

    mtor_gs_dict = {
        'kegg_msigdb': kegg_mtor_from_msigdb,
        'kegg': kegg_mtor_from_kegg,
        'biocarta': biocarta_mtor,
        'pid': pid_mtor
    }

    if mtor_source not in mtor_gs_dict:
        raise KeyError("Unsupported mTOR geneset source. Supported options are %s." % ','.join(mtor_gs_dict.keys()))

    # load MG/BMDM signatures
    tam_signature_fn = os.path.join(GIT_LFS_DATA_DIR, 'muller_2017_tam', '13059_2017_1362_MOESM5_ESM.xlsx')
    tam_signatures = pd.read_excel(tam_signature_fn, header=0, index_col=None)

    outdir = output.unique_output_dir()

    # Venn diagram showing various options and overlap
    fig, ax = plt.subplots()
    venn.venn_diagram(*mtor_gs_dict.values(), set_labels=mtor_gs_dict.keys(), ax=ax)
    fig.tight_layout()
    ax.set_facecolor('w')
    fig.savefig(os.path.join(outdir, "venn_mtor_genesets.png"), dpi=200)

    mtor_geneset = mtor_gs_dict[mtor_source]

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
    genesets = {
        'MG': tam_signatures['MG_Markers'].dropna().values,
        'BMDM': tam_signatures['Mac_Markers'].dropna().values,
        'mTOR': mtor_geneset
    }

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

    def scatter_plot_with_linregress(x, y, group_list, groups):
        res = pd.DataFrame(index=group_list, columns=['slope', 'intercept', 'rvalue', 'pvalue', 'stderr'])
        fig, axs = plt.subplots(2, int(np.ceil(len(group_list) * 0.5)), sharex=True, sharey=True)

        for i, sg in enumerate(group_list):
            sg_idx = (groups == sg)
            this_x = x.loc[sg_idx].values.astype(float)
            this_y = y.loc[sg_idx].values.astype(float)
            lr = stats.linregress(this_x, this_y)
            res.loc[sg] = lr

            ax = axs.flat[i]
            ols_plot(
                this_x,
                this_y,
                xlim=[-3.5, 3.5],
                ax=ax
            )
            rsq = lr.rvalue ** 2
            sl = lr.slope
            pval = lr.pvalue
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
            ax.set_title(sg)

        return {
            'fig': fig,
            'axs': axs,
            'linregress': res
        }


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

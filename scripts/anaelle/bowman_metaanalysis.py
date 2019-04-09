from rnaseq.gsea import ssgsea, run_one_ssgsea
import pandas as pd
from settings import DATA_DIR_NON_GIT, GIT_LFS_DATA_DIR
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
from utils.output import unique_output_dir


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


def get_de_tissue_tumour():

    # load tables S1A, S1B: lists of DE genes in healthy vs TAM
    fn_s1a = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1A.csv')
    fn_s1b = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1B.csv')

    s1a = pd.read_csv(fn_s1a, header=0, index_col=None)
    s1b = pd.read_csv(fn_s1b, header=0, index_col=None)

    # manual corrections
    s1a.replace('AI414108', 'Igsf9b', inplace=True)
    s1a.replace('Fam101a', 'Rflna', inplace=True)
    s1b.replace('Fam176b', 'Eva1b', inplace=True)
    s1b.replace('Gm14047', 'Il1bos', inplace=True)
    s1b.replace('Gpr114', 'Adgrg5', inplace=True)


    # get DE genes in (MG vs TAM-MG) and (monocytes vs TAM-BMDM)
    mg = s1a.MG.dropna()
    bmdm = s1b.BMDM.dropna()

    # convert gene symbols to ENS ID
    fn = os.path.join(GIT_LFS_DATA_DIR, 'ensembl', 'mouse', 'mart_export.txt.gz')
    ref = pd.read_csv(fn, sep='\t', index_col=None, header=0).set_index('Gene name')

    mg_ens = ref.loc[mg.values, 'Gene stable ID'].dropna().unique()
    bmdm_ens = ref.loc[bmdm.values, 'Gene stable ID'].dropna().unique()

    # now can paste these into DAVID / gProfile
    print '\n'
    print "DE in MG / TAM-MG"
    print '\n'.join(mg_ens)
    print '\n'
    print "DE in BMDM / TAM-BMDM"
    print '\n'.join(bmdm_ens)
    print '\n'


if __name__ == "__main__":

    outdir = unique_output_dir('bowman_meta')

    # rnaseq_type = 'counts'
    rnaseq_type = 'gliovis'
    # remove_idh1 = False
    remove_idh1 = True

    # cutoff for discarding genes
    fpkm_cutoff = 1.
    fpkm_min_samples = 10

    # which list should we use?
    list_name = 'S2'
    # list_name = 'S4'

    if list_name == 'S2':
        list_cols = ['MG', 'BMDM']
    elif list_name == 'S4':
        list_cols = ['TAM MG', 'TAM BMDM', 'Core MG', 'Core BMDM']
    else:
        raise NotImplementedError("Unrecognised list: %s" % list_name)

    ### load RNA-Seq data annotated by Brennan

    if rnaseq_type == 'counts':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'counts.csv')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')
    elif rnaseq_type == 'fpkm':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count_fpkm')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'fpkm.csv')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')
    elif rnaseq_type == 'gliovis':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'gliovis_tcga_gbmlgg_expression.csv')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'gliovis_tcga_gbmlgg_meta.csv')
    else:
        raise NotImplementedError("Unrecognised rnaseq data type")


    rnaseq_dat_raw = pd.read_csv(rnaseq_dat_fn, header=0, index_col=0)
    rnaseq_meta = pd.read_csv(rnaseq_meta_fn, header=0, index_col=0)

    if rnaseq_type == 'gliovis':
        # filter only GBM
        rnaseq_meta = rnaseq_meta.loc[rnaseq_meta.Histology == 'GBM']
        rnaseq_dat_raw = rnaseq_dat_raw.transpose().loc[:, rnaseq_meta.index]
        # add meta columns for compatibility
        idh1_status = pd.Series(data='Mut', index=rnaseq_meta.index, name='idh1_status')
        idh1_status.loc[rnaseq_meta.loc[rnaseq_meta.loc[:, 'IDH_codel.subtype'] == 'IDHwt'].index] = 'WT'
        rnaseq_meta.loc[:, 'idh1_status'] = idh1_status
        rnaseq_meta.loc[:, 'expression_subclass'] = rnaseq_meta.loc[:, 'Subtype.original']

    if remove_idh1:
        # filter IDH1 mutants
        idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')

        rnaseq_meta = rnaseq_meta.loc[idh1_wt]
        rnaseq_dat = rnaseq_dat_raw.loc[:, idh1_wt.values]
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

    # plot a histogram showing distribution of expression values
    xx = rnaseq_dat.values.flatten()
    if rnaseq_type == 'gliovis':
        xx = xx[xx != -1]
    else:
        xx = xx[xx != 0.]  # required to avoid MASSIVE spike at precisely 0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.distplot(np.log10(xx + 1.), bins=200, kde=False, ax=ax)
    if rnaseq_type == 'fpkm':
        ax.set_xlabel('log10(FPKM + 1)')
    elif rnaseq_type == 'counts':
        ax.set_xlabel('log10(CPM + 1)')
    ax.set_ylabel('Density')
    ax.axvline(np.log10(fpkm_cutoff + 1.), ls='--', c='k')

    fig.savefig(os.path.join(outdir, 'rnaseq_threshold_cutoff.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'rnaseq_threshold_cutoff.pdf'))

    if rnaseq_dat.index.duplicated().any():
        print "Warning: some gene symbols are duplicated."
        print ', '.join(rnaseq_dat.index[rnaseq_dat.index.duplicated()].tolist())

    # load signatures (Bowman et al.)
    fn_s1a = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1A.csv')
    fn_s1b = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1B.csv')
    fn_s2 = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S2.csv')
    fn_s4 = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S4.csv')

    s1a = pd.read_csv(fn_s1a, header=0, index_col=None)
    s1b = pd.read_csv(fn_s1b, header=0, index_col=None)
    s2 = pd.read_csv(fn_s2, header=0, index_col=None)
    s4 = pd.read_csv(fn_s4, header=0, index_col=None)

    # mouse signature in dictionary form
    genelist_mo = {}
    if list_name == 'S4':
        the_list_mo = s4
    elif list_name == 'S2':
        the_list_mo = s2
    else:
        raise NotImplementedError("Unrecognised list: %s" % list_name)

    for c in the_list_mo.columns:
        genelist_mo[c] = the_list_mo.loc[:, c].dropna().values.tolist()

    # generate list of orthologs of the relevant gene signatures
    from scripts.agdex_mouse_human_mb_microarray import generate_ortholog_table as got
    orth = got.homologs(got.mouse_tid, got.human_tid)
    orth.set_index('gene_symbol_%d' % got.mouse_tid, inplace=True)
    # convert to Series
    orth = orth.iloc[:, 0]

    # use this to generate human gene lists
    all_genes_in_set = set()

    the_list_hu = {}
    rna_list_hu = {}

    for c in the_list_mo.columns:
        l = orth.loc[genelist_mo[c]]
        n_matched = l.dropna().size
        print "Geneset %s. Found %d orthologous genes in human from a mouse list S4 of length %d. %d dropped." % (
            c, n_matched, l.size, l.isnull().sum()
        )
        the_list_hu[c] = l.dropna().values

    for c in the_list_hu:
        this_geneset = set(the_list_hu[c].tolist()).intersection(rnaseq_dat.index)
        removed = set(the_list_hu[c].tolist()).difference(rnaseq_dat.index)
        if len(removed):
            print "%d genes were removed from RNA-Seq geneset S4 %s as they are not found in the expression data. " \
                  "%d remaining of original %d.\n" % (
                      len(removed),
                      c,
                      len(this_geneset),
                      the_list_mo[c].dropna().size
            )
        rna_list_hu[c] = list(this_geneset)
        all_genes_in_set.update(this_geneset)

    # remove genes that have no appreciable expression level
    # >=10 samples must have FPKM >= 1
    to_keep = ((rnaseq_dat > fpkm_cutoff).sum(axis=1) > fpkm_min_samples) | (rnaseq_dat.index.isin(all_genes_in_set))
    print "Keeping %d / %d genes that are sufficiently abundant" % (to_keep.sum(), to_keep.size)
    rnaseq_dat = rnaseq_dat.loc[to_keep]

    # run ssGSEA
    # we'll also store the ecdfs
    rna_es = pd.DataFrame(index=the_list_hu.keys(), columns=rnaseq_dat.columns)
    rna_ecdf_in = dict()
    rna_ecdf_out = dict()
    for s_name in rnaseq_dat.columns:
        rna_ecdf_in[s_name] = dict()
        rna_ecdf_out[s_name] = dict()
        for g_name in rna_list_hu:
            rna_es.loc[g_name, s_name], rna_ecdf_in[s_name][g_name], rna_ecdf_out[s_name][g_name] = run_one_ssgsea(
                rnaseq_dat.loc[:, s_name],
                rna_list_hu[g_name],
                alpha=0.25,
                return_ecdf=True
            )

    # scale using the Z transform
    # z = es.subtract(es.mean(axis=1), axis=0).divide(es.std(axis=1), axis=0)
    rna_z = (rna_es - rna_es.values.flatten().mean()) / rna_es.values.flatten().std()

    fig = plt.figure(num="TCGA RNA-Seq")
    ax = fig.add_subplot(111)
    for g_name in the_list_hu:
        sns.kdeplot(rna_z.loc[g_name], ax=ax)
    ax.set_xlabel("Normalised ssGSEA score")
    ax.set_ylabel("Density")
    fig.savefig(os.path.join(outdir, 'rnaseq_ssgsea_score_tcga.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'rnaseq_ssgsea_score_tcga.pdf'))

    # now split by subgroup
    subgroup_order = [
        'Classical',
        'Mesenchymal',
        'Neural',
        'Proneural',
    ]
    subgroups = rnaseq_meta.groupby('expression_subclass').groups
    if remove_idh1:
        try:
            subgroups.pop('G-CIMP')
        except Exception:
            pass
    else:
        subgroup_order += ['G-CIMP']

    # boxplot by subgroup
    # for this purpose we need to normalise by gene set, not globally
    bplot = {}
    for g_name in rna_list_hu:
        the_data = rna_es.loc[g_name]
        the_data = (the_data - the_data.mean()) / the_data.std()
        bplot[g_name] = collections.OrderedDict()
        for sg in subgroup_order:
            # bplot[g_name][sg] = rna_z.loc[g_name, subgroups[sg]].values
            bplot[g_name][sg] = the_data.loc[subgroups[sg]].values


    for col in list_cols:
        lbl, tmp = zip(*bplot[col].items())
        tmp = [list(t) for t in tmp]
        fig = plt.figure(num=col)
        ax = fig.add_subplot(111)
        sns.boxplot(data=tmp, orient='v', ax=ax)
        ax.set_xticklabels(lbl, rotation=45)
        ax.set_ylabel("Normalised ssGSEA score")
        fig = ax.figure
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.png' % col), dpi=200)
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.pdf' % col))

    plt.show()

    # ITGA4, Tmem119 and P2ry12 markers vs ssGSEA score
    # checked in orth that these are replaced by the capitalized version in humans

    def signature_vs_gene(the_gene, geneset_name):
        the_expr = rnaseq_dat.loc[the_gene]
        # Z transform the signature scores for this gene set
        the_signature = rna_es.loc[geneset_name]
        the_signature = (the_signature - the_signature.mean()) / the_signature.std()
        # ensure the ordering is the same
        the_signature = the_signature.loc[the_expr.index]
        lr = stats.linregress(the_signature.astype(float), np.log2(the_expr + 1))
        x_lr = np.array([the_signature.min(), the_signature.max()])
        y_lr = lr.intercept + lr.slope * x_lr

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(the_signature, np.log2(the_expr + 1))
        ax.plot(x_lr, y_lr, 'k--')
        ax.set_xlabel('Normalised ssGSEA score')
        ax.set_ylabel('log2(%s)' % the_gene)
        fig.tight_layout()
        return ax

    for col in list_cols:
        for gene in ['ITGA4', 'TMEM119', 'P2RY12']:
            ax = signature_vs_gene(gene, col)
            fig = ax.figure
            fig.savefig(os.path.join(outdir, '%s_ssgsea_vs_%s.png' % (col, gene)), dpi=200)
            fig.savefig(os.path.join(outdir, '%s_ssgsea_vs_%s.pdf' % (col, gene)))

    # if mutants are present, recreate Fig 5F: TAM BMDM signature by IDH1 status
    if not remove_idh1:
        lbl = ['Mutant', 'WT']
        for g_name in list_cols:
            the_data = rna_es.loc[g_name]
            the_data = (the_data - the_data.mean()) / the_data.std()
            idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')
            tmp = [
                the_data.loc[~idh1_wt.values].values.tolist(),
                the_data.loc[idh1_wt.values].values.tolist(),
            ]
            fig = plt.figure(num="%s_vs_idh1" % g_name, figsize=(4.3, 6.7))
            ax = fig.add_subplot(111)
            ax.boxplot(tmp, widths=0.9)
            ax.set_xticklabels(lbl, rotation=45)
            ax.set_ylabel("Normalised ssGSEA score")
            fig = ax.figure
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, '%s_ssgsea_by_idh1_tcga.png' % g_name), dpi=200)
            fig.savefig(os.path.join(outdir, '%s_ssgsea_by_idh1_tcga.pdf' % g_name))

    yy = np.log2(rnaseq_dat.loc['ITGA4'].values + 1)

    fig = plt.figure(figsize=(3.4, 5.5))
    ax = fig.add_subplot(111)
    sns.boxplot(yy, orient='v', ax=ax)
    sns.stripplot(yy, jitter=True, orient='v', ax=ax, color=".3")
    ax.set_ylabel("log2(ITGA4 + 1)")
    ax.set_xlim([-2, 2])
    fig.savefig(os.path.join(outdir, 'tcga_itga4_boxplot.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'tcga_itga4_boxplot.pdf'))

    sns.boxplot()

    # define the mTOR signature(s)
    mtor_geneset_ad = [
        'EIF3H',
        'EIF4EBP1',
        'HIF1A',
        'PIK3R5',
        'PLD3',
        'PRKCA',
        'PRR5L',
        'RHOC',
        'RPS2',
        'RPS5',
        'RPS7',
        'RPS8',
        'RPS10',
        'RPS12',
        'RPS13',
        'RPS15',
        'RPS16',
        'RPS17',
        'RPS18',
        'RPS19',
        'RPS20',
        'RPS21',
        'RPS23',
        'RPS24',
        'RPS25',
        'RPS26',
        'RPS28',
        'RPS27A',
        'RPS27L',
        'RPS4Y1',
        'RPS6KA4',
        'RPTOR',
    ]

    # this from KEGG pathway hsa04150 (looked up via Entrez IDs)
    mtor_geneset_kegg = [
        'AKT1', 'AKT2', 'BRAF', 'EIF4B', 'EIF4E', 'EIF4EBP1', 'VEGFD', 'MTOR', 'HIF1A', 'IGF1', 'INS', 'PDPK1', 'PGF',
        'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PRKAA1', 'PRKAA2', 'MAPK1', 'MAPK3', 'RHEB',
        'RPS6', 'RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KB1', 'RPS6KB2', 'STK11', 'TSC2', 'VEGFA', 'VEGFB', 'VEGFC',
        'ULK1', 'PIK3R3', 'EIF4E2', 'ULK2', 'AKT3', 'PIK3R5', 'ULK3', 'RPS6KA6', 'CAB39', 'DDIT4', 'RPTOR', 'MLST8',
        'CAB39L', 'STRADA', 'RICTOR', 'EIF4E1B', 'TSC1'
    ]

    # remove any genes not in the data and report
    diff_ad = pd.Index(mtor_geneset_ad).difference(rnaseq_dat.index)
    if len(diff_ad):
        print "%d genes in the geneset mTOR (AD) are not in the data and will be removed: %s" % (
            len(diff_ad),
            ', '.join(diff_ad.tolist())
        )
        for t in diff_ad:
            mtor_geneset_ad.remove(t)

    diff_kegg = pd.Index(mtor_geneset_kegg).difference(rnaseq_dat.index)
    if len(diff_kegg):
        print "%d genes in the geneset mTOR (KEGG) are not in the data and will be removed: %s" % (
            len(diff_kegg),
            ', '.join(diff_kegg.tolist())
        )
        for t in diff_kegg:
            mtor_geneset_kegg.remove(t)

    # compute scores associated with each of these
    mtor_ad_scores = pd.Series([
        run_one_ssgsea(
            rnaseq_dat.loc[:, s_name],
            mtor_geneset_ad,
            alpha=0.25,
            return_ecdf=False
        ) for s_name in rnaseq_dat.columns
    ], index=rnaseq_dat.columns)

    mtor_kegg_scores = pd.Series([
        run_one_ssgsea(
            rnaseq_dat.loc[:, s_name],
            mtor_geneset_kegg,
            alpha=0.25,
            return_ecdf=False
        ) for s_name in rnaseq_dat.columns
    ], index=rnaseq_dat.columns)

    # establish an index to maintain ordering (unnecessary?)
    idx = rnaseq_dat.columns

    # extract three relevant scores (for all TCGA data)
    x = rna_es.loc[list_cols[0], idx]  # MG
    y = rna_es.loc[list_cols[1], idx]  # BMDM
    z = mtor_kegg_scores.loc[idx]

    # standardise each
    x = (x - x.mean()) / x.std()
    y = (y - y.mean()) / y.std()
    z = (z - z.mean()) / z.std()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colours = ['b', 'r', 'g', 'k']
    for sg, c in zip(subgroup_order, colours):
        sg_idx = (rnaseq_meta.expression_subclass == sg).loc[idx]
        ax.scatter(x.loc[sg_idx], y.loc[sg_idx], z.loc[sg_idx], c=c, marker='o', label=sg)
    ax.set_xlabel(list_cols[0])
    ax.set_ylabel(list_cols[1])
    ax.set_zlabel('mTOR KEGG')

    # run linear regression

    # is the correlation between MG / BMDM and mTOR higher in a given subgroup?
    fig, axs = plt.subplots(2, int(np.ceil(len(subgroup_order) * 0.5)), sharex=True, sharey=True)
    lr_mtor_mg = pd.DataFrame(index=subgroup_order, columns=['slope', 'intercept', 'rvalue', 'pvalue', 'stderr'])
    for i, sg in enumerate(subgroup_order):
        sg_idx = (rnaseq_meta.expression_subclass == sg).loc[idx]
        lr_mtor_mg.loc[sg] = stats.linregress(x.loc[sg_idx].values.tolist(), z.loc[sg_idx].values.tolist())
        ax = axs.flat[i]
        ols_plot(z.loc[sg_idx].values, x.loc[sg_idx].values.astype(float), xlim=[-3.5, 3.5], ax=ax)
        rsq = lr_mtor_mg.loc[sg].rvalue ** 2
        sl = lr_mtor_mg.loc[sg].slope
        pval = lr_mtor_mg.loc[sg].pvalue
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

    fig.savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.pdf"))

    # all
    sg = 'all'
    lr_mtor_mg.loc[sg] = stats.linregress(x.values.tolist(), z.values.tolist())
    _, ax = ols_plot(z.values, x.values.astype(float), xlim=[-3.5, 3.5])
    fig = ax.figure
    rsq = lr_mtor_mg.loc[sg].rvalue ** 2
    sl = lr_mtor_mg.loc[sg].slope
    pval = lr_mtor_mg.loc[sg].pvalue
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

    fig.savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_all.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_all.pdf"))



    fig, axs = plt.subplots(2, int(np.ceil(len(subgroup_order) * 0.5)), sharex=True)
    lr_mtor_bmdm = pd.DataFrame(index=subgroup_order, columns=['slope', 'intercept', 'rvalue', 'pvalue', 'stderr'])
    for i, sg in enumerate(subgroup_order):
        sg_idx = (rnaseq_meta.expression_subclass == sg).loc[idx]
        lr_mtor_bmdm.loc[sg] = stats.linregress(y.loc[sg_idx].values.tolist(), z.loc[sg_idx].values.tolist())
        ax = axs.flat[i]
        ols_plot(z.loc[sg_idx].values, y.loc[sg_idx].values.astype(float), xlim=[-3.5, 3.5], ax=ax)
        rsq = lr_mtor_bmdm.loc[sg].rvalue ** 2
        sl = lr_mtor_bmdm.loc[sg].slope
        pval = lr_mtor_bmdm.loc[sg].pvalue
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

    fig.savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.pdf"))

    # all
    sg = 'all'
    lr_mtor_bmdm.loc[sg] = stats.linregress(y.values.tolist(), z.values.tolist())
    _, ax = ols_plot(z.values, y.values.astype(float), xlim=[-3.5, 3.5])
    fig = ax.figure
    rsq = lr_mtor_bmdm.loc[sg].rvalue ** 2
    sl = lr_mtor_bmdm.loc[sg].slope
    pval = lr_mtor_bmdm.loc[sg].pvalue
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

    fig.savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_all.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_all.pdf"))

    # check for MG / BMDM correlation
    _, ax = ols_plot(y.values.astype(float), x.values.astype(float), xlim=(-3, 3))
    lr = stats.linregress(x.values.tolist(), y.values.tolist())
    fig = ax.figure
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
    fig.savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_all.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_all.pdf"))

    # export data for AD
    to_export = rnaseq_meta.copy()
    to_export.insert(0, 'mg_score_z', x.loc[to_export.index])
    to_export.insert(0, 'bmdm_score_z', y.loc[to_export.index])
    to_export.insert(0, 'mtor_score_z', z.loc[to_export.index])
    to_export.to_excel(os.path.join(outdir, 'bowman_data_with_scores.xlsx'))
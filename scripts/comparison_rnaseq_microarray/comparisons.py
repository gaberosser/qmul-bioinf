import load_references, load_illumina_data, load_rnaseq_data
from microarray.process import aggregate_by_probe_set
import numpy as np
from scipy import stats, cluster
import collections
import pandas as pd

MB_GROUPS = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EYS', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
)

REF_GROUPS = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
)


def annotate_by_MB_group(df):
    """
    Annotate the supplied dataframe by adding a new column, 'mb_group', that lists the MB group the gene is in
    or null if none.
    :param df:
    :return:
    """
    res = df.copy()
    res['mb_group'] = [np.nan] * len(df.index)

    for grp, arr in MB_GROUPS:
        res.loc[res.index.isin(arr), 'mb_group'] = grp

    return res


def plot_microarray_mb_gene_expression(mb_samples=('ICb1299-III', 'ICb1299-IV')):
    """
    Produce 2 publications:
    1) bar chart subplots showing the absolute normed intensity values for the MB-implicated genes in both
    healthy and MB samples.
    2) bar chart subplots showing the log2 fold change in those same genes
    :param mb_samples: Iterable with the sample names to use in the MB data
    :return:
    """
    from plotting import bar
    plt = bar.plt
    METHOD = 'median'
    HEALTHY_SAMPLE_NAMES = [
        'NT-1197',
        'NCb-1',
        'NCb-2',
        'A911105',
        'A508112',
        'A508285',
    ]
    HEALTHY_SAMPLE_NAMES += [t + '-R' for t in HEALTHY_SAMPLE_NAMES]

    # load full microarray data
    marray_data = load_illumina_data.load_normed_microarray_data(pval=0.01)
    # replace null with zero
    marray_data.fillna(value=0., inplace=True)
    # load probe set definitions
    probe_set = load_illumina_data.load_illumina_array_library()
    marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
    marray_by_gene = aggregate_by_probe_set(marray_ann, method=METHOD)
    mb_sample_names = list(mb_samples) + [t + '-R' for t in mb_samples]

    # pick out samples and aggregate
    mb = marray_by_gene.loc[:, mb_sample_names].mean(axis=1)
    he = marray_by_gene.loc[:, HEALTHY_SAMPLE_NAMES].mean(axis=1)

    # add MB group column
    marray_by_gene = annotate_by_MB_group(marray_by_gene)

    mb_grouped = dict(
        [(grp, mb.loc[marray_by_gene.mb_group == grp]) for grp, arr in MB_GROUPS]
    )
    he_grouped = dict(
        [(grp, he.loc[marray_by_gene.mb_group == grp]) for grp, arr in MB_GROUPS]
    )

    # figure 1: absolute TPM

    data = collections.OrderedDict([
        (grp, [he_grouped[grp], mb_grouped[grp]]) for grp, _ in REF_GROUPS
    ])
    fig, axs = bar.multi_grouped_bar_chart(data, xlabel_coords=(0.5, -.21))

    axs[-1].legend(['Healthy cerebellum', 'MB'])
    axs[0].set_ylabel('Normed intensity')
    ylim = list(axs[-1].get_ylim())
    ylim[0] = -1e-6
    axs[-1].set_ylim(ylim)
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.95, wspace=0.08, hspace=0.)

    # figure 2: log2 fold change

    LOG_MIN = -7
    LOG_MAX = 7
    log_fold_diff = {}
    for grp in mb_grouped:
        t = np.log2(mb_grouped[grp] / he_grouped[grp])
        log_fold_diff[grp] = t

    data = collections.OrderedDict([
        (grp, [log_fold_diff[grp]]) for grp, _ in REF_GROUPS
    ])
    fig, axs = bar.multi_grouped_bar_chart(data, xlabel_coords=(0.5, -.21), ylim=[LOG_MIN, LOG_MAX], colours=['gray'])

    axs[0].set_ylabel('Log2 fold change')
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.95, wspace=0.08, hspace=0.)


def plot_rna_seq_mb_gene_expression():
    """
    Produce 2 publications:
    1) bar chart subplots showing the absolute TPM values for RNA-Seq in the MB-implicated genes in both
    healthy and MB samples.
    2) bar chart subplots showing the log2 fold change in those same genes
    :return:
    """
    from plotting import bar
    plt = bar.plt
    # load 2015 RNA-Seq gene activity counts
    # ignore the confidence intervals, which are wrong
    rna_z = load_rnaseq_data.load_rnaseq_cufflinks_gene_count_data(unit='tpm')

    # load Allen RNA-Seq cerebellum activity TPM
    _, rna_a, meta = load_references.load_cerebellum_rnaseq_reference_data()

    # annotate both by MB group
    rna_z['mb_group'] = [np.nan] * len(rna_z.index)
    rna_a['mb_group'] = [np.nan] * len(rna_a.index)

    for grp, arr in MB_GROUPS:
        rna_z.loc[rna_z.index.isin(arr), 'mb_group'] = grp
        rna_a.loc[rna_a.index.isin(arr), 'mb_group'] = grp

    rna_z_mean = rna_z.mean(axis=1)
    rna_a_mean = rna_a.mean(axis=1)

    rna_a_grouped = dict(
        [(grp, rna_a_mean.loc[rna_a.mb_group == grp]) for grp, arr in MB_GROUPS]
    )
    rna_z_grouped = dict(
        [(grp, rna_z_mean.loc[rna_z.mb_group == grp]) for grp, arr in MB_GROUPS]
    )

    # figure 1: absolute TPM

    data = collections.OrderedDict([
        (grp, [rna_a_grouped[grp], rna_z_grouped[grp]]) for grp, _ in REF_GROUPS
    ])
    fig, axs = bar.multi_grouped_bar_chart(data, xlabel_coords=(0.5, -.21))

    axs[-1].legend(['Healthy cerebellum', 'MB'])
    axs[0].set_ylabel('TPM')
    ylim = list(axs[-1].get_ylim())
    ylim[0] = -1e-6
    axs[-1].set_ylim(ylim)
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.95, wspace=0.08, hspace=0.)

    # figure 2: log2 fold change

    LOG_MIN = -7
    LOG_MAX = 7
    log_fold_diff = {}
    for grp in rna_z_grouped:
        t = np.log2(rna_z_grouped[grp] / rna_a_grouped[grp])
        log_fold_diff[grp] = t

    data = collections.OrderedDict([
        (grp, [log_fold_diff[grp]]) for grp, _ in REF_GROUPS
    ])
    fig, axs = bar.multi_grouped_bar_chart(data, xlabel_coords=(0.5, -.21), ylim=[LOG_MIN, LOG_MAX], colours=['gray'])

    axs[0].set_ylabel('Log2 fold change')
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.95, wspace=0.08, hspace=0.)


def comparison_healthy_vs_mb_microarray_northcott_genes():
    from matplotlib import rc, pyplot as plt, gridspec as gridspec
    import seaborn as sns
    plt.interactive(True)
    sns.set_style('white')

    METHOD = 'median'
    HEALTHY_SAMPLE_NAMES = [
        'NT1197',
        'NCb1',
        'NCb2',
        'A911105',
        'A508112',
        'A508285',
    ]
    SAMPLE_GROUPS = (
        ('WNT', ('Pt1140', 'ICb1140-II', 'ICb1140-III', 'Pt1192', 'ICb1192-I', 'ICb1192-III', 'ICb1192-V')),
        ('SSH', ('Pt1338', 'ICb1338-I', 'ICb1338-III', 'ICb984-I', 'ICb984-III', 'ICb984-V')),
        ('Group C', (
            'ICb1197-I',
            'ICb1197-III',
            'Pt1494',
            'ICb1494-I',
            'ICb1494-III',
            'ICb1494-V',
            'Pt1572',
            'ICb1572-I',
            'ICb1572-III',
            'ICb1572-V',
            'Pt1595',
            'ICb1595-I',
            'ICb1595-III',
        )),
        ('Group D', (
            'Pt1078',
            'ICb1078-I',
            'ICb1078-III',
            'ICb1078-V',
            'Pt1299',
            'ICb1299-I',
            'ICb1299-III',
            'ICb1299-IV',
            'Pt1487',
            'ICb1487-I',
            'ICb1487-III',
        )),
    )

    all_northcott_patients = []
    for grp, arr in SAMPLE_GROUPS:
        all_northcott_patients.extend(arr)

    # load full microarray data
    marray_data = load_illumina_data.load_normed_microarray_data(pval=0.05)

    probe_set = load_illumina_data.load_illumina_array_library()
    marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
    marray_by_gene = load_illumina_data.aggregate_by_probe_set(marray_ann, method=METHOD)

    all_mb_genes = []
    for _, arr in REF_GROUPS:
        all_mb_genes.extend(arr)

    # standardised scores by gene
    marray_by_gene_stand = (
        marray_by_gene.subtract(marray_by_gene.mean(axis=1), axis=0)
            .divide(marray_by_gene.std(axis=1), axis=0)
    )




    # take mean over repeats
    for sn in load_illumina_data.SAMPLE_NAMES:
        marray_by_gene.loc[:, sn] = marray_by_gene.loc[:, [sn, sn + '-R']].mean(axis=1)
        marray_by_gene_stand.loc[:, sn] = marray_by_gene_stand.loc[:, [sn, sn + '-R']].mean(axis=1)
    marray_by_gene = marray_by_gene.loc[:, load_illumina_data.SAMPLE_NAMES]
    marray_by_gene_stand = marray_by_gene_stand.loc[:, load_illumina_data.SAMPLE_NAMES]

    VMAX = 15000

    # v1: absolute values

    fig = plt.figure(figsize=[5, 8])
    gs = gridspec.GridSpec(2, len(REF_GROUPS),
                           height_ratios=[1, 12],
                           width_ratios=[len(arr) for _, arr in REF_GROUPS])
    gs.update(
        left=0.2,
        right=0.95,
        top=0.95,
        bottom=0.15,
        wspace=0.,
        hspace=0.1)
    cbar_kws = {"orientation": "horizontal"}

    for i, (grp, arr) in enumerate(REF_GROUPS):
        ax = fig.add_subplot(gs[1:, i])
        if i == (len(REF_GROUPS) - 1):
            cbar = True
            cbar_ax = fig.add_subplot(gs[0, :])
        else:
            cbar = False
            cbar_ax = None
        sns.heatmap(
            marray_by_gene.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose(),
            ax=ax,
            vmin=0,
            vmax=VMAX,
            square=True,
            cmap='Reds',
            cbar=cbar,
            cbar_ax=cbar_ax,
            cbar_kws=cbar_kws
        )
        ax.set_xticklabels(arr, rotation=90)
        if i == 0:
            plt.yticks(rotation=0)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(.5, -.15)
    cbar_ax.set_title('$\log_2$(Normalised intensity)')

    fig.savefig("marray_all_samples_mb_gene_activity_heatmap.png", dpi=200)
    fig.savefig("marray_all_samples_mb_gene_activity_heatmap.pdf", dpi=200)

    # v2: log2(absolute) values

    fig = plt.figure(figsize=[5, 8])
    gs = gridspec.GridSpec(2, len(REF_GROUPS),
                           height_ratios=[1, 12],
                           width_ratios=[len(arr) for _, arr in REF_GROUPS])
    gs.update(
        left=0.2,
        right=0.95,
        top=0.95,
        bottom=0.15,
        wspace=0.,
        hspace=0.1)
    cbar_kws = {"orientation": "horizontal"}

    for i, (grp, arr) in enumerate(REF_GROUPS):
        ax = fig.add_subplot(gs[1:, i])
        if i == (len(REF_GROUPS) - 1):
            cbar = True
            cbar_ax = fig.add_subplot(gs[0, :])
        else:
            cbar = False
            cbar_ax = None
        sns.heatmap(
            np.log2(marray_by_gene.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose()),
            ax=ax,
            vmin=0,
            vmax=np.ceil(np.log2(VMAX)),
            square=True,
            cmap='Reds',
            cbar=cbar,
            cbar_ax=cbar_ax,
            cbar_kws=cbar_kws
        )
        ax.set_xticklabels(arr, rotation=90)
        if i == 0:
            plt.yticks(rotation=0)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(.5, -.15)
    cbar_ax.set_title('$\log_2$(Normalised intensity)')

    fig.savefig("marray_all_samples_mb_gene_log_activity_heatmap.png", dpi=200)
    fig.savefig("marray_all_samples_mb_gene_log_activity_heatmap.pdf", dpi=200)

    # v3: standardised by score values

    fig = plt.figure(figsize=[5, 8])
    gs = gridspec.GridSpec(2, len(REF_GROUPS),
                           height_ratios=[1, 12],
                           width_ratios=[len(arr) for _, arr in REF_GROUPS])
    gs.update(
        left=0.2,
        right=0.95,
        top=0.95,
        bottom=0.15,
        wspace=0.,
        hspace=0.1)
    cbar_kws = {"orientation": "horizontal"}

    for i, (grp, arr) in enumerate(REF_GROUPS):
        ax = fig.add_subplot(gs[1:, i])
        if i == (len(REF_GROUPS) - 1):
            cbar = True
            cbar_ax = fig.add_subplot(gs[0, :])
        else:
            cbar = False
            cbar_ax = None
        sns.heatmap(
            marray_by_gene_stand.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose(),
            ax=ax,
            vmin=-1,
            vmax=1.,
            square=True,
            cmap='RdBu_r',
            cbar=cbar,
            cbar_ax=cbar_ax,
            cbar_kws=cbar_kws
        )
        ax.set_xticklabels(arr, rotation=90)
        if i == 0:
            plt.yticks(rotation=0)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(.5, -.15)
    cbar_ax.set_title('Standardised score by gene')

    fig.savefig("marray_all_samples_mb_standardised_by_gene_activity_heatmap.png", dpi=200)
    fig.savefig("marray_all_samples_mb_standardised_by_gene_activity_heatmap.pdf", dpi=200)


if __name__ == '__main__':

    PVAL_MIN = 0.01
    from matplotlib import pyplot as plt
    plt.interactive(True)

    expr, meta = load_references.load_cerebellum_microarray_reference_data()
    # convert to by-gene data
    ref_gene_expr = load_references.microarray_gene_markers(expr)

    marray_data_normed = load_illumina_data.load_normed_microarray_data(pval=PVAL_MIN)
    cat = load_illumina_data.load_illumina_array_library()
    # convert to by-gene data
    m1 = load_illumina_data.convert_microarray_to_gene_activity(marray_data_normed, cat)

    # find intersecting gene set
    gene_set = m1.index.intersection(ref_gene_expr.index)

    m2 = m1.copy()
    m2 /= m2.sum(axis=0)

    # m2 = m1.loc[gene_set]
    # m2 /= m2.sum(axis=0)

    ref = ref_gene_expr.loc[gene_set]
    ref_mean = ref.mean(axis=1)
    ref_mean /= ref_mean.sum()

    # get groups of genes
    fig, axs = plt.subplots(ncols=len(MB_GROUPS), sharey=True, figsize=(10, 6))

    for i in range(len(MB_GROUPS)):
        grp, mb_arr = MB_GROUPS[i]
        _, ref_arr = REF_GROUPS[i]
        ax = axs[i]
        baseline = ref_mean[ref_mean.index.isin(ref_arr)]
        topline = m1[m1.index.isin(mb_arr)]

        ax = axs[i]
        topline.transpose().boxplot(ax=ax, rot=90)
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(0.5, -.21)
        if i == 0:
            ax.set_ylabel('Normalised intensity')

    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.99, wspace=0.05, hspace=0.)

    # # excise gene IDs
    # ex = expr.copy()
    # ex = ex.ix[:, 2:]



    # create triangular matrix of lin regression results
    # n = ex.columns.size
    # ex.columns = pd.Int64Index(range(n))

    # res = pd.DataFrame(columns=['i', 'j', 'linreg_intercept', 'linreg_slope', 'linreg_rsq', 'p_wilcox'])
    # for i in range(n):
    #     for j in range(i + 1, n):
    #         a = all_expr[[i, j]].dropna(axis=0)
    #         a0 = a[i]
    #         a1 = a[j]
    #         lr = stats.linregress(a0, a1)
    #         _, pwilcox = stats.wilcoxon(a0, a1)
    #         this_row = pd.Series(data=(i, j, lr.intercept, lr.slope, lr.rvalue ** 2, pwilcox), index=res.columns)
    #         res = res.append(this_row, ignore_index=True)
    #

    # with open('cerebellum_pairwise_comparison.pkl', 'rb') as f:
    #     res = dill.load(f)

    # cluster using Rsq as distance
    # Z = cluster.hierarchy.linkage(res.linreg_rsq)  # Z is a linkage matrix
    # dendro = cluster.hierarchy.dendrogram(Z, orientation='left')

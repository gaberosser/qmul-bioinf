import csv
import pandas as pd
import numpy as np
import collections
import re
import os
from scipy import stats
from microarray import aggregate_by_probe_set
from settings import DATA_DIR

RAW_MICROARRAY_TXT = os.path.join(DATA_DIR, 'microarray_GSE28192/raw/GSE28192_series_matrix.txt')
# RAW_MICROARRAY_TXT = 'data/microarray_GSE28192/raw/GSE28192_series_matrix.txt'
NORMED_MICROARRAY_DIR = os.path.join(DATA_DIR, 'microarray_GSE28192')
MICROARRAY_LIBRARY = os.path.join(DATA_DIR, 'microarray_GSE28192/probe_set/GPL6102-11574.txt')

SAMPLE_NAMES = [
    'A508112',
    'A508285',
    'A911105',
    'NCb1',
    'NCb2',
    'NT1197',
    'A602126',
    'A602127',
    'NFb',
    'Pt1078',
    'ICb1078-I',
    'ICb1078-III',
    'ICb1078-V',
    'Pt1140',
    'ICb1140-II',
    'ICb1140-III',
    'Pt1192',
    'ICb1192-V',
    'ICb1192-I',
    'ICb1192-III',
    'ICb1197-I',
    'ICb1197-III',
    'Pt1299',
    'ICb1299-I',
    'ICb1299-III',
    'ICb1299-IV',
    'Pt1338',
    'ICb1338-I',
    'ICb1338-III',
    'Pt1487',
    'ICb1487-I',
    'ICb1487-III',
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
    'ICb984-I',
    'ICb984-III',
    'ICb984-V',
]


def load_illumina_array_library():
    with open(MICROARRAY_LIBRARY, 'rb') as f:
        while True:
            line = f.readline()
            if re.match(r'ID', line):
                col_names = line.strip().split('\t')
                break
        c = csv.reader(f, delimiter='\t')
        all_data = list(c)
    return pd.DataFrame([t[1:] for t in all_data], columns=col_names[1:], index=[t[0] for t in all_data])


def load_normed_microarray_data(sample_names=None, pval=0.01, return_pvals=False):
    """
    Load normalised microarray data, optionally filtering by detection P-value
    These datasets were downloaded individually from the GEO site
    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28192
    They have been BG-subtracted and normed
    COLUMNS: VALUE, Avg_NBEADS, BEAD_STDERR, Detection Pval
    I'm assuming the VALUE is the key quantity
    :param sample_names: If supplied, this is an iterable containing the names of the samples to load
    :param return_pvals: If True, return a dataframe containing pvalues in addition to the intensities
    """

    if sample_names is None:
        sample_names = list(SAMPLE_NAMES)
        # add repetition runs
        for i in range(len(sample_names)):
            sample_names.append(sample_names[i] + '-R')

    res = pd.DataFrame()
    pvals = pd.DataFrame()

    for fn in sample_names:
        ff = os.path.join(NORMED_MICROARRAY_DIR, fn) + '.gz'
        this_df = pd.read_csv(
            ff,
            sep='\t',
            header=None,
            index_col=0,
            skiprows=6,
            usecols=[0, 1, 4],
            skip_blank_lines=True,
            skipfooter=1,
            engine='python',  # required to allow use of skipfooter
        )
        this_df.columns = [fn, 'pval']
        # mark non-sig results as null if requested
        if pval is not None:
            this_df[fn][this_df['pval'] > pval] = np.nan
        # remove pvalue column and append intensity
        res = pd.concat([res, this_df[fn]], axis=1)
        # also record the pvalue if required
        if return_pvals:
            this_pv = this_df.loc[:, 'pval']
            this_pv.name = fn
            pvals = pd.concat([pvals, this_pv], axis=1)

    if return_pvals:
        return res, pvals
    else:
        return res


def add_annotation_column(df, probe_set, from_col_name, to_col_name):
    """
    Return a copy of the supplied dataframe with an additional column called to_col_name.
    The index in the dataframe is ILM_ Illumina codes. Use the supplied catalogue to perform a lookup with from_col_name.
    :param df: Dataframe indexed by Illumina ID
    :param probe_set: The full Illumina probe set catalogue
    :param from_col_name:
    :param to_col_name:
    :return: Dataframe
    """
    res = df.copy()
    res.insert(0, to_col_name, probe_set[from_col_name])
    res.loc[res[to_col_name] == '', to_col_name] = np.nan
    return res


def add_gene_symbol_column(df, probe_set):
    """
    Return a copy of the supplied dataframe with an additional column called gene_symbol.
    The index in the dataframe is ILM_ Illumina codes. Use the supplied catalogue to perform a gene symbol lookup.
    :param df: Dataframe indexed by Illumina ID
    :param probe_set: The full Illumina probe set catalogue
    :return: Dataframe
    """
    return add_annotation_column(df, probe_set, 'Symbol', 'gene_symbol')


def add_entrez_column(df, probe_set):
    """
    Return a copy of the supplied dataframe with an additional column called entrez_id.
    The index in the dataframe is ILM_ Illumina codes. Use the supplied catalogue to perform ae Entrez ID lookup.
    :param df: Dataframe indexed by Illumina ID
    :param probe_set: The full Illumina probe set catalogue
    :return: Dataframe
    """
    return add_annotation_column(df, probe_set, 'Entrez_Gene_ID', 'entrez_id')


def load_raw_microarray_data():
    """
    Load the raw, unnormalised data.
    :return:
    """
    meta_map = {
        'Sample_title': 'title',
        'Sample_geo_accession': 'accession',
        'Sample_source_name_ch1': 'description',
        'Sample_characteristics_age': 'age',
        'Sample_characteristics_sex': 'sex',
        'Sample_characteristics_initial_diagnosis': 'initial_diagnosis',
        'Sample_characteristics_final_diagnosis': 'final_diagnosis',
    }
    # meta = collections.defaultdict(list)
    meta = pd.DataFrame()

    # raw = []
    # gene_ids = []
    with open(RAW_MICROARRAY_TXT, 'rb') as f:
        while True:
            line = f.readline().strip('\n')
            if len(line) == 0:
                continue
            header = re.match(r'^!(?P<hd>[^\t]*)', line).group('hd')
            if header in meta_map:
                meta[meta_map[header]] = [t.strip('"') for t in line.split('\t')[1:]]
            if line == '!series_matrix_table_begin':
                break
        # the intensity data start here
        # first line is GM numbers repeated
        raw = pd.read_csv(f, sep='\t', header=0, skip_footer=1, engine='python')

    # set index and convert accession -> titles
    tmp = meta.copy()
    tmp.set_index('accession', inplace=True)
    raw.set_index('ID_REF', inplace=True)
    raw.columns = tmp.loc[raw.columns, 'title']

    # set the index of meta
    meta.set_index('title', inplace=True)

    # standardise title formats
    titles = []
    for t in meta.index:
        tn = re.sub('^[iI][cC](?P<i>[0-9])', 'ICb\g<i>', t)
        tn = re.sub('^[iI][cC][bB]', 'ICb', tn)
        tn = re.sub('^[pP][tT]', 'Pt', tn)
        tn = re.sub(r'(?P<a>Pt|ICb|NT|NCb)-(?P<b>[0-9]*)', '\g<a>\g<b>', tn)
        tn = re.sub('(?P<a>-i+)', lambda x: x.group('a').upper(), tn)
        tn = re.sub(r'MB', '', tn)
        titles.append(tn)

    raw.columns = titles
    meta.index = titles

    return raw, meta


def convert_microarray_to_gene_activity(marray_data, ilm_cat, genes=None, method='median'):
    """
    Compute the gene-centric activity from microarray data. Many genes are represented on multiple probes.
    Following Zhao et al (PLoS One 2014), we take the most active gene whenever we encounter ambiguity (method='max').
    Also can take the sum, or mean, etc... specified by the method param.
    Optionally permit providing a list of genes of interest, otherwise all are computed.
    :param marray_data: pandas dataframe containing normalised intensity data
    :param ilm_cat: Illumina catalogue, used for conversion (pandas dataframe)
    :param genes: Optional list of genes to use, otherwise use all.
    :param method: 'max', 'mean', 'sum'
    :return:
    """
    marray_ann = add_gene_symbol_column(marray_data, ilm_cat)
    return aggregate_by_probe_set(marray_ann, method=method)


def microarray_activity_in_groups(marray_data_by_gene):
    groups = (
        ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
        ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
        ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
        ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
    )

    subdata = {}
    for grp, arr in groups:
        this_res = marray_data_by_gene.loc[arr, :]
        this_res[this_res.isnull()] = 0.
        subdata[grp] = this_res

    return subdata


def plot_microarray_activity_in_groups(marray_data_by_gene):
    MARRAY_SAMPLES = [
        'Pt1299',
        'ICb1299-I',
        'ICb1299-III',
        'ICb1299-IV',
    ]
    COLOURS = ['r', 'b', 'g', 'k']
    from matplotlib import pyplot as plt, rc
    plt.interactive(True)
    rc('text', usetex=True)
    width = 0.8

    marray_groups = microarray_activity_in_groups(marray_data_by_gene)

    fig, axes = plt.subplots(ncols=len(marray_groups), sharey=True)
    fig.subplots_adjust(wspace=0.05, right=0.99, top=0.99, bottom=0.15)

    first = True
    for ax, grp in zip(axes, marray_groups.keys()):
        # take mean of the two repeats in each case
        y = pd.DataFrame(columns=MARRAY_SAMPLES)
        logy = y.copy()
        for m in MARRAY_SAMPLES:
            t = marray_groups[grp].loc[:, [m, m + '-R']].mean(axis=1)
            t.name = m
            y.loc[:, m] = t
            logy.loc[:, m] = np.log2(t)

        ng, ns = y.shape
        b = []

        for i in range(ng):
            heights = logy.ix[i]
            lefts = np.linspace(i, i + width, ns + 1)[:-1]
            b.append(
                ax.bar(lefts, heights, width=width / float(ns), color=COLOURS)
            )

        ax.set(
            xticks=np.arange(ng) + width / 2.,
            xticklabels=y.index,
            xlabel=grp,
            xlim=[width - 1., ns],
        )
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=30, fontsize=12)
        # set fixed xlabel coords to make appearance uniform
        ax.xaxis.set_label_coords(0.5, -.1)

        if first:
            ax.set(ylabel='Log(normalised intensity)')
            first = False

    ax.legend(b[0][:], MARRAY_SAMPLES, loc='upper right', frameon=False)


def plot_corr_array_normal_cerebellum_samples(marray_data):
    from matplotlib import pyplot as plt, rc
    plt.interactive(True)
    import seaborn as sns
    sns.set_style('white')
    rc('text', usetex=True)

    sample_names = [
        # adult cerebellum
        'NT-1197',
        'NCb-1',
        'NCb-2',
        'A911105',
        'A508112',
        'A508285',
        # fetal brain
        'NFb',
        'A602126',
        'A602127',
    ]
    # get subdata
    subdata = pd.DataFrame()
    for s in sample_names:
        t = marray_data.loc[:, [s, s + '-R']].mean(axis=1)
        subdata[s] = t

    # generate correlation matrix
    c = subdata.corr()
    # zero upper triangle
    c.values[np.triu_indices_from(c, 1)] = np.nan

    # plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    h = ax.matshow(c, cmap='Reds', vmin=0.8, vmax=1., origin='lower')
    ax.set_xticks(range(subdata.columns.size))
    ax.set_yticks(range(subdata.columns.size))
    ax.set_xticklabels(subdata.columns, rotation=30)
    ax.set_yticklabels(subdata.columns)
    fig.colorbar(h)


def plot_scatter_array_microarray_passages(marray_data, norm=True):
    """
    The microarray data contains three different passages for the xenotransplant (I, III, IV).
    Produce a scatterplot grid to show how well these correlate
    :param marray_data: Normalised intensities loaded with load_normed_microarray_data
    :param norm: If True, divide by max intensity before running
    """

    from matplotlib import pyplot as plt, rc
    rc('text', usetex=True)
    MARRAY_SAMPLES = [
        'Pt1299',
        'ICb1299-I',
        'ICb1299-III',
        'ICb1299-IV',
    ]
    # reduce to max expression between repeats
    series = []
    for m in MARRAY_SAMPLES:
        s = marray_data.loc[:, [m, m + '-R']].mean(axis=1)
        s[s.isnull()] = 0.
        series.append(s)

    marray_sub = pd.DataFrame(series, index=MARRAY_SAMPLES).transpose()
    if norm:
        marray_sub /= marray_sub.max()

    xmin = marray_sub.min().min()
    xmax = marray_sub.max().max()
    x = np.linspace(xmin, xmax, 2)
    if norm:
        axlim = [0, 1]
    else:
        axlim = [0, np.ceil(xmax / 10000.) * 10000.]

    fig, axs = plt.subplots(len(MARRAY_SAMPLES), len(MARRAY_SAMPLES), sharex=True, sharey=True, figsize=(10, 10))
    for i in range(len(MARRAY_SAMPLES)):
        for j in range(len(MARRAY_SAMPLES)):
            m1 = MARRAY_SAMPLES[i]
            m2 = MARRAY_SAMPLES[j]
            d1 = marray_sub.loc[:, m1]
            d2 = marray_sub.loc[:, m2]
            lr = stats.linregress(d1, d2)
            axs[i, j].plot(x, lr.intercept + lr.slope * x, 'r--')
            axs[i, j].scatter(d1, d2, edgecolor='none', alpha=0.4)
            axs[i, j].set_aspect('equal')
            axs[i, j].set_xlim(axlim)
            axs[i, j].set_ylim(axlim)
            axs[i, j].legend(('$R^2 = %.3f$' % lr.rvalue ** 2,), loc='upper left', frameon=False)

    for j in range(len(MARRAY_SAMPLES)):
        labels = axs[2, j].get_xticklabels()
        plt.setp(labels, rotation=30)
        axs[j, 0].set_ylabel(MARRAY_SAMPLES[j], fontsize=14)
        axs[2, j].set_xlabel(MARRAY_SAMPLES[j], fontsize=14)

    plt.subplots_adjust(left=0.15, right=0.99, bottom=0.125, top=0.99, wspace=0., hspace=0.)


def plot_mb_group_activity_healthy_vs_mb(marray_data, cat):
    """
    Grouped bar chart showing the absolute expression levels of marker genes in the healthy samples and later passage MB
    xenografts.
    :param marray_data:
    :param cat:
    :return:
    """
    from matplotlib import pyplot as plt
    import seaborn as sns
    plt.interactive(True)
    sns.set_style('white')

    marray_ann = add_gene_symbol_column(marray_data, cat)

    REF_GROUPS = (
        ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
        ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
        ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
        ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
    )

    CEREBELLUM_SAMPLE_NAMES = [
        'NT-1197',
        'NCb-1',
        'NCb-2',
        'A911105',
        'A508112',
        'A508285',
    ]
    CEREBELLUM_SAMPLE_NAMES += [t + '-R' for t in CEREBELLUM_SAMPLE_NAMES]

    XENOGRAFT_SAMPLE_NAMES = [
        'ICb1299-III',
        'ICb1299-IV',
    ]
    XENOGRAFT_SAMPLE_NAMES += [t + '-R' for t in XENOGRAFT_SAMPLE_NAMES]

    # add new grouping column
    marray_ann['mb_group'] = [np.nan] * len(marray_ann.index)

    for grp, arr in REF_GROUPS:
        marray_ann.loc[marray_ann.gene_symbol.isin(arr), 'mb_group'] = grp

    # compare cerebellum sample median with passage III and IV
    grp_data = marray_ann.groupby(['gene_symbol', 'mb_group'])
    xeno_data = collections.defaultdict(dict)
    heal_data = collections.defaultdict(dict)
    for (gene, grp), arr in grp_data:
        xeno_data[grp][gene] = np.nanmedian(arr.loc[:, XENOGRAFT_SAMPLE_NAMES].values)
        heal_data[grp][gene] = np.nanmedian(arr.loc[:, CEREBELLUM_SAMPLE_NAMES].values)

    fig, axs = plt.subplots(ncols=len(xeno_data), sharey=True, figsize=(10, 5))
    WIDTH = 0.8
    MAX_N_X = max([len(b) for a, b in heal_data.items()])
    dw = WIDTH / 2.  # each slot has 2 bars
    for i in range(len(xeno_data)):
        grp = REF_GROUPS[i][0]
        x = xeno_data[grp]
        h = heal_data[grp]
        ax = axs[i]
        x0 = np.arange(len(h))
        l0, y0 = zip(*h.items())
        ax.bar(x0, y0, dw, facecolor='k', edgecolor='none')
        x1 = x0 + dw
        y1 = [x[l] for l in l0]
        ax.bar(x1, y1, dw, facecolor='b', edgecolor='none')

        ax.set_xticks(x1)
        ax.set_xticklabels(l0, rotation=90)
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(0.5, -.21)
        ax.set_xlim([WIDTH - 1, MAX_N_X])

    axs[-1].legend(['Healthy cerebellum', 'MB'])
    ylim = list(axs[-1].get_ylim())
    ylim[0] = -10
    axs[-1].set_ylim(ylim)
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.95, wspace=0.08, hspace=0.)

    return fig  # useful for programmatic saving
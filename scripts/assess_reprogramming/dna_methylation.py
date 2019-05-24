from methylation import loader, dmr, process
from plotting import clustering, common, venn, pca as pca_plotting
from stats import transformations
import pandas as pd
import numpy as np
import copy
import os
from utils import output, setops
import references
import collections
from scipy.stats import zscore, spearmanr
from scipy.cluster import hierarchy as hc
import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.decomposition import pca


def init_pool_shared_obj(obj):
    global the_obj
    the_obj = obj


def spearmanr_compute_one(i, j):
    return spearmanr(the_obj.iloc[:, i], the_obj.iloc[:, j])[0]


def spearman_rank_corr(dat, n_job=None):
    """
    For the supplied data, generate a square matrix with the Spearman rank corr coeff for each pairwise comparison.
    Parallelise this as it can be quite slow.
    :param dat: pd.DataFrame, samples in columns, features in rows.
    :return: Symmetric pd.DataFrame, containing the corr coeff for each comparison
    """
    pool = None
    jobs = {}
    m, n = dat.shape
    if (m * n) > 50000:
        if n_job is None:
            n_job = mp.cpu_count()
        if n_job > 1:
            pool = mp.Pool(processes=n_job, initializer=init_pool_shared_obj, initargs=(dat,))

    pdist = pd.DataFrame(1., index=dat.columns, columns=dat.columns)

    for i in range(dat.shape[1]):
        for j in range(i + 1, dat.shape[1]):
            if pool is None:
                pdist.iloc[i, j] = spearmanr(dat.iloc[:, i], dat.iloc[:, j])[0]
                pdist.iloc[j, i] = pdist.iloc[i, j]
            else:
                jobs[(i, j)] = pool.apply_async(
                    spearmanr_compute_one,
                    args=(i, j)
                )

    if pool is not None:
        pool.close()
        pool.join()
        for i, j in jobs:
            pdist.iloc[i, j] = pdist.iloc[j, i] = jobs[(i, j)].get()

    return pdist


def construct_colour_array_legend_studies(meta):
    studies = {}
    cc = pd.DataFrame('gray', index=meta.index, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'] == 'Our data'), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'] != 'Our data'), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'
    cc.loc[(meta['type'] == 'iNSC') & (meta['batch'] == 'Our data'), 'Cell type'] = '#9e3900'  # chestnut
    cc.loc[meta['type'] == 'iNSC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[meta['type'] == 'NPC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[meta['type'] == 'NSC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[(meta['type'] == 'NSC') & (meta.index.str.contains('fetal')), 'Cell type'] = '#ffaf47'  # light orange

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        studies[batches[i]] = c

    return cc, studies


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None, units='beta'):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = loader.load_by_patient(pids, norm_method=norm_method)
    anno = loader.load_illumina_methylationepic_annotation()

    # reference data
    if ref_names is not None:
        ref_obj = loader.load_reference(ref_names, norm_method=norm_method)
        if ref_name_filter is not None:
            ref_obj.filter_by_sample_name(ref_name_filter, exact=True)
        obj = loader.loader.MultipleBatchLoader([obj, ref_obj])

    me_data = obj.data.dropna()
    if units == 'm':
        me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    # dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def compute_dmr_clusters(anno, dmr_params):
    clusters = []
    cid = 0

    for cc in anno.CHR.unique():
        coords = anno.loc[anno.CHR == cc, 'MAPINFO'].sort_values()
        this_clust = dmr.identify_cluster(coords, dmr_params['n_min'], dmr_params['d_max'])

        for cl in this_clust.values():
            clusters.append(
                dmr.ProbeCluster(cl, anno, cluster_id=cid, chr=cc)
            )
            cid += 1
    return dmr.DmrResults(clusters=clusters, anno=anno)


def combine_data_meta(data_arr, meta_arr, units='beta'):
    if len(data_arr) != len(meta_arr):
        raise ValueError("data_arr and meta_arr must have the same size")

    # include all probes again
    dat = pd.concat(
        data_arr,
        axis=1,
        join='inner'
    )
    meta = pd.concat(
        meta_arr,
        axis=0,
        join='outer',
        sort=True
    )
    if units.lower() == 'm':
        # convert to M values
        dat = process.m_from_beta(dat)

    # drop any infinite valued probes (should be very few)
    inft = (~np.isfinite(dat)).sum(axis=1) > 0

    if inft.any():
        dat = dat.loc[~inft]
        print "Dropped %d probes with infinite M values" % inft.sum()

    return meta, dat


def aggregate_samples(search_repl_arr, data, meta):
    data = data.copy()
    meta = meta.copy()
    for srch, repl in search_repl_arr:
        idx = data.columns.str.contains(srch)
        if idx.sum() == 0:
            print "No columns matched the search string %s;  skipping" % srch
            continue

        new_col = data.loc[:, idx].mean(axis=1)

        data = data.loc[:, ~idx]
        meta_col = meta.loc[data.columns[0]]
        meta = meta.loc[data.columns]

        data.insert(data.shape[1], repl, new_col, allow_duplicates=True)
        meta.loc[repl] = meta_col
    return meta, data



if __name__ == "__main__":
    outdir = output.unique_output_dir("assess_reprogramming_methylation")
    pids = ['019', '030', '031', '050', '054']

    # these are the only two norming methods available in all data sets
    # norm_method = 'raw'
    norm_method = 'bmiq'
    # norm_method = 'pbc'

    min_val = 0.75  # value above which a probe is definitely methylated
    n_above_min = 3

    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    # our data
    me_obj, anno = load_methylation(pids, norm_method=norm_method)
    our_data = me_obj.data
    our_meta = me_obj.meta

    # discard unneeded samples
    our_meta = our_meta.loc[our_meta.type.isin(['iPSC', 'FB'])]
    our_data = our_data.loc[:, our_meta.index]
    our_meta.loc[:, 'batch'] = 'Our data'
    our_meta.insert(1, 'array_type', 'EPIC')

    # ref data

    # Encode EPIC data
    encode_epic_obj = loader.load_reference(
        'ENCODE_EPIC',
        norm_method=norm_method,
        samples=['H7 hESC', 'GM23248', 'GM23338', 'IMR-90']
    )
    encode_epic_obj.meta.insert(1, 'array_type', 'EPIC')
    encode_epic_obj.batch_id = 'Encode EPIC'
    # encode_epic_obj = loader.encode_epic(norm_method=norm_method, samples=['H7 hESC', 'GM23248', 'GM23338', 'IMR-90'])

    # Encode 450K data
    encode_450k_obj = loader.load_reference(
        'ENCODE_450K',
        norm_method=norm_method,
        samples=['H1 hESC']
    )
    encode_450k_obj.meta.insert(1, 'array_type', '450K')

    # Zimmerlin et al. (450K, 3 x hESC samples)
    zimmerlin_obj = loader.load_reference('GSE65214', norm_method=norm_method)
    zimmerlin_obj.batch_id = 'Zimmerlin et al.'
    zimmerlin_obj.meta.insert(1, 'array_type', '450K')

    # E-MTAB-6194
    e6194_obj = loader.load_reference(
        'E-MTAB-6194',
        norm_method=norm_method,
    )
    discard_ix = e6194_obj.meta.cell_line.isin([
        'NA07057',
        'HCT116',
        'HEL46.11',
    ])
    e6194_obj.meta = e6194_obj.meta.loc[~discard_ix]
    e6194_obj.data = e6194_obj.data[e6194_obj.meta.index]
    e6194_obj.meta.insert(1, 'array_type', 'EPIC')

    refs = [
        encode_epic_obj,
        encode_450k_obj,
        zimmerlin_obj,
        e6194_obj
    ]

    for r in refs:
        r.meta.index = ["%s_%s" % (t, r.batch_id) for t in r.meta.index]
        r.data.columns = r.meta.index

    # refs = [
    #     ('Kim et al.', loader.gse38216(norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2', 'H9 NPC 1', 'H9 NPC 2'])),
    #     ('Zimmerlin et al.', loader.gse65214(norm_method=norm_method)),
    #     ('Encode EPIC', loader.encode_epic(norm_method=norm_method)),
    #     ('Encode 450k', loader.encode_450k(norm_method=norm_method)),
    # ]
    # for bid, r in refs:
    #     r.batch_id = bid
    #     r.meta.index = ["%s_%s" % (t, bid) for t in r.meta.index]
    #     r.data.columns = r.meta.index
    # ref_obj = loader.loader.MultipleBatchLoader([t[1] for t in refs])
    # ref_meta = ref_obj.meta.loc[ref_obj.meta.index.str.contains('ESC')]
    # ref_data = ref_data.loc[:, ref_meta.index]

    # combine all loaders
    # this will reduce the probe list to the intersection (i.e. 450K)
    ref_obj = loader.loader.MultipleBatchLoader(refs)

    ref_meta = ref_obj.meta
    ref_data = ref_obj.data.dropna()

    # HipSci data
    hip_epic_ldr = loader.hipsci(norm_method=norm_method, n_sample=12, array_type='epic')
    hip_epic_meta = hip_epic_ldr.meta
    hip_epic_data = hip_epic_ldr.data

    # hip_450k_meta, hip_450k_data = loader.hipsci(norm_method=norm_method, n_sample=30, array_type='450k')

    hip_ldr = loader.hipsci(norm_method=norm_method, n_sample=12, array_type='all')
    hip_meta = hip_ldr.meta
    hip_data = hip_ldr.data
    hip_meta.batch = ["HipSci (%s)" % t for t in hip_meta.array_type]

    # clustering genome-wide
    # iPSC, FB, ESC

    # mix of HipSci samples by array_type
    meta, dat = combine_data_meta(
        (our_data, ref_data, hip_data),
        (our_meta, ref_meta, hip_meta)
    )

    # plot distribution of beta values (our data only)
    xi = np.linspace(0, 1, 101)
    xc = xi[:-1] + 0.5 * (xi[1] - xi[0])
    color_scalars = np.linspace(0, 1, our_data.shape[1])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color_scalars = np.linspace(0, 1, our_data.shape[1])
    for i, (c, ser) in enumerate(our_data.iteritems()):
        yh, _ = np.histogram(ser, xi, normed=True)
        ax.plot(xc, yh, label=c, color=plt.cm.jet(color_scalars[i]))
    ax.set_xlabel('Beta')
    ax.set_ylabel('Density')
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "edf_our_data_%s.png" % norm_method), dpi=200)

    # plot distribution of beta values (all data)
    leg_marker = set()
    plt_colours = {
        'IPSC': {'color': 'black'},
        'HPS': {'color': '#70a6ff', 'alpha': 0.6},
        'ESC': {'color': '#f77e7e'},
        'FB': {'color': '#a3d838'},
    }
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color_scalars = np.linspace(0, 1, dat.shape[1])
    for i, (c, ser) in enumerate(dat.iteritems()):
        kwds = {}
        for patt, d in plt_colours.items():
            if patt in c:
                kwds = d
                if patt not in leg_marker:
                    kwds['label'] = patt  # TODO: make this more flexible
                    leg_marker.add(patt)
                else:
                    kwds['label'] = None
                break
        yh, _ = np.histogram(ser, xi, normed=True)
        ax.plot(xc, yh, **kwds)
    ax.set_xlabel('Beta')
    ax.set_ylabel('Density')
    fig.tight_layout()
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(outdir, "edf_all_data_%s.png" % norm_method), dpi=200)

    # filter data (include only probes meth. in at least n samples)
    # dat = dat.loc[(dat > min_val).sum(axis=1) >= n_above_min]

    # two PCA plots, different groupings
    p = pca.PCA()
    pca_dat = p.fit_transform(dat.transpose())

    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10.3, 5.5))

    # 1. cell type
    groups = meta.type.copy()
    groups[groups.index.str.contains('IPSC')] = 'iPSC (this study)'
    # ordered dict to ensure our samples are always on top
    pca_colour_map = collections.OrderedDict([
        ('FB', '#fff89e'),
        ('iPSC', '#96ffb2'),
        ('ESC', 'green'),
        ('iPSC (this study)', 'blue'),
    ])
    ax = axs[0]

    grp_id, grp_nm = groups.factorize()
    for i, nm in enumerate(pca_colour_map.keys()):
        this_ix = groups == nm
        ax.scatter(
            pca_dat[this_ix, 0],
            pca_dat[this_ix, 1],
            c=pca_colour_map[nm],
            edgecolor='k',
            zorder=i+1,
            s=40,
            label=nm
        )
    ax.set_xlabel("Principle component 1")
    ax.set_ylabel("Principle component 2")
    ax.legend(frameon=True, edgecolor='k', facecolor='w', framealpha=0.5)

    # 2. array type
    groups = meta.batch.copy()
    pca_colour_map = dict(zip(meta.batch.unique(), common.COLOUR_BREWERS[len(meta.batch.unique())]))
    ax = axs[1]

    grp_id, grp_nm = groups.factorize()
    for i, nm in enumerate(pca_colour_map.keys()):
        this_ix = groups == nm
        ax.scatter(
            pca_dat[this_ix, 0],
            pca_dat[this_ix, 1],
            c=pca_colour_map[nm],
            edgecolor='k',
            zorder=i+1,
            s=40,
            label=nm
        )
    ax.set_xlabel("Principle component 1")
    ax.legend(frameon=True, edgecolor='k', facecolor='w', framealpha=0.5)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pca_all_data_two_label_sets.png"), dpi=200)

    # same again, zoomed for more detail on ESC/iPSC cluster
    ax.axis([-40, 5, -8, 8])
    fig.savefig(os.path.join(outdir, "pca_all_data_two_label_sets_zoom.png"), dpi=200)

    # Clustering / dendrogram
    # It is helpful to aggregate over some replicates for clarity here
    # NB before doing so, check the PCA (or clustering) plot(s) to be sure the samples are similar
    aggr_arr = [
        ('_HEL141.1_', 'HEL141.1_E-MTAB-6194'),
        ('_HEL140.1_', 'HEL140.1_E-MTAB-6194'),
        ('_HEL139.', 'HEL139.1_E-MTAB-6194'),
        ('_HFF_', 'HFF_E-MTAB-6194'),
        ('_H9_p50_', 'H9_E-MTAB-6194'),
    ]
    meta, dat = aggregate_samples(aggr_arr, dat, meta)


    cc, st = construct_colour_array_legend_studies(meta)
    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': st,
    }

    # MAD for ranking probes by variability across samples
    mad = transformations.median_absolute_deviation(dat).sort_values(ascending=False)

    # dendrogram

    # Spearman rank correlation distance
    # pdist = spearman_rank_corr(dat)
    # dist = hc.distance.squareform(1 - pdist.values)
    # lnk = hc.linkage(dist)
    # dend = clustering.dendrogram_with_colours(
    #     dat,
    #     cc,
    #     linkage=lnk,
    #     vertical=True,
    #     legend_labels=leg_dict,
    #     fig_kws={'figsize': [14, 6]}
    # )

    # Pearson correlation distance
    dend = clustering.dendrogram_with_colours(dat, cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})

    # Pearson with a limited number of probes
    # dend = clustering.dendrogram_with_colours(dat.loc[mad.index[:5000]], cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})

    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_all_probes.png"), dpi=200)

    # similar, but clustermap (dendrogram + heatmap)
    gc = clustering.plot_clustermap(
        dat.loc[mad.index[:5000]],
        cmap='RdBu_r',
        col_linkage=dend['linkage'],
        col_colors=cc
    )
    clustering.add_legend(leg_dict, gc.ax_heatmap, loc='right')
    gc.gs.update(bottom=0.2, right=0.82)

    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb_all_probes.png"), dpi=200)


from methylation import loader, process, dmr
from scripts.hgic_final import consts
import hgic_consts
from stats import transformations
from plotting import clustering, common, scatter, adjuster
from utils import output

import multiprocessing as mp
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import collections
import os


def hc_plot_dendrogram(
        data,
        row_colours,
        mad = None,
        n_ftr=3000,
        metric='correlation',
        **kwargs
):
    """
    For each value in n_gene_arr, plot a dendrogram showing the result of hierarchical clustering of the data using
    that many genes (selected in descending MAD order)
    :param data: Cols are samples, rows are genes (or similar)
    :param row_colours: As passed to dendrogram routine
    :param n_gene_arr: The values to test
    :return:
    """
    if 'fig_kws' not in kwargs:
        kwargs['fig_kws'] = {'figsize': (5.5, 10)}
    if mad is None:
        mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
    the_dat = data.loc[mad.index[:n_ftr]]
    fig_dict = clustering.dendrogram_with_colours(
        the_dat,
        row_colours,
        vertical=False,
        metric=metric,
        **kwargs
    )
    return fig_dict


def plot_pca(
        dat,
        colour_subgroups,
        p=None,
        components=(0, 1),
        marker_subgroups=None,
        ax=None,
        colour_map=None,
        marker_map=None,
        **kwargs
):
    if p is None:
        p = PCA()
        pca_data = p.fit_transform(dat.transpose())
    else:
        pca_data = p.transform(dat.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    ax = scatter.scatter_with_colour_and_markers(
        pca_data[:, components],
        colour_subgroups=colour_subgroups,
        colour_map=colour_map,
        marker_subgroups=marker_subgroups,
        marker_map=marker_map,
        ax=ax,
        **kwargs
    )

    ax.set_xlabel("PCA component %s (%.1f%%)" % (components[0] + 1, variance_explained[components[0]]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (components[1] + 1, variance_explained[components[1]]))

    return p, ax


def compute_one_fstat(row, samples):
    return stats.f_oneway(*[row[t] for t in samples.values()])

if __name__ == "__main__":
    """
    Here we assess the similarity of GIC cultures (in vitro) at both early and late passages and those established
    after xenografting (ex vivo).
    We also have bulk tumour samples from human patients (FFPE) and mouse xenografts (FFPE PDX).

    This is all with DNA methylation data (EPIC array).

    Except for FFPE, this is for one line (only): 019. May consider including the others to demonstrate how different
    lines/bulk vary.
    """
    # clustering parameters
    clust_n_ftr = 20000
    n_probe_to_show = 2000
    clustering_metric = 'euclidean'

    outdir = output.unique_output_dir()

    norm_method = 'swan'
    pdx_bulk_samples = [
        'SM18_108A_GBM019Luc_PDX1',
        'SM18_119A_GBM019Luc_PDX2'
    ]
    gic_late_samples = [
        'GBM019Luc_P12',
        'GBM019Luc_P3_PDX1',
        'GBM019Luc_P2_PDX2',
    ]

    # load all relevant data
    our_gic_obj = loader.load_by_patient(consts.PIDS, include_control=False, samples=consts.S1_METHYL_SAMPLES_GIC, norm_method=norm_method)
    our_ffpe_obj = loader.load_by_patient(consts.PIDS, type='ffpe', include_control=False, norm_method=norm_method)
    pdx_bulk = loader.load_reference('2018-12-14', norm_method=norm_method, samples=pdx_bulk_samples)
    gic_late = loader.load_reference('2018-12-06', norm_method=norm_method, samples=gic_late_samples)

    # add patient ID to samples
    our_gic_obj.meta.insert(
        0,
        'patient_id',
        our_gic_obj.meta.index.str.extract(r'GBM(?P<patient_id>[0-9]{3}).*')['patient_id'].values
    )
    our_ffpe_obj.meta.insert(0, 'patient_id', [hgic_consts.NH_ID_TO_PATIENT_ID_MAP[t] for t in our_ffpe_obj.meta.index])

    pdx_bulk.meta.insert(0, 'patient_id', '019')
    gic_late.meta.insert(0, 'patient_id', '019')

    # add descriptor
    our_gic_obj.meta.insert(1, 'descriptor', 'In vitro GIC')
    our_ffpe_obj.meta.insert(1, 'descriptor', 'Bulk GBM')
    pdx_bulk.meta.insert(1, 'descriptor', 'Bulk PDX')
    gic_late_desc = pd.Series('In vitro GIC', index=gic_late.meta.index)
    gic_late_desc.loc[gic_late_desc.index.str.contains('PDX')] = 'Ex vivo GIC'
    gic_late.meta.insert(1, 'descriptor', gic_late_desc)

    # combine into a single object with compatible probes
    obj = loader.loader.MultipleBatchLoader(
        [our_gic_obj, our_ffpe_obj, pdx_bulk, gic_late]
    )

    # colours for dendrogram / PCA
    n_desc = len(obj.meta.descriptor.unique())
    descriptor_cmap = common.continuous_cmap(0, n_desc - 1, cmap='Greys')
    descriptor_colours = dict(
        zip(obj.meta.descriptor.unique(), [descriptor_cmap(i) for i in range(n_desc)])
    )
    pid_colours = consts.PATIENT_COLOURS

    # convert to M values
    mdat = process.m_from_beta(obj.data)
    mdat_019 = mdat.loc[:, obj.meta.patient_id == '019']

    # extract meta for convenience
    meta = obj.meta
    meta_019 = obj.meta.loc[mdat_019.columns]

    row_colours_all = pd.DataFrame('white', index=mdat.columns, columns=['Type', 'Patient'])
    row_colours_all.loc[:, 'Type'] = obj.meta.descriptor.apply(descriptor_colours.get)
    row_colours_all.loc[:, 'Patient'] = obj.meta.patient_id.apply(pid_colours.get)

    # hierarchical clustering by M value, all samples
    plt_dict = hc_plot_dendrogram(mdat, row_colours_all, n_ftr=clust_n_ftr, metric=clustering_metric)

    leg_entry = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.,
    }

    lkg = plt_dict['linkage']
    leg_dict = collections.OrderedDict()
    leg_dict['Type'] = collections.OrderedDict()
    for k in sorted(descriptor_colours):
        leg_dict['Type'][k] = dict(leg_entry)
        leg_dict['Type'][k].update({'facecolor': descriptor_colours[k]})
    leg_dict['Patient'] = collections.OrderedDict()
    for k in sorted(pid_colours):
        leg_dict['Patient'][k] = dict(leg_entry)
        leg_dict['Patient'][k].update({'facecolor': pid_colours[k]})

    common.add_custom_legend(
        plt_dict['col_colour_ax'],
        leg_dict,
        loc_outside=True,
        loc_outside_horiz='right',
        loc_outside_vert='top' if clustering_metric == 'correlation' else 'bottom',
        bbox_to_anchor=(7.5, 1. if clustering_metric == 'correlation' else 0.),
        frameon=False
    )
    gs = plt_dict['gridspec']
    gs.update(right=0.5, bottom=0.1)
    plt_dict['dendrogram_ax'].set_frame_on(False)
    plt_dict['fig'].savefig(
        os.path.join(
            outdir,
            "all_samples_dendrogram_%d_probes_%s.png" % (clust_n_ftr, clustering_metric)
        ), dpi=200
    )

    # Again, again! But with only 019 samples this time
    row_colours = row_colours_all.loc[mdat_019.columns, ['Type']]
    plt_dict = hc_plot_dendrogram(
        mdat_019,
        row_colours,
        n_ftr=clust_n_ftr,
        metric=clustering_metric,
        fig_kws={'figsize': (4, 6)}
    )
    plt_dict['dendrogram_ax'].set_frame_on(False)

    common.add_custom_legend(
        plt_dict['col_colour_ax'],
        leg_dict['Type'],
        loc_outside=True,
        loc_outside_horiz='right',
        frameon=False,
        loc_outside_vert='top',
        bbox_to_anchor=(13.5, 1.)
    )
    gs = plt_dict['gridspec']
    gs.update(right=0.36, bottom=0.1)

    plt_dict['fig'].savefig(os.path.join(outdir, "dendrogram_019_only.png"), dpi=200)

    # PCA, all samples
    cmap = collections.OrderedDict([
        (pid, pid_colours[pid]) for pid in consts.PIDS
    ])

    mmap = pd.Series(
        common.FILLED_MARKERS[n_desc],
        index=meta.descriptor.unique()
    )

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    p, ax = plot_pca(
        mdat,
        meta.patient_id,
        colour_map=cmap,
        marker_subgroups=obj.meta.descriptor,
        marker_map=mmap,
        ax=ax
    )

    leg = ax.get_legend()
    leg.set_frame_on(False)
    ax.figure.subplots_adjust(left=0.12, right=0.75, top=0.98)

    ax.figure.savefig(os.path.join(outdir, "all_samples_pca.png"), dpi=200)

    # annotate a couple of samples
    to_annot = {
        'GBM019_P3n6': 'GBM019 P3 & 6',
        'GBM019_P4': 'GBM019 P4',
        'GBM019Luc_P12': 'GBM019 P12 (Luc)'
    }

    this_coords = p.transform(
        mdat.loc[:, to_annot.keys()].transpose()
    )
    anno = []
    for k, v in to_annot.items():
        this_coords = p.transform(mdat.loc[:, [k]].transpose())
        anno.append(ax.annotate(v, this_coords[0, :2]))
    ax.set_xlim([-800, 1400])
    ax.figure.savefig(os.path.join(outdir, "all_samples_pca_annotated.png"), dpi=200)

    # annotate all remaining GIC to determine whetehr passages are consistent
    [t.remove() for t in anno]
    anno = []
    for k in meta.index[meta.type == 'GBM']:
        this_str = meta.loc[k, 'passage']
        if not pd.isnull(this_str):
            this_coords = p.transform(mdat.loc[:, [k]].transpose())
            anno.append(ax.annotate(obj.meta.loc[k, 'passage'], this_coords[0, :2], fontsize=14))
    plt.draw()  # must have rendered plot before adjusting
    adjuster.adjust_text_radial_plus_repulsion(anno)
    ax.figure.savefig(os.path.join(outdir, "all_samples_pca_annotated_with_passage_numbers.png"), dpi=200)

    # same PCA definition, but only the 019 samples, so we can expand the view
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    p, ax = plot_pca(
        mdat_019,
        meta_019.patient_id,
        p=p,
        colour_map=cmap,
        marker_subgroups=meta_019.descriptor,
        marker_map=mmap,
        ax=ax
    )

    leg = ax.get_legend()
    leg.set_frame_on(False)
    ax.figure.subplots_adjust(left=0.12, right=0.75, top=0.98)

    ax.figure.savefig(os.path.join(outdir, "all_samples_pca_019_only.png"), dpi=200)

    anno = []
    for k in meta_019.index:
        this_coords = p.transform(mdat_019.loc[:, [k]].transpose())
        anno.append(ax.annotate(k, this_coords[0, :2]))

    ax.figure.savefig(os.path.join(outdir, "all_samples_pca_019_only_annotated.png"), dpi=200)

    # Again, again! But this time we compute the PCs with only 019 lines, then transform the remaining lines
    # using the same PCs

    p = PCA()
    pca_data = p.fit_transform(mdat_019.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax = scatter.scatter_with_colour_and_markers(
        pca_data[:, [0, 1]],
        colour_subgroups=obj.meta.loc[mdat_019.columns, 'patient_id'],
        colour_map=cmap,
        marker_subgroups=obj.meta.loc[mdat_019.columns, 'descriptor'],
        marker_map=mmap,
        ax=ax,
    )

    ax.set_xlabel("PCA component %s (%.1f%%)" % (1, variance_explained[0]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (2, variance_explained[1]))
    leg = ax.get_legend()
    leg.set_frame_on(False)
    fig.subplots_adjust(left=0.12, right=0.75, top=0.98)

    fig.savefig(os.path.join(outdir, "pca_019_only.png"), dpi=200)

    for i, col in enumerate(mdat_019.columns):
        this_coords = pca_data[i, :2]
        ax.annotate(col, this_coords)

    fig.savefig(os.path.join(outdir, "pca_019_only_annotated.png"), dpi=200)

    # overlay all data points, projected into the same 019-only transformation

    all_pca_data = p.transform(mdat.transpose())

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax = scatter.scatter_with_colour_and_markers(
        all_pca_data[:, [0, 1]],
        colour_subgroups=obj.meta.patient_id,
        colour_map=cmap,
        marker_subgroups=obj.meta.descriptor,
        marker_map=mmap,
        ax=ax,
    )

    ax.set_xlabel("PCA component %s (%.1f%%)" % (1, variance_explained[0]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (2, variance_explained[1]))
    leg = ax.get_legend()
    leg.set_frame_on(False)
    fig.subplots_adjust(left=0.12, right=0.75, top=0.98)
    fig.savefig(os.path.join(outdir, "pca_based_on_019_only.png"), dpi=200)

    for k, v in to_annot.items():
        this_coords = all_pca_data[mdat.columns == k, :2]
        ax.annotate(v, this_coords[0, :2])

    ax.figure.savefig(os.path.join(outdir, "pca_based_on_019_only_annotated.png"), dpi=200)

    # Identify 'passage aging' probes: most variable between earlier and later passages

    # Usable samples: drop 019 3+6 as this includes passages either side of 4
    this_samples = meta.index[meta.descriptor == 'In vitro GIC'].drop('GBM019_P3n6')
    this_passages = meta.loc[this_samples, 'passage'].str.split(';').apply(
        lambda x: np.mean([float(t) for t in x])
    )
    # for each probe, run a paired t test (WSRT is much slower)
    # to do this, separate each patient line into early and late (in the same order)
    g = this_passages.groupby(meta.loc[this_samples, 'patient_id'])
    early = g.idxmin().sort_index()
    late = g.idxmax()[early.index]

    mdat_e = mdat[early.values]
    mdat_l = mdat[late.values]

    ttest_paired = stats.ttest_rel(mdat_e, mdat_l, axis=1)
    # take top N probes by pvalue
    n_probes = 1000
    the_probes = mdat_e.index[np.argsort(ttest_paired.pvalue)[:n_probes]]

    # how do these vary between GIC019 late and GIC019 ex vivo?
    ## TODO

    # DMR analysis
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()
    anno = loader.load_illumina_methylationepic_annotation()
    anno = anno.loc[mdat.index]
    ffpe_gic_dmrs = dmr.DmrResults(anno=anno)
    ffpe_gic_dmrs.identify_clusters(**dmr_params)

    # use only GIC and FFPE here
    this_mdat = mdat.loc[:, obj.meta.descriptor.isin(['In vitro GIC', 'Bulk GBM'])]
    this_samples = this_mdat.columns.groupby(obj.meta.loc[this_mdat.columns, 'descriptor'])
    samples = [
        this_samples['In vitro GIC'],
        this_samples['Bulk GBM'],
    ]
    ffpe_gic_dmrs.test_clusters(
        this_mdat,
        samples=samples,
        n_jobs=dmr_params['n_jobs'],
        min_median_change=dmr_params['delta_m_min'],
        method=dmr_params['dmr_test_method'],
        alpha=dmr_params['alpha'],
        **dmr_params['test_kwargs']
    )

    # most variable DMRs between FFPE and GIC
    pa = pd.Series(dict([(k, v['padj']) for k, v in ffpe_gic_dmrs.results_significant.items()]))
    pa = pa.loc[pa.sort_values(ascending=True).index]
    # mc = pd.Series(dict([(k, v['median_change']) for k, v in ffpe_gic_dmrs.results_significant.items()]))
    # mc = mc.loc[mc.abs().sort_values(ascending=False).index]

    # colour by type
    colours = obj.meta.descriptor.apply(descriptor_colours.get)

    for j in [0, 10, 20]:
        fig, axs = plt.subplots(nrows=10, sharex=True, figsize=(8.0, 8.5))
        for i in range(len(axs)):
            mdat.loc[ffpe_gic_dmrs.clusters[pa.index[j+i]].pids].median(axis=0).plot.bar(
                color=colours,
                edgecolor='k',
                linewidth=1.0,
                ax=axs[i],
                width=0.9
            )
        fig.subplots_adjust(bottom=0.28, top=0.99, left=0.07, right=0.99, hspace=0.25)
        fig.savefig(os.path.join(outdir, "median_m_by_ffpe_gic_dmrs_%d-%d.png" % (j + 1, j + 10)))

    # most variable probes between patients (use F statistic)
    # use the same data
    this_samples = this_mdat.columns.groupby(
        obj.meta.loc[this_mdat.columns, 'patient_id']
    )

    pool = mp.Pool()

    jobs = {}
    for p, row in this_mdat.iterrows():
        jobs[p] = pool.apply_async(compute_one_fstat, args=(row, this_samples))

    f_stats = pd.Series(dict([(p, j.get().statistic) for p, j in jobs.items()]))
    f_stats.sort_values(ascending=False, inplace=True)

    colours = obj.meta.patient_id.apply(pid_colours.get)

    for j in [0, 10, 20]:
        fig, axs = plt.subplots(nrows=10, sharex=True, figsize=(8.0, 8.5))
        for i in range(len(axs)):
            mdat.loc[f_stats.index[j + i]].plot.bar(
                color=colours,
                edgecolor='k',
                linewidth=1.0,
                ax=axs[i],
                width=0.9
            )
        fig.subplots_adjust(bottom=0.28, top=0.99, left=0.07, right=0.99, hspace=0.25)
        fig.savefig(os.path.join(outdir, "m_by_patient_fstat_probes_%d-%d.png" % (j + 1, j + 10)))
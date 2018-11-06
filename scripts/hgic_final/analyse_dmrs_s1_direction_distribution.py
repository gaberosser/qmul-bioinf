from plotting import bar, common
from methylation import loader, dmr, process
import pandas as pd
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR
logger = log.get_console_logger()


def pair_dmr(me_meta, me_data, dmr_results_obj, pids, **dmr_params):
    dmr_res = {}

    for pid in pids:
        obj = dmr_results_obj.copy()
        the_idx1 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'GBM')
        the_idx2 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'iNSC')
        the_idx = the_idx1 | the_idx2
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups)
        the_samples = [the_samples['GBM'], the_samples['iNSC']]

        obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_res[pid] = obj

    return dmr.DmrResultCollection(**dmr_res)


def bed_file_from_probes(anno, probes, out_fn, probe_half_len=61):
    this_regions = {}
    this_anno = anno.loc[probes]
    for p, row in this_anno.iterrows():
        strand = '+' if row.Strand == 'F' else '-'
        # we'll prepend the chrom name with 'chr' to ensure compatibility with hg19 (built in to Homer)
        this_regions[p] = ["chr%s" % row.CHR, row.MAPINFO - probe_half_len, row.MAPINFO + probe_half_len, strand]

    genomics.write_bed_file(this_regions, out_fn)


def count_dm_by_direction(dmr_res, pids=consts.PIDS):
    """
    For each patient, count the regions/probes based on direction of differential methylation
    :param dmr_res: Dictionary, keyed by patient ID. Each value is another dictionary, keyed by probe or cluster ID.
    Each value of those contains an object with the attribute `median_change` that is used to assess direction.
    :return: pd.DataFrame, cols correspond to patients, rows correspond to hypo/hyper
    """
    dmr_direction = {}

    for pid in pids:
        this_res = pd.DataFrame.from_dict(dmr_res[pid]).transpose()
        dmr_direction[pid] = {
            'Hyper': (this_res['median_change'] > 0).sum(),
            'Hypo': (this_res['median_change'] < 0).sum(),
        }
    return pd.DataFrame.from_dict(dmr_direction)[pids].loc[['Hypo', 'Hyper']]


def direction_of_dm_bar_plot(
        dmr_res,
        pids=consts.PIDS,
        as_pct=True,
        legend=True,
        legend_loc='outside',
        **kwargs
):
    """
    Generate bar chart showing the split between hypo- and hypermethylation in the supplied differential methylation
    results.
    :param dmr_res: Dictionary, keyed by patient ID. Each value is another dictionary, keyed by probe or cluster ID.
    Each value of those contains an object with the attribute `median_change` that is used to assess direction.
    :param pids: Useful to specify the ordering of the PIDs.
    :param as_pct: If True (default), values are plotted as percentages.
    :param kwargs: Passed to the plotting function `bar.stacked_bar_chart()`. For example, might include `ax=...` to
    specify axis.
    :return:
    """

    for_plot = count_dm_by_direction(dmr_res, pids=pids)

    if as_pct:
        for_plot = for_plot.divide(for_plot.sum(), axis=1) * 100.
    colours = pd.Series({'Hyper': '#FF381F', 'Hypo': '#89CD61'})

    with sns.axes_style('whitegrid'):
        fig, ax = bar.stacked_bar_chart(for_plot, colours=colours, width=0.8, legend=legend, **kwargs)
    if as_pct:
        ax.set_ylim([0, 100])
    ax.set_xlim(np.array([0, len(pids)]) - 0.5)
    if legend:
        if legend_loc == 'outside':
            # shrink main axis and put legend on RHS
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
            common.legend_outside_axes(ax)
        else:
            ax.legend(loc=legend_loc)

    return fig, ax


def cpg_island_status(dmr_res, anno, clusters, all_probes=None, pids=consts.PIDS):
    """

    :param dmr_res:
    :param anno:
    :param clusters:
    :param all_probes:
    :param pids:
    :return:
    """
    if all_probes is None:
        # better to supply these, as we can then restrict ourselves to the probes actually represented in the data
        all_probes = anno.index

    # CpG island coverage
    pid_sets = {
        'background': {},
        'dmr': {},
        'hypo': {},
        'hyper': {},
    }
    for pid in pids:
        for k, v in pid_sets.items():
            if k == 'background':
                v[pid] = set(all_probes)
            else:
                v[pid] = set()
        for k in dmr_res[pid].keys():
            this_probe_ids = clusters[k].pids
            pid_sets['dmr'][pid].update(this_probe_ids)
            mc = dmr_res[pid][k]['median_change']
            if mc < 0:
                pid_sets['hypo'][pid].update(this_probe_ids)
            else:
                pid_sets['hyper'][pid].update(this_probe_ids)

    k_open_sea = 'open_sea'
    cats = {
        'N_Shore': 'shore',
        'S_Shore': 'shore',
        'Island': 'island',
        'N_Shelf': 'shelf',
        'S_Shelf': 'shelf',
        k_open_sea: 'open_sea',
    }
    empty_counts = {'open_sea': 0, 'island': 0, 'shore': 0, 'shelf': 0}
    island_counts = {}

    for pid_typ, pid_set in pid_sets.items():
        for pid in pids:
            p = pid_set[pid]
            this_cats = anno.loc[p, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea)
            this_counts = this_cats.value_counts().to_dict()
            island_counts.setdefault(pid_typ, {}).setdefault(pid, dict(empty_counts))
            for k, v in cats.items():
                island_counts[pid_typ][pid][v] += this_counts.get(k, 0)

    # sanity check
    for pid in pids:
        if not (pd.Series(island_counts['hyper'][pid]) + pd.Series(island_counts['hypo'][pid]) == pd.Series(island_counts['dmr'][pid])).all():
            raise ValueError("PID %s failed check # hypo + # hyper = # dmr" % pid)

    # save this in a 'nice' format for sharing
    cols = sorted(set(cats.values()))
    to_export = pd.DataFrame(
        index=pd.MultiIndex.from_product([pids, ['background', 'dmr', 'hypo', 'hyper']], names=['patient ID', 'probe list']),
        columns=cols
    )

    for pid_typ, pid_set in island_counts.items():
        for pid in pids:
            to_export.loc[(pid, pid_typ)] = pd.Series(pid_set[pid])[cols]

    return to_export


def dmr_to_dmp(dmr_res, clusters, me_data):
    dmp_ids = {}
    dmp_res = {}
    for pid in dmr_res:
        dmp_ids[pid] = set()
        dmp_res[pid] = {}
        gbm_samples = me_data.columns[me_data.columns.str.contains('GBM%s' % pid)]
        nsc_samples = me_data.columns[me_data.columns.str.contains('DURA%s' % pid)]
        for cl_id in dmr_res[pid]:
            dmp_ids[pid].update(clusters[cl_id].pids)
        dmp_ids[pid] = sorted(dmp_ids[pid])
        gbm_median = me_data.loc[dmp_ids[pid], gbm_samples].median(axis=1)
        nsc_median = me_data.loc[dmp_ids[pid], nsc_samples].median(axis=1)
        dmp_res[pid] = dict(
            [(k, {'median_change': v}) for k, v in (gbm_median - nsc_median).items()]
        )
    return {
        'dmp_ids': dmp_ids,
        'dmp_res': dmp_res,
    }


def dmr_to_dmp_specific(dmr_res, clusters, me_data):
    pids = dmr_res.keys()
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids
    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res[pid] for pid in pids])
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    dmp_ids_specific = {}
    for pid in pids:
        dmp_ids_specific[pid] = set()
        for cl_id in vs[pid]:
            dmp_ids_specific[pid].update(clusters[cl_id].pids)
        dmp_ids_specific[pid] = sorted(dmp_ids_specific[pid])

    venn_set_dmp, venn_ct_dmp = setops.venn_from_arrays(*[dmp_ids_specific[pid] for pid in pids])
    vs_dmp = dict([
        (p, venn_set_dmp[q]) for p, q in zip(pids, pu_sets)
    ])

    dmp_res_specific = {}
    for pid in pids:
        dmp_res_specific[pid] = {}
        gbm_samples = me_data.columns[me_data.columns.str.contains('GBM%s' % pid)]
        nsc_samples = me_data.columns[me_data.columns.str.contains('DURA%s' % pid)]
        gbm_median = me_data.loc[vs_dmp[pid], gbm_samples].median(axis=1)
        nsc_median = me_data.loc[vs_dmp[pid], nsc_samples].median(axis=1)
        dmp_res_specific[pid] = dict(
            [(k, {'median_change': v}) for k, v in (gbm_median - nsc_median).items()]
        )

    return {
        'dmp_ids': dmp_ids_specific,
        'dmp_res': dmp_res_specific
    }



def dm_probe_direction_panel_plot(
        dmr_res,
        dmp_res,
        pids=consts.PIDS,
):
    colours = {
        'dmr': '#689bed',
        'hypo': '#89CD61',
        'hyper': '#FF381F',
    }

    dmr_pct = {}
    for pid in pids:
        n_hypo = 0
        n_hyper = 0
        for r in dmr_res[pid].values():
            if r['median_change'] > 0:
                n_hyper += 1
            else:
                n_hypo += 1
        n_tot = float(n_hyper + n_hypo)
        dmr_pct[pid] = {
            'hypo': n_hypo / n_tot * 100.,
            'hyper': n_hyper / n_tot * 100.,
        }
    dmr_pct = pd.DataFrame(dmr_pct)[pids]

    probe_values = dict([
        (pid, np.array([v['median_change'] for v in dmp_res[pid].values()])) for pid in pids
    ])
    probe_pct = {}
    for pid, arr in probe_values.items():
        n_hypo = (arr < 0.).sum()
        n_hyper = (arr > 0.).sum()
        n_tot = float(n_hypo + n_hyper)
        probe_pct[pid] = {
            'hypo': n_hypo / n_tot * 100.,
            'hyper': n_hyper / n_tot * 100.,
        }
    probe_pct = pd.DataFrame(probe_pct)[pids]

    gs = plt.GridSpec(nrows=4, ncols=len(pids), height_ratios=[1, 1, 2, 1])
    fig = plt.figure(figsize=(10, 5))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res,
        as_pct=False,
        ax=fig.add_subplot(gs[2, :]),
        legend=False
    )
    ax.set_ylabel("Number DMRs")
    ax.set_xticklabels([])

    for i, pid in enumerate(pids):

        ax_kde = fig.add_subplot(gs[0, i])
        sns.kdeplot(
            probe_values[pid],
            ax=ax_kde,
            color='k',
        )
        # replot with two different shadings
        xx, yy = ax_kde.get_lines()[0].get_data()
        ax_kde.fill_between(xx[xx < 0], yy[xx < 0], facecolor=colours['hypo'], edgecolor='none', alpha=0.6)
        ax_kde.fill_between(xx[xx > 0], yy[xx > 0], facecolor=colours['hyper'], edgecolor='none', alpha=0.6)
        ax_kde.axvline(0., color='k', linestyle=':')
        ax_kde.yaxis.set_visible(False)
        if i == 0:
            ax_kde.set_ylabel('Probe density')
        ax_kde.xaxis.set_visible(False)
        ax_kde.set_xlim([-8, 12])
        ax_kde.set_title(pid)

        ax_pie = fig.add_subplot(gs[1, i])
        ax_pie.pie(
            probe_pct[pid].values,
            colors=[colours[t] for t in probe_pct.index],
            radius=1.,
            wedgeprops={'edgecolor': 'w', 'width': 0.5},
            startangle=90,
        )
        ax_pie.set_aspect('equal')

        ax_pie_dmr = fig.add_subplot(gs[3, i])
        ax_pie_dmr.pie(
            dmr_pct[pid].values,
            colors=[colours[t] for t in dmr_pct.index],
            radius=1.,
            wedgeprops={'edgecolor': 'w', 'width': 0.5},
            startangle=90,
        )
        ax_pie_dmr.set_aspect('equal')

        # if i == 0:
        #     ax_pie.set_ylabel('Probe direction')
        #     ax_pie_dmr.set_ylabel('DMR direction')
        #     ax_kde.set_ylabel('Probe density')

    fig.tight_layout()
    gs.update(wspace=0.02, hspace=0.05)

    return {
        'fig': fig,
        'gs': gs,
        'ax_main': ax
    }


def plot_panel_cpg_status(
        df,
        probe_values,
        pids=consts.PIDS
):
    probe_pct = {}
    for pid, arr in probe_values.items():
        n_hypo = (arr < 0.).sum()
        n_hyper = (arr > 0.).sum()
        n_tot = float(n_hypo + n_hyper)
        probe_pct[pid] = {
            'hypo': n_hypo / n_tot * 100.,
            'hyper': n_hyper / n_tot * 100.,
        }
    probe_pct = pd.DataFrame(probe_pct)[pids]

    colours = {
        'dmr': '#689bed',
        'hypo': '#89CD61',
        'hyper': '#FF381F',
    }

    nrows = df.loc[pids[0]].shape[0]
    nbar = df.loc[pids[0]].shape[1]
    buff = 0.6

    gs = plt.GridSpec(
        ncols=len(pids),
        nrows=nrows + 1,
        height_ratios=[1, 1] + [2] * (nrows - 1)
    )
    fig = plt.figure(figsize=(11.5, 5.5))

    label_map = {
        'dmr': 'All DMRs',
        'hypo': 'Hypo DMRs',
        'hyper': 'Hyper DMRs',
    }

    for i, pid in enumerate(pids):
        ax_kde = fig.add_subplot(gs[0, i])
        ax_pie = fig.add_subplot(gs[1, i])

        X = df.loc[pid].transpose().loc[['island', 'shore', 'shelf', 'open_sea']]
        X.index = ['Island', 'Shore', 'Shelf', 'Open sea']
        this_colours = [colours[t] for t in X.columns[1:]]  # first entry is background: not needed here

        # kde plot in top row
        sns.kdeplot(
            probe_values[pid],
            ax=ax_kde,
            color='k',
            shade='gray'
        )
        # replot with two different shadings
        xx, yy = ax_kde.get_lines()[0].get_data()
        ax_kde.fill_between(xx[xx < 0], yy[xx < 0], facecolor=colours['hypo'], edgecolor='none', alpha=0.6)
        ax_kde.fill_between(xx[xx > 0], yy[xx > 0], facecolor=colours['hyper'], edgecolor='none', alpha=0.6)

        ax_kde.axvline(0., linestyle=':', color='k')
        ax_kde.set_xlim([-9, 12])

        # pie plot in next row
        ax_pie.pie(
            probe_pct[pid].values,
            colors=[colours[t] for t in probe_pct.index],
            radius=1.,
            wedgeprops={'edgecolor': 'w', 'width': 0.5},
            startangle=90,
        )
        ax_pie.set_aspect('equal')

        ax_kde.xaxis.set_visible(False)
        ax_kde.yaxis.set_ticks([])
        ax_kde.set_title(pid)

        for j, c in enumerate(X.columns[1:]):
            ax_main = fig.add_subplot(gs[j + 2, i])

            ax_main.bar(
                range(X.shape[0]),
                X.loc[:, c].values,
                color=this_colours[j],
                edgecolor='k'
            )
            ax_main.set_xlim([-buff, nbar - 1. + buff])
            ax_main.set_ylim([0, 60.])
            ax_main.set_xticks(range(X.shape[0]))
            ax_main.set_xticklabels(X.index, rotation=90)

            # add points showing background values
            ax_main.scatter(
                range(X.shape[0]),
                X['background'],
                marker='o',
                edgecolor='k',
                facecolor='none',
                linewidths=1.,
                s=20,
                zorder=10,
            )

            if c != X.columns[-1]:
                ax_main.xaxis.set_ticklabels([])

            if i == 0:
                ax_main.set_ylabel(label_map[c])
            else:
                ax_main.yaxis.set_ticklabels([])

    fig.tight_layout()
    gs.update(hspace=0.2, wspace=0.1)
    return {
        'fig': fig,
        'gs': gs
    }


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tsgd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        dmr_res_s1 = tsgd.paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    # extract full (all significant) results
    dmr_res_all = dmr_res_s1.results_significant

    # patient-specific DMRs
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids
    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_all[pid] for pid in pids])
    vs = dict([
                  (p, venn_set[q]) for p, q in zip(pids, pu_sets)
                  ])
    dmr_res_specific = dict([
                                (
                                    pid,
                                    dict([(t, dmr_res_all[pid][t]) for t in vs[pid]])
                                ) for pid in pids
                                ])

    # full DMPs
    dd = dmr_to_dmp(
        dmr_res_all,
        dmr_res_s1.clusters,
        me_data
    )
    dmp_ids_all = dd['dmp_ids']
    dmp_res_all = dd['dmp_res']

    # patient-specific DMPs
    # define these as the probes in patient-specific DMRs
    dd = dmr_to_dmp_specific(
        dmr_res_all,
        dmr_res_s1.clusters,
        me_data
    )

    dmp_ids_specific = dd['dmp_ids']
    dmp_res_specific = dd['dmp_res']

    # 1) direction of methylation change for all DMRs
    # a) clusters: absolute and relative values
    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_all,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMRs")
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_all,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "all_dmr_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "all_dmr_direction.tiff"), dpi=200)

    # b) probes: otherwise as above
    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmp_res_all,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMPs")
    fig, ax = direction_of_dm_bar_plot(
        dmp_res_all,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "all_dmp_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "all_dmp_direction.tiff"), dpi=200)

    # c) probes with KDE and pie chart
    plot_dict = dm_probe_direction_panel_plot(
        dmr_res_all,
        dmp_res_all
    )

    plot_dict['fig'].savefig(os.path.join(outdir, "all_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "all_dmp_direction_panel.tiff"), dpi=200)

    # 2) direction of methylation change for patient-specific DMRs
    # a) clusters: absolute and relative values
    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_specific,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number patient-specific DMRs")
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_specific,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_specific_dmr_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_dmr_direction.tiff"), dpi=200)

    # b) probes: obtained by converting patient-specific DMRs to probes
    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmp_res_specific,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMPs")
    fig, ax = direction_of_dm_bar_plot(
        dmp_res_specific,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_specific_dmp_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_dmp_direction.tiff"), dpi=200)

    # c) probes with KDE and pie chart
    plot_dict = dm_probe_direction_panel_plot(
        dmr_res_specific,
        dmp_res_specific
    )

    plot_dict['fig'].savefig(os.path.join(outdir, "specific_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "specific_dmp_direction_panel.tiff"), dpi=200)

    # plot showing differing CpG distributions
    # 1) All DMRs

    cpg_island_counts = cpg_island_status(
        dmr_res_all,
        anno,
        dmr_res_s1.clusters,
        all_probes=me_data.index,
    )
    dist_as_pct = cpg_island_counts.divide(cpg_island_counts.sum(axis=1), axis=0) * 100.
    probe_values_all = dict([
                                (pid, np.array([v['median_change'] for v in dmp_res_all[pid].values()])) for pid in pids
                                ])

    plot_dict = plot_panel_cpg_status(
        dist_as_pct,
        probe_values_all
    )
    plot_dict['fig'].savefig(os.path.join(outdir, "probe_cpg_island_status_all.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "probe_cpg_island_status_all.tiff"), dpi=200)

    # statistical test: which distributions are significantly different from the background?
    cpg_island_pval_all = {}
    cpg_island_chisq_all = {}
    bg = cpg_island_counts.loc[(pids[0], 'background')]  # can use any PID here, they are all the same
    for pid in pids:
        cpg_island_chisq_all[pid] = {}
        cpg_island_pval_all[pid] = {}
        for typ in ['dmr', 'hypo', 'hyper']:
            cpg_island_chisq_all[pid][typ], cpg_island_pval_all[pid][typ] = stats.chisquare(
                cpg_island_counts.loc[(pid, typ)],
                bg
            )

    # 2) Patient-specific DMRs
    probe_values_specific = dict([
                                     (pid, np.array([v['median_change'] for v in dmp_res_specific[pid].values()])) for pid in pids
                                     ])
    cpg_island_counts = cpg_island_status(
        dmr_res_specific,
        anno,
        dmr_res_s1.clusters,
        all_probes=me_data.index,
    )
    dist_as_pct = cpg_island_counts.divide(cpg_island_counts.sum(axis=1), axis=0) * 100.

    plot_dict = plot_panel_cpg_status(
        dist_as_pct,
        probe_values_specific
    )
    plot_dict['fig'].savefig(os.path.join(outdir, "probe_cpg_island_status_specific.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "probe_cpg_island_status_specific.tiff"), dpi=200)

    # statistical test: which distributions are significantly different from the background?
    cpg_island_pval_specific = {}
    cpg_island_chisq_specific = {}
    bg = cpg_island_counts.loc[(pids[0], 'background')]  # can use any PID here, they are all the same
    for pid in pids:
        cpg_island_chisq_specific[pid] = {}
        cpg_island_pval_specific[pid] = {}
        for typ in ['dmr', 'hypo', 'hyper']:
            cpg_island_chisq_specific[pid][typ], cpg_island_pval_specific[pid][typ] = stats.chisquare(
                cpg_island_counts.loc[(pid, typ)],
                bg
            )


    # export the probe regions for each patient to a BED file for motif enrichment
    # probes are CpG +/- 60 bp: https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_cpg_loci_identification.pdf

    island_status_map = {
        'N_Shore': 'shore',
        'S_Shore': 'shore',
        'Island': 'island',
        'N_Shelf': 'shelf',
        'S_Shelf': 'shelf',
        'open_sea': 'open_sea',
    }
    probe_half_len = 61

    # 1) All probes
    out_subdir = os.path.join(outdir, "bedfiles_all")
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    for pid in pids:
        this_mc = dmp_res_all[pid]
        this_ps_dict = {
            'all': this_mc.keys(),
            'hypo': [k for k, v in this_mc.items() if v['median_change'] < 0],
            'hyper': [k for k, v in this_mc.items() if v['median_change'] > 0],
        }
        for typ, ps in this_ps_dict.items():
            bed_fn = os.path.join(out_subdir, "%s_%s_oligo_mappings.bed" % (pid, typ))
            bed_file_from_probes(anno, ps, bed_fn, probe_half_len=probe_half_len)

            # repeat with even more granularity: include probe types
            island_statuses = anno.loc[ps, 'Relation_to_UCSC_CpG_Island'].fillna('open_sea')
            island_statuses = island_statuses.apply(island_status_map.get)
            for cs in ['island', 'shore', 'shelf', 'open_sea']:
                ps_subtype = island_statuses.index[island_statuses == cs]
                bed_fn = os.path.join(out_subdir, "%s_%s_%s_oligo_mappings.bed" % (pid, typ, cs))
                bed_file_from_probes(anno, ps, bed_fn)

    # 2) Patient-specific probes
    out_subdir = os.path.join(outdir, "bedfiles_patient_specific")
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    for pid in pids:
        this_mc = dmp_res_specific[pid]
        this_ps_dict = {
            'all': this_mc.keys(),
            'hypo': [k for k, v in this_mc.items() if v['median_change'] < 0],
            'hyper': [k for k, v in this_mc.items() if v['median_change'] > 0],
        }
        for typ, ps in this_ps_dict.items():
            bed_fn = os.path.join(out_subdir, "%s_%s_oligo_mappings.bed" % (pid, typ))
            bed_file_from_probes(anno, ps, bed_fn, probe_half_len=probe_half_len)

            # repeat with even more granularity: include probe types
            island_statuses = anno.loc[ps, 'Relation_to_UCSC_CpG_Island'].fillna('open_sea')
            island_statuses = island_statuses.apply(island_status_map.get)
            for cs in ['island', 'shore', 'shelf', 'open_sea']:
                ps_subtype = island_statuses.index[island_statuses == cs]
                bed_fn = os.path.join(out_subdir, "%s_%s_%s_oligo_mappings.bed" % (pid, typ, cs))
                bed_file_from_probes(anno, ps, bed_fn)

    # repeat the whole thing with S2 data

    norm_method_s2 = 'pbc'
    external_ref_names_dm = ['gse38216']
    external_ref_samples_dm = ['H9 NPC 1', 'H9 NPC 2']
    external_refs_dm = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]
    external_refs_dm_labels = [t[0] for t in external_refs_dm]

    # load methylation with external references
    me_obj_with_ref, anno = tsgd.load_methylation(
        pids,
        ref_names=external_ref_names_dm,
        ref_name_filter=external_ref_samples_dm,
        norm_method=norm_method_s2,
        patient_samples=consts.ALL_METHYL_SAMPLES  # technically contains H9 too, but these will be ignored
    )

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s2

    the_hash = tsgd.dmr_results_hash(me_obj_with_ref.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        dmr_res_s2 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s2 = dmr.compute_cross_dmr(
            me_obj_with_ref.data,
            me_obj_with_ref.meta,
            anno,
            pids,
            dmr_params,
            external_references=external_refs_dm
        )
        # Save DMR results to disk
        dmr_res_s2.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    ## Now move to the full (syngeneic + reference) dataset

    # recreate the S1 plots with the (reduced) S2 data

    dmr_res_s2_syn = {}
    for pid in pids:
        dmr_res_s2_syn[pid] = dmr_res_s2[pid][pid].results_significant

    dd = dmr_to_dmp(
        dmr_res_s2_syn,
        dmr_res_s2.clusters,
        me_obj_with_ref.data
    )
    dmp_res_s2_syn = dd['dmp_res']

    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_s2_syn[pid] for pid in pids])
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])
    dmr_res_s2_syn_specific = dict([
        (
            pid,
            dict([(t, dmr_res_s2_syn[pid][t]) for t in vs[pid]])
        ) for pid in pids
    ])

    dd = dmr_to_dmp_specific(
        dmr_res_s2_syn,
        dmr_res_s2.clusters,
        me_obj_with_ref.data
    )
    dmp_res_s2_syn_specific = dd['dmp_res']

    plot_dict = dm_probe_direction_panel_plot(
        dmr_res_s2_syn,
        dmp_res_s2_syn
    )
    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_all_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_all_dmp_direction_panel.tiff"), dpi=200)

    plot_dict = dm_probe_direction_panel_plot(
        dmr_res_s2_syn_specific,
        dmp_res_s2_syn_specific
    )

    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_specific_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_specific_dmp_direction_panel.tiff"), dpi=200)

    # check similar panel in reference comparisons
    # do this separately for the two references
    dmr_res_s2_ref = {}
    for r in external_refs_dm_labels:
        dmr_res_s2_ref[r] = {}
        for pid in pids:
            dmr_res_s2_ref[r][pid] = dmr_res_s2[pid][r].results_significant

    for r in external_refs_dm_labels:
        this_dmr_res = dmr_res_s2_ref[r]
        dd = dmr_to_dmp(
            this_dmr_res,
            dmr_res_s2.clusters,
            me_obj_with_ref.data
        )
        this_dmp_res = dd['dmp_res']

        venn_set, venn_ct = setops.venn_from_arrays(*[this_dmr_res[pid] for pid in pids])
        vs = dict([
            (p, venn_set[q]) for p, q in zip(pids, pu_sets)
        ])
        this_dmr_res_specific = dict([
            (
                pid,
                dict([(t, this_dmr_res[pid][t]) for t in vs[pid]])
            ) for pid in pids
        ])

        dd = dmr_to_dmp_specific(
            this_dmr_res,
            dmr_res_s2.clusters,
            me_obj_with_ref.data
        )
        this_dmp_res_specific = dd['dmp_res']

        plot_dict = dm_probe_direction_panel_plot(
            this_dmr_res,
            this_dmp_res
        )
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_all_dmp_direction_panel.png" % r), dpi=200)
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_all_dmp_direction_panel.tiff" % r), dpi=200)

        plot_dict = dm_probe_direction_panel_plot(
            this_dmr_res_specific,
            this_dmp_res_specific
        )
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_specific_dmp_direction_panel.png" % r), dpi=200)
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_specific_dmp_direction_panel.tiff" % r), dpi=200)

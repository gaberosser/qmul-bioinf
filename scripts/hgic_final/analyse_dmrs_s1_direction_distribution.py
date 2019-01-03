from plotting import bar, common
from methylation import loader, dmr, process
import pandas as pd
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import collections
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt, colors, gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from sklearn.neighbors import KernelDensity
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR
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
        stacked=True,
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
    colours = pd.Series(
        {
            'Hyper': consts.METHYLATION_DIRECTION_COLOURS['hyper'],
            'Hypo': consts.METHYLATION_DIRECTION_COLOURS['hypo']
        }
    )

    with sns.axes_style('whitegrid'):
        if stacked:
            fig, ax = bar.stacked_bar_chart(for_plot, colours=colours, width=0.8, legend=legend, **kwargs)
        else:
            ax = for_plot.transpose().plot.bar(
                color=[colours[t] for t in for_plot.index],
                legend=legend,
                **kwargs
            )
            fig = ax.figure
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


def direction_kde(dat, xi=None, separated=True):
    """
    Generate a KDE of the supplied DMR direction data.
    :param dat:
    :param xi: If supplied, this is a 1D array containing lookup points at which to evaluate the density.
    :param separated: If requested, generate a separate KDE for positive and negative results.
    :return: Array of KDE evaluated results
    """
    if xi is None:
        xi = np.linspace(
            np.floor(dat.min()),
            np.ceil(dat.max()),
            128
        )
    if separated:
        k_minus = stats.gaussian_kde(dat[dat <= 0])
        y_minus = k_minus.evaluate(xi[xi <= 0])
        k_plus = stats.gaussian_kde(dat[dat > 0])
        y_plus = k_plus.evaluate(xi[xi > 0])
        y = np.concatenate((y_minus, y_plus))
    else:
        k = stats.gaussian_kde(dat)
        y = k.evaluate(xi)
    return xi, y


def dm_probe_direction_panel_plot(
        dmr_res,
        dmp_res,
        me_data,
        cell_type,
        patient_id,
        pids=consts.PIDS,
        stacked_bar=False,
        log_residual=True,
        ref_name=None
):
    colours = {
        'dmr': '#689bed',
        'hypo': consts.METHYLATION_DIRECTION_COLOURS['hypo'],
        'hyper': consts.METHYLATION_DIRECTION_COLOURS['hyper'],
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

    gs = plt.GridSpec(
        nrows=6,
        ncols=len(pids),
        # height_ratios=[1, 1, 1, 2, 1]
    )

    fig = plt.figure(figsize=(10, 5))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res,
        as_pct=False,
        ax=fig.add_subplot(gs[4, :]),
        legend=False,
        stacked=stacked_bar
    )
    ax.set_ylabel("Number DMRs")
    ax.set_xticklabels([])

    for i, pid in enumerate(pids):
        ax_diff_hist = fig.add_subplot(gs[0, i])
        ax_diff_pie = fig.add_subplot(gs[1, i])

        if ref_name is None:
            this_dat = me_data.loc[:, patient_id == pid]
        else:
            this_dat = me_data.loc[:, (patient_id == pid) | (me_data.columns.str.contains(ref_name))]
        this_ct = cell_type.loc[this_dat.columns]

        quantify_dual_probe_beta_plot(
            this_dat,
            this_ct,
            hist_ax=ax_diff_hist,
            pie_ax=ax_diff_pie,
            edgecolor='none',
            nbin=40,
            log_scale=log_residual
        )
        # ax_diff_hist.set_xlim([-1, 1])
        ax_diff_hist.xaxis.set_visible(False)
        ax_diff_hist.yaxis.set_visible(False)

        ax_diff_hist.set_title(pid)

        ax_kde = fig.add_subplot(gs[2, i])
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

        ax_pie = fig.add_subplot(gs[3, i])
        ax_pie.pie(
            probe_pct[pid].values,
            colors=[colours[t] for t in probe_pct.index],
            radius=1.,
            wedgeprops={'edgecolor': 'w', 'width': 0.5},
            startangle=90,
        )
        ax_pie.set_aspect('equal')

        ax_pie_dmr = fig.add_subplot(gs[5, i])
        ax_pie_dmr.pie(
            dmr_pct[pid].values,
            colors=[colours[t] for t in dmr_pct.index],
            radius=1.,
            wedgeprops={'edgecolor': 'w', 'width': 0.5},
            startangle=90,
        )
        ax_pie_dmr.set_aspect('equal')

    fig.tight_layout()
    gs.update(wspace=0.02, hspace=0.05)

    return {
        'fig': fig,
        'gs': gs,
        'ax_main': ax,
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


def shaded_pie_chart_direction(
    dat,
    edges,
    ax=None,
    cmap_offset=0.2
):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    props = [(dat < 0).sum(), (dat > 0).sum()]

    centres = edges[:-1] + 0.5 * (edges[1] - edges[0])
    counts, bins= np.histogram(dat, edges)

    # inner pie: 2 categories
    ax.pie(
        props,
        colors=[consts.METHYLATION_DIRECTION_COLOURS[t] for t in ['hypo', 'hyper']],
        radius=.5,
        wedgeprops={'edgecolor': 'k'},
        startangle=90,
    )

    hypo_cmap = plt.get_cmap('Greens')
    hyper_cmap = plt.get_cmap('Reds')

    # only use a portion of the full cmap range, to avoid having entirely transparent (white) regions
    grad_colours = [hypo_cmap(-x + cmap_offset) if x < 0 else hyper_cmap(x + cmap_offset) for x in centres]

    ax.pie(
        counts,
        colors=grad_colours,
        radius=1.,
        wedgeprops={'edgecolor': 'none', 'width': 0.5},
        startangle=90
    )
    circ = plt.Circle([0, 0], radius=1., edgecolor='k', linewidth=1., facecolor='none')
    ax.add_patch(circ)

    ax.set_aspect('equal')

    return ax


def shaded_histogram_direction(dat, ax=None, nbin=100, **kwargs):
    # define some defaults, which will be applied where they have not been otherwise specified
    default_kwargs = {
        'linewidth': 1.,
        'edgecolor': 'k',
        'normed': True,
        'facecolor': 'none'
    }

    # shaded colours to indicate extent of beta difference
    cmap_offset = 0.2
    hypo_cmap = plt.get_cmap('Greens')
    hyper_cmap = plt.get_cmap('Reds')

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    for k, v in default_kwargs.items():
        if k not in kwargs:
            kwargs[k] = v

    edges = np.linspace(-1, 1, num=nbin + (1 - nbin % 2))  # need to ensure an EVEN number of edges (so zero is an edge)
    centres = edges[:-1] + 0.5 * (edges[1] - edges[0])
    counts, bins, patches = ax.hist(dat, edges, **kwargs)

    for ptch, bd in zip(patches, centres):
        if bd < 0:
            ptch.set_facecolor(hypo_cmap(-bd + cmap_offset))
        else:
            ptch.set_facecolor(hyper_cmap(bd + cmap_offset))

    return {
        'ax': ax,
        'edges': edges,
        'centres': centres,
        'counts': counts
    }


def quantify_dual_probe_beta_plot(m_dat, cell_type, hist_ax=None, pie_ax=None, nbin=100, log_scale=False, **kwargs):
    """

    :param m_dat:
    :param cell_type:
    :param hist_ax:
    :param pie_ax: Axes to use for the shaded pie chart. If not supplied, add an inset axis (oooooh!).
    :param nbin:
    :param kwargs:
    :return:
    """
    ix, ct = pd.Index(cell_type).factorize()
    if len(ct) != 2:
        raise ValueError("Expecting two cell types, found %d." % len(ct))

    b_dat = dmr.beta_from_m(m_dat)
    b_dat = b_dat.groupby(by=cell_type, axis=1).mean()
    b_diff = b_dat[ct[0]] - b_dat[ct[1]]

    # res = shaded_histogram_direction(b_diff, ax=hist_ax, nbin=nbin, **kwargs)
    res = beta_difference_trace(m_dat, cell_type, ax=hist_ax, nbin=nbin, log_scale=log_scale)

    hist_ax = res['ax']

    # need to ensure an EVEN number of edges (so zero is an edge)
    edges = np.linspace(-1, 1, num=nbin + (1 - nbin % 2))

    # inset pie chart
    if pie_ax is None:
        pie_ax = inset_axes(hist_ax, width="50%", height="70%", loc=2)

    shaded_pie_chart_direction(b_diff, edges, ax=pie_ax)

    return hist_ax, pie_ax


def beta_difference_trace(
    m_dat,
    cell_type,
    nbin=50,
    ax=None,
    cmap_offset=0.2,
    log_scale=False
):
    ix, ct = pd.Index(cell_type).factorize()
    if len(ct) != 2:
        raise ValueError("Expecting two cell types, found %d." % len(ct))

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    edges = np.linspace(0, 1, nbin + 1)
    de = edges[1] - edges[0]
    centres = edges[:-1] + 0.5 * de

    b_dat = dmr.beta_from_m(m_dat)
    b_dat = b_dat.groupby(by=cell_type, axis=1).mean()
    b_diff = b_dat[ct[0]] - b_dat[ct[1]]
    b_diff_abs = b_diff.abs()

    residuals = []
    for i in range(len(edges) - 1):
        ix = (edges[i] <= b_diff_abs) & (b_diff_abs < edges[i + 1])
        residuals.append(np.sign(b_diff[ix]).sum())

    hypo_cmap = plt.get_cmap('Greens')
    hyper_cmap = plt.get_cmap('Reds')

    cs = [
        hyper_cmap(ctr + cmap_offset) if t > 0 else hypo_cmap(ctr + cmap_offset)
        # consts.METHYLATION_DIRECTION_COLOURS['hyper'] if t > 0
        # else consts.METHYLATION_DIRECTION_COLOURS['hypo']
        for t, ctr in zip(residuals, centres)
    ]

    if log_scale:
        # split into positive and negative components
        residuals = np.array(residuals)

        pos_res = np.log10(residuals[residuals > 0])
        pos_c = np.array(cs)[residuals > 0]
        pos_e = edges[:-1][residuals > 0]

        neg_res = -np.log10(-residuals[residuals < 0])
        neg_c = np.array(cs)[residuals < 0]
        neg_e = edges[:-1][residuals < 0]

        ax.bar(pos_e, pos_res, color=pos_c, width=de, align='edge')
        ax.bar(neg_e, neg_res, color=neg_c, width=de, align='edge')
    else:
        ax.bar(edges[:-1], residuals, color=cs, width=de, align='edge')
    ax.axhline(0, color='k', linestyle='--', linewidth=1.)
    ax.set_xlim([0, 1])

    return {
        'ax': ax,
        'edges': edges,
        'residuals': residuals
    }


def get_binned_dmr_locations(
    dmr_res,
    clusters,
    chrom_lengths,
    window_size=20000,
    split_by_direction=False,
    coord_summary_method='first'
):
    """

    :param dmr_res:
    :param clusters:
    :param chrom_lengths:
    :param window_size:
    :param split_by_direction:
    :param coord_summary_method: Method used to reduce the list of CpG coordinates to a single one.
    Default is 'first', meaning take the 5'-most coordinate (first in the list). Other options: 'last', 'median', 'mean'
    :return:
    """
    if split_by_direction:
        dmr_loci_hypo = {}
        dmr_loci_binned_hypo = {}
        dmr_loci_hyper = {}
        dmr_loci_binned_hyper = {}
    else:
        dmr_loci = {}
        dmr_loci_binned = {}

    for pid in dmr_res:
        if split_by_direction:
            dmr_loci_binned_hypo[pid] = {}
            this_loci_hypo = collections.defaultdict(list)
            dmr_loci_binned_hyper[pid] = {}
            this_loci_hyper = collections.defaultdict(list)
        else:
            dmr_loci_binned[pid] = {}
            this_loci = collections.defaultdict(list)
        # this_loci = dict([(chrom, []) for chrom in chroms])
        for cluster_id, cl in dmr_res[pid].items():
            # get the chrom and locus
            pc = clusters[cluster_id]
            # use the requested method to get a representative coordinate from the list
            if coord_summary_method == 'first':
                the_coord = pc.coord_list[0]
            elif coord_summary_method == 'last':
                the_coord = pc.coord_list[-1]
            elif coord_summary_method == 'median':
                the_coord = np.median(pc.coord_list)
            elif coord_summary_method == 'mean':
                the_coord = np.mean(pc.coord_list)
            else:
                raise ValueError("Unsupported coordinate summary method '%s'." % coord_summary_method)

            if split_by_direction:
                if cl['median_change'] > 0:
                    this_loci_hyper[pc.chr].append(the_coord)
                else:
                    this_loci_hypo[pc.chr].append(the_coord)
            else:
                this_loci[pc.chr].append(the_coord)
        if split_by_direction:
            dmr_loci_hyper[pid] = this_loci_hyper
            dmr_loci_hypo[pid] = this_loci_hypo
        else:
            dmr_loci[pid] = this_loci
        # run the histogram process on each chrom
        if split_by_direction:
            for chrom, arr in this_loci_hyper.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned_hyper[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
            for chrom, arr in this_loci_hypo.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned_hypo[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
        else:
            for chrom, arr in this_loci.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
    if split_by_direction:
        return {
            'dmr_loci_hyper': dmr_loci_hyper,
            'dmr_loci_hypo': dmr_loci_hypo,
            'dmr_binned_hyper': dmr_loci_binned_hyper,
            'dmr_binned_hypo': dmr_loci_binned_hypo,
        }
    else:
        return {
            'dmr_loci': dmr_loci,
            'dmr_binned': dmr_loci_binned
        }


from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def dmr_location_plot(
        dmr_loci,
        chrom_length,
        cg_density=None,
        unmapped_density=None,
        window_size=20000,
        unmapped_threshold_pct=10,
        max_n_per_row=15,
        cg_fmax=0.95,
        cg_fmin=0.,
        plot_border=True
):
    chroms = chrom_length.keys()
    dmr_loci_hyper, dmr_loci_hypo = dmr_loci

    xmax = max(chrom_length.values())

    if cg_density is not None:
        # get the mean CG density so that we can centre the colour scale there
        all_cg_densities = []
        for chrom in chroms:
            this_cg = cg_density[chrom]
            this_cg_pct = this_cg / float(window_size) * 100.
            all_cg_densities.extend(this_cg_pct.values)
        all_cg_densities = sorted(all_cg_densities)
        cg_mean = np.mean(all_cg_densities)
        cg_vmin = all_cg_densities[int(cg_fmin * len(all_cg_densities))]
        cg_vmax = all_cg_densities[int(cg_fmax * len(all_cg_densities))]
        if cg_vmin > cg_mean:
            raise ValueError("The vmin value calculated for the background density is greater than the mean.")
        if cg_vmax < cg_mean:
            raise ValueError("The vmax value calculated for the background density is less than the mean.")

    fig = plt.figure(figsize=(10, 8))
    ncols = int(np.ceil(len(chrom_length) / float(max_n_per_row)))
    nrows = len(chrom_length) / ncols + len(chrom_length) % ncols
    ncols = 2
    gs_main = plt.GridSpec(
        nrows=nrows,
        ncols=ncols,
        left=0.01,
        right=.99,
        bottom=0.01,
        top=.99,
        hspace=0.03,
        wspace=0.03,
    )

    # we want to traverse by row then column, so create an array of axes beforehand
    # order='C' would give us traversal by column then row
    main_ax_arr = np.array([gs_main[i] for i in range(len(chroms))]).reshape((nrows, ncols)).flatten(order='F')
    ax_gc_dict = {}
    ax_hypo_dict = {}
    ax_hyper_dict = {}

    for i, chrom in enumerate(chroms):
        gs = gridspec.GridSpecFromSubplotSpec(
            3,
            1,
            main_ax_arr[i],
            height_ratios=[5, 1, 5],
            hspace=0.
        )

        ax_gc = fig.add_subplot(gs[1])
        ax_hyper = fig.add_subplot(gs[0], sharex=ax_gc)
        ax_hypo = fig.add_subplot(gs[2], sharex=ax_gc)

        ax_gc_dict[chrom] = ax_gc
        ax_hyper_dict[chrom] = ax_hyper
        ax_hypo_dict[chrom] = ax_hypo

        if cg_density is not None:
            this_cg = cg_density[chrom]
            this_cg_pct = this_cg / float(window_size) * 100.

            xx = np.array([this_cg.index.tolist() + [chrom_length[chrom]]] * 2)
            yy = np.zeros_like(xx); yy[1] = 1.
            cc = np.ma.masked_less([this_cg_pct.values], cg_vmin)
            # cc = np.array([this_cg_pct.values])

            norm = MidpointNormalize(midpoint=cg_mean, vmin=cg_vmin, vmax=cg_vmax)
            ax_gc.pcolor(xx, yy, cc, cmap='PuOr', norm=norm)
            ax_gc.set_xlim([0, xmax])

        if unmapped_density is not None:

            this_unmapped = unmapped_density[chrom]
            this_unmapped_pct = this_unmapped / float(window_size) * 100.

            uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmapped_threshold_pct)
            # since we don't care about the extent of unmapping, replace all values with a single one
            # uu[~uu.mask] = 0.3
            uu[~uu.mask] = 0.8

            ax_gc.pcolor(xx, yy, uu, cmap='Greys', vmax=1., vmin=0.)
            ax_gc.set_xlim([0, xmax])

        if plot_border:
            # draw a border around the extent of the chromosome
            border = plt.Rectangle(
                [0, 0.01],
                chrom_length[chrom],
                .98,
                edgecolor='k',
                facecolor='none',
                linewidth=1.,
                zorder=100.
            )
            ax_gc.add_patch(border)

        this_hypo = dmr_loci_hypo[chrom]
        this_hyper = dmr_loci_hyper[chrom]

        # KDE estimation gives us a nice representation of the DMR location distribution
        # NB this library expects a 2D array and complains otherwise!
        k_hypo = KernelDensity(bandwidth=window_size, kernel='gaussian')
        k_hypo.fit(np.array(this_hypo)[:, None])  # this increases the dim of the 1D array
        logd_hypo = k_hypo.score_samples(xx[0, None].transpose())

        # blank out baseline (for plotting purposes)
        logd_hypo[logd_hypo < -100] = -np.inf

        d_hypo = np.ma.masked_equal(np.exp(logd_hypo), 0.)
        hypo_max = d_hypo.max()

        k_hyper = KernelDensity(bandwidth=window_size, kernel='gaussian')
        k_hyper.fit(np.array(this_hyper)[:, None])  # this increases the dim of the 1D array
        logd_hyper = k_hyper.score_samples(xx[0, None].transpose())

        # blank out baseline (for plotting purposes)
        logd_hyper[logd_hyper < -100] = -np.inf

        d_hyper = np.ma.masked_equal(np.exp(logd_hyper), 0.)
        hyper_max = d_hyper.max()

        # plot the KDEs
        ax_hyper.fill_between(xx[0], d_hyper, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hyper'])
        ax_hyper.set_ylim([0, hyper_max * 1.02])

        # hypo: needs to be plotted upside down
        ax_hypo.invert_yaxis()
        ax_hypo.fill_between(xx[0], d_hypo, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hypo'])
        ax_hypo.set_ylim([hypo_max * 1.02, 0.])

        # plotting tick marks is a nice idea, but in practice gets messy - may work for smaller datasets?
        # ax_hyper.plot(xx[0], np.full_like(xx[0], -0.1 * d_hyper.max()), '|k', markeredgewidth=1)

        for ax in [ax_gc, ax_hypo, ax_hyper]:
            ax.set_facecolor('w')
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

    gs_main.update(right=1.2)

    return {
        'gs_main': gs_main,
        'main_ax_arr': main_ax_arr,
        'ax_gc_dict': ax_gc_dict,
        'ax_hyper_dict': ax_hyper_dict,
        'ax_hypo_dict': ax_hypo_dict,
    }


def fit_kde_dmr_location(locs, xi, bandwidth, normed=True):
    # the KDE library expects a 2D array, even if data are not multidimensional
    locs = np.array(locs)
    if len(locs.shape) == 1:
        locs = locs[:, None]
    elif len(locs.shape) == 2:
        if locs.shape[0] == 1 and locs.shape[1] > 1:
            locs = locs.transpose()
        elif locs.shape[1] == 1 and locs.shape[0] > 1:
            pass
        else:
            raise ValueError("locations array must be 1D or (quasi-)2D")
    else:
        raise ValueError("locations array must be 1D or (quasi-)2D")

    xi = np.array(xi)
    if len(xi.shape) == 1:
        xi = xi[:, None]
    elif len(xi.shape) == 2:
        if xi.shape[0] == 1 and xi.shape[1] > 1:
            xi = xi.transpose()
        elif xi.shape[1] == 1 and xi.shape[0] > 1:
            pass
        else:
            raise ValueError("xi array must be 1D or (quasi-)2D")
    else:
        raise ValueError("xi array must be 1D or (quasi-)2D")

    k = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
    k.fit(locs)  # this increases the dim of the 1D array
    logd_score = k.score_samples(xi)
    # blank out baseline (for plotting purposes)
    logd_score[logd_score < -100] = -np.inf
    score = np.ma.masked_equal(np.exp(logd_score), 0.)
    if not normed:
        # convert to unnormalised density
        score *= locs.size
    return score


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    subgroups = consts.SUBGROUPS

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    # set this to True if output bed files are required (this is quite slow due to the large number of combinations)
    write_bed_files = False

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

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
    ss_sets = {}
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets[grp] = k
    ss_part_sets = tsgd.partial_subgroup_specific(pids, subgroup_ind)

    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_all[pid] for pid in pids])
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    # patient specific DMRs
    dmr_res_specific = dict([
        (
            pid,
            dict([(t, dmr_res_all[pid][t]) for t in vs[pid]])
        ) for pid in pids
    ])

    # (full) subgroup specific DMRs
    dmr_res_subgroup_specific_full = dict([(k, venn_set[v]) for k, v in ss_sets.items()])



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

    # Example plots, useful for presentation slides
    cell_type_colours = {
        'GBM': 'r',
        'iNSC': 'k'
    }
    pid = '018'
    patient_id = me_meta.patient_id
    cell_type = me_meta.type

    fig = plt.figure()
    ax = fig.add_subplot(111)
    b_dat = process.beta_from_m(me_data.loc[:, patient_id == pid])
    ix, ct = cell_type.loc[b_dat.columns].factorize()
    for ii, jj in enumerate(ix):
        c = ct[jj]
        col = b_dat.columns[ii]
        sns.kdeplot(b_dat[col], color=cell_type_colours[c], ax=ax)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel('Density')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_%s_beta_kde.png" % pid), dpi=200)

    b_diff = b_dat.groupby(me_meta.type, axis=1).mean()
    b_diff = b_diff[ct[0]] - b_diff[ct[1]]
    res = shaded_histogram_direction(b_diff)
    ax = res['ax']
    fig = ax.figure
    ax.set_xlabel(r'$\Delta\beta$')
    ax.set_ylabel('Frequency')
    fig.tight_layout()

    iax = inset_axes(ax, width="50%", height="70%", loc=2)

    shaded_pie_chart_direction(b_diff, res['edges'], ax=iax)
    fig.savefig(os.path.join(outdir, "patient_%s_delta_beta_quantification.png" % pid), dpi=200)

    res = beta_difference_trace(
        me_data.loc[:, patient_id == pid],
        cell_type.loc[patient_id == pid],
        log_scale=False
    )
    ax = res['ax']
    ax.set_ylabel(r'Residual $\Delta\beta$')
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_%s_residual_delta_beta_quantification.png" % pid), dpi=200)

    res = beta_difference_trace(
        me_data.loc[:, patient_id == pid],
        cell_type.loc[patient_id == pid],
        log_scale=True
    )
    ax = res['ax']
    ax.set_ylabel(r'$\log_{10}$ residual $\Delta\beta$')
    ax.set_yticks([])
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_%s_log10_residual_delta_beta_quantification.png" % pid), dpi=200)



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
        dmp_res_all,
        me_data,
        me_meta.type,
        me_meta.patient_id,
        log_residual=True
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
        dmp_res_specific,
        me_data,
        me_meta.type,
        me_meta.patient_id
    )

    plot_dict['fig'].savefig(os.path.join(outdir, "specific_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "specific_dmp_direction_panel.tiff"), dpi=200)

    # investigate (genomic) distribution of DMRs

    # extract cluster attributes
    dmr_s1_clusters = dmr_res_s1[pids[0]].clusters

    # Get the distribution of CpGs in bins across the entire genome
    window_size = int(2e5)
    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    cg_density, chrom_length = genomics.cg_content_windowed(fa_fn, features=chroms, window_size=window_size)
    # reorder
    chrom_length = collections.OrderedDict([(k, chrom_length[k]) for k in chroms])

    # find unmapped regions so we can mask them on the plot
    unmapped_density, _ = genomics.cg_content_windowed(fa_fn, features=chroms, window_size=window_size, motif='N')

    # get all DMRs (whether significant or not) to use as a null hypothesis
    # the maximum DMR size is 7203 bp (mean ~ 1000), i.e. much smaller than the window size
    # TODO: think this justifies binning by the number of DMRs within a window (?)
    coord_summary_method = 'first'
    tmp = get_binned_dmr_locations(
        dmr_res_s1.results,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=False,
        coord_summary_method=coord_summary_method
    )
    dmr_loci_all = tmp['dmr_loci'][pids[0]]
    dmr_all_binned = tmp['dmr_binned'][pids[0]]

    #1 ) All DMRs
    tmp = get_binned_dmr_locations(
        dmr_res_all,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )
    dmr_loci_hyper = tmp['dmr_loci_hyper']
    dmr_loci_hypo = tmp['dmr_loci_hypo']

    tmp = get_binned_dmr_locations(
        dmr_res_specific,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )

    specific_dmr_loci_hyper = tmp['dmr_loci_hyper']
    specific_dmr_loci_hypo = tmp['dmr_loci_hypo']

    pid = pids[0]
    tmp = dmr_location_plot(
        [dmr_loci_hyper[pid], dmr_loci_hypo[pid]],
        chrom_length,
        # cg_density=cg_density,
        cg_density=dmr_all_binned,
        # unmapped_density=unmapped_density,
        unmapped_density=None,
        window_size=window_size
    )

    unmap_threshold_pct = 10. # unmapped % above this value will be masked
    xmax = max(chrom_length.values())

    # TODO: ongoing work on a polar plot

    all_cg_densities = []
    for chrom in chroms:
        this_cg = cg_density[chrom]
        this_cg_pct = this_cg / float(window_size) * 100.
        all_cg_densities.extend(this_cg_pct.values)
    cg_pct_mean = np.mean(all_cg_densities)

    gap_radians = 2 * np.pi / 200.
    sum_chrom_length = float(sum(chrom_length.values()))
    radians_per_bp = (2 * np.pi - gap_radians * len(chroms)) / sum_chrom_length
    inner_r = 2.
    outer_r = 2.1
    kde_width = 0.5

    # generate all KDEs first, so we can rescale them correctly at plot time
    hypo_kdes = {}
    hyper_kdes = {}

    hypo_kdes_specific = {}
    hyper_kdes_specific = {}



    for chrom in chroms:
        this_cg = cg_density[chrom]
        xx = np.array(this_cg.index.tolist() + [chrom_length[chrom]])
        dummy_res = np.ma.masked_all(xx.shape)
        this_hypo = dmr_loci_hypo[pid][chrom]
        this_hyper = dmr_loci_hyper[pid][chrom]
        this_hypo_specific = specific_dmr_loci_hypo[pid][chrom]
        this_hyper_specific = specific_dmr_loci_hyper[pid][chrom]

        if len(this_hypo):
            hypo_kdes[chrom] = fit_kde_dmr_location(this_hypo, xx, window_size, normed=False)
        else:
            hypo_kdes[chrom] = dummy_res

        if len(this_hyper):
            hyper_kdes[chrom] = fit_kde_dmr_location(this_hyper, xx, window_size, normed=False)
        else:
            hyper_kdes[chrom] = dummy_res

        if len(this_hypo_specific):
            hypo_kdes_specific[chrom] = fit_kde_dmr_location(this_hypo_specific, xx, window_size, normed=False)
        else:
            hypo_kdes_specific[chrom] = dummy_res

        if len(this_hyper_specific):
            hyper_kdes_specific[chrom] = fit_kde_dmr_location(this_hyper_specific, xx, window_size, normed=False)
        else:
            hyper_kdes_specific[chrom] = dummy_res

    # rescale KDEs
    # we want the density to reflect the number of DMRs

    hypo_kdes_rescaled = {}
    hyper_kdes_rescaled = {}
    hypo_kdes_specific_rescaled = {}
    hyper_kdes_specific_rescaled = {}

    hypo_max = max([t.max() for t in hypo_kdes.values()])
    hyper_max = max([t.max() for t in hyper_kdes.values()])

    for chrom in chroms:

        d_hypo = hypo_kdes[chrom]
        hypo_kdes_rescaled[chrom] = d_hypo / hypo_max * kde_width

        d_hyper = hyper_kdes[chrom]
        hyper_kdes_rescaled[chrom] = d_hyper / hyper_max * kde_width

        d_hypo_specific = hypo_kdes_specific[chrom]
        hypo_kdes_specific_rescaled[chrom] = d_hypo_specific / hypo_max * kde_width

        d_hyper_specific = hyper_kdes_specific[chrom]
        hyper_kdes_specific_rescaled[chrom] = d_hyper_specific / hyper_max * kde_width

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')

    curr_theta = 0
    for chrom in chroms:
        this_cg = cg_density[chrom]
        this_cg_pct = this_cg / float(window_size) * 100.

        d_hypo = hypo_kdes_rescaled[chrom]
        d_hyper = hyper_kdes_rescaled[chrom]
        d_hypo_specific = hypo_kdes_specific_rescaled[chrom]
        d_hyper_specific = hyper_kdes_specific_rescaled[chrom]

        xx = np.array([this_cg.index.tolist() + [chrom_length[chrom]]] * 2)
        th = curr_theta + np.array(this_cg.index.tolist() + [chrom_length[chrom]]) * radians_per_bp
        tt = np.array([th, th])
        rr = np.zeros_like(tt) + inner_r; rr[1] = outer_r

        cc = np.array([this_cg_pct.values])

        norm = MidpointNormalize(midpoint=cg_pct_mean, vmin=0, vmax=3.)
        ax.pcolor(tt, rr, cc, cmap='RdYlBu_r')

        this_unmapped = unmapped_density[chrom]
        this_unmapped_pct = this_unmapped / float(window_size) * 100.

        uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmap_threshold_pct)
        # since we don't care about the extent of unmapping, replace all values with a single one
        uu[~uu.mask] = 0.8
        ax.pcolor(tt, rr, uu, cmap='Greys', vmax=1., vmin=0.)

        # plot the KDEs
        ax.fill_between(th, y1=d_hyper + outer_r, y2=outer_r, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hyper'])
        ax.fill_between(th, y1=inner_r, y2=inner_r - d_hypo, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hypo'])

        ax.set_facecolor('w')

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        print "Chrom %s started at angle %.2f and ended at %.2f" % (chrom, curr_theta, th[-1] + gap_radians)
        curr_theta = th[-1] + gap_radians
        
    import ipdb; ipdb.set_trace()

    # the old code (for reference)
    if False:

        fig = plt.figure(figsize=(10, 8))
        nrows = 11
        ncols = 2
        gs_main = plt.GridSpec(
            nrows=nrows,
            ncols=ncols,
            left=0.01,
            right=.99,
            bottom=0.01,
            top=.99,
            hspace=0.03,
            wspace=0.03,
        )

        # we want to traverse by row then column, so create an array of axes beforehand
        # order='C' would give us traversal by column then row
        main_ax_arr = np.array([gs_main[i] for i in range(len(chroms))]).reshape((nrows, ncols)).flatten(order='F')

        # TODO: this will be in a loop once it's finalised
        for i, chrom in enumerate(chroms):
            gs = gridspec.GridSpecFromSubplotSpec(
                3,
                1,
                main_ax_arr[i],
                height_ratios=[5, 1, 5],
                hspace=0.
            )

            # gs = plt.GridSpec(nrows=3, ncols=1, height_ratios=[5, 1, 5])
            ax_gc = fig.add_subplot(gs[1])
            ax_hyper = fig.add_subplot(gs[0], sharex=ax_gc)
            ax_hypo = fig.add_subplot(gs[2], sharex=ax_gc)

            this_cg = cg_density[chrom]
            this_cg_pct = this_cg / float(window_size) * 100.

            xx = np.array([this_cg.index.tolist() + [chrom_length[chrom]]] * 2)
            yy = np.zeros_like(xx); yy[1] = 1.
            cc = np.array([this_cg_pct.values])

            this_unmapped = unmapped_density[chrom]
            this_unmapped_pct = this_unmapped / float(window_size) * 100.

            uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmap_threshold_pct)
            # since we don't care about the extent of unmapping, replace all values with a single one
            uu[~uu.mask] = 0.3

            ax_gc.pcolor(xx, yy, cc, cmap='copper_r', vmax=5., vmin=0.)
            ax_gc.pcolor(xx, yy, uu, cmap='Greys', vmax=1., vmin=0.)
            ax_gc.set_xlim([0, xmax])

            this_hypo = dmr_loci_hypo[pid][chrom]
            this_hyper = dmr_loci_hyper[pid][chrom]

            # KDE estimation gives us a nice representation of the DMR location distribution
            # NB this library expects a 2D array and complains otherwise!
            k_hypo = KernelDensity(bandwidth=window_size, kernel='gaussian')
            k_hypo.fit(np.array(this_hypo)[:, None])  # this increases the dim of the 1D array
            logd_hypo = k_hypo.score_samples(xx[0, None].transpose())

            # blank out baseline (for plotting purposes)
            logd_hypo[logd_hypo < -100] = -np.inf

            d_hypo = np.ma.masked_equal(np.exp(logd_hypo), 0.)
            hypo_max = d_hypo.max()

            k_hyper = KernelDensity(bandwidth=window_size, kernel='gaussian')
            k_hyper.fit(np.array(this_hyper)[:, None])  # this increases the dim of the 1D array
            logd_hyper = k_hyper.score_samples(xx[0, None].transpose())

            # blank out baseline (for plotting purposes)
            logd_hyper[logd_hyper < -100] = -np.inf

            d_hyper = np.ma.masked_equal(np.exp(logd_hyper), 0.)
            hyper_max = d_hyper.max()

            # plot the KDEs
            ax_hyper.fill_between(xx[0], d_hyper, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hyper'])
            ax_hyper.set_ylim([0, hyper_max * 1.02])

            # hypo: needs to be plotted upside down
            ax_hypo.invert_yaxis()
            ax_hypo.fill_between(xx[0], d_hypo, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hypo'])
            ax_hypo.set_ylim([hypo_max * 1.02, 0.])

            # plotting tick marks is a nice idea, but in practice gets messy - may work for smaller datasets?
            # ax_hyper.plot(xx[0], np.full_like(xx[0], -0.1 * d_hyper.max()), '|k', markeredgewidth=1)

            for ax in [ax_gc, ax_hypo, ax_hyper]:
                ax.set_facecolor('w')
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)

        gs_main.update(right=1.2)
        fig.savefig(os.path.join(outdir, "full_dmr_location_plot_%s.png" % pid), dpi=200)
        

    # 2) Patient-specific DMRs
    tmp = get_binned_dmr_locations(
        dmr_res_specific,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True
    )

    dmr_loci_hyper = tmp['dmr_loci_hyper']
    dmr_loci_hypo = tmp['dmr_loci_hypo']
    dmr_binned_hyper = tmp['dmr_binned_hyper']
    dmr_binned_hypo = tmp['dmr_binned_hypo']

    unmap_threshold_pct = 10. # unmapped % above this value will be masked
    xmax = max(chrom_length.values())
    pid = pids[0]

    fig = plt.figure(figsize=(10, 8))
    nrows = 11
    ncols = 2
    gs_main = plt.GridSpec(
        nrows=nrows,
        ncols=ncols,
        left=0.01,
        right=.99,
        bottom=0.01,
        top=.99,
        hspace=0.03,
        wspace=0.03,
    )

    # we want to traverse by row then column, so create an array of axes beforehand
    # order='C' would give us traversal by column then row
    main_ax_arr = np.array([gs_main[i] for i in range(len(chroms))]).reshape((nrows, ncols)).flatten(order='F')

    # TODO: this will be in a loop once it's finalised
    for i, chrom in enumerate(chroms):
        gs = gridspec.GridSpecFromSubplotSpec(
            3,
            1,
            main_ax_arr[i],
            height_ratios=[5, 1, 5],
            hspace=0.
        )

        # gs = plt.GridSpec(nrows=3, ncols=1, height_ratios=[5, 1, 5])
        ax_gc = fig.add_subplot(gs[1])
        ax_hyper = fig.add_subplot(gs[0], sharex=ax_gc)
        ax_hypo = fig.add_subplot(gs[2], sharex=ax_gc)

        this_cg = cg_density[chrom]
        this_cg_pct = this_cg / float(window_size) * 100.

        xx = np.array([this_cg.index.tolist() + [chrom_length[chrom]]] * 2)
        yy = np.zeros_like(xx); yy[1] = 1.
        cc = np.array([this_cg_pct.values])

        this_unmapped = unmapped_density[chrom]
        this_unmapped_pct = this_unmapped / float(window_size) * 100.

        uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmap_threshold_pct)
        # since we don't care about the extent of unmapping, replace all values with a single one
        uu[~uu.mask] = 0.3

        ax_gc.pcolor(xx, yy, cc, cmap='copper_r', vmax=5., vmin=0.)
        ax_gc.pcolor(xx, yy, uu, cmap='Greys', vmax=1., vmin=0.)
        ax_gc.set_xlim([0, xmax])

        this_hypo = dmr_loci_hypo[pid][chrom]
        this_hyper = dmr_loci_hyper[pid][chrom]

        # KDE estimation gives us a nice representation of the DMR location distribution
        # NB this library expects a 2D array and complains otherwise!
        if len(this_hypo) > 0:
            k_hypo = KernelDensity(bandwidth=window_size, kernel='gaussian')
            k_hypo.fit(np.array(this_hypo)[:, None])  # this increases the dim of the 1D array
            logd_hypo = k_hypo.score_samples(xx[0, None].transpose())
        else:
            logd_hypo = np.ones_like(xx[0]) * -np.inf

        # blank out baseline (for plotting purposes)
        logd_hypo[logd_hypo < -100] = -np.inf

        d_hypo = np.ma.masked_equal(np.exp(logd_hypo), 0.)
        hypo_max = d_hypo.max()

        if len(this_hyper) > 0:
            k_hyper = KernelDensity(bandwidth=window_size, kernel='gaussian')
            k_hyper.fit(np.array(this_hyper)[:, None])  # this increases the dim of the 1D array
            logd_hyper = k_hyper.score_samples(xx[0, None].transpose())
        else:
            logd_hyper = np.ones_like(xx[0]) * -np.inf

        # blank out baseline (for plotting purposes)
        logd_hyper[logd_hyper < -100] = -np.inf

        d_hyper = np.ma.masked_equal(np.exp(logd_hyper), 0.)
        hyper_max = d_hyper.max()

        # plot the KDEs
        ax_hyper.fill_between(xx[0], d_hyper, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hyper'])
        ax_hyper.set_ylim([0, hyper_max * 1.02])

        # hypo: needs to be plotted upside down
        ax_hypo.invert_yaxis()
        ax_hypo.fill_between(xx[0], d_hypo, alpha=0.9, color=consts.METHYLATION_DIRECTION_COLOURS['hypo'])
        ax_hypo.set_ylim([hypo_max * 1.02, 0.])

        # plotting tick marks is a nice idea, but in practice gets messy - may work for smaller datasets?
        # ax_hyper.plot(xx[0], np.full_like(xx[0], -0.1 * d_hyper.max()), '|k', markeredgewidth=1)

        for ax in [ax_gc, ax_hypo, ax_hyper]:
            ax.set_facecolor('w')
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

    gs_main.update(right=1.2)
    fig.savefig(os.path.join(outdir, "specific_dmr_location_plot_%s.png" % pid), dpi=200)

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

    if write_bed_files:
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
    # we may only use the Gibco line here, to avoid issues associated with changing the array type / probe selection
    s2_850_only = True

    if s2_850_only:
        norm_method_s2 = 'swan'
        external_ref_names_dm = None
        external_ref_samples_dm = None
        external_refs_dm = [
            ('GIBCO', 'NSC'),
        ]
    else:
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
        patient_samples=consts.ALL_METHYL_SAMPLES  # any samples not found in patient data will be ignored (H9 NSC)
    )
    me_obj_with_ref.meta.insert(
        0,
        'patient_id',
        me_obj_with_ref.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>')
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
        dmp_res_s2_syn,
        me_obj_with_ref.data,
        me_obj_with_ref.meta.type,
        me_obj_with_ref.meta.patient_id,
    )
    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_all_dmp_direction_panel.png"), dpi=200)
    plot_dict['fig'].savefig(os.path.join(outdir, "s2_syngeneic_all_dmp_direction_panel.tiff"), dpi=200)

    plot_dict = dm_probe_direction_panel_plot(
        dmr_res_s2_syn_specific,
        dmp_res_s2_syn_specific,
        me_obj_with_ref.data,
        me_obj_with_ref.meta.type,
        me_obj_with_ref.meta.patient_id,
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

        # remove iNSC samples before plotting panel
        me_dat_for_plot = me_obj_with_ref.data.loc[:,
            (me_obj_with_ref.meta.type == 'GBM') | me_obj_with_ref.data.columns.str.contains(r)
        ]
        me_meta_for_plot = me_obj_with_ref.meta.loc[me_dat_for_plot.columns]

        plot_dict = dm_probe_direction_panel_plot(
            this_dmr_res,
            this_dmp_res,
            me_dat_for_plot,
            me_meta_for_plot.type,
            me_meta_for_plot.patient_id,
            ref_name=r
        )
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_all_dmp_direction_panel.png" % r), dpi=200)
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_all_dmp_direction_panel.tiff" % r), dpi=200)

        plot_dict = dm_probe_direction_panel_plot(
            this_dmr_res_specific,
            this_dmp_res_specific,
            me_obj_with_ref.data,
            me_obj_with_ref.meta.type,
            me_obj_with_ref.meta.patient_id,
        )
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_specific_dmp_direction_panel.png" % r), dpi=200)
        plot_dict['fig'].savefig(os.path.join(outdir, "s2_%s_specific_dmp_direction_panel.tiff" % r), dpi=200)

from plotting import bar, common
from methylation import loader, dmr, process
import pandas as pd
from utils import output, setops, genomics
import multiprocessing as mp
import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR


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


def bed_file_from_probes(anno, probes, out_fn):
    this_regions = {}
    this_anno = anno.loc[probes]
    for p, row in this_anno.iterrows():
        strand = '+' if row.Strand == 'F' else '-'
        # we'll prepend the chrom name with 'chr' to ensure compatibility with hg19 (built in to Homer)
        this_regions[p] = ["chr%s" % row.CHR, row.MAPINFO - probe_half_len, row.MAPINFO + probe_half_len, strand]

    genomics.write_bed_file(this_regions, out_fn)


def direction_of_dm_bar_plot(
        dmr_res,
        venn_set=None,
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
    :param venn_set: If supplied, this is a dictionary keyed by pids with iterable values giving cluster/probe IDs.
    Limit the clusters/probes included to those IDs.
    :param pids: Useful to specify the ordering of the PIDs.
    :param as_pct: If True (default), values are plotted as percentages.
    :param kwargs: Passed to the plotting function `bar.stacked_bar_chart()`. For example, might include `ax=...` to
    specify axis.
    :return:
    """
    dmr_direction = {}

    for pid in pids:
        if venn_set is None:
            keys = dmr_res[pid].keys()
        else:
            keys = venn_set[pid]
        this_res = pd.DataFrame.from_dict([dmr_res[pid][k] for k in keys])
        dmr_direction[pid] = {
            'Hyper': (this_res['median_change'] > 0).sum(),
            'Hypo': (this_res['median_change'] < 0).sum(),
        }
    for_plot = pd.DataFrame.from_dict(dmr_direction)[pids].loc[['Hypo', 'Hyper']]
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


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS

    outdir = output.unique_output_dir("dmr_without_classes")
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
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s1 = tsgd.paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    # 1) direction of methylation change for all DMRs
    # a) clusters: absolute and relative values
    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_sign_s1,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMRs")
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_sign_s1,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "all_dmr_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "all_dmr_direction.tiff"), dpi=200)

    # b) probes: otherwise as above
    dmp_ids = {}
    dmp_res = {}
    for pid in pids:
        dmp_ids[pid] = set()
        dmp_res[pid] = {}
        gbm_samples = me_data.columns[me_data.columns.str.contains('GBM%s' % pid)]
        nsc_samples = me_data.columns[me_data.columns.str.contains('DURA%s' % pid)]
        for cl_id in dmr_res_sign_s1[pid]:
            dmp_ids[pid].update(dmr_res_s1.clusters[cl_id].pids)
        dmp_ids[pid] = sorted(dmp_ids[pid])
        gbm_median = me_data.loc[dmp_ids[pid], gbm_samples].median(axis=1)
        nsc_median = me_data.loc[dmp_ids[pid], nsc_samples].median(axis=1)
        dmp_res[pid] = dict(
            [(k, {'median_change': v}) for k, v in (gbm_median - nsc_median).items()]
        )

    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMPs")
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "all_dmp_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "all_dmp_direction.tiff"), dpi=200)

    # 2) direction of methylation change for patient-specific DMRs
    # a) clusters: absolute and relative values
    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_sign_s1[pid].keys() for pid in pids])
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_sign_s1,
        venn_set=vs,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number patient-specific DMRs")
    fig, ax = direction_of_dm_bar_plot(
        dmr_res_sign_s1,
        venn_set=vs,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_specific_dmr_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_dmr_direction.tiff"), dpi=200)

    # b-i) probes: obtained by converting regions -> probes, then getting those specific to each patient
    # arguably, this is a bit of a strange way to get 'patient-specific probes' (see b-ii)
    dmp_ids = {}
    for pid in pids:
        dmp_ids[pid] = set()
        for cl_id in dmr_res_sign_s1[pid]:
            dmp_ids[pid].update(dmr_res_s1.clusters[cl_id].pids)
        dmp_ids[pid] = sorted(dmp_ids[pid])

    venn_set, venn_ct = setops.venn_from_arrays(*[dmp_ids[pid] for pid in pids])
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    dmp_res = {}
    for pid in pids:
        dmp_res[pid] = {}
        gbm_samples = me_data.columns[me_data.columns.str.contains('GBM%s' % pid)]
        nsc_samples = me_data.columns[me_data.columns.str.contains('DURA%s' % pid)]
        gbm_median = me_data.loc[vs[pid], gbm_samples].median(axis=1)
        nsc_median = me_data.loc[vs[pid], nsc_samples].median(axis=1)
        dmp_res[pid] = dict(
            [(k, {'median_change': v}) for k, v in (gbm_median - nsc_median).items()]
        )

    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMPs")
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_specific_direct_dmp_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_direct_dmp_direction.tiff"), dpi=200)

    # b-ii) probes: obtained by converting patient-specific DMRs to probes
    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_sign_s1[pid].keys() for pid in pids])
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids
    vs_dmrs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    dmp_ids = {}
    for pid in pids:
        dmp_ids[pid] = set()
        for cl_id in vs_dmrs[pid]:
            dmp_ids[pid].update(dmr_res_s1.clusters[cl_id].pids)
        dmp_ids[pid] = sorted(dmp_ids[pid])

    venn_set, venn_ct = setops.venn_from_arrays(*[dmp_ids[pid] for pid in pids])
    vs = dict([
        (p, venn_set[q]) for p, q in zip(pids, pu_sets)
    ])

    dmp_res = {}
    for pid in pids:
        dmp_res[pid] = {}
        gbm_samples = me_data.columns[me_data.columns.str.contains('GBM%s' % pid)]
        nsc_samples = me_data.columns[me_data.columns.str.contains('DURA%s' % pid)]
        gbm_median = me_data.loc[vs[pid], gbm_samples].median(axis=1)
        nsc_median = me_data.loc[vs[pid], nsc_samples].median(axis=1)
        dmp_res[pid] = dict(
            [(k, {'median_change': v}) for k, v in (gbm_median - nsc_median).items()]
        )

    gs = plt.GridSpec(2, 1, height_ratios=[1, 4])
    fig = plt.figure(figsize=(6, 7))
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        as_pct=False,
        ax=fig.add_subplot(gs[1]),
        legend_loc='upper right'
    )
    ax.set_ylabel("Number DMPs")
    fig, ax = direction_of_dm_bar_plot(
        dmp_res,
        ax=fig.add_subplot(gs[0]),
        legend=False
    )
    ax.xaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "patient_specific_dmp_via_dmr_direction.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_dmp_via_dmr_direction.tiff"), dpi=200)


    ## to move to function
    dmr_res = dmr_res_sign_s1
    venn_set = None

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
                v[pid] = set(me_obj.data.index)
            else:
                v[pid] = set()

        this_res = pd.DataFrame.from_dict(dmr_res[pid])
        for k in dmr_res[pid].keys():
            this_probe_ids = dmr_res_s1.clusters[k].pids
            if venn_set is not None:
                this_probe_ids = set(this_probe_ids).intersection(venn_set[pid])
            pid_sets['dmr'].update(this_probe_ids)
            mc = dmr_res[pid][k]['median_change']
            if mc < 0:
                pid_sets['hypo'].update(this_probe_ids)
            else:
                pid_sets['hyper'].update(this_probe_ids)





    # for pid in pids:
    #     if venn_set is None:
    #         keys = dmr_res[pid].keys()
    #     else:
    #         keys = venn_set[pid]
    #     this_res = pd.DataFrame.from_dict([dmr_res[pid][k] for k in keys])
    #     dmr_direction[pid] = {
    #         'Hyper': (this_res['median_change'] > 0).sum(),
    #         'Hypo': (this_res['median_change'] < 0).sum(),
    #     }


    raise StopIteration

    ## island status

    pid_sets = {
        'background': {},
        'dmr': {},
        'hypo': {},
        'hyper': {},
    }

    for pid, s in zip(pids, pu_sets):
        for k, v in pid_sets.items():
            if k == 'background':
                v[pid] = set(me_obj.data.index)
            else:
                v[pid] = set()

        the_dmr_res = dmr_res[pid]
        cids = venn_set[s]
        this_res = pd.DataFrame.from_dict([the_dmr_res.results[t] for t in cids])

        for i, c in enumerate(cids):
            the_probe_ids = the_dmr_res.clusters[c].pids

            pid_sets['background'][pid].difference_update(the_probe_ids)
            pid_sets['dmr'][pid].update(the_probe_ids)

            if the_dmr_res.results[c]['median_change'] > 0:
                # hypermethylation
                pid_sets['hyper'][pid].update(the_probe_ids)
            else:
                # hypomethylation
                pid_sets['hypo'][pid].update(the_probe_ids)

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

    pid_sets_by_island_status = {}
    for pid_typ in ['dmr', 'hypo', 'hyper']:
        pid_sets_by_island_status[pid_typ] =  {}
        for pid in pids:
            pid_sets_by_island_status[pid_typ][pid] = {}
            probes = pid_sets[pid_typ][pid]
            this_cats = anno.loc[probes, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea)
            for k, v in cats.items():
                this_cat_idx = this_cats.loc[probes] == k
                pid_sets_by_island_status[pid_typ][pid].setdefault(v, set()).update(this_cat_idx.index[this_cat_idx])

    # sanity check
    for pid in pids:
        if not (pd.Series(island_counts['hyper'][pid]) + pd.Series(island_counts['hypo'][pid]) == pd.Series(island_counts['dmr'][pid])).all():
            raise ValueError("PID %s failed check # hypo + # hyper = # dmr" % pid)

    # all
    this_counts = anno.loc[:, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea).value_counts().to_dict()
    island_counts_all = dict([
        (v, this_counts.get(k, 0)) for k, v in cats.items()
    ])

    # save this in a 'nice' format for sharing
    cols = sorted(set(cats.values()))
    to_export = pd.DataFrame(
        index=pd.MultiIndex.from_product([pids, ['background', 'dmr', 'hypo', 'hyper']], names=['patient ID', 'probe list']),
        columns=cols
    )

    for pid_typ, pid_set in island_counts.items():
        for pid in pids:
            to_export.loc[(pid, pid_typ)] = pd.Series(pid_set[pid])[cols]

    to_export.loc[('all', 'all'), cols] = pd.Series(island_counts_all)[cols]

    # nice plot showing these differing distributions
    dist_as_pct = to_export.divide(to_export.sum(axis=1), axis=0) * 100.

    fig, axs = plt.subplots(ncols=len(pids), sharex=True, sharey=True, figsize=(11.5, 5.5))

    for i, pid in enumerate(['017', '019', '030', '031', '018', '050', '054', '061', '026', '052']):
        ax = axs.flat[i]
        X = dist_as_pct.loc[pid].transpose()
        X.columns = ['B/G', 'All DMRs', 'Hypo', 'Hyper']
        bar.stacked_bar_chart(X, legend=False, ax=ax)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        ax.set_title(pid)
        ax.set_ylim([0, 100])
    axs[0].set_ylabel('Percentage of probes')
    axs[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.subplots_adjust(left=0.05, right=0.89, bottom=0.15, top=0.95, wspace=0.05)
    fig.savefig(os.path.join(outdir, "probe_cpg_dist.png"), dpi=200)

    # export the probe regions for each patient to a BED file for motif enrichment
    # probes are CpG +/- 60 bp: https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_cpg_loci_identification.pdf
    probe_half_len = 61
    for pid in pids:
        for typ in ['dmr', 'hyper', 'hypo']:
            ps = pid_sets[typ][pid]
            bed_fn = os.path.join(outdir, "%s_%s_oligo_mappings.bed" % (pid, typ))
            bed_file_from_probes(anno, ps, bed_fn)

    # repeat with even more granularity: include probe types
    for pid in pids:
        for typ in ['dmr', 'hyper', 'hypo']:
            for cs in cats.values():
                ps = pid_sets_by_island_status[typ][pid][cs]
                bed_fn = os.path.join(outdir, "%s_%s_%s_oligo_mappings.bed" % (pid, typ, cs))
                bed_file_from_probes(anno, ps, bed_fn)


    # for each patient, repeat this process but with the full set of DMRs
    n_dmr_by_direction_full = {}
    pid_sets_full = {}

    for p in pids:
        the_dmr_res = dmr_res[p]
        cids, attrs = zip(*the_dmr_res.results_significant.items())
        n_dmr_by_direction_full[p] = {
            'Hyper': len([t for t in attrs if t['median_change'] > 0]),
            'Hypo': len([t for t in attrs if t['median_change'] < 0]),
        }
        # probes
        the_probe_ids_dmr = set()
        the_probe_ids_hypo = set()
        the_probe_ids_hyper = set()
        for c, a in zip(cids, attrs):
            the_probe_ids_dmr.update(the_dmr_res.clusters[c].pids)
            if a['median_change'] > 0:
                the_probe_ids_hyper.update(the_dmr_res.clusters[c].pids)
            else:
                the_probe_ids_hypo.update(the_dmr_res.clusters[c].pids)
        pid_sets_full.setdefault('dmr', {})[p] = the_probe_ids_dmr
        pid_sets_full.setdefault('hyper', {})[p] = the_probe_ids_hyper
        pid_sets_full.setdefault('hypo', {})[p] = the_probe_ids_hypo

    pid_sets_full_by_island_status = {}
    for pid_typ in ['dmr', 'hypo', 'hyper']:
        pid_sets_full_by_island_status[pid_typ] =  {}
        for pid in pids:
            pid_sets_full_by_island_status[pid_typ][pid] = {}
            probes = pid_sets_full[pid_typ][pid]
            this_cats = anno.loc[probes, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea)
            for k, v in cats.items():
                this_cat_idx = this_cats.loc[probes] == k
                pid_sets_full_by_island_status[pid_typ][pid].setdefault(v, set()).update(this_cat_idx.index[this_cat_idx])

    # another plot of hypo vs hyper cluster counts, this time on the full lists

    for_plot_full = pd.DataFrame.from_dict(n_dmr_by_direction_full)[[
        '017', '019', '030', '031', '018', '050', '054', '061', '026', '052'
    ]].loc[['Hypo', 'Hyper']]
    for_plot_pct_full = for_plot_full.divide(for_plot_full.sum(), axis=1) * 100.
    colours = pd.Series({'Hyper': '#FF381F', 'Hypo': '#89CD61'})
    fig, ax = bar.stacked_bar_chart(for_plot_pct_full, colours=colours)
    ax.set_ylabel("Percentage of clusters")
    ax.set_ylim([0, 100])
    # shrink main axis and put legend on RHS
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(os.path.join(outdir, "pct_clusters_by_dm_direction_full.png"), dpi=200)

    # export bed regions for the full lists
    for pid in pids:
        for typ in ['dmr', 'hyper', 'hypo']:
            ps = pid_sets_full[typ][pid]
            bed_fn = os.path.join(outdir, "%s_%s_oligo_mappings_full.bed" % (pid, typ))
            bed_file_from_probes(anno, ps, bed_fn)

    # output to bed file
    for pid in pids:
        for typ in ['dmr', 'hyper', 'hypo']:
            for cs in cats.values():
                ps = pid_sets_full_by_island_status[typ][pid][cs]
                bed_fn = os.path.join(outdir, "%s_%s_%s_oligo_mappings_full.bed" % (pid, typ, cs))
                bed_file_from_probes(anno, ps, bed_fn)
from plotting import bar, venn
from methylation import loader, dmr, process
import pandas as pd
from utils import output, setops, genomics
import multiprocessing as mp
import os
from matplotlib import pyplot as plt
import seaborn as sns


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None):
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
    me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    # dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


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


def bed_file_from_regions():
    pass


if __name__ == "__main__":
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    norm_method_s1 = 'swan'
    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }
    outdir = output.unique_output_dir("dmr_without_classes")

    me_obj, anno = load_methylation(pids, norm_method=norm_method_s1)
    clusters = {}
    for cc in anno.CHR.unique():
        coords = anno.loc[anno.CHR == cc, 'MAPINFO'].sort_values()
        clusters[cc] = dmr.identify_cluster(coords, dmr_params['n_min'], dmr_params['d_max'])

    cl = []
    ix = 0
    for _, the_cl in clusters.items():
        for cid, x in the_cl.items():
            cl.append(
                dmr.ProbeCluster(x, anno, cluster_id=ix)
            )
            ix += 1

    obj = dmr.DmrResults(clusters=cl, anno=anno)

    dmr_res = pair_dmr(me_obj.meta, me_obj.data, obj, pids, **dmr_params)

    # identify patient-specific DMRs
    dmr_by_member = [dmr_res[pid].clusters_significant.keys() for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids

    ## balance of direction of change for patient-specific clusters

    dmr_direction = {}
    for pid, s in zip(pids, pu_sets):
        this_res = pd.DataFrame.from_dict([dmr_res[pid].results[t] for t in venn_set[s]])
        dmr_direction[pid] = {
            'Hyper': (this_res.median_change > 0).sum(),
            'Hypo': (this_res.median_change < 0).sum(),
        }
    for_plot = pd.DataFrame.from_dict(dmr_direction)[[
        '017', '019', '030', '031', '018', '050', '054', '061', '026', '052'
    ]].loc[['Hypo', 'Hyper']]
    for_plot_pct = for_plot.divide(for_plot.sum(), axis=1) * 100.
    colours = pd.Series({'Hyper': '#FF381F', 'Hypo': '#89CD61'})

    fig, ax = bar.stacked_bar_chart(for_plot_pct, colours=colours)
    ax.set_ylabel("Percentage of clusters")
    ax.set_ylim([0, 100])
    # shrink main axis and put legend on RHS
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(os.path.join(outdir, "pct_clusters_by_dm_direction.png"), dpi=200)

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
            this_counts = anno.loc[p, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea).value_counts().to_dict()
            island_counts.setdefault(pid_typ, {}).setdefault(pid, dict(empty_counts))
            for k, v in cats.items():
                island_counts[pid_typ][pid][v] += this_counts.get(k, 0)

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
            this_regions = {}
            ps = pid_sets[typ][pid]
            this_anno = anno.loc[ps]
            for p, row in this_anno.iterrows():
                strand = '+' if row.Strand == 'F' else '-'
                # we'll prepend the chrom name with 'chr' to ensure compatibility with hg19 (built in to Homer)
                this_regions[p] = ["chr%s" % row.CHR, row.MAPINFO - probe_half_len, row.MAPINFO + probe_half_len, strand]

            bed_fn = os.path.join(outdir, "%s_%s_oligo_mappings.bed" % (pid, typ))
            genomics.write_bed_file(this_regions, bed_fn)

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
            this_regions = {}
            ps = pid_sets_full[typ][pid]
            this_anno = anno.loc[ps]
            for p, row in this_anno.iterrows():
                strand = '+' if row.Strand == 'F' else '-'
                # we'll prepend the chrom name with 'chr' to ensure compatibility with hg19 (built in to Homer)
                this_regions[p] = ["chr%s" % row.CHR, row.MAPINFO - probe_half_len, row.MAPINFO + probe_half_len, strand]

            bed_fn = os.path.join(outdir, "%s_%s_oligo_mappings_full.bed" % (pid, typ))
            genomics.write_bed_file(this_regions, bed_fn)


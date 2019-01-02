from plotting import bar, common
from methylation import loader, dmr, process
import pandas as pd
from statsmodels.sandbox.stats import multicomp
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import collections
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR
logger = log.get_console_logger()

"""
Here we analyse the direction and genomic locations of DMRs in the individual patients.
We note a bias in the direction, which is amplified considerably when we consider patient-specific DMRs.
"""

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


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

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

    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    window_size = int(2e5)
    coord_summary_method = 'first'
    cg_locs, chrom_length = genomics.motif_locations(fa_fn, features=chroms)

    # reorder
    chrom_length = collections.OrderedDict([(k, chrom_length[k]) for k in chroms])

    dmr_s1_clusters = dmr_res_s1[pids[0]].clusters

    unmapped_density, _ = genomics.cg_content_windowed(fa_fn, features=chroms, window_size=window_size, motif='N')

    # get all DMRs (whether significant or not) to use as a null hypothesis
    tmp = get_binned_dmr_locations(
        dmr_res_s1.results,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=False,
        coord_summary_method=coord_summary_method
    )
    dmr_loci_all = tmp['dmr_loci'][pids[0]]

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


    def ks_test_dmr_locations(
        dmr_locs,
        null_locs,
        mht_method='fdr_bh',
    ):
        """
        Run the 2 sample KS test on DMR locations, using the supplied null locations as the reference sample.
        :param dmr_locs: Nested dictionary keyed by PID then feature (chromosome name)
        :param null_locs: Dictionary keyed by feature (chromosome name)
        :param mht_method: The method to use for MHT correction. Passed to multicomp.multipletests.
        If None, skip MHT and use the unadjusted pvalues.
        :return:
        """
        ks_pval_flat = []
        ks_stat_flat = []
        ks_flat_map = []
        for pid, val in dmr_locs.items():
            for feat in val:
                tmp = stats.ks_2samp(dmr_locs[pid][feat], null_locs[feat])
                ks_stat_flat.append(tmp[0])
                ks_pval_flat.append(tmp[1])
                ks_flat_map.append((pid, feat))

        # correct for MHT
        if mht_method is None:
            padj = ks_pval_flat
        else:
            _, padj, _, _ = multicomp.multipletests(ks_pval_flat, alpha=0.05, method=mht_method)


        res_fdr = dict([(pid, {}) for pid in dmr_locs.keys()])
        res_stat = dict([(pid, {}) for pid in dmr_locs.keys()])

        i = 0
        for pid, feat in ks_flat_map:
            res_fdr[pid][feat] = padj[i]
            res_stat[pid][feat] = ks_stat_flat[i]
            i += 1

        return res_fdr, res_stat

    hyper_all_fdr, hyper_all_stat = ks_test_dmr_locations(
        dmr_loci_hyper,
        dmr_loci_all
    )
    hyper_all_fdr = pd.DataFrame(hyper_all_fdr).loc[chroms, pids]

    hypo_all_fdr, hypo_all_stat = ks_test_dmr_locations(
        dmr_loci_hypo,
        dmr_loci_all
    )
    hypo_all_fdr = pd.DataFrame(hypo_all_fdr).loc[chroms, pids]

    norm = MidpointNormalize(midpoint=1., vmin=0, vmax=3.)
    gs = plt.GridSpec(1, 3, width_ratios=[9, 9, 1])
    fig = plt.figure(figsize=(9.5, 5.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
    ax_cbar = fig.add_subplot(gs[2])

    sns.heatmap(-np.log10(hyper_all_fdr), cmap='RdBu_r', norm=norm, vmin=0., vmax=3., ax=ax_hyper, cbar=False)
    sns.heatmap(-np.log10(hypo_all_fdr), cmap='RdBu_r', norm=norm, vmin=0., vmax=3., ax=ax_hypo, cbar_ax=ax_cbar)
    plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
    ax_hyper.set_ylabel('Chromosome')
    ax_hyper.set_xlabel('Patient')
    ax_hyper.set_title('Hypermethylated DMRs')
    ax_hypo.set_xlabel('Patient')
    ax_hypo.set_title('Hypomethylated DMRs')
    ax_cbar.set_ylabel('$-\log_{10}(p)$')
    plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

    gs.update(left=0.06, right=0.94, top=0.95, wspace=0.1)

    # repeat for patient-specific DMRs

    tmp = get_binned_dmr_locations(
        dmr_res_specific,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )

    dmr_loci_hyper_specific = tmp['dmr_loci_hyper']
    dmr_loci_hypo_specific = tmp['dmr_loci_hypo']

    hyper_specific_fdr, hyper_specific_stat = ks_test_dmr_locations(
        dmr_loci_hyper_specific,
        dmr_loci_all
    )
    hyper_specific_fdr = pd.DataFrame(hyper_specific_fdr).loc[chroms, pids]

    hypo_specific_fdr, hypo_specific_stat = ks_test_dmr_locations(
        dmr_loci_hypo_specific,
        dmr_loci_all
    )
    hypo_specific_fdr = pd.DataFrame(hypo_specific_fdr).loc[chroms, pids]

    norm = MidpointNormalize(midpoint=1., vmin=0, vmax=3.)
    gs = plt.GridSpec(1, 3, width_ratios=[9, 9, 1])
    fig = plt.figure(figsize=(9.5, 5.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
    ax_cbar = fig.add_subplot(gs[2])

    sns.heatmap(-np.log10(hyper_specific_fdr), cmap='RdBu_r', norm=norm, vmin=0., vmax=3., ax=ax_hyper, cbar=False)
    sns.heatmap(-np.log10(hypo_specific_fdr), cmap='RdBu_r', norm=norm, vmin=0., vmax=3., ax=ax_hypo, cbar_ax=ax_cbar)
    plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
    ax_hyper.set_ylabel('Chromosome')
    ax_hyper.set_xlabel('Patient')
    ax_hyper.set_title('Hypermethylated DMRs')
    ax_hypo.set_xlabel('Patient')
    ax_hypo.set_title('Hypomethylated DMRs')
    ax_cbar.set_ylabel('$-\log_{10}(p)$')
    plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

    gs.update(left=0.06, right=0.94, top=0.95, wspace=0.1)

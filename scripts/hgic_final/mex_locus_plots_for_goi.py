import os
import collections
import numpy as np
import pickle
from matplotlib import pyplot as plt, colors
import pandas as pd

from utils import log, output, setops
from plotting import genomics as plt_genomics
from scripts.hgic_final import two_strategies_combine_de_dmr as tscdd
from methylation import dmr
from rnaseq import loader as rnaseq_loader
from integrator import rnaseq_methylationarray
import consts
from settings import INTERMEDIATE_DIR

logger = log.get_console_logger()

"""
Here we apply a similar plotting approach to that in two_strategies_combine_de_dmr, but this time we 
do so with an explicitly specified list of genes of interest (rather than an automated selection criteria).
"""
if __name__ == '__main__':
    GOI = [
        'PTGER4',
        'NTRK2',
        'ALDH3B1'
    ]

    pids = consts.PIDS

    de_params = consts.DE_PARAMS
    dmr_params = consts.DMR_PARAMS
    norm_method_s1 = 'swan'

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

    outdir = output.unique_output_dir()

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # some quantities relating to set membership

    # load methylation
    # For S1, we just need the paired comparison. This is important - any other samples lead to a change in the final
    # probe list (any NA rows are dropped from ALL samples).
    me_obj, anno = tscdd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tscdd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    # get samples used in each comparison
    dmr_comparison_groups = collections.OrderedDict([(pid, {}) for pid in consts.PIDS])
    gg = me_data.columns.groupby(zip(me_meta.patient_id, me_meta.type))
    for (pid, typ), samples in gg.items():
        dmr_comparison_groups[pid][typ] = samples


    rnaseq_obj = rnaseq_loader.load_by_patient(pids)  # quicker than tscdd method that loads refs too

    # rnaseq_obj = tscdd.load_rnaseq(
    #     pids,
    #     external_ref_names_de,
    #     strandedness=external_ref_strandedness_de,
    # )

    rnaseq_obj.filter_by_sample_name(consts.ALL_RNASEQ_SAMPLES)

    # only keep the required syngeneic samples for this analysis
    dat_s1 = rnaseq_obj.data.loc[
             :, rnaseq_obj.meta.index.isin(consts.S1_RNASEQ_SAMPLES)
             ]
    meta_s1 = rnaseq_obj.meta.loc[dat_s1.columns]

    the_hash = tscdd.de_results_hash(meta_s1.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise AttributeError("Unable to load pre-computed DE results, expected at %s" % fn)

    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

    # get the joint table
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    # link specified GOIs with DMRs
    de_dm_list = []
    for g in GOI:
        # find all DMRs linked to this
        the_clusters = [k for k, v in dmr_res_s1.clusters.items() if any([t[0] == g for t in v.genes])]
        if len(the_clusters) == 0:
            logger.warning("Gene %s has no associated DMRs. Skipping", g)
        else:
            logger.info("Gene %s: found %d associated DMRs.", g, len(the_clusters))

        # plot everything
        de_dm_list.extend([(c, g) for c in the_clusters])

    # set plotting attributes
    plot_colours = {'GBM': '#de8100', 'iNSC': '#1f8bff'}
    plot_markers = {'GBM': 'o', 'iNSC': '^'}
    plot_zorder = {'GBM': 21, 'iNSC': 20}
    plot_alpha = {'GBM': 0.5, 'iNSC': 0.7}

    reds = plt.cm.Reds(np.linspace(.1, .7, 128))
    greens = plt.cm.Greens_r(np.linspace(0.3, .9, 128))
    colours_comb = np.vstack((greens, reds))
    cmap = colors.LinearSegmentedColormap.from_list("mex_cmap", colours_comb)

    # debug / test
    # plt_obj = plt_genomics.MexLocusPlotter()
    # gene_list = sorted(ps_de_dm_list)
    # plt_obj.set_plot_parameters(
    #     colours=plot_colours,
    #     markers=plot_markers,
    #     zorder=plot_zorder,
    #     alpha=plot_alpha,
    #     de_direction_colours=cmap,
    #     dm_direction_colours=cmap,
    #     de_vmin=-5,
    #     de_vmax=5,
    #     dm_vmin=-8,
    #     dm_vmax=8
    # )
    # plt_obj.set_mvalues(me_data)
    # plt_obj.set_dmr_res(dmr_res_s1, dmr_comparison_groups)
    # plt_obj.set_de_res(de_res_full_s1)
    # g = gene_list[0][1]
    # the_obj = plt_obj.plot_gene(g)
    # the_obj.set_title(g)

    # shortlist
    tscdd.multipage_pdf_mex_plots(
        os.path.join(outdir, "de_dmr_mex_locus_plots.pdf"),
        sorted(de_dm_list),
        me_data,
        dmr_res_s1,
        dmr_comparison_groups,
        de_res_full_s1,
        plot_colours=plot_colours,
        plot_markers=plot_markers,
        plot_zorder=plot_zorder,
        plot_alpha=plot_alpha,
        direction_cmap=cmap
    )
import collections
import os
import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler

from plotting import common, heatmap
from rnaseq import loader as rnaseq_loader
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from settings import HGIC_LOCAL_DIR, INTERMEDIATE_DIR
from utils import output, setops, log, reference_genomes

logger = log.get_console_logger()


if __name__ == '__main__':
    """
    Extract some specific pathways from the IPA results and plot the gene expression and logFC to show how the
    constitutent parts behave across patient lines.
    """

    # pathways of interest: allow combining multiple
    poi = {
        'mtor': ('mTOR Signaling', 'EIF2 Signaling'),
        'gpcr': ('G-Protein Coupled Receptor Signaling', 'cAMP-mediated signaling')
    }

    # set a minimum pval for pathways to be used
    alpha = 0.01
    plogalpha = -np.log10(alpha)

    pids = consts.PIDS
    de_params = consts.DE_PARAMS

    patient_colours = {
        '018': '#ccffcc',
        '019': '#4dff4d',
        '030': '#00cc00',
        '031': '#004d00',
        '017': '#ffcccc',
        '050': '#ff4d4d',
        '054': '#cc0000',
        '061': '#660000',
        '026': '#ff80ff',
        '052': '#800080'
    }

    de_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/ipa/pathways')
    dm_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/methylation/s0_individual_patients_direct_comparison/ipa/pathways')
    dedm_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq_methylation_combined/s0_individual_patients_direct_comparison/ipa/pathways')

    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

    outdir = output.unique_output_dir()

    #######################################################
    # DE
    #######################################################
    # data
    rnaseq_obj = rnaseq_loader.load_by_patient(pids, include_control=False, source='star')
    rnaseq_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    dat_s1 = rnaseq_obj.data
    meta_s1 = rnaseq_obj.meta

    cpm = dat_s1.divide(dat_s1.sum(axis=0), axis=1) * 1e6

    # DE results
    the_hash = tsgd.de_results_hash(meta_s1.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise AttributeError("Unable to find pre-computed S1 comparison results.")

    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])
    vs, vc = setops.venn_from_arrays(*[de_res_s1[pid].index for pid in pids])
    de_res_wide = setops.venn_set_to_wide_dataframe(de_res_s1, vs, pids, cols_to_include=['logFC', 'FDR'])

    ipa_de_res = collections.OrderedDict()
    for pid in pids:
        fn = os.path.join(de_indir, "full_de_patient{pid}.xls".format(pid=pid))
        this_df = pd.read_excel(fn, skiprows=1, header=0, index_col=0)
        this_df.columns = ['-logp', 'ratio', 'z', 'genes']
        this_df.insert(3, 'n_gene', this_df.genes.str.split(',').apply(len))
        # filter to include only relevant pathways
        ipa_de_res[pid] = this_df.loc[this_df['-logp'] >= plogalpha]

    # for plotting
    groups = [
        (pid, dat_s1.columns[dat_s1.columns.str.contains(pid)]) for pid in pids
    ]

    for name, p_arr in poi.items():
        this_ipa = dict([
            (pid, ipa_de_res[pid].reindex(list(p_arr))) for pid in pids
        ])

        # get all genes involved
        union_genes = set()
        for v in this_ipa.values():
            for arr in v.dropna().genes.str.split(',').values:
                union_genes.update(arr)

        union_genes_ens = reference_genomes.gene_symbol_to_ensembl(union_genes)
        union_genes_comp = union_genes_ens.reset_index().set_index('Ensembl Gene ID').squeeze()

        # plot TPM
        this_cpm = cpm.loc[union_genes_ens]

        # sort by mean logFC across patients
        this_de = de_res_wide.loc[union_genes_ens]
        this_logfc = this_de.loc[:, this_de.columns.str.contains('logFC')]
        this_logfc.columns = pids
        sort_ix = (
            this_logfc.sum(axis=1) / float(len(pids))
        ).sort_values(ascending=False).index
        this_logfc = this_logfc.loc[sort_ix]
        this_cpm = this_cpm.loc[sort_ix]

        this_logfc.index = union_genes_comp.loc[sort_ix]
        this_cpm.index = union_genes_comp.loc[sort_ix]


        # standardise: mean centred data required for sensible decomposition
        # standardisation occurs along the FEATURES axis, which is dim 1
        scaler = StandardScaler(with_std=True)

        # features on the ROWS, mean centre by gene
        scaler = scaler.fit(this_cpm.transpose())
        this_cpm_z = pd.DataFrame(
            scaler.transform(this_cpm.transpose()).transpose(),
            index=this_cpm.index,
            columns=this_cpm.columns
        )

        fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            groups,
            this_cpm_z.transpose(),
            fig_kwargs={'figsize': (7, 10)},
            heatmap_kwargs={'square': False}
        )

        gs.update(left=0.15, bottom=0.2, right=0.95, top=0.95)
        if this_cpm.shape[0] > 50:
            plt.setp(axs[0].yaxis.get_ticklabels(), fontsize=6)
        for ax in axs:
            ax.xaxis.label.set_visible(False)

        fig.savefig(os.path.join(outdir, "%s_standardised_cpm_heatmap.png" % name), dpi=200)

        # now repeat with similar plot, but showing logFC (only)
        fig = plt.figure(figsize=(4.5, 9.5))
        ax = fig.add_subplot(111)
        h = ax.pcolor(
            np.ma.masked_where(this_logfc.isnull(), this_logfc.astype(float))[::-1],
            cmap='RdBu_r',
            edgecolor='k',
            linewidth=0.5,
            vmin=-6,
            vmax=6
        )
        cbar = fig.colorbar(h, shrink=0.6)
        cbar.set_label('logFC (GIC - iNSC)')
        ax.set_xlim([-.1, len(pids) + .1])
        ax.set_ylim([-.1, this_logfc.shape[0] + .1])
        ax.set_xticks(np.arange(len(pids)) + 0.5)
        ax.set_xticklabels(pids, rotation=90)
        ax.set_yticks(np.arange(this_logfc.shape[0]) + 0.5)
        ax.set_yticklabels(this_logfc.index[::-1], rotation=0)

        if this_cpm.shape[0] > 50:
            plt.setp(ax.yaxis.get_ticklabels(), fontsize=6)

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "%s_logfc_heatmap.png" % name), dpi=200)

    # Plot bar chart showing z values (where available) across different pathways (and comparisons?)
    # start with S1 data (no references)
    pois = [
        'Chondroitin Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis (Late Stages)',
        'Chondroitin Sulfate Biosynthesis (Late Stages)',
        'mTOR Signaling',
        'EIF2 Signaling',
        'G-Protein Coupled Receptor Signaling',
        'cAMP-mediated signaling'
    ]

    tmp = {}
    i = 0
    for pid in pids:
        for pw in pois:
            tmp[i] = {'patient_id': pid, 'pathway': pw, 'z': ipa_de_res[pid].reindex([pw]).squeeze()['z']}
            i += 1
    df = pd.DataFrame(tmp).transpose()
    pw_to_keep = df.groupby('pathway').apply(lambda x: ~x.isnull().all()).loc[pois, 'z']
    pw_to_keep = pw_to_keep.index[pw_to_keep]

    fig, axs = plt.subplots(nrows=len(pw_to_keep), sharex=True, sharey=True, figsize=(5.5, 6.5))
    big_ax = common.add_big_ax_to_subplot_fig(fig)
    for i, pw in enumerate(pw_to_keep):
        ax = axs[i]
        ax.bar(
            range(len(pids)),
            df.loc[df.pathway == pw, 'z'],
            color=[patient_colours[pid] for pid in pids],
            edgecolor='k',
            linewidth=0.8
        )
        ax.set_title(pw)
        ax.set_ylim([-6, 6])
        ax.axhline(0, c='k', linestyle='--', alpha=0.6, linewidth=0.8)
        ax.xaxis.set_ticks(range(len(pids)))

    axs[-1].xaxis.set_ticklabels(pids, rotation=90)
    axs[-1].set_xlim([-.5, len(pids) - .5])
    # add manual legend
    legend_dict = collections.OrderedDict([
        (pid, {'class': 'patch', 'facecolor': patient_colours[pid], 'edgecolor': 'k', 'linewidth': 0.5}) for pid in pids
    ])
    common.add_custom_legend(axs[0], {'': legend_dict}, loc_outside=True, loc_outside_vert='top')
    big_ax.set_ylabel('IPA inferred z value')
    fig.subplots_adjust(bottom=0.06, left=0.1, right=0.8, top=0.95, hspace=0.3)
    fig.savefig(os.path.join(outdir, "selected_pathways_s1_ipa_z_values.png"), dpi=200)



import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from utils import output, ipa, setops
from settings import GIT_LFS_DATA_DIR
from scripts.hgic_final import ipa_results_s1_s2 as irss
from scripts.hgic_final import consts


if __name__ == '__main__':
    # set a minimum pval for pathways to be used
    alpha = 0.005
    plogalpha = -np.log10(alpha)
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    pids = consts.PIDS

    outdir = output.unique_output_dir()

    indir = os.path.join(GIT_LFS_DATA_DIR, 'ipa_from_biplots')

    # single components
    for q in [99, 995]:
        ipa_pathways_single = {}
        for dim in range(1, 4):
            fn = os.path.join(indir, "%d_%d.txt" % (dim, q))
            this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
            this.columns = ['-logp', 'ratio', 'z', 'genes']
            # add ngenes column
            this.insert(3, 'n_gene', this.genes.str.split(',').apply(len))
            this.index = [x.decode('utf-8') for x in this.index]
            ipa_pathways_single[dim] = this

        ipa_single = pd.DataFrame(index=sorted(setops.reduce_union(*[ipa_pathways_single[i].index for i in range(1, 4)])))
        [ipa_single.insert(0, i, ipa_pathways_single[i]['-logp']) for i in range(1, 4)[::-1]]
        ipa_single.fillna(0., inplace=True)
        # drop rows with no significant results
        ipa_single = ipa_single.loc[(ipa_single > -np.log10(0.05)).sum(axis=1) > 0]
        p_order = ipa_single.sum(axis=1).sort_values(ascending=False).index
        fig = plt.figure(figsize=(7., 9.8))
        ax = fig.add_subplot(111)
        sns.heatmap(
            ipa_single.loc[p_order],
            mask=ipa_single.loc[p_order] == 0,
            cmap='YlOrRd',
            linewidths=.2,
            linecolor='w',
            cbar_kws={"orientation": 'vertical', "shrink": 0.6},
            ax=ax
        )
        plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90, fontsize=10)
        fig.subplots_adjust(left=0.6, right=0.95, bottom=0.05, top=.99)
        fig.savefig(os.path.join(outdir, "single_pc_ipa_pval_heatmap_q%d.png" % q), dpi=200)
        fig.savefig(os.path.join(outdir, "single_pc_ipa_pval_heatmap_q%d.pdf" % q), dpi=200)
        fig.savefig(os.path.join(outdir, "single_pc_ipa_pval_heatmap_q%d.tiff" % q), dpi=200)

    # PC (1, 2) [mean logFC]
    ipa_pc1_2_mean = {}
    for q in ['0.99', '0.995']:
        fn = os.path.join(indir, 'pc_1_2_q%s.txt' % q)
        this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
        this.columns = ['-logp', 'ratio', 'z', 'genes']
        # add ngenes column
        this.insert(3, 'n_gene', this.genes.str.split(',').apply(len))
        this.index = [x.decode('utf-8') for x in this.index]

        ipa_pc1_2_mean[q] = this

    ipa_pc1_2 = pd.DataFrame(index=sorted(setops.reduce_union(*[v.index for v in ipa_pc1_2_mean.values()])))
    [ipa_pc1_2.insert(0, q, ipa_pc1_2_mean[q]['-logp']) for q in ['0.99', '0.995']]
    ipa_pc1_2.fillna(0., inplace=True)
    ipa_pc1_2 = ipa_pc1_2[['0.99', '0.995']]
    # drop rows with no significant results
    ipa_pc1_2 = ipa_pc1_2.loc[(ipa_pc1_2> -np.log10(0.05)).sum(axis=1) > 0]


    p_order = ipa_pc1_2.sum(axis=1).sort_values(ascending=False).index
    fig = plt.figure(figsize=(7., 9.8))
    ax = fig.add_subplot(111)
    sns.heatmap(
        ipa_pc1_2.loc[p_order],
        mask=ipa_pc1_2.loc[p_order] == 0,
        cmap='YlOrRd',
        linewidths=.2,
        linecolor='w',
        cbar_kws={"orientation": 'vertical', "shrink": 0.6},
        ax=ax
    )
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90, fontsize=10)
    fig.subplots_adjust(left=0.6, right=0.95, bottom=0.05, top=.99)
    fig.savefig(os.path.join(outdir, "pc1-2_ipa_pval_heatmap.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "pc1-2_ipa_pval_heatmap.pdf"), dpi=200)
    fig.savefig(os.path.join(outdir, "pc1-2_ipa_pval_heatmap.tiff"), dpi=200)

    ipa_res = {}
    ipa_res_wide = {}

    for first_dim in [2, 3]:

        this_ipa_res = ipa.load_raw_reports(indir, "pc_%d_%d_{0}_q0.99.txt" % (first_dim, first_dim + 1), pids)
        for k, v in this_ipa_res.items():
            rele_ix = v.index[v['-logp'] >= plogalpha_relevant]
            this_ipa_res[k] = v.loc[rele_ix]
        ipa_res[(first_dim, first_dim + 1)] = this_ipa_res

        # wideform version of this (i.e. 30 blocks)
        ipa_res_wide[(first_dim, first_dim + 1)] = irss.ipa_results_to_wideform(this_ipa_res, plogalpha)

        dd = irss.generate_plotting_structures(this_ipa_res, ipa_res_wide[(first_dim, first_dim + 1)].index, plogalpha_relevant)
        all_p = dd['p']
        all_z = dd['z']
        all_in = dd['in']

        p_order = all_p.sum(axis=1).sort_values(ascending=False).index

        fig = plt.figure(figsize=(7., 9.8))
        ax = fig.add_subplot(111)
        sns.heatmap(
            all_p.loc[p_order, pids],
            mask=all_p.loc[p_order, pids].isnull(),
            cmap='YlOrRd',
            linewidths=.2,
            linecolor='w',
            cbar_kws={"orientation": 'vertical', "shrink": 0.6},
            ax=ax
        )
        plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90, fontsize=10)
        fig.subplots_adjust(left=0.6, right=0.95, bottom=0.05, top=.99)
        fig.savefig(os.path.join(outdir, "pc_%d_%d_ipa_pval_heatmap.png" % (first_dim, first_dim + 1)), dpi=200)
        fig.savefig(os.path.join(outdir, "pc_%d_%d_ipa_pval_heatmap.pdf" % (first_dim, first_dim + 1)), dpi=200)
        fig.savefig(os.path.join(outdir, "pc_%d_%d_ipa_pval_heatmap.tiff" % (first_dim, first_dim + 1)), dpi=200)


import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from utils import output, ipa
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

    de_comps = {
        'syngeneic': 'syngeneic',
        'h9': 'H9',
        'gibco': 'GIBCO'
    }

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

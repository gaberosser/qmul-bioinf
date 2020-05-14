"""
Key idea: We observe a large number of patient-specific DEs / DMRs / DE+DMRs. Is this to be expected from a random
selection model?
"""
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from scripts.hgic_final import consts
from settings import HGIC_LOCAL_DIR

from utils import log, output, setops
from plotting import common

def run_patient_specific_permutations(n_tot, n_all, n_perm=1000):
    # perms
    n_spec = {pid: [] for pid in n_tot}

    for i in range(n_perm):
        this = []
        for pid in n_tot:
            ix = np.random.permutation(n_all)[:n_tot[pid]]
            this.append(ix)
        this_spec = setops.specific_features(*this)
        for pid, s in zip(n_tot.keys(), this_spec):
            n_spec[pid].append(len(s))

    return n_spec


def plot_perms_kde_vs_obs(n_spec_perm, n_spec_obs, order=None, colours=consts.PATIENT_COLOURS, xlabel=None):
    if order is None:
        order = n_spec_perm.keys()
    with sns.axes_style('whitegrid'):
        fig, axs = plt.subplots(nrows=len(n_spec_perm), sharex=True, sharey=True, figsize=(5, 5.5))
        # big_ax = common.add_big_ax_to_subplot_fig(fig)
        for i, pid in enumerate(order):
            ax = axs[i]
            sns.kdeplot(np.array(n_spec_perm[pid]), shade=True, color=colours[pid], ax=ax, alpha=0.7)
            ax.yaxis.set_ticks([])
            ax.axvline(n_spec_obs[pid], c='k', lw=1.5)
            ax.set_ylabel(pid)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        fig.tight_layout()

    return fig, axs


if __name__ == '__main__':
    n_perm = 1000

    # DE

    pids = consts.PIDS
    outdir = output.unique_output_dir()
    # load previously generated DE results
    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq',
        'full_de_syngeneic_only.xlsx'
    )
    de_res = pd.read_excel(fn, header=0, index_col=0)
    all_ens = de_res.index[(de_res[pids] == 'Y').any(axis=1)]

    de_per_pat = {pid: de_res.index[de_res[pid] == 'Y'] for pid in pids}
    n_tot = {pid: de_per_pat[pid].size for pid in pids}

    vs, vc = setops.venn_from_arrays(*[de_per_pat[pid] for pid in pids])
    pp = setops.specific_sets(pids)
    n_ps = {pid: vc[pp[pid]] for pid in pids}

    # perms
    n_all = len(all_ens)
    n_spec = run_patient_specific_permutations(n_tot, n_all, n_perm=n_perm)

    fig, axs = plot_perms_kde_vs_obs(n_spec, n_ps, xlabel='Number of patient-specific DE genes', order=pids)
    fig.savefig(os.path.join(outdir, "patient_specific_de.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_de.tiff"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_de.pdf"))

    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'methylation',
        'full_dmr_syngeneic_only.xlsx'
    )
    dm_res = pd.read_excel(fn)
    # total number of (DM cluster, gene) pairs
    # n_all_de_dm = dm_res.shape[0]
    # reported elsewhere, so use that here (it's approximate)
    n_all_de_dm = 17000

    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq_methylation_combined',
        'de_dmr_concordant_syngeneic_only.xlsx'
    )
    de_dm_res = pd.read_excel(fn)
    # number of patient specific DE/DMRs
    de_dm_per_pat = {pid: de_dm_res.index[de_dm_res[pid] == 'Y'] for pid in pids}
    n_tot_de_dm = {pid: de_dm_per_pat[pid].size for pid in pids}
    tt = setops.specific_features(*[de_dm_per_pat[pid] for pid in pids])
    n_ps_de_dmr = {pid: len(t) for pid, t in zip(pids, tt)}

    n_spec_perm_de_dm = run_patient_specific_permutations(n_tot_de_dm, n_all_de_dm, n_perm=n_perm)

    fig, axs = plot_perms_kde_vs_obs(n_spec_perm_de_dm, n_ps_de_dmr, xlabel='Number of patient-specific DE/DMRs', order=pids)
    fig.savefig(os.path.join(outdir, "patient_specific_de_dmr.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_de_dmr.tiff"), dpi=200)
    fig.savefig(os.path.join(outdir, "patient_specific_de_dmr.pdf"))
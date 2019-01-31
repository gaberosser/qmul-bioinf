import pandas as pd
import numpy as np
import os
import multiprocessing as mp
from bisect import bisect_left
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from methylation import dmr, loader
from utils import setops, output, log
from plotting import common, venn
from matplotlib import pyplot as plt
import seaborn as sns
logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    n_iter = 1000

    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tsgd.dmr_results_hash(consts.S1_METHYL_SAMPLES, dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        anno = loader.load_illumina_methylationepic_annotation(split_genes=False)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise IOError("We require a pre-computed file, %s, which could not be found." % fn)

    # extract full (all significant) results
    dmr_res_all = dmr_res_s1.results_significant

    clusters = dmr_res_s1[pids[0]].clusters

    n_by_patient = dict([
        (pid, len(dmr_res_all[pid])) for pid in pids
    ])

    specific_dmrs = dict(zip(pids, setops.specific_features(*[dmr_res_all[pid] for pid in pids])))
    n_by_patient_specific = dict([
        (pid, len(specific_dmrs[pid])) for pid in pids
    ])

    ntot = sum(n_by_patient.values())

    # 1) Null: DMRs are picked uniformly randomly from the pool with variable marginal totals for each patient.
    # Marginal totals are given by the (real) number of DMRs in each patient.
    rvs = dict([
        (
            pid,
            [np.random.choice(range(ntot), replace=False, size=n_by_patient[pid]) for i in range(n_iter)]
        ) for pid in pids
    ])

    inters_1 = [[len(x) for x in setops.specific_features(*[rvs[pid][i] for pid in pids])] for i in range(n_iter)]
    inters_1 = dict(zip(pids, zip(*inters_1)))

    clist = common.get_best_cmap(len(pids))
    patient_colours = dict(zip(pids, clist))

    fig, axs = plt.subplots(nrows=len(pids), sharex=True, sharey=True)
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(top='off', bottom='off', left='off', right='off', labelcolor='none')
    big_ax.grid(False)
    big_ax.set_ylabel('Density (a.u.)')

    for pid, ax in zip(pids, axs):
        sns.kdeplot(np.array(inters_1[pid]), color=patient_colours[pid], shade=True, ax=ax)
        ax.axvline(n_by_patient_specific[pid], color='k')
        ax.yaxis.set_ticks([])
        ax.set_ylabel(pid)

    axs[-1].set_xlim([0, np.array(inters_1.values()).max() * 1.1])
    axs[-1].set_xlabel('Number of DMRs')
    fig.subplots_adjust(hspace=0.05, left=0.1, right=0.98, bottom=0.1, top=0.98)

    fig.savefig(os.path.join(outdir, "null_one_number_specific_dmrs.png"), dpi=200)

    # 2) Null: DMRs are picked with a propensity given by the number of patients who share them and with marginal totals
    # as before.
    # This means that only the DMRs seen in our patients have a non-zero weighting??
    weights = np.array(setops.feature_membership_count(*[dmr_res_all[pid] for pid in pids]).values())
    weights = weights / float(weights.sum())

    ntot = len(weights)
    rvs = dict([
        (
            pid,
            np.array([np.random.choice(range(ntot), replace=False, size=n_by_patient[pid], p=weights) for i in range(n_iter)])
        ) for pid in pids
    ])

    inters_2 = [[len(x) for x in setops.specific_features(*[rvs[pid][i] for pid in pids])] for i in range(n_iter)]
    inters_2 = dict(zip(pids, zip(*inters_2)))

    fig, axs = plt.subplots(nrows=len(pids), sharex=True, sharey=True)
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(top='off', bottom='off', left='off', right='off', labelcolor='none')
    big_ax.grid(False)
    big_ax.set_ylabel('Density (a.u.)')
    xmax = 500.

    for pid, ax in zip(pids, axs):
        sns.kdeplot(np.array(inters_2[pid]), color=patient_colours[pid], shade=True, ax=ax)
        ax.axvline(n_by_patient_specific[pid], color='k')
        ax.yaxis.set_ticks([])
        ax.set_ylabel(pid)
        # aax = ax.twinx()
        # aax.grid(False)
        # aax.set_yticks([])
        ix = bisect_left(sorted(inters_2[pid]), n_by_patient_specific[pid])
        # two tailed test (?) - need to check for direction
        pval = min(ix / float(n_iter), 1 - ix / float(n_iter))
        if pval == 0:
            lbl = 'p < %.3f' % (1 / float(n_iter))
        else:
            lbl = 'p = %.3f' % pval
        ax.text(
            0.99,
            0.5,
            lbl,
            horizontalalignment='right',
            verticalalignment='center',
            transform=ax.transAxes
        )

        # aax.set_ylabel(lbl, rotation=0, horizontalalignment='right')

    axs[-1].set_xlim([0, xmax])
    axs[-1].set_xlabel('Number of DMRs')
    fig.subplots_adjust(hspace=0.05, left=0.1, right=0.98, bottom=0.1, top=0.98)

    fig.savefig(os.path.join(outdir, "null_two_number_specific_dmrs.png"), dpi=200)

    # to illustrate the differences, create an UpSet plot
    subgroup_set_colours = {
        'RTK I full': '#0d680f',
        'RTK II full': '#820505',
        'MES full': '#7900ad',
        'RTK I partial': '#6ecc70',
        'RTK II partial': '#d67373',
        'MES partial': '#cc88ea',
        'Expanded core': '#4C72B0',
        'Specific': '#f4e842',
    }
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

    for_plot = [rvs[pid][0] for pid in pids]
    upset = venn.upset_plot_with_groups(
        for_plot,
        pids,
        subgroup_ind,
        subgroup_set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset['figure'].savefig(os.path.join(outdir, "upset_de.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir, "upset_de.tiff"), dpi=200)

    # quantify this: distribution for each set
    # this is slow when n_iter is large, so we'll apply multiprocessing
    pool = mp.Pool()
    jobs = {}

    simulated_set_sizes = {}
    for i in range(n_iter):
        jobs[i] = pool.apply_async(
            setops.venn_from_arrays,
            args=tuple(rvs[pid][i] for pid in pids)
        )

    pool.close()
    pool.join()

    for i, j in jobs.items():
        _, this_vc = j.get()
        for k, v in this_vc.items():
            simulated_set_sizes.setdefault(k, []).append(v)

    df = pd.DataFrame(simulated_set_sizes).transpose()
    set_size = [sum([int(t) for t in x]) for x in df.index]
    df.insert(0, 'set_size', set_size)
    df.sort_values(by='set_size', inplace=True)

    # plot showing simulated range and our value
    _, vc = setops.venn_from_arrays(*[dmr_res_all[pid] for pid in pids])
    df_min = df.drop('set_size', axis=1).min(axis=1)
    df_range = df.drop('set_size', axis=1).max(axis=1) - df_min

    num = np.array([vc[k] for k in df.index])
    # Z transform (with offset to avoid infty)
    mu = df.drop('set_size', axis=1).mean(axis=1)
    den = df.drop('set_size', axis=1).std(axis=1) + 1e-6
    zz = (num - mu) / den
    # take absolute and sign separately
    logz_abs = np.log10(zz.abs())
    # any entries where observed and simulated are zero should be masked
    logz_abs.loc[np.isinf(logz_abs)] = np.nan
    z_sign = np.sign(zz)

    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111)

    cmap = plt.cm.Vega10.colors
    cols = [cmap[i-1] for i in df.set_size]

    ax.bar(
        range(df.shape[0]),
        df_range,
        width=.9,
        bottom=df_min,
        color=cols,
        edgecolor='none',
        zorder=10,
    )
    ax.scatter(
        range(df.shape[0]),
        num,
        facecolor=cols,
        edgecolor='k',
        linewidth=.5,
        marker='o',
        s=10,
        zorder=9,
    )
    ax.set_xticklabels([])
    ax.set_xlim([0, df.shape[0] + 1])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "num_dmr_in_set_obs_sim.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(
        range(df.shape[0]),
        logz_abs * z_sign,
        facecolor=cols,
        edgecolor='k',
        linewidth=.5,
        marker='o',
        s=10,
        zorder=9
    )
    ax.axhline(0, color='k', linestyle='--', zorder=8)
    ax.set_ylabel(r"$\log_{10}(z)$")
    fig.tight_layout()
    # remove minus signs
    yticklabels = [t.get_text().replace(u'\u2212', '') for t in ax.yaxis.get_ticklabels()]
    ax.yaxis.set_ticklabels(yticklabels)
    ax.xaxis.set_ticklabels([])
    fig.savefig(os.path.join(outdir, "num_dmr_in_set_logz.png"), dpi=200)

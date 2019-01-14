import numpy as np
import os
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from methylation import dmr, loader
from utils import setops, output, log
from plotting import common
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

    for pid, ax in zip(pids, axs):
        sns.kdeplot(np.array(inters_2[pid]), color=patient_colours[pid], shade=True, ax=ax)
        ax.axvline(n_by_patient_specific[pid], color='k')
        ax.yaxis.set_ticks([])
        ax.set_ylabel(pid)
        aax = ax.twinx()
        aax.grid(False)
        aax.set_yticks([])
        aax.set_ylabel('p = ')

    axs[-1].set_xlim([0, 500.])
    axs[-1].set_xlabel('Number of DMRs')
    fig.subplots_adjust(hspace=0.05, left=0.1, right=0.98, bottom=0.1, top=0.98)

    fig.savefig(os.path.join(outdir, "null_two_number_specific_dmrs.png"), dpi=200)






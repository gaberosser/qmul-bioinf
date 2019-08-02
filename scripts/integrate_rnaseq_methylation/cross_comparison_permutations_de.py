import itertools
import multiprocessing as mp
import operator
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from load_data import rnaseq_data
from rnaseq import filter, differential_expression
from utils import output, setops, reference_genomes


def add_gene_symbols(df):
    """
    Add gene symbols to the DataFrame df which is indexed by Ensembl IDs
    """
    gs = reference_genomes.ensembl_to_gene_symbol(df.index)
    # resolve any duplicates arbitrarily (these should be rare)
    gs = gs.loc[~gs.index.duplicated()]
    df.insert(0, 'Gene Symbol', gs)


def add_fc_direction(df):
    direction = pd.Series(index=df.index, name='Direction')
    direction.loc[df.logFC < 0] = 'down'
    direction.loc[df.logFC > 0] = 'up'
    df.insert(df.shape[1], 'Direction', direction)


def run_one_de(the_data, the_groups, the_comparison, lfc=1, fdr=0.01, method='QLGLM'):
    if method == 'QLGLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmqlfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
    elif method == 'GLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
    elif method == 'exact':
        res = differential_expression.edger_exacttest(the_data, the_groups, pair=the_comparison[::-1], lfc=lfc, fdr=fdr)

    add_gene_symbols(res)
    add_fc_direction(res)

    return res


def compute_cross_de(rnaseq_obj, main_pids, additional_pids, external_references=(('GIBCO', 'NSC'),), lfc=1, fdr=0.01, method='QLGLM'):
    """
    Compute DE between every patient GBM sample and every _other_ healthy patient sample, in addition to paired DE.
    We can also include one or more external references (e.g. Gibco, the default).
    :param rnaseq_obj:
    :param main_pids: The PIDs of the patients in the study
    :param additional_pids: The PIDs of the additional iNSC lines to include. These will be treated as additional
    reference lines.
    :param external_references:
    :param lfc:
    :param fdr:
    :param method:
    :param njob:
    :return:
    """
    if method not in {'QLGLM', 'GLM', 'exact'}:
        raise NotImplementedError("Unsupported method.")
    de = {}

    for pid in main_pids:

        # cross comparison
        for pid2 in list(main_pids) + list(additional_pids):
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(pid2) & (rnaseq_obj.meta.loc[:, 'type'] == 'iNSC'))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_data = filter.filter_by_cpm(the_data, min_n_samples=1)
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', 'iNSC']

            de[(pid, pid2)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)

        # external reference comparison
        for er, er_type in external_references:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(er) & (rnaseq_obj.meta.loc[:, 'type'] == er_type))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_data = filter.filter_by_cpm(the_data, min_n_samples=1)
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', er_type]
            de[(pid, er)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)


    return de


def ncr(n, r):
    r = min(r, n - r)
    if r == 0: return 1
    numer = reduce(operator.mul, xrange(n, n - r, -1))
    denom = reduce(operator.mul, xrange(1, r + 1))
    return numer // denom


if __name__ == "__main__":
    # if this is specified, we load the DMR results from a JSON rather than recomputing them to save time
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'integrate_rnaseq_methylation')

    outdir = output.unique_output_dir("cross_validate_de_permutations", reuse_empty=True)
    ref_name = 'GIBCONSC_P4'
    # all n=2 samples and RTK II samples
    pids = ['017', '019', '030', '031', '050', '054']
    additional_pids = ['018', '049', '061', '026', '044', '052']
    cmap = 'RdYlGn_r'

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    intersecter = lambda x, y: set(x).intersection(y)
    unioner = lambda x, y: set(x).union(y)

    # Load RNA-Seq from STAR
    rnaseq_obj = rnaseq_data.load_by_patient(pids + additional_pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('IPSC')]
    rnaseq_obj.data = rnaseq_obj.data.loc[:, rnaseq_obj.meta.index]

    de_res = compute_cross_de(rnaseq_obj, pids, additional_pids, **de_params)

    # counts of DE genes
    de_counts = pd.DataFrame(index=pids, columns=pids + additional_pids + ['GIBCO'])
    for pid in pids:
        for pid2 in pids + additional_pids + ['GIBCO']:
            de_counts.loc[pid, pid2] = de_res[(pid, pid2)].shape[0]

    # now we need to compare the paired results with every other result (Gibco and other iNSC)
    pair_only = pd.DataFrame(index=pids, columns=pids + additional_pids + ['GIBCO'])
    ref_only = pd.DataFrame(index=pids, columns=pids + additional_pids + ['GIBCO'])
    pair_and_ref_concordant = pd.DataFrame(index=pids, columns=pids + additional_pids + ['GIBCO'])
    pair_and_ref_discordant = pd.DataFrame(index=pids, columns=pids + additional_pids + ['GIBCO'])

    # loop over GBM samples
    for pid in pids:
        # syngeneic comparison
        the_pair = de_res[(pid, pid)]

        # loop over (i)NSC samples
        # when this is the same as the syngeneic comparison, there will (obviously) be no 'pair only' or 'ref only'
        # genes!
        for pid2 in pids +additional_pids + ['GIBCO']:
            the_ref = de_res[(pid, pid2)]
            the_sets, _ = setops.venn_from_arrays(the_pair.index, the_ref.index)
            pair_only.loc[pid, pid2] = the_sets['10']
            ref_only.loc[pid, pid2] = the_sets['01']
            # for overlapping genes: separate based on direction (matching or non matching)
            the_conc_idx = (the_pair.loc[the_sets['11']].Direction == the_ref.loc[the_sets['11']].Direction)
            pair_and_ref_concordant.loc[pid, pid2] = the_pair.loc[the_sets['11']].loc[the_conc_idx].index
            pair_and_ref_discordant.loc[pid, pid2] = the_pair.loc[the_sets['11']].loc[~the_conc_idx].index

    # can get counts like this
    po_counts = pair_only.applymap(len)
    ro_counts = ref_only.applymap(len)

    # the permutation part
    # here we find the genes that are core (for a given N) over every possible permutation of references
    po_core_genes = {}

    for pid in pids:
        po_core_genes[pid] = {}
        the_row = pair_only.loc[pid]
        the_refs = the_row.index.difference([pid])
        # all possible combinations of references
        for ref_selection in itertools.combinations(the_refs, len(pids)):
            this_portion = the_row.loc[list(ref_selection)]
            po_core_genes[pid][ref_selection] = pd.Series(index=range(1, len(pids) + 1))
            gl, cts = setops.venn_from_arrays(*this_portion.values)
            for N in range(1, len(pids) + 1):
                po_core_genes[pid][ref_selection].loc[N] = reduce(
                    unioner, (gl[k] for k in setops.binary_combinations_sum_gte(len(pids), N))
                )

    # this strangeness is to remind me that there is currently only 1 ref (but that could change)
    n_perm = ncr(len(pids) + len(additional_pids) + 1 - 1, len(pids))

    # run over results and compute
    # 1) similarity score: the number of genes in each list divided by the number in the union over all permutations
    # 2) union of genes over all permutations
    # 3) intersection of genes over all permutations

    similarities = {}
    isct = {}
    unn = {}
    for N in range(1, len(pids) + 1):
        similarities[N] = pd.DataFrame(index=pids, columns=range(n_perm))
        isct[N] = pd.Series(index=pids)
        unn[N] = pd.Series(index=pids)
        for pid in pids:
            unn[N].loc[pid] = reduce(
                unioner, (t.loc[N] for t in po_core_genes[pid].values())
            )
            isct[N].loc[pid] = reduce(
                intersecter, (t.loc[N] for t in po_core_genes[pid].values())
            )
            norm = float(len(unn[N].loc[pid]))
            for i, t in enumerate(po_core_genes[pid].values()):
                similarities[N].loc[pid, i] = len(t.loc[N]) / norm

    isct_over_unn = pd.DataFrame(
        [isct[i].apply(len) / unn[i].apply(len) for i in range(1, len(pids) + 1)],
        index=range(1, len(pids) + 1)
    )
    ax = isct_over_unn.plot.line()
    ax.set_ylabel('Consistent / union')
    ax.set_xlabel('N')
    ax.figure.savefig(os.path.join(outdir, "intersection_over_union.png"), dpi=200)

    # compute the pairwise overlap for every possible pair
    pool = mp.Pool()
    jobs = {}

    for N in range(1, len(pids) + 1):
        for pid in pids:
            arr = [t.loc[N] for t in po_core_genes[pid].values()]
            jobs[(N, pid)] = pool.apply_async(setops.pairwise_similarity, args=(arr,), kwds={'method': 'union'})

    pool.close()
    pool.join()
    pwise_similarities = {}
    for N in range(1, len(pids) + 1):
        pwise_similarities[N] = {}
        for pid in pids:
            pwise_similarities[N][pid] = jobs[(N, pid)].get(1e6)

    # plot overall similarity
    # per PID
    sim_bins = np.linspace(0, 1, 200)
    for pid in pids:
        fig, axs = plt.subplots(nrows=len(pids), num=pid, figsize=(6, 7.5))
        big_ax = fig.add_subplot(111, frameon=False)
        big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        big_ax.grid(False)
        big_ax.set_ylabel("Frequency")
        big_ax.set_xlabel("Similarity score")
        for i, ax in enumerate(axs):
            ax.hist(similarities[i + 1].loc[pid].values, bins=sim_bins)
            if i != (len(pids) - 1):
                ax.xaxis.set_ticklabels([])
            ax.text(0.05, 0.95, "N=%d" % (i + 1), transform=ax.transAxes, fontsize=14, va='top')

        fig.tight_layout(pad=0.1, h_pad=0.1)
        fig.savefig(os.path.join(outdir, "similarity_score_%s.png" % pid), dpi=200)

    fig, axs = plt.subplots(nrows=len(pids), figsize=(6, 7.5))
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)
    big_ax.set_ylabel("Frequency")
    big_ax.set_xlabel("Overall similarity score")
    for i, ax in enumerate(axs):
        ax.hist(similarities[i + 1].values.flatten(), bins=sim_bins)
        if i != (len(pids) - 1):
            ax.xaxis.set_ticklabels([])
        ax.text(0.05, 0.95, "N=%d" % (i + 1), transform=ax.transAxes, fontsize=14, va='top')

    fig.tight_layout(pad=0.1, h_pad=0.1)
    fig.savefig(os.path.join(outdir, "combined_overall_similarity_score.png"), dpi=200)


    # plot pairwise similarity
    for pid in pids:
        fig, axs = plt.subplots(nrows=len(pids), num=pid, figsize=(6, 7.5))
        big_ax = fig.add_subplot(111, frameon=False)
        big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        big_ax.grid(False)
        big_ax.set_ylabel("Frequency")
        big_ax.set_xlabel("Pairwise similarity score")
        for i, ax in enumerate(axs):
            ax.hist(pwise_similarities[i + 1][pid], bins=sim_bins)
            if i != (len(pids) - 1):
                ax.xaxis.set_ticklabels([])
            ax.text(0.05, 0.95, "N=%d" % (i + 1), transform=ax.transAxes, fontsize=14, va='top')

        fig.tight_layout(pad=0.1, h_pad=0.1)
        fig.savefig(os.path.join(outdir, "pwise_similarity_score_%s.png" % pid), dpi=200)

    # combine PIDs to generate one giant list for each N
    combined_pwise_sim = dict([
        (
            N,
            reduce(lambda x, y: x + y, pwise_similarities[N].values())
        ) for N in range(1, len(pids) + 1)
    ])

    # plot pairwise similarity
    fig, axs = plt.subplots(nrows=len(pids), figsize=(6, 7.5))
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)
    big_ax.set_ylabel("Frequency")
    big_ax.set_xlabel("Pairwise similarity score")
    for i, ax in enumerate(axs):
        ax.hist(combined_pwise_sim[i + 1], bins=sim_bins)
        if i != (len(pids) - 1):
            ax.xaxis.set_ticklabels([])
        ax.text(0.05, 0.95, "N=%d" % (i + 1), transform=ax.transAxes, fontsize=14, va='top')

    fig.tight_layout(pad=0.1, h_pad=0.1)
    fig.savefig(os.path.join(outdir, "combined_pwise_similarity_score.png"), dpi=200)
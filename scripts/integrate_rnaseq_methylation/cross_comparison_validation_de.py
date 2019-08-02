import collections
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from load_data import rnaseq_data
from rnaseq import differential_expression
from settings import LOCAL_DATA_DIR
from utils import output, setops, excel, ipa, reference_genomes


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


def compute_cross_de(
        rnaseq_obj,
        pids,
        external_references=(('GIBCO', 'NSC'),),
        lfc=1,
        fdr=0.01,
        method='QLGLM',
        njob=None,
):
    """
    Compute DE between every patient GBM sample and every _other_ healthy patient sample, in addition to paired DE.
    We can also include one or more external references (e.g. Gibco, the default).
    :param rnaseq_obj:
    :param pids:
    :param external_references:
    :param lfc:
    :param fdr:
    :param method:
    :param njob: Number of multiprocessing workers to use (default: None - automatically selected)
    :return:
    """
    if method not in {'QLGLM', 'GLM', 'exact'}:
        raise NotImplementedError("Unsupported method.")
    de = {}

    pool = mp.Pool()
    jobs = {}

    for pid in pids:

        # cross comparison
        for pid2 in pids:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(pid2) & (rnaseq_obj.meta.loc[:, 'type'] == 'iNSC'))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', 'iNSC']
            if njob == 1:
                de[(pid, pid2)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)
            else:
                jobs[(pid, pid2)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method)
                )

        # external reference comparison
        for er, er_type in external_references:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(er) & (rnaseq_obj.meta.loc[:, 'type'] == er_type))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', er_type]
            if njob == 1:
                de[(pid, er)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)
            else:
                jobs[(pid, er)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method)
                )

    if njob != 1:
        pool.close()
        pool.join()
        for k, j in jobs.items():
            de[k] = j.get(1e6)

    return de


if __name__ == "__main__":

    outdir = output.unique_output_dir("cross_validate_de", reuse_empty=True)
    # all n=2 samples and RTK II samples
    pids = ['017', '019', '030', '031', '050', '054']
    cmap = 'RdYlGn_r'

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    subgroups = {
        'RTK I': ['019', '030', '031'],
        'RTK II': ['017', '050', '054'],
    }

    intersecter = lambda x, y: set(x).intersection(y)
    unioner = lambda x, y: set(x).union(y)

    # Load RNA-Seq from STAR
    rnaseq_obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')

    # load additional references if required
    h9_obj = rnaseq_data.gse61794(annotate_by='Ensembl Gene ID')
    h1_obj = rnaseq_data.gse38993(annotate_by='Ensembl Gene ID')
    rnaseq_obj = rnaseq_data.MultipleBatchLoader([rnaseq_obj, h1_obj, h9_obj])

    # discard unmapped, etc
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('PSC')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('fibroblast')]
    rnaseq_obj.data = rnaseq_obj.data.loc[:, rnaseq_obj.meta.index]

    # load RNA-Seq from Salmon (for normalised comparison)
    # disabled for now
    if False:
        salmon_dat = rnaseq_data.load_salmon_by_patient_id(pids)
        idx = salmon_dat.index.str.replace(r'.[0-9]+$', '')
        salmon_dat.index = idx
        fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90',
                          'gene_to_transcript.txt')
        gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

        # aggregate to gene level
        genes = gene_transcript.loc[salmon_dat.index, 'Gene stable ID']
        salmon_dat = salmon_dat.groupby(genes).sum()

        # discard unmapped, etc
        salmon_dat = salmon_dat.loc[:, ~salmon_dat.columns.str.contains('IPSC')]

    external_refs = [
        ('GIBCO', 'NSC'),
        # ('H9', 'NSC'),
        # ('H1', 'NSC'),
    ]
    external_ref_labels = [t[0] for t in external_refs]
    ref_samples = reduce(
        lambda x, y: x+y,
        [rnaseq_obj.meta.index[rnaseq_obj.meta.index.str.contains(t[0])].tolist() for t in external_refs]
    )
    cols = pids + [t[0] for t in external_refs]
    de_res = compute_cross_de(rnaseq_obj, pids, external_references=external_refs, **de_params)

    # counts of DE genes
    de_counts = pd.DataFrame(index=pids, columns=cols)
    for pid in pids:
        for pid2 in cols:
            de_counts.loc[pid, pid2] = de_res[(pid, pid2)].shape[0]

    # now we need to compare the paired results with every other result (Gibco and other iNSC)
    pair_only = pd.DataFrame(index=pids, columns=cols)
    ref_only = pd.DataFrame(index=pids, columns=cols)
    pair_and_ref_concordant = pd.DataFrame(index=pids, columns=cols)
    pair_and_ref_discordant = pd.DataFrame(index=pids, columns=cols)
    # loop over GBM samples
    for pid in pids:
        # syngeneic comparison
        the_pair = de_res[(pid, pid)]

        # loop over (i)NSC samples
        # when this is the same as the syngeneic comparison, there will (obviously) be no 'pair only' or 'ref only'
        # genes!
        for pid2 in cols:
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

    ## genes that are pair-only in every possible ref comparison
    po_each = [
        sorted(
            reduce(intersecter, pair_only.loc[pid, ~pair_only.columns.str.contains(pid)])
        ) for pid in pids
    ]
    po_each = pd.Series(po_each, index=pids)

    # export gene lists here
    po_export = {}
    for pid in pids:
        po_export["GBM%s_pair_only" % pid] = de_res[(pid, pid)].loc[po_each.loc[pid]]
    excel.pandas_to_excel(po_export, os.path.join(outdir, "pair_only_all_consistent.xlsx"))
    subdir = os.path.join(outdir, "ipa_all_consistent")
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    ipa.results_to_ipa_format(po_export, outdir=subdir)

    # now relax this requirement: which genes would be included if we require their inclusion in N of the cells
    # (rather than all)?
    possible_counts = range(1, pair_only.shape[1])
    po_each_threshold = pd.DataFrame(index=pids, columns=possible_counts)
    for pid in pids:
        this_counter = collections.Counter()
        # iterate over each column
        # we can include the empty diagonal cell, since it will not affect the counting
        for col in pair_only.columns:
            for e in pair_only.loc[pid, col]:
                this_counter[e] += 1
        # progressively filter the gene list based on counts
        the_genes = this_counter.keys()
        for i in possible_counts:
            the_genes = [k for k in this_counter if this_counter[k] >= i]
            po_each_threshold.loc[pid, i] = the_genes

    # export gene lists here - only consider N-1 and N-2
    for i in [1, 2]:
        j = len(pids) - i
        the_list = po_each_threshold.loc[:, j]
        po_export = {}
        for pid in pids:
            po_export["GBM%s_pair_only" % pid] = de_res[(pid, pid)].loc[the_list.loc[pid]]
        excel.pandas_to_excel(po_export, os.path.join(outdir, "pair_only_%d_consistent.xlsx" % j))
        subdir = os.path.join(outdir, "ipa_%d_consistent" % j)
        if not os.path.isdir(subdir):
            os.makedirs(subdir)
        ipa.results_to_ipa_format(po_export, outdir=subdir)

    # ...how many of these are shared between patients?
    # consider K to K-3 (inclusive)
    for i in possible_counts:
        _, cts = setops.venn_from_arrays(*po_each_threshold.loc[:, i].values)
        this_tally = []
        K = len(pids)
        print "N = %d" % i
        for j in [K, K-1, K-2, K-3]:
            this_ct = sum([cts[k] for k in setops.binary_combinations_sum_gte(K, j)])
            print "%d genes shared by >=%d patients" % (this_ct, j)

        # also look at the overlap within the subgroups
        for grp_name, grp_members in subgroups.items():
            # get the group member results
            this_po_each_threshold = po_each_threshold.loc[grp_members]
            _, cts = setops.venn_from_arrays(*this_po_each_threshold.loc[:, i].values)
            the_idx = ''.join(['1'] * len(grp_members))
            print "%d DE genes shared by all patients in subgroup %s" % (cts[the_idx], grp_name)

    # for reference: what do these numbers look like in the Gibco comparison (only)?
    po_gibco_common_counts = pd.Series(index=possible_counts, dtype=int)
    _, cts = setops.venn_from_arrays(*pair_only.loc[:, 'GIBCO'].values)
    for j in possible_counts:
        po_gibco_common_counts.loc[j] = sum([cts[k] for k in setops.binary_combinations_sum_gte(K, j)])
        print "%d DE genes shared by >=%d patients in the pair-only Gibco comparison" % (
            int(po_gibco_common_counts.loc[j]),
            j
        )

    # What is present in X vs Y_i that isn't in X vs any other Y?
    po_diff = pd.DataFrame(index=pair_only.index, columns=pair_only.columns)
    for pid in pids:
        for pid2 in pair_only.columns:
            the_ref = pair_only.loc[pid, pid2]
            all_else = pair_only.loc[pid, pair_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            po_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    # find DE genes that are always PO when a (non-matching) iNSC reference is used, but NOT when an external reference
    # is used.
    po_intersection_insc = pd.DataFrame(index=pids, columns=external_ref_labels + ['any'])
    for pid in pids:
        # first, find the genes that are always PO when an iNSC reference is used
        tmp = reduce(intersecter, pair_only.loc[pid, pair_only.index[pair_only.index != pid]])
        for er in external_ref_labels:
            # second, find the union of genes that are PO when this external reference is used
            tmp2 = pair_only.loc[pid, er]
            # we want anything in the first part that is NOT in the second part
            po_intersection_insc.loc[pid, er] = tmp.difference(tmp2)
        # now, find the union of genes that are PO when ANY of the external references is used
        tmp2 = reduce(unioner, pair_only.loc[pid, external_ref_labels])
        po_intersection_insc.loc[pid, 'any'] = tmp.difference(tmp2)

    # find DE genes
    po_specific_to_reference = [
        sorted(
            reduce(intersecter, po_diff.loc[~po_diff.index.str.contains(pid), pid])
        ) for pid in cols
    ]
    po_specific_to_reference = pd.Series(po_specific_to_reference, index=cols)

    # get the genes that consistently differ in the pair comparison only and NOT in Gibco (across all patients)
    # these will have an expression pattern in Gibco similar to GBM, so that they do NOT appear
    po_gibco_diff = po_specific_to_reference.loc['GIBCO']
    po_gibco_diff_gs = reference_genomes.ensembl_to_gene_symbol(po_gibco_diff)
    po_gibco_diff_gs = po_gibco_diff_gs.where(~po_gibco_diff_gs.isnull(), po_gibco_diff)

    po_dat = rnaseq_obj.data.loc[po_gibco_diff]
    po_dat.index = po_gibco_diff_gs
    po_dat = np.log2(po_dat + 1)

    # po_dat = salmon_dat.loc[po_gibco_diff]
    # po_dat.index = po_gibco_diff_gs
    # # dropna() here loses one gene - LINC01090 / ENSG00000231689
    # # all others are present
    # po_dat = np.log2(po_dat.dropna() + 0.01)

    # rearrange columns
    the_cols = (
        po_dat.columns[po_dat.columns.str.contains('GBM')].tolist() +
        ref_samples +
        po_dat.columns[po_dat.columns.str.contains('DURA')].tolist()
    )
    spacing1 = po_dat.columns.str.contains('GBM').sum()
    spacing2 = spacing1 + len(ref_samples) + 1  # +1 required as we will already have added a space to the left of this
    po_dat = po_dat.loc[:, the_cols]

    # insert spacing columns
    po_dat.insert(spacing1, '', np.nan)
    po_dat.insert(spacing2, ' ', np.nan)

    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(po_dat, cmap=cmap, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "consistently_in_pair_only.png"), dpi=200)

    # export those same genes to a file, adding gene symbols
    for_export = rnaseq_obj.data.loc[po_gibco_diff, the_cols]
    gs = reference_genomes.ensembl_to_gene_symbol(for_export.index)
    for_export.insert(0, 'gene_symbol', gs)
    for_export.to_excel(os.path.join(outdir, 'consistently_in_pair_only.xlsx'))

    ## TODO: fix from here

    # get the genes that consistently appear in the Gibco reference comparison only and NOT in any other reference
    # these will have a different expression pattern in Gibco to GBM (while the level in iNSC will not differ from GBM)
    ro_diff = pd.DataFrame(index=ref_only.index, columns=ref_only.columns)
    for pid in pids:
        for pid2 in ref_only.columns:
            the_ref = ref_only.loc[pid, pid2]
            all_else = ref_only.loc[pid, ref_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            ro_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    # get counts like this
    ro_diff.applymap(len)

    # this computes the set of genes that ALWAYS appears in the ref_only comparison for each possible ref
    ro_each = [
        sorted(
            reduce(lambda x, y: set(x).intersection(y), ro_diff.loc[~ro_diff.index.str.contains(pid), pid])
        ) for pid in cols
    ]
    ro_each = pd.Series(ro_each, index=cols)

    # ro_gibco_diff = sorted(reduce(lambda x, y: set(y).intersection(x), ro_diff.loc[:, 'GIBCO']))
    ro_gibco_diff = ro_each.loc['GIBCO']
    ro_gibco_diff_gs = reference_genomes.ensembl_to_gene_symbol(ro_gibco_diff)
    # the lincRNA symbols are missing, so keep ENSG for those
    ro_gibco_diff_gs = ro_gibco_diff_gs.where(~ro_gibco_diff_gs.isnull(), other=ro_gibco_diff)

    ro_dat = rnaseq_obj.data.loc[ro_gibco_diff]
    ro_dat.index = ro_gibco_diff_gs
    ro_dat = np.log2(ro_dat + 1)

    # ro_dat = salmon_dat.loc[ro_gibco_diff]
    # ro_dat.index = ro_gibco_diff_gs
    # # dropna() here
    # # all others are present
    # ro_dat = np.log2(ro_dat.dropna() + 0.01)

    # rearrange columns
    cols = (
        ro_dat.columns[ro_dat.columns.str.contains('GBM')].tolist() +
        ['GIBCO_NSC_P4'] +
        ro_dat.columns[ro_dat.columns.str.contains('DURA')].tolist()
    )
    ro_dat = ro_dat.loc[:, cols]

    # insert spacing columns
    idx = np.where(ro_dat.columns.str.contains('GIBCO'))[0][0]
    ro_dat.insert(idx, '', np.nan)
    ro_dat.insert(idx + 2, ' ', np.nan)

    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(ro_dat, cmap=cmap, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "consistently_in_ref_only.png"), dpi=200)

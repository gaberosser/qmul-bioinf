import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from load_data import rnaseq_data
from rnaseq import differential_expression, general
from settings import LOCAL_DATA_DIR
from utils import output, setops, excel, ipa, reference_genomes

if __name__ == "__main__":

    outdir = output.unique_output_dir("cross_validate_de_multiple_refs", reuse_empty=True)
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
        ('H9', 'NSC'),
        ('H1', 'NSC'),
    ]
    external_ref_labels = [t[0] for t in external_refs]
    ref_samples = reduce(
        lambda x, y: x+y,
        [rnaseq_obj.meta.index[rnaseq_obj.meta.index.str.contains(t[0])].tolist() for t in external_refs]
    )
    cols_combined_ref = ['ref_intersect', 'ref_intersect2', 'ref_union']
    cols = pids + [t[0] for t in external_refs] + cols_combined_ref
    cols_comb_only = pids + cols_combined_ref
    de_res = differential_expression.compute_cross_de(rnaseq_obj, pids, external_references=external_refs, **de_params)

    # add the combined DE results for the refs combined
    for pid in pids:
        # complete intersection
        the_idx = sorted(reduce(intersecter, [de_res[(pid, t)].index for t in external_ref_labels]))
        one_cols = de_res[(pid, pid)].columns
        tups = reduce(lambda x, y: x + y, [zip([t] * one_cols.size, one_cols.tolist()) for t in external_ref_labels])
        the_cols = pd.MultiIndex.from_tuples(tups, names=['ref', 'field'])
        the_block = pd.DataFrame(index=the_idx, columns=the_cols)
        for t in external_ref_labels:
            the_block.loc[the_idx, t] = de_res[(pid, t)].loc[the_idx].values
        de_res[(pid, 'ref_intersect')] = the_block

        # intersect 2
        this_venn, _ = setops.venn_from_arrays(*[de_res[(pid, t)].index for t in external_ref_labels])
        the_idx = reduce(unioner, [this_venn[k] for k in setops.binary_combinations_sum_gte(len(external_refs), 2)])
        the_block = pd.DataFrame(index=the_idx, columns=the_cols)
        for t in external_ref_labels:
            try:
                the_block.loc[the_idx, t] = de_res[(pid, t)].loc[the_idx].values
            except KeyError:
                # no matches for this ref - no problem
                pass
        de_res[(pid, 'ref_intersect2')] = the_block

        # union
        the_idx = sorted(reduce(unioner, [de_res[(pid, t)].index for t in external_ref_labels]))
        the_block = pd.DataFrame(index=the_idx, columns=the_cols)
        for t in external_ref_labels:
            try:
                the_block.loc[the_idx, t] = de_res[(pid, t)].loc[the_idx].values
            except KeyError:
                # no matches for this ref - no problem
                pass
        de_res[(pid, 'ref_union')] = the_block

    # counts of DE genes
    de_counts = pd.DataFrame(index=pids, columns=cols)
    for pid in pids:
        for pid2 in cols:
            de_counts.loc[pid, pid2] = de_res[(pid, pid2)].shape[0]

    # we can integrate the multiple refs in a number of ways
    # loop over these now
    for c in cols_combined_ref:

        # now we need to compare the paired results with every other result (Gibco and other iNSC)
        cols = pids + [c]
        pair_only = pd.DataFrame(index=pids, columns=cols)
        ref_only = pd.DataFrame(index=pids, columns=cols)

        # loop over GBM samples
        for pid in pids:
            # syngeneic comparison
            the_pair = de_res[(pid, pid)]

            # loop over refs
            # when this is the same as the syngeneic comparison, there will (obviously) be no 'pair only' or 'ref only'
            # genes!
            for pid2 in cols:
                the_ref = de_res[(pid, pid2)]
                the_sets, _ = setops.venn_from_arrays(the_pair.index, the_ref.index)
                pair_only.loc[pid, pid2] = the_sets['10']
                ref_only.loc[pid, pid2] = the_sets['01']

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
        po_intersection_insc = pd.Series(index=pids)
        for pid in pids:
            # first, find the genes that are always PO when an iNSC reference is used
            tmp = reduce(intersecter, pair_only.loc[pid, pair_only.index[pair_only.index != pid]])
            # second, find the union of genes that are PO when this external reference is used
            tmp2 = pair_only.loc[pid, c]
            # we want anything in the first part that is NOT in the second part
            po_intersection_insc.loc[pid] = tmp.difference(tmp2)

        # find DE genes that are always unique to a given reference (regardless of the GBM)
        po_specific_to_reference = [
            sorted(
                reduce(intersecter, po_diff.loc[~po_diff.index.str.contains(pid), pid])
            ) for pid in cols
        ]
        po_specific_to_reference = pd.Series(po_specific_to_reference, index=cols)

        # get the genes that consistently differ in the pair comparison only and NOT in Gibco (across all patients)
        # these will have an expression pattern in Gibco similar to GBM, so that they do NOT appear
        po_ref_diff = po_specific_to_reference.loc[c]
        po_ref_diff_gs = reference_genomes.ensembl_to_gene_symbol(po_ref_diff)
        po_ref_diff_gs = po_ref_diff_gs.where(~po_ref_diff_gs.isnull(), po_ref_diff)

        po_dat = rnaseq_obj.data.loc[po_ref_diff]
        po_dat.index = po_ref_diff_gs
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

        fig = plt.figure(figsize=(7, 3 + 10 / 50. * len(po_ref_diff)))
        ax = fig.add_subplot(111)
        ax = sns.heatmap(po_dat, cmap=cmap, ax=ax)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "consistently_in_pair_only_%s.png" % c), dpi=200)

        # export those same genes to a file, adding gene symbols
        for_export = rnaseq_obj.data.loc[po_ref_diff, the_cols]
        general.add_gene_symbols_to_ensembl_data(for_export)
        for_export.to_excel(os.path.join(outdir, 'consistently_in_pair_only_%s.xlsx' % c))

        # export the pair_only genes (GBM vs iNSC paired -- GBM vs ref) to a file
        print "Combination type: %s" % c
        po_de_export = {}
        for pid in pids:
            # get PO DE genes
            the_idx = pair_only.loc[pid, c]
            # subtract the correction genes
            the_corr = po_ref_diff
            print "Subtracting %d correction genes from %d PO DE genes to leave %d PO DE genes." % (
                len(the_corr),
                len(the_idx),
                len(set(the_idx).difference(the_corr))
            )
            the_idx = list(set(the_idx).difference(the_corr))
            po_de_export[pid] = de_res[(pid, pid)].loc[the_idx]
        excel.pandas_to_excel(po_de_export, os.path.join(outdir, 'pair_only_de_lists_corrected_%s.xlsx' % c))




    ## TODO: refactor from here for ref_only

    # # get the genes that consistently appear in the Gibco reference comparison only and NOT in any other reference
    # # these will have a different expression pattern in Gibco to GBM (while the level in iNSC will not differ from GBM)
    # ro_diff = pd.DataFrame(index=ref_only.index, columns=ref_only.columns)
    # for pid in pids:
    #     for pid2 in ref_only.columns:
    #         the_ref = ref_only.loc[pid, pid2]
    #         all_else = ref_only.loc[pid, ref_only.columns != pid2]
    #         union_all_else = reduce(set.union, all_else, set())
    #         ro_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))
    #
    # # get counts like this
    # ro_diff.applymap(len)
    #
    # # this computes the set of genes that ALWAYS appears in the ref_only comparison for each possible ref
    # ro_each = [
    #     sorted(
    #         reduce(lambda x, y: set(x).intersection(y), ro_diff.loc[~ro_diff.index.str.contains(pid), pid])
    #     ) for pid in cols
    # ]
    # ro_each = pd.Series(ro_each, index=cols)
    #
    # # ro_gibco_diff = sorted(reduce(lambda x, y: set(y).intersection(x), ro_diff.loc[:, 'GIBCO']))
    # ro_gibco_diff = ro_each.loc['GIBCO']
    # ro_gibco_diff_gs = references.ensembl_to_gene_symbol(ro_gibco_diff)
    # # the lincRNA symbols are missing, so keep ENSG for those
    # ro_gibco_diff_gs = ro_gibco_diff_gs.where(~ro_gibco_diff_gs.isnull(), other=ro_gibco_diff)
    #
    # ro_dat = rnaseq_obj.data.loc[ro_gibco_diff]
    # ro_dat.index = ro_gibco_diff_gs
    # ro_dat = np.log2(ro_dat + 1)
    #
    # # ro_dat = salmon_dat.loc[ro_gibco_diff]
    # # ro_dat.index = ro_gibco_diff_gs
    # # # dropna() here
    # # # all others are present
    # # ro_dat = np.log2(ro_dat.dropna() + 0.01)
    #
    # # rearrange columns
    # cols = (
    #     ro_dat.columns[ro_dat.columns.str.contains('GBM')].tolist() +
    #     ['GIBCO_NSC_P4'] +
    #     ro_dat.columns[ro_dat.columns.str.contains('DURA')].tolist()
    # )
    # ro_dat = ro_dat.loc[:, cols]
    #
    # # insert spacing columns
    # idx = np.where(ro_dat.columns.str.contains('GIBCO'))[0][0]
    # ro_dat.insert(idx, '', np.nan)
    # ro_dat.insert(idx + 2, ' ', np.nan)
    #
    # fig = plt.figure(figsize=(7, 10))
    # ax = fig.add_subplot(111)
    # ax = sns.heatmap(ro_dat, cmap=cmap, ax=ax)
    # plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    # plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "consistently_in_ref_only.png"), dpi=200)

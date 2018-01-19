import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import venn
import references
from rnaseq import differential_expression, general
from utils import output, setops, excel, ipa
from load_data import rnaseq_data


if __name__ == "__main__":

    outdir = output.unique_output_dir("cruk_trial_2", reuse_empty=True)
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
    cols = pids + [t[0] for t in external_refs]
    de_res = differential_expression.compute_cross_de(rnaseq_obj, pids, external_references=external_refs, **de_params)

    pair_only = pd.DataFrame(index=pids, columns=cols)
    # ref_only = pd.DataFrame(index=pids, columns=cols)

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
            # ref_only.loc[pid, pid2] = the_sets['01']
    po_counts = pair_only.applymap(len)

    # For each reference, get the DE genes that are pair only in that reference and not in any of the iNSC
    po_diff = pd.DataFrame(index=pair_only.index, columns=external_ref_labels)
    for pid in pids:
        for pid2 in external_ref_labels:
            the_ref = pair_only.loc[pid, pid2]
            all_else = pair_only.loc[pid, pids]
            union_all_else = reduce(set.union, all_else, set())
            po_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    # Intersection down the columns gives us a correction list for each reference
    po_specific_to_ref = po_diff.apply(lambda x: reduce(intersecter, x))

    # Intersection across the references gives us a final list that need correcting
    po_specific_to_all_refs = sorted(reduce(intersecter, po_specific_to_ref))

    # get the genes that consistently differ in the pair comparison only and NOT in Gibco (across all patients)
    # these will have an expression pattern in Gibco similar to GBM, so that they do NOT appear
    po_specific_to_all_refs_gs = references.ensembl_to_gene_symbol(po_specific_to_all_refs)
    po_specific_to_all_refs_gs = po_specific_to_all_refs_gs.where(~po_specific_to_all_refs_gs.isnull(), po_specific_to_all_refs)

    po_dat = rnaseq_obj.data.loc[po_specific_to_all_refs]
    po_dat.index = po_specific_to_all_refs_gs
    po_dat = np.log2(po_dat + 1)

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

    fig = plt.figure(figsize=(7, 3 + 10 / 50. * len(po_specific_to_all_refs)))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(po_dat, cmap=cmap, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "consistently_in_pair_only_across_all_refs.png"), dpi=200)

    # export those same genes to a file, adding gene symbols
    for_export = rnaseq_obj.data.loc[po_specific_to_all_refs, the_cols]
    general.add_gene_symbols_to_ensembl_data(for_export)
    for_export.to_excel(os.path.join(outdir, 'consistently_in_pair_only_across_all_refs.xlsx'))

    # correct the reference PO lists, take the intersection, then export to a file
    po_de_export = {}
    for pid in pids:
        this_row = pair_only.loc[pid, external_ref_labels]
        this_genes_pre = reduce(intersecter, this_row)
        this_genes = sorted(this_genes_pre.difference(po_specific_to_all_refs))
        print "PID %s. Subtracted %d correction genes from the %d PO intersection genes to leave %d PO genes" % (
            pid, len(po_specific_to_all_refs), len(this_genes_pre), len(this_genes)
        )
        po_de_export[pid] = de_res[(pid, pid)].loc[this_genes]

    excel.pandas_to_excel(po_de_export, os.path.join(outdir, 'pair_only_de_lists_corrected.xlsx'))

    # plot: how many DE genes are present in each reference comparison?

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            de_res[(pid, r)].index for r in external_ref_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn3(cts, set_labels=external_ref_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'number_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'number_de_multiple_references.tiff'), dpi=200)

    # plot: how many DE genes are in the pair only comparison when each reference is used?
    # NB apply correction

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_3 = pd.DataFrame(0, index=pids, columns=external_ref_labels)
    n_pair_only_2 = pd.DataFrame(0, index=pids, columns=external_ref_labels)

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            set(pair_only.loc[pid, r]).difference(po_specific_to_all_refs) for r in external_ref_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn3(cts, set_labels=external_ref_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)

        for i, r in enumerate(external_ref_labels):
            # this will fail for anything other than 3 refs
            n_pair_only_3.loc[pid, r] = cts['111']
            n_pair_only_2.loc[pid, r] = cts['111']
            for k in setops.binary_combinations_sum_eq(len(external_ref_labels), 2):
                if k[i] == '1':
                    n_pair_only_2.loc[pid, r] += cts[k]
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'number_po_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'number_po_de_multiple_references.tiff'), dpi=200)

    # plot: overlap between individual references in terms of PO genes shared
    pct_pair_only_3 = n_pair_only_3 / po_counts.loc[:, external_ref_labels] * 100.
    pct_pair_only_2 = n_pair_only_2 / po_counts.loc[:, external_ref_labels] * 100.

    ax = pct_pair_only_2.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
    ax.set_xlabel('Patient')
    ax.set_ylabel("% DE genes shared")
    ax.set_ylim([0, 100])
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_2_of_3.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_2_of_3.tiff"), dpi=200)

    ax = pct_pair_only_3.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
    ax.set_xlabel('Patient')
    ax.set_ylabel("% DE genes shared")
    ax.set_ylim([0, 100])
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_3_of_3.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_3_of_3.tiff"), dpi=200)

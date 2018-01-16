import multiprocessing as mp
import os
import re
import collections
import numpy as np
from plotting import venn
import pandas as pd
from matplotlib import pyplot as plt, gridspec
import seaborn as sns
import references
from rnaseq import filter, differential_expression
from settings import LOCAL_DATA_DIR
from utils import output, setops, excel, ipa
from load_data import rnaseq_data


def add_gene_symbols(df):
    """
    Add gene symbols to the DataFrame df which is indexed by Ensembl IDs
    """
    gs = references.ensembl_to_gene_symbol(df.index)
    # resolve any duplicates arbitrarily (these should be rare)
    gs = gs.loc[~gs.index.duplicated()]
    df.insert(0, 'Gene Symbol', gs)


def add_fc_direction(df):
    direction = pd.Series(index=df.index, name='Direction')
    direction.loc[df.logFC < 0] = 'down'
    direction.loc[df.logFC > 0] = 'up'
    df.insert(df.shape[1], 'Direction', direction)


def run_one_de(the_data, the_groups, the_comparison, lfc=1, fdr=0.01, method='QLGLM', return_full=False):
    if method == 'QLGLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmqlfit(
            the_data,
            the_groups,
            the_contrast,
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )
    elif method == 'GLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmfit(
            the_data,
            the_groups,
            the_contrast,
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )
    elif method == 'exact':
        res = differential_expression.edger_exacttest(
            the_data,
            the_groups,
            pair=the_comparison[::-1],
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )

    add_gene_symbols(res)
    add_fc_direction(res)

    return res


def venn_set_to_dataframe(data, venn_set, set_labels, include_sets=None, full_data=None):
    """
    Given the input DE data and Venn sets, generate a long format dataframe containing all the data, one column
    per patient and one row per gene.
    Optionally filter the sets to include only a subset.
    Optionally include non-significant results too.
    :param data: Dict containing DE results, keyed by the entries of set_labels
    :param venn_set:
    :param set_labels:
    :param include_sets:
    :param full_data: If supplied, this has the same format as `data`, but the lists are complete so that even non-
    significant results can be accessed.
    :return:
    """
    if include_sets is not None:
        venn_set = dict([
            (k, v) for k, v in venn_set.items() if k in include_sets
        ])

    res = []
    for k in venn_set:
        the_genes = venn_set[k]
        # populate with individual patient results
        blocks = []
        consistency_check = []
        for i, t in enumerate(k):
            pid = set_labels[i]
            this_datum = pd.DataFrame(index=the_genes, columns=[pid, "%s_logFC" % pid, "%s_FDR" % pid])
            if t == '1':
                this_datum.loc[the_genes, pid] = 'Y'
                this_datum.loc[the_genes, "%s_logFC" % pid] = data[pid].loc[the_genes, 'logFC']
                this_datum.loc[the_genes, "%s_FDR" % pid] = data[pid].loc[the_genes, 'FDR']
                cc = data[pid].loc[the_genes, 'Direction']
                cc.name = pid
                consistency_check.append(cc)
            else:
                this_datum.loc[the_genes, pid] = 'N'
                if full_data is not None:
                    this_datum.loc[the_genes, "%s_logFC" % pid] = full_data[pid].loc[the_genes, 'logFC']
                    this_datum.loc[the_genes, "%s_FDR" % pid] = full_data[pid].loc[the_genes, 'FDR']

            blocks.append(this_datum)

        core_block = pd.concat(blocks, axis=1)
        # assess consistency of DE direction
        consist = pd.Series(index=the_genes)

        if len(consistency_check) > 0:
            consistency_check = pd.concat(consistency_check, axis=1)
            idx = consistency_check.apply(lambda col: col == consistency_check.iloc[:, 0]).all(axis=1)
            consist.loc[idx] = 'Y'
            consist.loc[~idx] = 'N'

        core_block.insert(core_block.shape[1], 'consistent', consist)
        res.append(core_block)

    # check: no genes should be in more than one data entry
    for i, k in enumerate(venn_set):
        for j, k2 in enumerate(venn_set):
            if k == k2: continue
            bb = len(res[i].index.intersection(res[j].index))
            if bb > 0:
                raise AttributeError("Identified %d genes that are in BOTH %s and %s" % (bb, k, k2))

    res = pd.concat(res, axis=0)

    # add gene symbols
    add_gene_symbols(res)

    return res


if __name__ == "__main__":
    outdir = output.unique_output_dir("compare_paired_de", reuse_empty=True)

    # all n=2 samples and RTK II samples
    pids = ['019', '030', '031', '017', '050', '054']
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

    # only keep gene counts
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]

    # compute DE between hGIC and paired iNSC
    de_res = {}
    de_res_full = {}
    for pid in pids:
        hgic_samples = rnaseq_obj.meta.index[rnaseq_obj.meta.index.str.contains(pid)]
        the_data = rnaseq_obj.data.loc[:, hgic_samples]
        the_groups = rnaseq_obj.meta.loc[hgic_samples, 'type']
        the_comparison = ['GBM', 'iNSC']
        de_res[pid] = run_one_de(the_data, the_groups, the_comparison, **de_params)
        de_res_full[pid] = run_one_de(the_data, the_groups, the_comparison, return_full=True, **de_params)
        print "GBM %s paired comparison, %d DE genes" % (pid, de_res[pid].shape[0])

    venn_set, venn_ct = setops.venn_from_arrays(*[de_res[pid].index for pid in pids])

    # add null set manually
    de_genes_all = reduce(lambda x, y: set(x).union(y), venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(de_res_full[pids[0]].index.difference(de_genes_all))
    venn_ct[k_null] = len(venn_set[k_null])

    # check direction is the same
    venn_set_consistent = {}
    venn_set_inconsistent = {}
    for k in venn_set:
        the_genes = venn_set[k]
        the_pids = [pids[i] for i, t in enumerate(k) if t == '1']
        the_de_direction = pd.DataFrame(dict([(pid, de_res[pid].loc[the_genes, 'Direction']) for pid in the_pids]))
        # check consistency
        idx = the_de_direction.apply(lambda col: col == the_de_direction.iloc[:, 0]).all(axis=1)
        venn_set_consistent[k] = the_de_direction.loc[idx].index
        venn_set_inconsistent[k] = the_de_direction.loc[~idx].index

    # generate list and save to Excel file
    data = venn_set_to_dataframe(de_res, venn_set, pids, full_data=de_res_full)

    # save
    data.to_excel(os.path.join(outdir, 'de_list_compare_patients.xlsx'))

    # generate an expanded core gene set, defined as genes that are DE in both RTK I and RTK II (any number of patients)
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])
    expanded_core_sets = []
    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        nmatch = 0
        for grp, grp_idx in subgroup_ind.items():
            if this_k[grp_idx].any():
                nmatch += 1
        if nmatch > 1:
            # add to the expanded core set
            expanded_core_sets.append(k)



    expanded_core = venn_set_to_dataframe(de_res, venn_set, pids, include_sets=expanded_core_sets)

    # save
    expanded_core.to_excel(os.path.join(outdir, 'expanded_core_gene_list.xlsx'))

    # subgroup-specific lists (just for completeness)
    sg_specific_sets = {}
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        sg_specific_sets[grp] = k

    subgroup_specific = venn_set_to_dataframe(de_res, venn_set, pids, include_sets=sg_specific_sets.values())
    subgroup_specific.to_excel(os.path.join(outdir, 'subgroup_specific_gene_list.xlsx'))


    # UpsetR attribute plots
    set_labels = ['019', '030', '031', '017', '050', '054']
    data_for_upset = [de_res[pid].index for pid in set_labels]  # this will be supplied to the function

    # 1. Descending order
    upset1 = venn.upset_set_size_plot(
        data_for_upset,
        set_labels,
        min_size=10,
        n_plot=30
    )
    upset1['figure'].savefig(os.path.join(outdir, "upset_descending.png"), dpi=200)

    # 2. Grouped by number of participants in the sets
    upset2 = venn.upset_set_size_plot(
        data_for_upset,
        set_labels,
        order_by_n_members=True,
        min_size=30,
        n_plot=50
    )
    upset2['figure'].savefig(os.path.join(outdir, "upset_grouped_descending.png"), dpi=200)


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


if __name__ == "__main__":
    outdir = output.unique_output_dir("compare_paired_de", reuse_empty=True)

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

    # only keep gene counts
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]

    # compute DE between hGIC and paired iNSC
    de_res = {}
    for pid in pids:
        hgic_samples = rnaseq_obj.meta.index[rnaseq_obj.meta.index.str.contains(pid)]
        the_data = rnaseq_obj.data.loc[:, hgic_samples]
        the_groups = rnaseq_obj.meta.loc[hgic_samples, 'type']
        the_comparison = ['GBM', 'iNSC']
        de_res[pid] = run_one_de(the_data, the_groups, the_comparison, **de_params)
        print "GBM %s paired comparison, %d DE genes" % (pid, de_res[pid].shape[0])

    venn_set, venn_ct = setops.venn_from_arrays(*[de_res[pid].index for pid in pids])

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
    idx = []
    cols = reduce(
        lambda x, y: x + y,
        [
            ["%s" % pid, "%s_logFC" % pid, "%s_FDR" % pid] for pid in pids
        ]
    )
    data = []
    for k in venn_set:
        the_genes = venn_set[k]
        # populate with individual patient results
        blocks = []
        for i, t in enumerate(k):
            pid = pids[i]
            this_datum = pd.DataFrame(index=the_genes, columns=[pid, "%s_logFC" % pid, "%s_FDR" % pid])
            if t == '1':
                this_datum.loc[the_genes, pid] = 'Y'
                this_datum.loc[the_genes, "%s_logFC" % pid] = de_res[pid].loc[the_genes, 'logFC']
                this_datum.loc[the_genes, "%s_FDR" % pid] = de_res[pid].loc[the_genes, 'FDR']
            else:
                this_datum.loc[the_genes, pid] = 'N'
                # this_datum.loc[the_genes, "%s_logFC" % pid] = de_res[pid].loc[the_genes, 'logFC']
                # this_datum.loc[the_genes, "%s_FDR" % pid] = de_res[pid].loc[the_genes, 'FDR']
            blocks.append(this_datum)
        core_block = pd.concat(blocks, axis=1)
        consist = pd.Series(index=the_genes)
        consist.loc[venn_set_consistent[k]] = 'Y'
        consist.loc[venn_set_inconsistent[k]] = 'N'
        core_block.insert(core_block.shape[1], 'consistent', consist)
        data.append(core_block)

    # check: no genes should be in more than one data entry
    for i, k in enumerate(venn_set):
        for j, k2 in enumerate(venn_set):
            if k == k2: continue
            bb = len(data[i].index.intersection(data[j].index))
            if bb > 0: print "%s, %s, %d" % (k, k2, bb)

    data = pd.concat(data, axis=0)

    # add gene symbols
    add_gene_symbols(data)

    # save
    data.to_excel(os.path.join(outdir, 'de_list_compare_patients.xlsx'))

    # plot UpsetR basic plot
    # will move this to plotting.venn when complete
    set_labels = pids  # this will be supplied to the function
    args = [de_res[pid].index for pid in set_labels]  # this will be supplied to the function
    venn_set, venn_ct = setops.venn_from_arrays(*args)
    include_singletons = False
    n_to_plot = 30
    width = 0.9
    ms = 10  # marker size in pt
    lightgrey = '#cecece'

    # could alternatively sort by no. intersecting sets?
    if include_singletons:
        ordered_counts = sorted(venn_ct.items(), key=lambda x: x[1], reverse=True)
    else:
        # exclude any results with only one set
        ordered_counts = sorted(
            [t for t in venn_ct.items() if len(t[0].replace('0', '')) > 1],
            key=lambda x: x[1],
            reverse=True
        )
    ordered_counts = ordered_counts[:n_to_plot]

    gs_kw = dict(
        left=0.05,
        right=0.99,
        top=0.99,
        bottom=0.1,
        wspace=0.1,
        hspace=0.01,
        height_ratios=[6, 3],
        width_ratios=[3, 6],
    )

    gs = gridspec.GridSpec(nrows=2, ncols=2, **gs_kw)
    fig = plt.figure(figsize=(9, 6))
    ax_tl = fig.add_subplot(gs[0, 0])
    ax_set_size = fig.add_subplot(gs[1, 0])
    ax_intersect = fig.add_subplot(gs[1, 1], sharey=ax_set_size)
    ax_main = fig.add_subplot(gs[0, 1], sharex=ax_intersect)

    # hide some things
    ax_tl.set_visible(False)
    plt.setp(ax_intersect.get_yticklabels(), visible=False)
    plt.setp(ax_main.get_xticklabels(), visible=False)
    plt.setp(ax_intersect.get_xticklabels(), visible=False)

    # data
    x_arr = np.arange(len(ordered_counts)) + 0.5
    y_arr = np.arange(len(set_labels))

    # main bar chart
    ax_main.bar(x_arr, [t[1] for t in ordered_counts], width=width)
    ax_main.set_ylabel('Number of DE genes in set')

    # bottom right set intersections
    # grey markers everywhere
    for y in y_arr:
        ax_intersect.plot(x_arr, np.ones_like(x_arr) * y, marker='o', mfc=lightgrey, mec='none', ms=ms, ls='none')
    # black markers only on sets that are included
    for i, (k, v) in enumerate(ordered_counts):
        x = x_arr[i]
        y = [j for j, u in enumerate(k) if u == '1']
        ax_intersect.plot(x * np.ones(len(y)), y, marker='o', mfc='k', mec='k', ms=ms, ls='none')

    # bottom left : singleton set size
    singleton_sizes = [len(t) for t in args]
    ax_set_size.barh(y_arr + 0.5, singleton_sizes, -width, align='edge')
    ax_set_size.invert_xaxis()
    ax_set_size.set_ylim([-.5, len(set_labels) - .5])
    ax_set_size.yaxis.tick_right()
    ax_set_size.set_yticks(y_arr)
    ax_set_size.set_yticklabels(set_labels)
    ax_set_size.set_xlabel("Number of DE genes in single comparison")






from load_data import rnaseq_data
from rnaseq.differential_expression import edger
from rnaseq.filter import filter_by_cpm
import pandas as pd
import numpy as np
from scipy import stats
import references
import os
from utils import output, setops
from matplotlib import pyplot as plt
import seaborn as sns


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


if __name__ == '__main__':
    outdir = output.unique_output_dir("paired_rnaseq")
    pids = ['017', '050', '054', '061']
    # pids = [t for t in rnaseq_data.PATIENT_LOOKUP_STAR if t != 'GIBCO']
    obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    obj.data = obj.data.loc[obj.data.index.str.contains('ENSG')]

    dat_filt = filter_by_cpm(obj.data, min_n_samples=2)

    # look at the distribution of counts (ECDF)

    q75 = {}

    colours = plt.cm.jet(np.linspace(0, 1., obj.meta.shape[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    focus = 'DURA061_NSC_N6_P4'
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)
        q75[sname] = the_data.iloc[np.where(pct >= 0.75)[0][0]]
        cdf = the_data.cumsum() / float(the_data.sum())
        alpha = None if sname == focus else 0.3
        ax.plot(np.log2(the_data + 1), pct, c=colours[i], label=sname, alpha=alpha)
        # ax.plot(the_data, cdf, c=colours[i], label=sname)
    ax.legend(loc='lower right')
    ax.set_xlabel('Log2(count + 1)')
    ax.set_ylabel('Percentile')

    fig.savefig(os.path.join(outdir, 'percentile_vs_log_count.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_log_count.pdf'))

    # repeat and zoom in on the top 1% most expressed genes
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)

        ax.plot(the_data.loc[pct>0.99], pct[pct>0.99], c=colours[i], label=sname)
    ax.legend(loc='lower right')
    ax.set_xlabel('Count')
    ax.set_ylabel('Percentile')
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct.pdf'))

    # highlight DURA061 samples in the same plot
    colours = plt.cm.prism(np.linspace(0, 1., obj.meta.shape[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)
        c = colours[i] if 'DURA061' in sname else 'gray'
        alpha = None if 'DURA061' in sname else 0.3
        label = sname if 'DURA061' in sname else None
        ax.plot(the_data.loc[pct > 0.99], pct[pct > 0.99], c=c, label=label, alpha=alpha)
    ax.legend(loc='lower right')
    ax.set_xlabel('Count')
    ax.set_ylabel('Percentile')
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct_focus061NSC.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct_focus061NSC.pdf'))

    # scatterplots where replicates exist (all but NSC017)
    fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True)
    for i, pid in enumerate(pids):
        for j, typ in enumerate(('GBM', 'DURA')):
            the_idx = dat_filt.columns.str.contains("%s%s" % (typ, pid))
            if the_idx.sum() != 2:
                axs[j][i].set_title("%s%s (n=1)" % (typ, pid))
                continue
            the_data = np.log2(dat_filt.loc[:, the_idx] + 1.)
            axs[j][i].scatter(the_data.iloc[:, 0], the_data.iloc[:, 1])

            axs[j][i].set_title("%s%s r=%.3f" % (
                typ, pid, stats.linregress(the_data.iloc[:, 0], the_data.iloc[:, 1]).rvalue
            ))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'replicate_logcorr.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'replicate_logcorr.pdf'))


    """
    So it seems we have a problem with 061: the two NSC replicates differ more than other replicates
    This is the underlying reason for the small number of DE genes (which we will now see...)
    """

    # Here I tried to 'compress' the data by removing the most/least variable genes
    # it doesn't lead to more DE genes in 061

    data_compressed = pd.DataFrame(index=dat_filt.index, columns=dat_filt.columns)
    q_lower = 0.05
    q_upper = 0.95
    for col, arr in dat_filt.iteritems():
        arr = pd.DataFrame(arr, copy=True)
        # identify cutoffs and insert nan
        pct = stats.rankdata(arr, method='max') / float(arr.size)
        arr.loc[pct < q_lower] = np.nan
        arr.loc[pct > q_upper] = np.nan
        data_compressed.loc[:, col] = arr
    data_compressed.dropna(axis=0, inplace=True)

    # Here I just remove the top 5 genes in the 061 NSC P4 sample
    # doesn't help

    focus = 'DURA061_NSC_N6_P4'
    data_capped = pd.DataFrame(dat_filt, copy=True)
    to_drop = data_capped.loc[:, focus].sort_values(ascending=False).index[:5]
    data_capped.drop(to_drop, inplace=True)

    de = {}
    de_up = {}
    de_down = {}
    de_counts = {}
    de_counts_up = {}
    de_counts_down = {}
    de_matched = {}
    de_gibco = {}

    source = obj.data
    # source = data_capped
    # source = data_compressed


    xl_writer = pd.ExcelWriter(os.path.join(outdir, "individual_gene_lists.xlsx"))

    for pid in pids:
        try:
            the_idx = obj.meta.index.str.contains(pid)
            the_data = obj.data.loc[:, the_idx]

            the_data = filter_by_cpm(the_data, min_n_samples=1)
            the_genes = the_data.index

            the_groups = obj.meta.loc[the_idx, 'type'].values
            the_contrast = "GBM - iNSC"

            de_matched[pid] = edger(the_data, the_groups, the_contrast)

            # repeat with gibco reference
            # use the same genes, rather than filtering again
            the_idx = (obj.meta.index.str.contains(pid) & (obj.meta.type == 'GBM')) | obj.meta.index.str.contains('GIBCO')
            the_data = obj.data.loc[the_genes, the_idx]
            the_groups = obj.meta.loc[the_idx, 'type'].values
            the_contrast = "GBM - NSC"

            de_gibco[pid] = edger(the_data, the_groups, the_contrast)

            # Separate into sets
            # all
            de[pid], de_counts[pid] = setops.venn_from_arrays(de_matched[pid].index, de_gibco[pid].index)

            # up only
            idx_up_match = de_matched[pid].loc[de_matched[pid].logFC > 0].index
            idx_up_ref = de_gibco[pid].loc[de_gibco[pid].logFC > 0].index
            de_up[pid], de_counts_up[pid] = setops.venn_from_arrays(idx_up_match, idx_up_ref)

            # down only
            idx_down_match = de_matched[pid].loc[de_matched[pid].logFC < 0].index
            idx_down_ref = de_gibco[pid].loc[de_gibco[pid].logFC < 0].index
            de_down[pid], de_counts_down[pid] = setops.venn_from_arrays(idx_down_match, idx_down_ref)

            # write to files, one worksheet per list (5 per individual)
            # paired comparison (all)
            block = de_matched[pid].sort_values('PValue')
            add_gene_symbols(block)
            add_fc_direction(block)
            block.to_excel(xl_writer, 'GBM%s_pair_all' % pid)

            # reference comparison (all)
            block = de_gibco[pid].sort_values('PValue')
            add_gene_symbols(block)
            add_fc_direction(block)
            block.to_excel(xl_writer, 'GBM%s_ref_all' % pid)

            # paired comparison only
            block = de_matched[pid].loc[de[pid]['10']].sort_values('PValue')
            add_gene_symbols(block)
            add_fc_direction(block)
            block.to_excel(xl_writer, 'GBM%s_pair_only' % pid)

            # paired and ref comparison
            # main block is taken from the matched comparison
            block_match = de_matched[pid].loc[de[pid]['11']].sort_values('PValue')
            add_fc_direction(block_match)
            # add second block of values for the ref comparison, changing colnames
            block_ref = de_gibco[pid].loc[block_match.index]
            add_fc_direction(block_ref)
            block_ref.columns = ["%s_REF" % t for t in block_ref.columns]

            block = pd.concat((block_match, block_ref), axis=1)

            # add column with agreement
            agree = (np.sign(block.logFC) == np.sign(block.logFC_REF))
            block.insert(block.shape[1], 'agreement', agree.astype(str))

            if (~agree).any():
                print "GBM%s: %d / %d match and ref DE genes have discordant direction" % (
                    pid, (~agree).sum(), block.shape[0]
                )

            add_gene_symbols(block)
            block.to_excel(xl_writer, 'GBM%s_pair_and_ref' % pid)


            # ref comparison only
            block = de_gibco[pid].loc[de[pid]['01']].sort_values('PValue')
            add_gene_symbols(block)
            add_fc_direction(block)
            block.to_excel(xl_writer, 'GBM%s_ref_only' % pid)

        except Exception as exc:
            print "Patient %s failed." % pid
            print repr(exc)

    xl_writer.save()

    # Generate a series of Venn diagrams
    from plotting import venn
    fig, axs = plt.subplots(ncols=len(pids), nrows=3, figsize=[12, 8])
    set_labels = ('Isogenic', 'Reference')
    for i, pid in enumerate(pids):
        venn.venn2(de_counts[pid], ax=axs[0, i], set_labels=set_labels)
        axs[0, i].set_title("GBM%s all" % pid)
        venn.venn2(de_counts_up[pid], ax=axs[1, i], set_labels=set_labels)
        axs[1, i].set_title("Up in GBM%s vs NSC" % pid)
        venn.venn2(de_counts_down[pid], ax=axs[2, i], set_labels=set_labels)
        axs[2, i].set_title("Down in GBM%s vs NSC" % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_gene_counts.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "de_gene_counts.pdf"))

    # onwards! let's bring in the methylation data

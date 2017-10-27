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
    # RTK II samples
    # pids = ['017', '050', '054', '061']
    # all n=2 samples
    pids = ['018', '044', '049', '050', '052', '054', '061']
    # all samples
    # pids = [t for t in rnaseq_data.PATIENT_LOOKUP_STAR if t != 'GIBCO']
    obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    obj.data = obj.data.loc[obj.data.index.str.contains('ENSG')]

    dat_filt = filter_by_cpm(obj.data, min_n_samples=2)

    de = {}
    de_up = {}
    de_down = {}
    de_counts = {}
    de_counts_up = {}
    de_counts_down = {}
    de_matched = {}
    de_gibco = {}

    # all gene lists combined in one file (one sheet per list)
    xl_writer = pd.ExcelWriter(os.path.join(outdir, "individual_gene_lists.xlsx"))

    for pid in pids:
        # try:
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
        block.to_excel(os.path.join(outdir, 'GBM%s_pair_all.xlsx' % pid))

        # reference comparison (all)
        block = de_gibco[pid].sort_values('PValue')
        add_gene_symbols(block)
        add_fc_direction(block)
        block.to_excel(xl_writer, 'GBM%s_ref_all' % pid)
        block.to_excel(os.path.join(outdir, 'GBM%s_ref_all.xlsx' % pid))

        # paired comparison only
        block = de_matched[pid].loc[de[pid]['10']].sort_values('PValue')
        add_gene_symbols(block)
        add_fc_direction(block)
        block.to_excel(xl_writer, 'GBM%s_pair_only' % pid)
        block.to_excel(os.path.join(outdir, 'GBM%s_pair_only.xlsx' % pid))

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
        block.to_excel(os.path.join(outdir, 'GBM%s_pair_and_ref.xlsx' % pid))

        # ref comparison only
        block = de_gibco[pid].loc[de[pid]['01']].sort_values('PValue')
        add_gene_symbols(block)
        add_fc_direction(block)
        block.to_excel(xl_writer, 'GBM%s_ref_only' % pid)
        block.to_excel(os.path.join(outdir, 'GBM%s_ref_only.xlsx' % pid))

        # except Exception as exc:
        #     print "Patient %s failed." % pid
        #     print repr(exc)

    xl_writer.save()

    # Generate a series of Venn diagrams
    from plotting import venn
    figwidth = 3 * len(pids) ** 0.9
    fig, axs = plt.subplots(ncols=len(pids), nrows=3, figsize=[figwidth, 8])
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

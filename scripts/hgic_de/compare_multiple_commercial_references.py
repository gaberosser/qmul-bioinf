import multiprocessing as mp
import os
import numpy as np
from plotting import venn
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from rnaseq import filter, differential_expression
from utils import output, setops
from load_data import rnaseq_data


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_trial_2", reuse_empty=True)

    # all n=2 samples and RTK II samples
    pids = ['017', '019', '030', '031', '050', '054']
    cmap = 'RdYlGn_r'
    njob = mp.cpu_count()

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

    if njob != 1:
        pool = mp.Pool(njob)

    # Load RNA-Seq from STAR
    rnaseq_obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')

    # load additional references if required
    refs = [
        ('H1', rnaseq_data.gse38993(annotate_by='Ensembl Gene ID')),
        ('H9', rnaseq_data.gse61794(annotate_by='Ensembl Gene ID', collapse_replicates=False))
    ]
    all_refs = [t[0] for t in refs] + ['GIBCO']

    rnaseq_obj = rnaseq_data.MultipleBatchLoader([rnaseq_obj] + [t[1] for t in refs])

    # only keep gene counts
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]

    # discard iPSC
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('PSC')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('NHF1-hTERT')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('fibroblast')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('ESC')]

    rnaseq_obj.data = rnaseq_obj.data.loc[:, rnaseq_obj.meta.index]

    # compute DE between hGIC and references
    jobs = {}
    de_res = {}
    for pid in pids:
        hgic_samples = rnaseq_obj.meta.index[
            rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')
        ]
        # GBM vs matched iNSC
        the_ref = 'Paired'
        ref_type = 'iNSC'
        ref_samples = rnaseq_obj.meta.index[
            rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'iNSC')
        ]
        the_data = rnaseq_obj.data.loc[:, hgic_samples | ref_samples]
        the_groups = rnaseq_obj.meta.loc[hgic_samples | ref_samples, 'type']
        the_comparison = ['GBM', ref_type]
        if njob == 1:
            de_res[(pid, the_ref)] = differential_expression.run_one_de(the_data, the_groups, the_comparison, **de_params)
            print "GBM %s paired comparison, %d DE genes" % (pid, de_res[(pid, the_ref)].shape[0])
        else:
            jobs[(pid, the_ref)] = pool.apply_async(
                differential_expression.run_one_de,
                args=(the_data, the_groups, the_comparison),
                kwds=de_params
            )

        # GBM vs reference
        for the_ref in all_refs:
            ref_samples = rnaseq_obj.meta.index[rnaseq_obj.meta.index.str.contains(the_ref)]
            ref_type = rnaseq_obj.meta.loc[ref_samples, 'type'][0]
            the_data = rnaseq_obj.data.loc[:, hgic_samples | ref_samples]
            the_groups = rnaseq_obj.meta.loc[hgic_samples | ref_samples, 'type']
            the_comparison = ['GBM', ref_type]
            if njob == 1:
                de_res[(pid, the_ref)] = differential_expression.run_one_de(the_data, the_groups, the_comparison, **de_params)
                print "GBM %s, ref %s, %d DE genes" % (pid, the_ref, de_res[(pid, the_ref)].shape[0])
            else:
                jobs[(pid, the_ref)] = pool.apply_async(
                    differential_expression.run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=de_params
                )

    if njob != 1:
        pool.close()
        pool.join()
        for k, j in jobs.items():
            de_res[k] = j.get(1e6)
            print "GBM %s, ref %s, %d DE genes" % (k[0], k[1], de_res[k].shape[0])

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
            de_res[(pid, r)] for r in all_refs
        ]
        venn_sets, cts = setops.venn_from_arrays(*[t.index for t in the_lists])
        venn.venn3(cts, set_labels=all_refs, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'number_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'number_de_multiple_references.tiff'), dpi=200)

    # plot: how many DE genes are in the pair only comparison when each reference is used?
    pair_only = pd.DataFrame(index=pids, columns=all_refs)
    pair_only_core = dict([
        (k, pair_only.copy()) for k in range(2, len(all_refs) + 1)
    ])

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists_ref = [
            de_res[(pid, r)] for r in all_refs
        ]
        the_pair = de_res[(pid, 'Paired')]
        for ref in all_refs:
            the_ref = de_res[(pid, ref)]
            the_sets, _ = setops.venn_from_arrays(the_pair.index, the_ref.index)
            pair_only.loc[pid, ref] = the_sets['10']
        venn_sets, cts = setops.venn_from_arrays(*pair_only.loc[pid].values)
        for k in range(2, len(all_refs) + 1):
            for ref in all_refs:
                the_idx = all_refs.index(ref)
                pair_only_core[k].loc[pid, ref] = reduce(
                    lambda x, y: x + y,
                    [venn_sets[t] for t in setops.binary_combinations_sum_gte(len(all_refs), k) if t[the_idx] == '1']
                )
        venn.venn3(cts, set_labels=all_refs, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'number_po_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'number_po_de_multiple_references.tiff'), dpi=200)

    # proportion of each pair only DE count that is shared by all
    for k in range(2, len(all_refs) + 1):
        ax = (pair_only_core[k].applymap(len) / pair_only.applymap(len) * 100).plot.bar()
        ax.set_xlabel('Patient')
        # ax.set_title('Percentage of pair only DE genes that are in %d / %d reference comparisons' % (k, len(all_refs)))
        ax.set_ylim([0, 100])
        ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_%d_of_%d.png" % (k, len(all_refs))), dpi=200)
        ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_%d_of_%d.tiff" % (k, len(all_refs))), dpi=200)


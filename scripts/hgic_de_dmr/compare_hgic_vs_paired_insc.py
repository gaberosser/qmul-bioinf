import os

import numpy as np
import pandas as pd

from load_data import rnaseq_data
from plotting import venn
from rnaseq import differential_expression
from utils import output, setops

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
        de_res[pid] = differential_expression.run_one_de(the_data, the_groups, the_comparison, **de_params)
        # de_res_full[pid] = differential_expression.run_one_de(the_data, the_groups, the_comparison, return_full=True, **de_params)
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
    data = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, full_data=de_res_full)

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

    expanded_core = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=expanded_core_sets)

    # save
    expanded_core.to_excel(os.path.join(outdir, 'expanded_core_gene_list.xlsx'))

    # subgroup-specific lists (just for completeness)
    sg_specific_sets = {}
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        sg_specific_sets[grp] = k

    subgroup_specific = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=sg_specific_sets.values())
    subgroup_specific.to_excel(os.path.join(outdir, 'subgroup_specific_gene_list.xlsx'))


    # UpsetR attribute plots
    set_labels = ['019', '030', '031', '017', '050', '054']
    data_for_upset = [de_res[pid].index for pid in set_labels]  # this will be supplied to the function



    # set colours for UpsetR plot
    sets_full = {}
    sets_partial = {}
    sets_unique = []

    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        if this_k.sum() == 1:
            sets_unique.append(k)
        elif this_k.sum() == 2:
            for grp, grp_idx in subgroup_ind.items():
                if this_k[grp_idx].sum() == this_k.sum():
                    sets_partial.setdefault(grp, []).append(k)
        elif this_k.sum() == 3:
            for grp, grp_idx in subgroup_ind.items():
                if this_k[grp_idx].sum() == this_k.sum():
                    sets_full.setdefault(grp, []).append(k)

    set_colours = [
        ('RTK I full', {'sets': sets_full['RTK I'], 'colour': '#0d680f'}),
        ('RTK I partial', {'sets': sets_partial['RTK I'], 'colour': '#6ecc70'}),
        ('RTK II full', {'sets': sets_full['RTK II'], 'colour': '#820505'}),
        ('RTK II partial', {'sets': sets_partial['RTK II'], 'colour': '#d67373'}),
        ('Expanded core', {'sets': expanded_core_sets, 'colour': '#4C72B0'}),
        ('Unique', {'sets': sets_unique, 'colour': '#c3d15c'})
    ]

    # 1. Descending order
    upset1 = venn.upset_set_size_plot(
        data_for_upset,
        set_labels,
        set_colours=set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset1['figure'].savefig(os.path.join(outdir, "upset_descending.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir, "upset_descending.tiff"), dpi=200)

    # 2. Grouped by number of participants in the sets
    upset2 = venn.upset_set_size_plot(
        data_for_upset,
        set_labels,
        set_colours=set_colours,
        order_by_n_members=True,
        min_size=30,
        n_plot=50,
        default_colour = 'gray'
    )
    upset2['figure'].savefig(os.path.join(outdir, "upset_grouped_descending.png"), dpi=200)
    upset2['figure'].savefig(os.path.join(outdir, "upset_grouped_descending.tiff"), dpi=200)

    # expression heatmap to accompany it, just showing GBM (collapsed replicates)
    from plotting import clustering

    dat_gbm_aggr = pd.DataFrame(index=rnaseq_obj.data.index, columns=pids)

    # aggregate pairs
    for pid in pids:
        the_data = rnaseq_obj.data.loc[:, rnaseq_obj.data.columns.str.contains("GBM%s" % pid)]
        the_aggr = the_data.mean(axis=1)
        dat_gbm_aggr.loc[:, pid] = the_aggr

    # select genes
    g = set(venn_set['111111'])
    for x in sets_full.values() + sets_partial.values():
        for k in x:
            g.update(venn_set[k])

    the_data = dat_gbm_aggr.loc[g]

    # remove any rows that have no variation
    the_data = the_data.loc[~(the_data.diff(axis=1).iloc[:, 1:] == 0).all(axis=1)]
    the_data = np.log2(the_data + 1)

    cg = clustering.plot_clustermap(the_data, cmap='RdBu_r', vmax=12., figsize=(3.6, 8))
    cg.gs.update(bottom=0.1)
    cg.savefig(os.path.join(outdir, "clustermap_gbm_by_subgroup_gene_sets.png"), dpi=200)
    cg.savefig(os.path.join(outdir, "clustermap_gbm_by_subgroup_gene_sets.tiff"), dpi=200)


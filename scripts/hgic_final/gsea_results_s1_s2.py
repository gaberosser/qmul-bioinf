import os
from settings import HGIC_LOCAL_DIR
import pandas as pd
from rnaseq import gsea
from utils import output, setops, excel
from plotting import common
import consts
from matplotlib import pyplot as plt, patches, collections, gridspec
import seaborn as sns
import numpy as np
import collections


if __name__ == "__main__":
    # Carry out a very similar analysis to that run in ipa_results_s1_s2, including the same plot.

    # set a minimum pval for pathways to be used
    alpha = 0.001
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.1
    # small offset to avoid zero FDR
    # set this slightly lower than the smallest non-zero FDR value
    eps = 1e-6

    outdir = output.unique_output_dir()

    indir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/gsea/results/raw'
    )
    msigdb_c5_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/gsea/msigdb/c5.all.v6.2.symbols.gmt'
    )
    pids = consts.PIDS
    comparison_names = collections.OrderedDict([
        ('', 'syngeneic'),
        ('_h9_nsc', 'h9',),
        ('_gibco_nsc', 'gibco')
    ])

    # load C5 (all GO term pathways) for filtering and network analysis
    c5_gmt = gsea.read_gmt_file(msigdb_c5_fn)

    keep_pathways = c5_gmt.keys()

    res = collections.OrderedDict()
    res_full = collections.OrderedDict()
    for pid in pids:
        for c in comparison_names:
            fn = os.path.join(indir, "%s%s.csv" % (pid, c))
            this = pd.read_csv(fn, sep='\t', header=0, index_col=0, usecols=[0, 3, 5, 7])
            this.columns = ['n_gene', 'nes', 'fdr']
            this = this.reindex(keep_pathways).dropna(how='all')
            res_full["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha_relevant]
            res["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha]

    pathways_sign = sorted(setops.reduce_union(*[t.index for t in res.values()]))
    pathways_rele = sorted(setops.reduce_union(*[t.index for t in res_full.values()]))

    excel.pandas_to_excel(res, os.path.join(outdir, "gsea_results_significant_by_patient.xlsx"))

    # use this list to export a second wideform Excel file with the top list of pathways
    for_export = pd.DataFrame(index=pathways_sign, columns=['n_gene'])
    nes_columns = []
    fdr_columns = []
    for k, v in res.items():
        for_export.loc[v.index, 'n_gene'] = v.n_gene
        this_yn = pd.Series('N', index=pathways_sign)
        this_yn.loc[v.index] = 'Y'
        for_export.insert(
            for_export.shape[1],
            k,
            this_yn
        )
        for_export.insert(
            for_export.shape[1],
            "%s_nes" % k,
            res_full[k].reindex(pathways_sign)['nes']
        )
        nes_columns.append("%s_nes" % k)
        for_export.insert(
            for_export.shape[1],
            "%s_fdr" % k,
            res_full[k].reindex(pathways_sign)['fdr']
        )
        fdr_columns.append("%s_fdr" % k)
    for_export.to_excel(os.path.join(outdir, "gsea_results_significant.xlsx"))

    # extract fdr data only for plotting
    fdr_dat = for_export[fdr_columns]
    fdr_dat = fdr_dat[fdr_dat < alpha]
    fdr_dat.columns = res.keys()  # dict is ordered, so this is OK
    all_in = ~fdr_dat.isnull()
    all_in.columns = res.keys()  # dict is ordered, so this is OK
    log_fdr_dat = np.log10(fdr_dat + 1e-6) * -1


    # number syngen. only, ref. only and intersection
    n_set = pd.DataFrame(0, index=pathways_sign, columns=['Syngen. only', 'Ref. only', 'Intersect.'], dtype=int)
    so = {}
    ro = {}
    inters = {}
    for pid in pids:
        s = all_in.index[all_in.loc[:, "%s_syngeneic" % pid]]
        r = all_in.index[all_in.loc[:, ["%s_%s" % (pid, t) for t in comparison_names.values()[1:]]].any(axis=1)]
        vs, _ = setops.venn_from_arrays(s, r)
        n_set.loc[vs['10'], 'Syngen. only'] += 1
        n_set.loc[vs['01'], 'Ref. only'] += 1
        n_set.loc[vs['11'], 'Intersect.'] += 1
    n_set = n_set.fillna(0)

    from ipa_results_s1_s2 import pathway_involvement_heatmap_by_p
    comparison_dict = {
        'syngeneic': 'Syngen.',
        'h9': 'H9',
        'gibco': 'Gibco'
    }

    # plot 1) P values, ordered by sum of -log(P)
    p_order = log_fdr_dat.sum(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        log_fdr_dat,
        n_set,
        p_order,
        pids,
        comparison_dict
    )
    ax = plot_dict['axs'][0]
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=4)
    eax = plot_dict['eax']
    plt.setp([t for t in eax.get_children() if isinstance(t, plt.Text)], fontsize=4)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.pdf"), dpi=200)

    # plot 2) P values, ordered by mean of -log(P)
    p_order = log_fdr_dat.mean(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        log_fdr_dat,
        n_set,
        p_order,
        pids,
        comparison_dict
    )
    ax = plot_dict['axs'][0]
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=4)
    eax = plot_dict['eax']
    plt.setp([t for t in eax.get_children() if isinstance(t, plt.Text)], fontsize=4)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.pdf"), dpi=200)

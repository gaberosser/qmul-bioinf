import collections
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import MultiComparison

from hgic_consts import NH_ID_TO_PATIENT_ID_MAP
from plotting import common
from rnaseq import loader, gsva
from scripts.anaelle import tcga_tam_vs_mtor as ttm
from scripts.hgic_final import consts
from utils import output, log, reference_genomes

logger = log.get_console_logger()


def nh_id_to_patient_id(arr):
    nh_id = pd.Index(arr).str.replace(r'(_?)(DEF|SP).*', '')
    return [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]


if __name__ == "__main__":
    """
    Use the BTR bulk tumour RNA-Seq data to validate a link between the mTOR pathway and the proportion of microglial
    and macrophage immune infiltrate in the bulk samples.

    mTOR is assessed using a known set of genes.

    Tumour-associated bone marrow-derived macrophages (TAM-BMDM) and microglia (TAM-MG) are distinguished using a
    human signature from Muller et al. (Genome Biol 2017) or a converted mouse signature from Bowman et al. (???).
    """
    # tam_signature_source = 'bowman'
    tam_signature_source = 'muller'
    mtor_source = 'kegg_msigdb'  # ('kegg', 'pid', 'biocarta')

    # significance cutoff
    alpha = 0.05

    # load mTOR signatures
    mtor_gs_dict = ttm.mtor_signature_dict()

    # load MG/BMDM signatures
    tam_gs_dict = ttm.tam_signature_dict()

    mtor_geneset = mtor_gs_dict[mtor_source]
    tam_genesets = tam_gs_dict[tam_signature_source]

    genesets = dict(tam_genesets)
    genesets['mTOR'] = mtor_geneset

    subgroups = consts.SUBGROUPS

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    outdir = output.unique_output_dir()

    obj = loader.load_by_patient(consts.PIDS, type='ffpe', source='salmon', include_control=False)
    obj.filter_by_sample_name(consts.FFPE_RNASEQ_SAMPLES)
    obj.meta.insert(0, 'patient_id', nh_id_to_patient_id(obj.meta.index))
    obj.meta.insert(1, 'subgroup', [subgroups_lookup[pid] for pid in obj.meta.patient_id])

    rnaseq_dat = obj.data.copy()
    # use gene symbol identifiers
    gs = reference_genomes.ensembl_to_gene_symbol(rnaseq_dat.index).dropna()
    rnaseq_dat = rnaseq_dat.loc[gs.index]
    rnaseq_dat.index = gs.values

    groups = obj.meta.subgroup
    group_list = groups.unique()

    # check that signature genes are all found in the data
    for k, v in genesets.items():
        for i, t in enumerate(v):
            if t in ttm.manual_gene_name_correction:
                v[i] = ttm.manual_gene_name_correction[t]
        g_in = rnaseq_dat.index.intersection(v)
        if set(g_in) != set(v):
            missing = set(v).difference(rnaseq_dat.index)
            logger.warn(
                "%d genes in the %s signature do not match with the data index and will be dropped: %s.",
                len(missing),
                k,
                ', '.join(missing)
            )
            genesets[k] = g_in

    # run ssGSEA then Z transform the results
    es = gsva.ssgsea(rnaseq_dat, genesets)
    es_z = ttm.z_transform(es, axis=1)

    # export
    for_export = es_z.transpose()
    for_export.insert(for_export.shape[1], 'Patient ID', obj.meta.patient_id)
    for_export.insert(for_export.shape[1], 'Subgroup', obj.meta.subgroup)
    for_export.to_excel(os.path.join(outdir, "signature_scores_and_subgroups.xlsx"))


    # boxplots showing signature scores in the different niches
    bplot = {}
    anova_res = {}
    tukey_res = {}
    for k in genesets:
        the_data = es_z.loc[k]
        bplot[k] = collections.OrderedDict()
        for sg in group_list:
            bplot[k][sg] = the_data.loc[groups.fillna('').str.contains(sg)].values

        anova_res[k] = stats.f_oneway(*bplot[k].values())
        mc = MultiComparison(the_data, groups, group_order=group_list)
        tukey_res[k] = mc.tukeyhsd(alpha=alpha)

        lbl, tmp = zip(*bplot[k].items())
        tmp = [list(t) for t in tmp]
        fig = plt.figure(num=k, figsize=(5, 4))
        ax = fig.add_subplot(111)
        sns.boxplot(data=tmp, orient='v', ax=ax, color='0.5')
        ax.set_xticklabels(lbl, rotation=45)
        ax.set_ylabel("Normalised ssGSEA score")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.png' % k.lower()), dpi=200)
        fig.savefig(os.path.join(outdir, '%s_ssgsea_by_subgroup_tcga.pdf' % k.lower()))

    # can annotate these manually based on statistics?

    # scatterplots showing correlation in the different niches
    dict_mg = ttm.scatter_plot_with_linregress(es_z.loc['mTOR'], es_z.loc['MG'], group_list, groups)
    dict_mg['fig'].set_size_inches((11, 5.5))
    dict_mg['fig'].subplots_adjust(left=0.05, right=0.98, top=0.95, hspace=0.25)
    dict_mg['fig'].savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_mg['fig'].savefig(os.path.join(outdir, "mtor_vs_mg_correlation_by_tcga_subgroup.pdf"))

    dict_bmdm = ttm.scatter_plot_with_linregress(es_z.loc['mTOR'], es_z.loc['BMDM'], group_list, groups)
    dict_bmdm['fig'].set_size_inches((11, 5.5))
    dict_bmdm['fig'].subplots_adjust(left=0.05, right=0.98, top=0.95, hspace=0.25)
    dict_bmdm['fig'].savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_bmdm['fig'].savefig(os.path.join(outdir, "mtor_vs_bmdm_correlation_by_tcga_subgroup.pdf"))

    # check for MG / BMDM correlation
    dict_both = ttm.scatter_plot_with_linregress(es_z.loc['MG'], es_z.loc['BMDM'], group_list, groups)
    dict_both['fig'].set_size_inches((11, 5.5))
    dict_both['fig'].subplots_adjust(left=0.05, right=0.98, top=0.95, hspace=0.25)
    dict_both['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_subgroup.png"), dpi=300)
    dict_both['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation_by_tcga_subgroup.pdf"))

    # again but across groups
    dict_both_uniform = ttm.scatter_plot_with_linregress(es_z.loc['MG'], es_z.loc['BMDM'])
    dict_both_uniform['fig'].set_size_inches([6, 4])
    dict_both_uniform['fig'].tight_layout()
    dict_both_uniform['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation.png"), dpi=200)
    dict_both_uniform['fig'].savefig(os.path.join(outdir, "mg_vs_bmdm_correlation.pdf"))

    # summary plot with all information
    slope_cmap = plt.get_cmap('RdBu_r')
    slope_norm = common.MidpointNormalize(vmin=-.5, vmax=1.2, midpoint=0.)
    slope_sm = plt.cm.ScalarMappable(cmap=slope_cmap, norm=slope_norm)

    group_list_extended = list(group_list) + ['All']

    for_plot = collections.OrderedDict([
        # ('MG-BMDM', dict_both),
        ('mTOR-BMDM', dict_bmdm),
        ('mTOR-MG', dict_mg)
    ])

    s, p = ttm.get_slope_and_pval(
        for_plot,
        col_order=group_list_extended,
    )

    # p_to_size = lambda t: min(150., 45 - 12 * np.log10(t))
    p_to_size = lambda t: min(500., 100 + 20 * np.log10(t) ** 2)

    x = range(len(group_list_extended))
    y_fun = lambda t: [t] * len(group_list_extended)

    fig, ax = plt.subplots(figsize=(6, 2.4))
    for i, k in enumerate(s.index):
        ax.scatter(
            x,
            y_fun(i),
            color=[slope_sm.to_rgba(t) for t in s.loc[k]],
            s=[p_to_size(t) for t in p.loc[k]],
            edgecolor='k',
            linewidth=[.5 if t > alpha else 1.5 for t in p.loc[k]]
        )
    ax.grid('off')
    ax.set_facecolor('w')
    ax.set_xticks(x)
    ax.set_xticklabels(p.columns)
    ax.set_yticks(range(p.shape[0]))
    ax.set_yticklabels(p.index)
    ax.set_ylim([-0.5, p.shape[0] - 0.5])

    slope_sm.set_array(s.values)
    cbar = fig.colorbar(slope_sm)
    cbar.set_label('Slope')
    fig.tight_layout()

    fig.savefig(os.path.join(outdir, "correlation_summary.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "correlation_summary.pdf"))

    plt.close('all')

    # alternative: line plot
    plt_dict = ttm.line_plot_pvalues_slope(p, s, alpha=alpha)
    plt_dict['fig'].set_size_inches([6., 2.5])
    plt_dict['gs'].update(bottom=0.23, top=0.9, hspace=0.4, right=0.75)
    plt_dict['fig'].savefig(os.path.join(outdir, "correlation_summary_line.png"), dpi=200)
    plt_dict['fig'].savefig(os.path.join(outdir, "correlation_summary_line.pdf"))

    # bespoke analysis using individual genes from the mTOR pathway rather than the whole signature
    mtor_genes = [
        'RPS6',
        'EIF4EBP1',
        'VEGFA',
        'VEGFB',
        'VEGFC',
        'VEGFD',
        'IL10',
        'IL6',
    ]

    ix = reduce(lambda x, y: x + y, [['%s_%s' % (t, g) for g in mtor_genes] for t in ['MG', 'BMDM']])
    s_g = pd.DataFrame(index=ix, columns=group_list_extended)
    p_g = s_g.copy()

    for g in mtor_genes:
        for t in ['MG', 'BMDM']:
            the_dict = ttm.scatter_plot_with_linregress(
                es_z.loc[t],
                np.log10(rnaseq_dat.loc[g] + 0.01),
                group_list,
                groups
            )
            the_dict['fig'].set_size_inches((11, 5.5))
            the_dict['fig'].tight_layout()
            the_dict['fig'].savefig(os.path.join(outdir, "%s_%s_scatterplots.png" % (t, g)), dpi=200)

            s_g.loc['%s_%s' % (t, g)] = [the_dict['statsmodels'][k].params[-1] for k in group_list_extended]
            p_g.loc['%s_%s' % (t, g)] = [the_dict['statsmodels'][k].f_pvalue for k in group_list_extended]

    x = range(len(group_list_extended))
    y_fun = lambda t: [t] * len(group_list_extended)

    for k in ['MG', 'BMDM']:
        fig, ax = plt.subplots(figsize=(6, 3.4))
        for i, g in enumerate(mtor_genes):
            the_key = '%s_%s' % (k, g)
            ax.scatter(
                x,
                y_fun(i),
                color=[slope_sm.to_rgba(t) for t in s_g.loc[the_key]],
                s=[p_to_size(t) for t in p_g.loc[the_key]],
                edgecolor='k',
                linewidth=[.5 if t > alpha else 1.5 for t in p_g.loc[the_key]]
            )

        ax.grid('off')
        ax.set_facecolor('w')
        ax.set_xticks(x)
        ax.set_xticklabels(group_list_extended)

        ix = ['%s_%s' % (k, g) for g in mtor_genes]

        ax.set_yticks(range(len(ix)))
        ax.set_yticklabels(ix)

        slope_sm.set_array(s_g.values)
        cbar = fig.colorbar(slope_sm)
        cbar.set_label('Slope')
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, "correlation_summary_individual_genes_%s.png" % k), dpi=200)
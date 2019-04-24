from rnaseq import gsva, loader
import pandas as pd
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR, DATA_DIR_NON_GIT
from plotting import venn

import os
import re
import references
import datetime
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.multicomp import MultiComparison

import seaborn as sns
import numpy as np
import collections
from scipy import stats
from plotting import common
from utils import setops, output, log
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP

from scripts.anaelle import tcga_tam_vs_mtor as ttm

logger = log.get_console_logger()


if __name__ == "__main__":
    """
    Use the IVY GAP databaseto validate a link between the mTOR pathway and the proportion of microglial and macrophage
    immune infiltrate in the bulk samples.

    This includes FFPE bulk RNA-Seq data of multi-site-sampled GBM tumours. The tumour niche is recorded in all cases.

    mTOR is assessed using a known set of genes.

    Tumour-associated bone marrow-derived macrophages (TAM-BMDM) and microglia (TAM-MG) are distinguished using a
    human signature from Muller et al. (Genome Biol 2017) or a converted mouse signature from Bowman et al. (???).
    """
    rnaseq_units = 'fpkm_norm'  # ('fpkm', 'tpm')
    # tam_signature_source = 'bowman'
    tam_signature_source = 'muller'
    mtor_source = 'kegg_msigdb'  # ('kegg', 'pid', 'biocarta')

    # significance cutoff
    alpha = 0.01

    # load mTOR signatures
    mtor_gs_dict = ttm.mtor_signature_dict()

    # load MG/BMDM signatures
    tam_gs_dict = ttm.tam_signature_dict()

    mtor_geneset = mtor_gs_dict[mtor_source]
    tam_genesets = tam_gs_dict[tam_signature_source]

    genesets = dict(tam_genesets)
    genesets['mTOR'] = mtor_geneset

    outdir = output.unique_output_dir()

    obj = loader.ivygap(units=rnaseq_units)
    rnaseq_dat = obj.data

    # define tumour niches
    groups = obj.meta.structure_abbreviation.replace(re.compile(r'-.*'), '')
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

    s = pd.DataFrame(index=['MG-BMDM', 'mTOR-BMDM', 'mTOR-MG'], columns=group_list)
    p = s.copy()

    s.loc['MG-BMDM'] = [dict_both['statsmodels'][k].params[-1] for k in group_list]
    s.loc['mTOR-BMDM'] = [dict_bmdm['statsmodels'][k].params[-1] for k in group_list]
    s.loc['mTOR-MG'] = [dict_mg['statsmodels'][k].params[-1] for k in group_list]

    p.loc['MG-BMDM'] = [dict_both['statsmodels'][k].f_pvalue for k in group_list]
    p.loc['mTOR-BMDM'] = [dict_bmdm['statsmodels'][k].f_pvalue for k in group_list]
    p.loc['mTOR-MG'] = [dict_mg['statsmodels'][k].f_pvalue for k in group_list]

    p_to_size = lambda t: min(150., 45 - 12 * np.log10(t))

    x = range(len(group_list))
    y_fun = lambda t: [t] * len(group_list)

    fig, ax = plt.subplots(figsize=(6, 2.4))
    ax.scatter(
        x,
        y_fun(0),
        color=[slope_sm.to_rgba(t) for t in s.loc['MG-BMDM']],
        s=[p_to_size(t) for t in p.loc['MG-BMDM']],
        edgecolor='k',
        linewidth=[.5 if t > alpha else 1.5 for t in p.loc['MG-BMDM']]
    )
    ax.scatter(
        x,
        y_fun(1),
        color=[slope_sm.to_rgba(t) for t in s.loc['mTOR-BMDM']],
        s=[p_to_size(t) for t in p.loc['mTOR-BMDM']],
        edgecolor='k',
        linewidth=[.5 if t > alpha else 1.5 for t in p.loc['mTOR-BMDM']]
    )
    ax.scatter(
        x,
        y_fun(2),
        color=[slope_sm.to_rgba(t) for t in s.loc['mTOR-MG']],
        s=[p_to_size(t) for t in p.loc['mTOR-MG']],
        edgecolor='k',
        linewidth=[.5 if t > alpha else 1.5 for t in p.loc['mTOR-MG']]
    )
    ax.grid('off')
    ax.set_facecolor('w')
    ax.set_xticks(x)
    ax.set_xticklabels(group_list)
    ax.set_yticks(range(p.shape[0]))
    ax.set_yticklabels(p.index)

    slope_sm.set_array(s.values)
    cbar = fig.colorbar(slope_sm)
    cbar.set_label('Slope')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "correlation_summary.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "correlation_summary.pdf"))
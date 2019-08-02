import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from plotting import common
from rnaseq import loader
from scripts.hgic_final import consts
from settings import HGIC_LOCAL_DIR
from utils import output, setops, dictionary, reference_genomes

if __name__ == "__main__":
    pids = consts.PIDS
    eps = 0.01

    outdir = output.unique_output_dir()
    obj = loader.load_by_patient(consts.PIDS, source='salmon', include_control=False)
    obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    meta = obj.meta
    meta.insert(0, 'patient_id', meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # load genes in GAG-related IPA pathways
    pathways = [
        'Chondroitin Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis (Late Stages)',
        'Chondroitin Sulfate Biosynthesis (Late Stages)',
    ]

    # pathways by DE
    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq',
        's0_individual_patients_direct_comparison',
        'ipa',
        'pathways',
        'full_de_all.xls'
    )
    # result is a dict keyed by PID
    res = pd.read_excel(fn, sheet_name=None, index_col=0)

    # get the list of genes for each pathway
    pathway_genes = {}
    for p in pathways:
        pathway_genes[p] = set()
        for this in res.values():
            if p in this.index:
                pathway_genes[p].update(this.loc[p, 'genes'].split(','))

    goi = sorted(setops.reduce_union(*pathway_genes.values()))

    # goi = [
    #     "SULT1E1",
    #     "CHST1",
    #     "CHST2",
    #     "CHST3",
    #     "CHST4",
    #     "CHST6",
    #     "CHST10",
    #     "CHST11",
    #     "CHST14",
    #     "NDST2",
    #     "HAS1",
    #     "HAS2",
    #     "HAS3"
    # ]

    ee = reference_genomes.gene_symbol_to_ensembl(goi)

    dat = obj.data.loc[ee]
    dat.index = ee.index
    dat_mean = np.log10(dat.groupby(by=[meta.patient_id, meta.type], axis=1).mean() + eps)
    # dat_mean = dat.groupby(by=[meta.patient_id, meta.type], axis=1).mean()

    # the indexing with a MultiIndex is an unbelievable hassle
    # let's just extract two separate dfs
    dat_mean_gic = dat_mean.loc[:, (pids, 'GBM')]
    dat_mean_gic.columns = dat_mean_gic.columns.levels[0]
    dat_mean_gic = dat_mean_gic[pids]

    dat_mean_insc = dat_mean.loc[:, (pids, 'iNSC')]
    dat_mean_insc.columns = dat_mean_insc.columns.levels[0]
    dat_mean_insc = dat_mean_insc[pids]

    the_dat = pd.concat((dat_mean_insc, dat_mean_gic), axis=1)
    col_colours = pd.Series(['0.8'] * len(pids) + ['0.2'] * len(pids), index=the_dat.columns)

    cg = sns.clustermap(
        the_dat,
        z_score=0,
        col_cluster=False,
        figsize=(6.5, 6.5),
        col_colors=['0.9'] * len(pids) + ['0.6'] * len(pids),
    )
    ax = cg.ax_heatmap
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    cg.cax.set_title('Row-normalised\nZ value', horizontalalignment='left')
    cg.gs.update(left=0.06, top=0.92, right=0.82)

    cg.ax_col_colors.text(0.25, 0.5, 'iNSC',
         horizontalalignment='center',
         verticalalignment='center',
         transform=cg.ax_col_colors.transAxes)

    cg.ax_col_colors.text(0.75, 0.5, 'GIC',
         horizontalalignment='center',
         verticalalignment='center',
         transform=cg.ax_col_colors.transAxes)

    cg.savefig(os.path.join(outdir, "gag_genes_logtpm_clustermap.png"), dpi=200)
    cg.savefig(os.path.join(outdir, "gag_genes_logtpm_clustermap.tiff"), dpi=200)

    # same again, buut disable row clustering and provide our own function-based groups
    function_to_gene = {
        'Heparan sulphate sulfotransferase': the_dat.index[the_dat.index.str.contains(r'^HS[36]')].tolist(),
        'Bifunctional heparan deacetylase / sulfotransferase': the_dat.index[the_dat.index.str.contains(r'^NDST')].tolist(),
        'Sulfotransferase': the_dat.index[the_dat.index.str.contains(r'^SULT')].tolist() + ['UST', 'CHSY3'],
        'Glucuronosyltransferase': ['B3GAT1', 'B3GAT2', 'CHPF', 'CHPF2', 'CSGALNACT1'],
        'Carbohydrate transferase': the_dat.index[the_dat.index.str.contains(r'^CHST')].tolist(),
        'Dermatan sulfate epimerase': the_dat.index[the_dat.index.str.contains(r'^DSE')].tolist(),
    }
    gene_to_function = dictionary.complement_dictionary_of_iterables(function_to_gene, squeeze=True)
    function_colours = dict(zip(function_to_gene.keys(), common.get_best_cmap(len(function_to_gene))))

    # reorder data
    the_dat = the_dat.loc[reduce(lambda x, y: x + y, function_to_gene.values())]
    row_colours = pd.DataFrame(
        [function_colours[gene_to_function[t]] for t in the_dat.index],
        index=the_dat.index,
        columns=['Function']
    )

    # standardise (Z)
    z = the_dat.subtract(the_dat.mean(axis=1), axis=0).divide(the_dat.std(axis=1), axis=0)

    cg = sns.clustermap(
        z,
        col_cluster=False,
        row_cluster=False,
        figsize=(6.5, 6.5),
        col_colors=['0.9'] * len(pids) + ['0.6'] * len(pids),
        row_colors=row_colours
    )
    ax = cg.ax_heatmap
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    cg.cax.set_title('Row-normalised\nZ value', horizontalalignment='left')
    cg.gs.update(left=0.06, top=0.92, right=0.82)

    cg.ax_col_colors.text(0.25, 0.5, 'iNSC',
         horizontalalignment='center',
         verticalalignment='center',
         transform=cg.ax_col_colors.transAxes)

    cg.ax_col_colors.text(0.75, 0.5, 'GIC',
         horizontalalignment='center',
         verticalalignment='center',
         transform=cg.ax_col_colors.transAxes)

    # add custom legend
    leg_ax = cg.ax_col_dendrogram
    leg_dict = dict([
        (k, {'class': 'patch', 'facecolor': function_colours[k], 'edgecolor': 'k', 'linewidth': 1.}) for k in function_to_gene
    ])
    leg_dict = {'': leg_dict}
    leg = common.add_custom_legend(leg_ax, legend_dict=leg_dict)
    leg_ax.legend(loc='upper center', bbox_to_anchor=(0.6, 1.7), handles=leg)

    cg.savefig(os.path.join(outdir, "gag_genes_logtpm_clustermap_row_by_function.png"), dpi=200)
    cg.savefig(os.path.join(outdir, "gag_genes_logtpm_clustermap_row_by_function.tiff"), dpi=200)


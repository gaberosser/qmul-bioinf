from rnaseq import loader
from scripts.hgic_final import consts
from matplotlib import pyplot as plt, gridspec
import seaborn as sns
import numpy as np
import os
import pandas as pd

from settings import HGIC_LOCAL_DIR
from utils import output, setops
import references

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

    ee = references.gene_symbol_to_ensembl(goi)

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

    gs = gridspec.GridSpec(5, 3, width_ratios=[10, 10, 1])
    fig = plt.figure(figsize=(5.5, 5.))
    ax_insc = fig.add_subplot(gs[:, 0])
    ax_gic = fig.add_subplot(gs[:, 1])

    cax = fig.add_subplot(gs[1:-1, 2])

    # we need to fix the colour scale in advance
    vmin = -2
    # vmin = None
    vmax = 2
    # vmax = 50

    sns.heatmap(
        dat_mean_insc,
        vmin=vmin,
        vmax=vmax,
        cmap='Reds',
        cbar=False,
        ax=ax_insc,
        yticklabels=True
    )

    sns.heatmap(
        dat_mean_gic,
        vmin=vmin,
        vmax=vmax,
        cmap='Reds',
        cbar=True,
        cbar_ax=cax,
        ax=ax_gic,
        yticklabels=False
    )

    plt.setp(ax_insc.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax_insc.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax_gic.xaxis.get_ticklabels(), rotation=90)

    ax_gic.yaxis.label.set_visible(False)
    ax_gic.xaxis.label.set_visible(False)
    ax_insc.yaxis.label.set_visible(False)
    ax_insc.xaxis.label.set_visible(False)

    ax_insc.set_title('iNSC')
    ax_gic.set_title('GIC')
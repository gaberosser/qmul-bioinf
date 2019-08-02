import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from rnaseq import loader, general
from scripts.hgic_final import consts
from utils import output


def chunker(arr, chunk_size):
    arr = np.array(arr)
    for i in range(0, len(arr), chunk_size):
        yield arr[i:i + chunk_size]


if __name__ == "__main__":
    """
    Generate a series of boxplots for genes of interest that have been identified from the boxplots
    Do this for both syngeneic and reference comparisons.
    """
    pids = consts.PIDS

    # genes of interest
    gois = [
        'C1orf61',
        'RAMP1',
        'SCG2',
        'PMP2',
        'METTL7B',
        'MLC1',
        'ST8SIA5',
        'CRABP1',
        'MSX2',
        'KRT19',
        'KRT8',
        'KRT18',
        'NFIX'
    ]

    row_per_fig = 5

    # load TPM data
    cols_syn_gic = consts.S1_RNASEQ_SAMPLES_GIC
    cols_syn_insc = consts.S1_RNASEQ_SAMPLES_INSC
    cols_ref_nsc = [
        'GIBCO_NSC_P4',
        'H9_NSC_1',
        'H9_NSC_2'
    ]

    outdir = output.unique_output_dir()

    obj1 = loader.load_by_patient(pids, source='salmon', include_control=True)
    obj2 = loader.load_references('GSE61794', source='salmon')
    obj = loader.MultipleBatchLoader([obj1, obj2])
    obj.filter_by_sample_name(consts.ALL_RNASEQ_SAMPLES)

    dat = obj.data.copy()
    general.add_gene_symbols_to_ensembl_data(dat)

    # check that all GOIs are present
    vc = dat['Gene Symbol'].value_counts()
    for g in gois:
        if g not in vc:
            raise KeyError("Gene %s was not found." % g)
        if vc[g] != 1:
            raise AttributeError("Gene %s has %d hits." % (g, vc[g]))

    # create plots in multiple figures
    # syngeneic
    the_dat = dat.loc[dat['Gene Symbol'].isin(gois), cols_syn_gic + cols_syn_insc].transpose()
    the_dat.columns = dat.loc[the_dat.columns, 'Gene Symbol']
    the_dat = np.log10(the_dat + 0.01)

    the_type = obj.meta.loc[cols_syn_gic + cols_syn_insc, 'type']
    the_dat.insert(0, 'Cell type', the_type)

    fig, axs = plt.subplots(nrows=len(gois), ncols=1, sharex=True, sharey=True, figsize=(4.5, 8.))
    for i, g in enumerate(gois):
        ax = axs[i]
        this_dat = pd.melt(the_dat.loc[:, ['Cell type', g]], id_vars='Cell type')
        sns.boxplot(x='value', y='Cell type', data=this_dat, orient='h', ax=ax, fliersize=0)
        sns.swarmplot(x='value', y='Cell type', data=this_dat, orient='h', ax=ax, edgecolor='k', linewidth=1., alpha=0.5)
        ax.set_ylabel(g, fontsize=9.)
        if g != gois[-1]:
            ax.xaxis.label.set_visible(False)
        else:
            ax.set_xlabel(r'$\log_{10}(\mathrm{TPM} + 0.01)$')
    axs[-1].set_ylim([-.6, 1.6])
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.07, top=0.98, hspace=0.08)
    fig.savefig(os.path.join(outdir, "syngeneic_boxplots.png"), dpi=200)

    # reference
    palette = {
        'GIBCO': 'y',
        'H9': 'forestgreen',
        '': 'royalblue'
    }
    the_dat = dat.loc[dat['Gene Symbol'].isin(gois), cols_syn_gic + cols_ref_nsc].transpose()
    the_dat.columns = dat.loc[the_dat.columns, 'Gene Symbol']
    the_dat = np.log10(the_dat + 0.01)

    the_type = obj.meta.loc[cols_syn_gic + cols_ref_nsc, 'type']
    the_batch = pd.Series('', index=the_type.index)
    the_batch.loc[the_batch.index.str.contains('GIBCO')] = 'GIBCO'
    the_batch.loc[the_batch.index.str.contains('H9')] = 'H9'

    the_dat.insert(0, 'Cell type', the_type)
    the_dat.insert(0, 'Batch', the_batch)

    fig, axs = plt.subplots(nrows=len(gois), ncols=1, sharex=True, sharey=True, figsize=(4.5, 8.))
    for i, g in enumerate(gois):
        ax = axs[i]
        this_dat = pd.melt(the_dat.loc[:, ['Cell type', 'Batch', g]], id_vars=['Cell type', 'Batch'])
        sns.boxplot(x='value', y='Cell type', data=this_dat, orient='h', ax=ax, fliersize=0)
        sns.swarmplot(x='value', y='Cell type', data=this_dat, hue='Batch', orient='h', ax=ax, edgecolor='k', linewidth=1., alpha=0.5, palette=palette)
        ax.set_ylabel(g, fontsize=9.)
        if g != gois[-1]:
            ax.xaxis.label.set_visible(False)
        else:
            ax.set_xlabel(r'$\log_{10}(\mathrm{TPM} + 0.01)$')
        ax.get_legend().set_visible(False)
    axs[-1].set_ylim([-.6, 1.6])
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.07, top=0.98, hspace=0.08)
    fig.savefig(os.path.join(outdir, "reference_boxplots.png"), dpi=200)
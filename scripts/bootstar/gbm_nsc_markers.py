import os
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import ticker
import seaborn as sns
from scipy.cluster import hierarchy

from load_data import rnaseq_data
from microarray import process
from utils.output import unique_output_dir
from plotting import clustering, bar, heatmap

from scripts.rnaseq import gtf_reader

import references


if __name__ == '__main__':
    gene_lengths = {
        'PDGFRA': 6576,
        'SLC1A3': 4170,
    }

    OUTDIR = unique_output_dir("jb.marker_levels", reuse_empty=True)

    # GSE73721 (reference astrocytes, oligos, ...)
    obj73721 = rnaseq_data.gse73721(source='star', annotate_by='Ensembl Gene ID')

    # remove unneeded samples
    to_keep73721 = (
        obj73721.data.columns.str.contains('yo ctx astro')
        | obj73721.data.columns.str.contains('Hippocampus astro')
        | obj73721.data.columns.str.contains('oligo')
    )

    # GSE61794 (H9-derived NSC x 2)
    obj61794 = rnaseq_data.gse61794(source='star', annotate_by='Ensembl Gene ID')
    # combining replicates
    rc = obj61794.meta.read_count.sum()
    obj61794.meta = pd.DataFrame(
        data={
            'cell_type': 'NSC',
            'srr': 'SRR1586371-2',
            'read_count': rc,
            'sample': 'H9 NSC',
        },
        index=['SRR1586371-2']
    )
    obj61794.data = pd.DataFrame(obj61794.data.sum(axis=1), columns=['H9 NSC'])

    # WTCHG ALL samples
    objwtchg_all = rnaseq_data.all_hgic_loader(annotate_by='Ensembl Gene ID', include_derived=True)
    to_keep_wtchg = (
        'GIBCO_NSC_P4',
        'DURA018_NSC_N2_P6',
        'DURA018_NSC_N4_P4',
        'DURA019_NSC_N8C_P2',
        'DURA030_NSC_N16B6_P1',
        'DURA031_NSC_N44B_P2'
    )


    # rRNA gene IDs
    rrna_ensg = set(gtf_reader.get_rrna())

    # MT gene_ids
    mt_ensg = set(gtf_reader.get_mitochondrial())

    # combine the data
    data = pd.concat((
        obj73721.data.loc[:, to_keep73721],
        obj61794.data,
        objwtchg_all.data.loc[:, to_keep_wtchg]
    ), axis=1)

    # combine the metadata
    meta = pd.concat((
        obj73721.meta.loc[to_keep73721],
        obj61794.meta,
        objwtchg_all.meta.loc[to_keep_wtchg, :],
    ), axis=0)

    data = data.loc[data.index.str.contains('ENSG')]

    # remove rRNA
    data = data.loc[~data.index.isin(rrna_ensg)]

    # remove MT RNA
    data = data.loc[~data.index.isin(mt_ensg)]

    # normalise by read counts: per million
    data_n = data.divide(data.sum(axis=0), axis=1)

    # extract genes of interest
    genes = references.gene_symbol_to_ensembl(['SLC1A3', 'PDGFRA'])

    this_fpkm = data_n.loc[genes] * 1e6
    this_fpkm.index = genes.index

    # normalise by gene length: per million per kilobase
    for g in genes.index:
        this_fpkm.loc[g] = this_fpkm.loc[g] / gene_lengths[g] * 1e3

    ax = this_fpkm.transpose().plot.bar()
    ax.set_ylabel('FPKM')
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(OUTDIR, 'all_markers_fpkm.png'), dpi=200)

    # now just our own samples and the reference
    this_fpkm_lim = this_fpkm.loc[:, this_fpkm.columns.str.contains('NSC')]
    ax = this_fpkm_lim.transpose().plot.bar()
    ax.set_ylabel('FPKM')
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(OUTDIR, 'all_markers_fpkm2.png'), dpi=200)
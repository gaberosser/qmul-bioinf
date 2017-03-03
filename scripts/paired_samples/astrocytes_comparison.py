import os
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

from load_data import rnaseq_data
from microarray import process
from utils.output import unique_output_dir
from plotting import clustering

if __name__ == "__main__":
    N_GENES = 1500
    OUTDIR = unique_output_dir("astrocytes", reuse_empty=True)

    # GSE73721 (reference astrocytes, oligos, ...)
    obj73721 = rnaseq_data.gse73721(source='star', annotate_by='Ensembl Gene ID')

    # GSE61794 (NSC x 2)
    obj61794 = rnaseq_data.gse61794(source='star', annotate_by='Ensembl Gene ID')

    # GBM paired samples
    objwtchg = rnaseq_data.gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='Ensembl Gene ID')

    # combine the data
    data = pd.concat((obj73721.data, obj61794.data, objwtchg.data), axis=1)
    data = data.loc[data.index.str.contains('ENSG')]

    mad = process.median_absolute_deviation(data, axis=1).sort_values(ascending=False)
    top_mad = mad.iloc[:N_GENES].index

    clustering.plot_clustermap(data.loc[top_mad], z_score=0, show_gene_labels=False)

    # for reference, we can also load the original 'pre-processed' GSE73721 data
    # but these are indexed by mouse gene??
    fpkm, meta = rnaseq_data.brainrnaseq_preprocessed()


    # find common genes based on lowercase match
    common_genes = a.index.str.lower().intersection(b.index.str.lower())

    a = a.loc[a.index.str.lower().isin(common_genes)]
    b = b.loc[b.index.str.lower().isin(common_genes)]

    # meet on common gene name grounds
    a_idx = np.array(a.index.str.capitalize())
    a_idx[a_idx == 'Spata2l'] = 'Spata2L'
    a.index = a_idx

    # combine
    c = pd.concat((a, b), axis=1)
    genes = [
        'Slc1a3',
        'Gja1',
        'Aldh1l1',
        'S100b',
        'Aqp4',
        'Gfap',
        'Cd44',
        'Aldoc',
        'Gfap',
        'Gapdh',
        'B2m'
    ]
    for g in genes:
        fig = plt.figure(figsize=(4.5, 7))
        ax = fig.add_subplot(111)
        ax.set_position([0.5, 0.08, 0.48, 0.9])
        c.loc[g].plot.barh(width=0.8, ax=ax)
        ax.set_xlabel(g)
        fig.savefig(os.path.join(OUTDIR, '%s.png' % g), dpi=200)
        fig.savefig(os.path.join(OUTDIR, '%s.pdf' % g))
        # plt.close(fig)

    # TODO: if required, make a single chart with all early-stage markers side by side

    # produce some bar charts

    # could normalise here by sample?
    # d = c.divide(c.sum(axis=0), axis=1)

    # find N most variable genes in the BrainRnaSeq dataset
    mad = process.median_absolute_deviation(b).sort_values(ascending=False)
    idx = mad.index[:N_GENES]

    d = c.loc[idx]

    col_colors = pd.DataFrame(index=d.columns, columns=['group'])
    cmap = {
        'tumor': '#a6cee3',
        'Hippocampus': '#1f78b4',
        'Fetal': '#b2df8a',
        re.compile(r'ctx [0-9]* *astro', flags=re.IGNORECASE): '#33a02c',
        'ctx neuron': '#fb9a99',
        'oligo': '#e31a1c',
        'myeloid': '#fdbf6f',
        'ctx endo': '#ff7f00',
        'whole cortex': '#cab2d6',
        'DURA': 'black',
    }
    for t, col in cmap.items():
        col_colors.loc[col_colors.index.str.contains(t)] = col

    z = hierarchy.linkage(d.transpose(), method='average', metric='correlation')

    cg = clustering.plot_clustermap(d, show_gene_labels=False,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z,
        z_score=0,
        xticklabels=True
    )
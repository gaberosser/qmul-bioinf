import os
import re

import pandas as pd

from astrocytes_comparison import plot_all_clustermaps
from load_data import rnaseq_data
from scripts.rnaseq import gtf_reader
from utils.output import unique_output_dir

N_GENES = 500


if __name__ == "__main__":

    SHOW_GENE_LABELS = False
    OUTDIR = unique_output_dir("btr_review", reuse_empty=True)
    INCLUDE_ALL_NSC = True
    COMBINE_REPLICATES = True  # merges the H9 replicates
    SIDETRACK_1 = False

    # GSE73721 (reference astrocytes, oligos, ...)
    obj73721 = rnaseq_data.gse73721(source='star', annotate_by='Ensembl Gene ID')

    # remove unneeded samples
    to_keep73721 = (
        ~obj73721.data.columns.str.contains('whole cortex')
        & ~obj73721.data.columns.str.contains('endo')
        & ~obj73721.data.columns.str.contains('tumor')
        & ~obj73721.data.columns.str.contains('myeloid')
    )

    # GSE61794 (H9-derived NSC x 2)
    obj61794 = rnaseq_data.gse61794(source='star', annotate_by='Ensembl Gene ID')
    if COMBINE_REPLICATES:
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


    # GBM paired samples
    objwtchg_paired = rnaseq_data.gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='Ensembl Gene ID')

    # WTCHG ALL samples
    objwtchg_all = rnaseq_data.all_wtchg_loader(source='star', annotate_by='Ensembl Gene ID')
    to_keepwtchg = (
        objwtchg_all.data.columns.str.contains('NSC')
        | objwtchg_all.data.columns.str.contains('ASTRO')
    )
    objwtchg_all.data = objwtchg_all.data.loc[:, to_keepwtchg]
    objwtchg_all.meta = objwtchg_all.meta.loc[to_keepwtchg, :]

    # Pollard (NSC x 2)
    objpollard = rnaseq_data.pollard_nsc(source='star', annotate_by='Ensembl Gene ID')
    objpollard.meta.iloc[:, -1] = ['Fetal NSC (1)', 'Fetal NSC (2)']
    objpollard.data.columns = ['Fetal NSC (1)', 'Fetal NSC (2)']

    # H9 references
    objh9_1 = rnaseq_data.gse24399_merged_loader(annotate_by='Ensembl Gene ID')
    objh9_1.meta.iloc[0, -1] = 'H9 ESC (1)'
    objh9_1.data.columns = ['H9 ESC (1)']

    objh9_2 = rnaseq_data.gse77920_loader(annotate_by='Ensembl Gene ID')
    objh9_2.meta.iloc[0, -1] = 'H9 ESC (2)'
    objh9_2.data.columns = ['H9 ESC (2)']

    # rRNA gene IDs
    rrna_ensg = set(gtf_reader.get_rrna())

    # MT gene_ids
    mt_ensg = set(gtf_reader.get_mitochondrial())

    # combine the data
    if INCLUDE_ALL_NSC:
        data = pd.concat((
            obj73721.data.loc[:, to_keep73721],
            obj61794.data,
            objwtchg_all.data,
            objpollard.data,
            objh9_1.data,
            objh9_2.data,
        ), axis=1)

        # combine the metadata
        meta = pd.concat((
            obj73721.meta.loc[to_keep73721],
            obj61794.meta,
            objwtchg_all.meta,
            objpollard.meta,
            objh9_1.meta,
            objh9_2.meta,
        ), axis=0)

    else:
        data = pd.concat((obj73721.data.loc[:, to_keep73721], obj61794.data, objwtchg_paired.data, objpollard.data), axis=1)

        # combine the metadata
        meta = pd.concat((
            obj73721.meta.loc[to_keep73721],
            obj61794.meta,
            objwtchg_paired.meta,
            objpollard.meta
        ), axis=0)

    data = data.loc[data.index.str.contains('ENSG')]

    col_colors = pd.DataFrame(index=data.columns, columns=['group'])

    sample_groups = {
        'Hippocampal astrocytes': {'regex': 'Hippocampus', 'colour' : '#7fc97f'},
        'Fetal cortical astrocytes': {'regex': 'Fetal ctx', 'colour' : '#beaed4'},
        'Cortical astrocytes': {'regex': re.compile(r'[0-9]*yo *ctx [0-9]* *astro', flags=re.IGNORECASE), 'colour' : '#fdc086'},
        'Cortical neurons': {'regex': 'ctx neuron', 'colour': '#ffff99'},
        'Oligodendrocytes': {'regex': 'oligo', 'colour': '#386cb0'},
        'Our iNSC': {'regex': '_NSC', 'colour': '#777777'},
        'H9 iNSC': {'regex': 'H9 NSC', 'colour': '#cccccc'},
        'Pollard NSC': {'regex': 'Fetal NSC', 'colour': '#cccccc'},
        'Our induced astrocytes': {'regex': 'ASTRO', 'colour': 'black'},
        'ESC': {'regex': 'H9 ESC', 'colour': '#ff3333'},
    }

    for d in sample_groups.values():
        col_colors.loc[col_colors.index.str.contains(d['regex'])] = d['colour']

    # remove rRNA
    data_rr = data.loc[~data.index.isin(rrna_ensg)]

    # remove MT in addition to rRNA
    data_rr_mt = data_rr.loc[~data_rr.index.isin(mt_ensg)]

    # all data except our astrocytes
    data_no_iastro = data_rr_mt.loc[:, ~data_rr_mt.columns.str.contains('_ASTRO')]

    filestem = os.path.join(OUTDIR, 'clustermap_no_iastro')
    col_order = plot_all_clustermaps(data_no_iastro, filestem, col_colors=col_colors)

    # our NSC, ref NSC, H9 ESC
    data_nsc_esc = data_rr_mt.loc[:, data_rr_mt.columns.str.contains('NSC') | data_rr_mt.columns.str.contains('ESC')]
    filestem = os.path.join(OUTDIR, 'clustermap_nsc_esc')
    col_order = plot_all_clustermaps(data_nsc_esc, filestem, col_colors=col_colors)

    filestem = os.path.join(OUTDIR, 'clustermap_nsc_esc1')
    col_order = plot_all_clustermaps(
        data_nsc_esc.loc[:, data_nsc_esc.columns != "H9 ESC (2)"],
        filestem,
        col_colors=col_colors
    )

    filestem = os.path.join(OUTDIR, 'clustermap_nsc_esc2')
    col_order = plot_all_clustermaps(
        data_nsc_esc.loc[:, data_nsc_esc.columns != "H9 ESC (1)"],
        filestem,
        col_colors=col_colors
    )



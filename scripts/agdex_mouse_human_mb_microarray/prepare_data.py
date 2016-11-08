import os
import numpy as np
from microarray import illumina, process, annotation
from scripts.agdex_mouse_human_mb_microarray import load_data
from settings import DATA_DIR


def load_hkg_list():
    fn = os.path.join(DATA_DIR, 'mouse_housekeeping/nature05453_supptable1.csv')
    hkg_all = pd.read_csv(fn)
    hkg = hkg_all.loc[:, 'Near Ubiquitious']
    return hkg[hkg.notnull()].values


if __name__ == '__main__':
    from scripts.agdex_mouse_human_mb_microarray import generate_ortholog_table
    from scripts.comparison_rnaseq_microarray import load_references
    from plotting import corr
    import pandas as pd
    from matplotlib import pyplot as plt
    import seaborn as sns

    FIGSIZE = (10.25, 8.)
    fig_kwargs = {'figsize': FIGSIZE}

    # load healthy mouse cerebellum data
    mo_he = load_data.load_annotated_microarray_gse54650()  # indexed by Entrez gene ID

    # load mouse MB data
    mo_mb, chd7 = load_data.load_annotated_microarray_sb_data() # indexed by Entrez gene ID

    # reduce to common genes
    common_genes = mo_he.index.intersection(mo_mb.index)

    mo_mb = mo_mb.loc[common_genes]
    mo_he = mo_he.loc[common_genes]

    # combine
    mo_all = pd.concat((mo_he, mo_mb), axis=1)

    # YuGene
    yg_mo_all = process.yugene_transform(mo_all)

    # apply YuGene transform
    # yg_mo_he = process.yugene_transform(mo_he)

    # apply YuGene
    # yg_mb = process.yugene_transform(mo_mb)

    # yg_mb = yg_mb.loc[common_genes]
    # yg_mo_he = yg_mo_he.loc[common_genes]

    # plot correlation using ALL matching genes
    # if True:
    if False:
        corr.plot_correlation_coefficient_array(yg_mo_all, vmin=0.8, fig_kwargs=fig_kwargs)
        plt.tight_layout()

    # HKG - use a few 'standard' ones
    # from http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-000A.html
    # https://www.qiagen.com/us/spotlight-pages/newsletters-and-magazines/articles/endogenous-controls/
    mo_hkg = [
        'Rn18s',
        'Actb',
        'Gapdh',  # missing?
        'Rpl13a',
        'B2m',
        'Hmbs',
        'Pgk1',
        'Hsp90ab1',
        # 'Hprt',  # this one is a bit noisy?
        'Pgk1',
    ]

    # use the list of genes published in a Nature paper
    # hkg = load_hkg_list()
    # TODO: some of the names used in the HKG list seem to match on SYNONYMS
    # e.g. 0610008A10Rik is a synonym of Aph1b
    # Initially 66 / 121 HKGs are matched, but we can probably do better
    # NB: These performed very poorly in terms of improving correlation between the two mouse datasets

    mo_annot = annotation.load_from_r_format('mogene10sttranscriptcluster.db')
    # translate from gene symbol to Entrez ID
    mo_hkg_translation = mo_annot.loc[mo_annot.loc[:, 'gene_symbol'].isin(mo_hkg), ['gene_symbol', 'entrez_id']]
    mo_hkg_entrez_id = mo_hkg_translation.loc[:, 'entrez_id']  # we silently lose hkg (e.g. Gapdh) here

    mo_hkg_entrez_id = yg_mo_all.index.intersection(mo_hkg_entrez_id)
    yg_mo_hkg_all = yg_mo_all.loc[mo_hkg_entrez_id]

    # if True:
    if False:
        corr.plot_correlation_coefficient_array(yg_mo_hkg_all, vmin=0.8, fig_kwargs=fig_kwargs)
        plt.tight_layout()

    # load human samples and avg over repeats
    hu_mb, hu_mb_meta = load_data.load_annotated_microarray_gse37382()


    hu_he, hu_he_meta = load_references.load_cerebellum_microarray_reference_data()
    hu_he = load_references.microarray_entrez_markers(hu_he, method='median')

    hu_all = pd.concat((hu_he, hu_mb), axis=1).dropna(axis=0, how='any')
    yg_hu_all = process.yugene_transform(hu_all)

    # write to txt for R
    mo_all.to_csv('mo.txt.gz', compression='gzip')
    yg_mo_all.to_csv('mo_yg.txt.gz', compression='gzip')
    hu_all.to_csv('hu.txt.gz', compression='gzip')
    yg_hu_all.to_csv('hu_yg.txt.gz', compression='gzip')

    mouse_tid = 10090
    human_tid = 9606
    orth_en = generate_ortholog_table.homologs(human_tid, mouse_tid, field='gene_id')  # gene_id is Entrez ID
    orth_gs = generate_ortholog_table.homologs(human_tid, mouse_tid, field='gene_symbol')
    # use the mouse gene ID as global index (anything would do)
    orth_tot = pd.concat((orth_gs, orth_en), axis=1).set_index('gene_symbol_%s' % human_tid)


    # write to txt for R
    orth_tot.to_csv('homolog_mapping.txt')

    # group names
    hu_he_names = hu_all.columns[:9]
    hu_mb_names = hu_all.columns[9:]
    mo_he_names = mo_he.columns
    mo_mb_names = mo_mb.columns
    mo_mb_grp1_names = mo_mb.columns[:3]
    mo_mb_grp2_names = mo_mb.columns[3:]

    # write phenoData to txt
    # hu_meta = load_illumina_data.load_sample_metadata()
    hu_pdata = pd.DataFrame(columns=['type', 'subgroup'], index=hu_all.columns)
    hu_pdata.loc[hu_he.columns, :] = ['hu.control', 'hu.control']
    hu_subgroups = hu_mb_meta.loc[:, 'subgroup'].replace('Group 3', 'C').replace('Group 4', 'D')
    hu_pdata.loc[hu_mb.columns, 'subgroup'] = hu_subgroups

    hu_pdata.loc[hu_mb.columns, 'type'] = 'hu.mb'
    hu_pdata.to_csv('hu_pdata.txt')

    mo_pdata = pd.DataFrame(columns=['type', 'chd7'], index=mo_all.columns)
    mo_pdata.loc[mo_he_names, :] = ['mo.control', 'mo.control']
    mo_pdata.loc[mo_mb_names, 'type'] = 'mo.mb'
    mo_pdata.loc[mo_mb_names, 'chd7'] = ['mo.chd7' if t else 'mo.no_chd7' for t in chd7]
    mo_pdata.to_csv('mo_pdata.txt')

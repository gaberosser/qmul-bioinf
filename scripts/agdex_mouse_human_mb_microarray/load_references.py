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
    from scripts.agdex_mouse_human_mb_microarray import load_sb_mouse_data
    from scripts.comparison_rnaseq_microarray import load_illumina_data
    from plotting import corr
    import pandas as pd
    from matplotlib import pyplot as plt
    import seaborn as sns

    FIGSIZE = (10.25, 8.)
    fig_kwargs = {'figsize': FIGSIZE}

    # load healthy mouse cerebellum data
    healthy = load_data.load_annotated_microarray_gse54650()  # indexed by Entrez gene ID

    # load mouse MB data
    mb, chd7 = load_data.load_annotated_microarray_sb_data() # indexed by Entrez gene ID

    # apply YuGene transform
    yg_healthy = process.yugene_transform(healthy)
    # yg_healthy_log = process.yugene_transform(np.log10(healthy))

    # apply YuGene
    yg_mb = process.yugene_transform(mb)

    # reduce to common genes
    common_genes = healthy.index.intersection(mb.index)

    mb = mb.loc[common_genes]
    healthy = healthy.loc[common_genes]
    # healthy_log = np.log10(healthy)

    yg_mb = yg_mb.loc[common_genes]
    yg_healthy = yg_healthy.loc[common_genes]
    # yg_healthy_log = yg_healthy_log.loc[common_genes]

    # this shows the hists of the raw data and the YG transformed data
    # Most importantly, we can see that a log transformation is definitely required on the healthy data
    # I guess this must have already been carried out on the SB samples.
    if False:
        fig, axs = plt.subplots(nrows=2, ncols=2)
        axs[0, 0].hist(mb.values.flatten(), 100, normed=True)
        axs[0, 1].hist(healthy.values.flatten(), 100, normed=True)
        # axs[0, 2].hist(healthy_log.values.flatten(), 100, normed=True)
        axs[1, 0].hist(yg_mb.values.flatten(), 100, normed=True)
        axs[1, 1].hist(yg_healthy.values.flatten(), 100, normed=True)
        # axs[1, 2].hist(yg_healthy_log.values.flatten(), 100, normed=True)

        axs[0, 0].set_title('SB')
        axs[0, 1].set_title('Healthy')
        # axs[0, 2].set_title('log10(Healthy)')
        axs[1, 0].set_title('YuGene SB')
        axs[1, 1].set_title('YuGene Healthy')
        # axs[1, 2].set_title('YuGene log10(Healthy)')

    yg_all = pd.concat((yg_healthy, yg_mb), axis=1)

    # plot correlation using ALL matching genes
    if True:
    # if False:
        corr.plot_correlation_coefficient_array(yg_all, vmin=0.8, fig_kwargs=fig_kwargs)
        plt.tight_layout()

    # HKG - use a few 'standard' ones
    # from http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-000A.html
    # https://www.qiagen.com/us/spotlight-pages/newsletters-and-magazines/articles/endogenous-controls/
    hkg = [
        'Rn18s',
        'Actb',
        'Gapdh',
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

    annot = annotation.load_from_r_format('mogene10sttranscriptcluster.db')
    hkg_entrez_id = annot.loc[annot.loc[:, 'gene_symbol'].isin(hkg), 'entrez_id']  # we silently lose hkg (e.g. Gapdh) here

    common_hkg = yg_all.index.intersection(hkg_entrez_id)

    yg_hkg_all = yg_all.loc[common_hkg]

    if True:
    # if False:
        corr.plot_correlation_coefficient_array(yg_hkg_all, vmin=0.8, fig_kwargs=fig_kwargs)
        plt.tight_layout()

    # load human samples
    hu_all = load_illumina_data.load_normed_microarray_data(pval=0.01)
    ilm_probes = load_illumina_data.load_illumina_array_library()
    hu_all = load_illumina_data.convert_microarray_to_gene_activity(hu_all, ilm_probes)
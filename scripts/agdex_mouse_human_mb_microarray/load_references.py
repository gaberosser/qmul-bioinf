import os
import numpy as np
from microarray import illumina, process
from settings import DATA_DIR


def load_annotated_microarray_gse54650():
    directory = os.path.join(DATA_DIR, 'microarray_GSE54650')
    probe_fn = os.path.join(directory, 'probe_set', 'GPL6246-18741.txt')
    sample_names = [
        'Cer_CT18',
        'Cer_CT20',
        'Cer_CT22',
        'Cer_CT24',
        'Cer_CT26',
        'Cer_CT28',
        'Cer_CT30',
        'Cer_CT32',
        'Cer_CT34',
        'Cer_CT36',
        'Cer_CT38',
        'Cer_CT40',
        'Cer_CT42',
        'Cer_CT44',
        'Cer_CT46',
        'Cer_CT48',
        'Cer_CT50',
        'Cer_CT52',
        'Cer_CT54',
        'Cer_CT56',
        'Cer_CT58',
        'Cer_CT60',
        'Cer_CT62',
        'Cer_CT64',
    ]
    arr_data = illumina.load_data(directory, sample_names, extension='.gz')
    probes = illumina.load_full_probeset_definitions(probe_fn, format='txt')
    probes = probes.loc[probes.loc[:, 'gene_assignment'].notnull(), 'gene_assignment']
    probes = probes.apply(lambda x: x.split('//')[1].strip())  # gene symbol

    # attach gene symbol and aggregate
    arr_data.insert(0, 'gene_symbol', probes)
    arr_data = arr_data.loc[arr_data.loc[:, 'gene_symbol'].notnull()]
    arr_data = process.aggregate_by_probe_set(arr_data)

    return arr_data


def load_hkg_list():
    fn = os.path.join(DATA_DIR, 'mouse_housekeeping/nature05453_supptable1.csv')
    hkg_all = pd.read_csv(fn)
    hkg = hkg_all.loc[:, 'Near Ubiquitious']
    return hkg[hkg.notnull()].values


if __name__ == '__main__':
    from scripts.agdex_mouse_human_mb_microarray import load_sb_mouse_data
    from plotting import corr
    import pandas as pd
    from matplotlib import pyplot as plt
    import seaborn as sns

    healthy = load_annotated_microarray_gse54650()

    # apply YuGene transform
    yg_healthy = process.yugene_transform(healthy)
    yg_healthy_log = process.yugene_transform(np.log10(healthy))

    # load mouse MB data
    mb, chd7 = load_sb_mouse_data.load_annotated_array_data()

    # apply YuGene
    yg_mb = process.yugene_transform(mb)

    # reduce to common genes
    common_genes = healthy.index.intersection(mb.index)

    mb = mb.loc[common_genes]
    healthy = healthy.loc[common_genes]
    healthy_log = np.log10(healthy)

    yg_mb = yg_mb.loc[common_genes]
    yg_healthy = yg_healthy.loc[common_genes]
    yg_healthy_log = yg_healthy_log.loc[common_genes]

    # this shows the hists of the raw data and the YG transformed data
    # Most importantly, we can see that a log transformation is definitely required on the healthy data
    # I guess this must have already been carried out on the SB samples.
    if False:
        fig, axs = plt.subplots(nrows=2, ncols=3)
        axs[0, 0].hist(mb.values.flatten(), 100, normed=True)
        axs[0, 1].hist(healthy.values.flatten(), 100, normed=True)
        axs[0, 2].hist(healthy_log.values.flatten(), 100, normed=True)
        axs[1, 0].hist(yg_mb.values.flatten(), 100, normed=True)
        axs[1, 1].hist(yg_healthy.values.flatten(), 100, normed=True)
        axs[1, 2].hist(yg_healthy_log.values.flatten(), 100, normed=True)

        axs[0, 0].set_title('SB')
        axs[0, 1].set_title('Healthy')
        axs[0, 2].set_title('log10(Healthy)')
        axs[1, 0].set_title('YuGene SB')
        axs[1, 1].set_title('YuGene Healthy')
        axs[1, 2].set_title('YuGene log10(Healthy)')

    yg_all = pd.concat((yg_healthy_log, yg_mb), axis=1)

    # plot correlation using ALL matching genes
    corr.plot_correlation_coefficient_array(yg_all, vmin=0.8)

    # HKG
    hkg = load_hkg_list()

    # TODO: some of the names used in the HJG list seem to match on SYNONYMS
    # e.g. 0610008A10Rik is a synonym of Aph1b
    # Initially 66 / 121 HKGs are matched, but we can probably do better

    common_hkg = yg_all.index.intersection(hkg)

    yg_hkg_all = yg_all.loc[common_hkg]

    corr.plot_correlation_coefficient_array(yg_hkg_all, vmin=0.8)
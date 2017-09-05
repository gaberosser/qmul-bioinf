import os
from load_data import rnaseq_data
from plotting import corr, clustering
from matplotlib import pyplot as plt
from settings import DATA_DIR_NON_GIT
import pandas as pd
import seaborn as sns


if __name__ == "__main__":
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')

    lanedirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
    ]

    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'mouse', 'star_alignment') for d in lanedirs]
    samples_ensc_med = [u'eNSC3med', u'eNSC5med', u'eNSC6med']
    samples_ensc_mouse = [u'eNSC3mouse', u'eNSC5mouse', u'eNSC6mouse',]
    samples_insc_mouse = [u'mDura3N1mouse', u'mDura5N24Amouse', u'mDura6N6mouse']
    samples_insc_human = [u'mDura3N1human', u'mDura5N24Ahuman', u'mDura6N6human']

    # samples = samples_ensc_med + samples_ensc_mouse + samples_insc_mouse + samples_insc_human
    samples = list(samples_ensc_med)
    # optional
    samples += ['CRL3034_METASTATICMB_P4']

    obj = rnaseq_data.all_samples_multilane_loader(
        countdirs, metafiles, source='star', annotate_by='Ensembl Gene ID', samples=samples
    )
    data = obj.data.loc[obj.data.index.str.contains('ENS')]
    meta = obj.meta

    # load replacements
    indir_rep = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170506', '170829_K00150_0236_AHL5YHBBXX')
    obj_rep = rnaseq_data.StarCountLoader(
        count_dir=os.path.join(indir_rep, 'mouse', 'star_alignment'),
        strandedness='r',
        samples=samples_ensc_med,
        meta_fn=os.path.join(indir_rep, 'sources.csv'))
    data_rep = obj_rep.data.loc[obj.data.index.str.contains('ENS')]
    meta_rep = obj_rep.meta
    # change sample names
    data_rep.columns = [t + '_rep' for t in samples_ensc_med]
    meta_rep.index = data_rep.columns

    # add new data  (optional)
    data = pd.concat((data, data_rep), axis=1)
    meta = pd.concat((meta, meta_rep), axis=0)

    cpm = data.divide(meta.loc[:, 'read_count'].values, axis=1) * 1e6
    keep = (cpm > 1).sum(axis=1) > 5

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    # ax = corr.plot_correlation_coefficient_array(data.loc[keep], vmin=0.6, ax=ax)
    ax = corr.plot_correlation_coefficient_array(data.loc[keep], vmin=0.98, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    fig.tight_layout()

    # clustering.plot_correlation_clustermap(data.loc[keep])


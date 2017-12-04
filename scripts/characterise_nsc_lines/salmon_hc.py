from load_data import rnaseq_data
from plotting import clustering
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
import os
from settings import LOCAL_DATA_DIR


def hist_logvalues(data, thresholds=None, eps=1e-6):
    all_vals = data.values.flatten().astype(float)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(np.log10(all_vals + 1e-6), 100)
    if thresholds:
        for x in thresholds:
            ax.axvline(np.log10(x), c='k', ls='--')
    return ax


def cluster_logdata_with_threshold(data, min_val, n, mad=None, min_over=2, eps=1e-2):
    func = lambda x: np.log2(data + eps)
    return cluster_data_with_threshold(data, min_val, n, mad=mad, min_over=min_over, transform=func)


def cluster_data_with_threshold(data, min_tpm, n, mad=None, min_over=2, transform=None):
    idx = (data > min_tpm).sum(axis=1) > min_over
    dat = data.loc[idx]
    if transform is not None:
        dat = transform(dat)
    if mad is None:
        mad = transformations.median_absolute_deviation(dat).sort_values(ascending=False)
    else:
        if len(mad.index.intersection(data.loc[idx].index)) != idx.sum():
            raise AttributeError("If a pre-computed MAD is supplied, it must contain all required entries")
    this_mad = mad.loc[dat.loc[idx].index]
    the_dat = dat.loc[this_mad.index[:n]]

    cm = clustering.plot_clustermap(the_dat, cmap='RdBu_r', metric='correlation')
    cm.gs.update(bottom=0.2)
    return cm, mad


def cluster_data_with_threshold(data, min_tpm, n, mad=None, min_over=2, eps=1e-2):
    idx = (data > min_tpm).sum(axis=1) > min_over
    p = data.loc[idx]
    if mad is None:
        mad = transformations.median_absolute_deviation(p).sort_values(ascending=False)
    else:
        if len(mad.index.intersection(p.index)) != idx.sum():
            raise AttributeError("If a pre-computed MAD is supplied, it must contain all required entries")
    this_mad = mad.loc[data.loc[idx].index]
    the_dat = p.loc[this_mad.index[:n]]

    cm = clustering.plot_clustermap(the_dat, cmap='RdBu_r', metric='correlation')
    cm.gs.update(bottom=0.2)
    return cm, mad


if __name__ == "__main__":
    units = 'estimated_counts'
    outdir = unique_output_dir("salmon_insc_characterisation")
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']
    if units == 'tpm':
        min_val = 1
        min_n = 4
        min_vals = [.1, 1, 10]
    elif units == 'estimated_counts':
        min_val = 10
        min_n = 4
        min_vals = [1, 5, 10]

    patient_data = rnaseq_data.load_salmon_by_patient_id(pids, units=units)

    # discard GBM
    patient_data = patient_data.loc[:, ~patient_data.columns.str.contains('GBM')]

    # update index to remove accession version
    idx = patient_data.index.str.replace(r'.[0-9]+$', '')
    patient_data.index = idx



    # ax = hist_logvalues(patient_data, thresholds=min_vals)
    #
    # ax.figure.savefig(os.path.join(outdir, "log10_intensities_with_min_tpm_threshold.png"), dpi=200)
    # ax.figure.savefig(os.path.join(outdir, "log10_intensities_with_min_tpm_threshold.pdf"))
    #
    # mad = None
    # for n_t in [1000, 2000, 5000][::-1]:
    #     for min_val in min_vals:
    #         cm, mad = cluster_logdata_with_threshold(patient_data, min_val, n_t, mad=mad)
    #         fname = "clustering_corr_log10_top%d_by_mad_mintpm_%.1f.{ext}" % (n_t, min_val)
    #         cm.savefig(os.path.join(outdir, fname).format(ext='pdf'))
    #         cm.savefig(os.path.join(outdir, fname).format(ext='png'), dpi=200)

    # now aggregate to gene level and repeat
    # TODO: move to rnaseq module or similar

    fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90', 'gene_to_transcript.txt')
    gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

    # shouldn't be necessary, but remove transcripts that have no translation
    to_keep = patient_data.index.intersection(gene_transcript.index)
    if len(to_keep) != patient_data.shape[0]:
        to_drop = patient_data.index.difference(gene_transcript.loc[:, 'Transcript stable ID'])
        print "Discarding %d transcripts that have no associated gene: %s" % (
            len(to_drop), ', '.join(to_drop)
        )
        patient_data = patient_data.loc[to_keep]

    # gene list in same order as data
    genes = gene_transcript.loc[patient_data.index, 'Gene stable ID']

    patient_data_by_gene = patient_data.groupby(genes).sum()
    # discard genes expressed at low values
    idx = (patient_data_by_gene > min_val).sum(axis=1) > min_n

    pdbg_norm = patient_data_by_gene.loc[idx]
    pdbg_norm = pdbg_norm.divide(pdbg_norm.sum(axis=0), axis=1)

    ax = hist_logvalues(patient_data_by_gene, thresholds=min_vals)
    ax.figure.savefig(os.path.join(outdir, "log10_intensities_by_gene_with_min_tpm_threshold.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "log10_intensities_by_gene_with_min_tpm_threshold.pdf"))

    mad_by_gene = None
    for n_t in [1000, 2000, 5000][::-1]:
        # no m in val as we already applied it
        cm, mad_by_gene = cluster_logdata_with_threshold(pdbg_norm, 0, n_t, mad=mad_by_gene, eps=1e-12)
        fname = "clustering_by_gene_corr_log10_top%d_by_mad.{ext}" % n_t
        # cm.savefig(os.path.join(outdir, fname).format(ext='pdf'))
        cm.savefig(os.path.join(outdir, fname).format(ext='png'), dpi=200)

    # bring in reference data
    ref_dats = [
        rnaseq_data.gse73721_salmon(units=units),
        rnaseq_data.pollard_salmon(units=units),
        rnaseq_data.gse80732_salmon(units=units),
        rnaseq_data.gse64882_salmon(units=units),
        rnaseq_data.gse84166_salmon(units=units)
    ]
    ref = pd.concat(ref_dats, axis=1)
    ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    #discard fetal ctx astro (Barres)
    ref = ref.loc[:, ~ref.columns.str.contains('Fetal ctx')]

    ref_by_gene = ref.groupby(genes).sum()

    # now let's try clustering everything together
    altogether_by_gene = pd.concat((patient_data_by_gene, ref_by_gene), axis=1)

    # discard genes expressed at low values
    idx = (altogether_by_gene > min_val).sum(axis=1) > min_n

    abg_norm = altogether_by_gene.loc[idx]
    abg_norm = abg_norm.divide(abg_norm.sum(axis=0), axis=1)

    # this looks bad...
    cm, mad_all = cluster_logdata_with_threshold(abg_norm, 0, n=5000, eps=1e-12)

    abg_vsd = transformations.variance_stabilizing_transform(altogether_by_gene)



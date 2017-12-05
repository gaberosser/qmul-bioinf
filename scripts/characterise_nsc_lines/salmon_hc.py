from load_data import rnaseq_data
from plotting import clustering
from scripts.rnaseq import gtf_reader
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


def cluster_logdata_with_threshold(data, min_val=None, n=None, mad=None, min_over=2, eps=1e-2, **kwargs):
    func = lambda x: np.log2(data + eps)
    return cluster_data_with_threshold(data, min_val, n, mad=mad, min_over=min_over, transform=func, **kwargs)


def cluster_data_with_threshold(data, min_val=None, n=None, mad=None, min_over=2, transform=None, **kwargs):
    if min_val is not None and min_over is not None:
        idx = (data > min_val).sum(axis=1) > min_over
        data = data.loc[idx]

    if transform is not None:
        data = transform(data)

    if n is not None:
        if mad is None:
            mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
        else:
            mad = mad.sort_values(ascending=False)
            if len(mad.index.intersection(data.index)) != data.shape[0]:
                raise AttributeError("If a pre-computed MAD is supplied, it must contain all required entries")

        data = data.loc[mad.index[:n]]

    cm = clustering.plot_clustermap(data, cmap='RdBu_r', metric='correlation', **kwargs)
    cm.gs.update(bottom=0.2)
    return cm, mad


if __name__ == "__main__":
    # units = 'estimated_counts'
    units = 'tpm'
    # remove_mt = True
    remove_mt = False
    outdir = unique_output_dir("salmon_insc_characterisation")
    n_gene_try = [1000, 2000, 3000, 5000][::-1]  # largest first, so we can reuse the MAD array
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']
    if units == 'tpm':
        min_val = 1
        min_n = 4
        eps = .01
    elif units == 'estimated_counts':
        min_val = 10
        min_n = 4
        eps = .01

    if remove_mt:
        mt_ensg = set(gtf_reader.get_mitochondrial())

    patient_data = rnaseq_data.load_salmon_by_patient_id(pids, units=units)

    # discard GBM
    patient_data = patient_data.loc[:, ~patient_data.columns.str.contains('GBM')]

    # update index to remove accession version
    idx = patient_data.index.str.replace(r'.[0-9]+$', '')
    patient_data.index = idx

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

    # discard mitochondrial genes
    if remove_mt:
        idx = ~patient_data_by_gene.index.isin(mt_ensg)
        pdbg = patient_data_by_gene.loc[idx]
        # renorm
        if units == 'tpm':
            pdbg = pdbg.divide(pdbg.sum(), axis=1) * 1e6
    else:
        pdbg = patient_data_by_gene

    # discard genes expressed at low values
    idx = (pdbg > min_val).sum(axis=1) > min_n
    pdbg = pdbg.loc[idx]

    if units == 'estimated_counts':
        # here we can normalise by library size if desired
        pass

    ax = hist_logvalues(patient_data_by_gene, thresholds=[min_val])
    ax.figure.savefig(os.path.join(outdir, "log2_intensities_by_gene_with_min_tpm_threshold.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "log2_intensities_by_gene_with_min_tpm_threshold.pdf"))

    mad = transformations.median_absolute_deviation(pdbg).sort_values(ascending=False)
    pdbg_log = np.log2(pdbg + eps)
    mad_log = transformations.median_absolute_deviation(pdbg_log).sort_values(ascending=False)
    row_colours = pd.DataFrame('gray', index=pdbg_log.columns, columns=[''])
    row_colours.loc[row_colours.index.str.contains('IPSC')] = '#fdc086'
    row_colours.loc[row_colours.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'
    row_colours.loc[row_colours.index.str.contains('GIBCO')] = '#96daff'


    for n_t in n_gene_try:
        fname = "clustering_by_gene_corr_top%d_by_mad.{ext}" % n_t
        fname_log = "clustering_by_gene_corr_log_top%d_by_mad.{ext}" % n_t

        d = clustering.dendrogram_with_colours(pdbg_log.loc[mad_log.index[:n_t]], row_colours, fig_kws={'figsize': (10, 5.5)})
        d['fig'].savefig(os.path.join(outdir, fname_log.format(ext='png')), dpi=200)

    # bring in reference data
    ref_dats = [
        rnaseq_data.gse73721_salmon(units=units),
        rnaseq_data.pollard_salmon(units=units),
        rnaseq_data.gse80732_salmon(units=units),
        rnaseq_data.gse64882_salmon(units=units),
        rnaseq_data.gse84166_salmon(units=units),
        rnaseq_data.gse61794_salmon(units=units),
    ]

    ref = pd.concat(ref_dats, axis=1)
    ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    # discard Barres irrelevant samples
    ref = ref.loc[:, ~ref.columns.str.contains('Fetal ctx')]
    ref = ref.loc[:, ~ref.columns.str.contains('tumor')]
    ref = ref.loc[:, ~ref.columns.str.contains('ctx endo')]
    ref = ref.loc[:, ~ref.columns.str.contains('whole cortex')]
    ref = ref.loc[:, ~ref.columns.str.contains('myeloid')]
    # discard immortalised cell line
    ref = ref.loc[:, ~ref.columns.str.contains('LM-NSC')]
    # discard fibroblasts
    ref = ref.loc[:, ~ref.columns.str.contains(r'^fibroblast')]
    #

    # Take mean over various replicates (for simplicity)
    cols = ref.loc[:, ref.columns.str.contains('Pollard NSC')].columns
    ref.insert(0, 'Fetal NSC (Caren et al.)', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'H1-[12]')].columns
    ref.insert(0, 'H1', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'H1-EPS')].columns
    ref.insert(0, 'H1-EPS', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'ES1-EPS')].columns
    ref.insert(0, 'ES1-EPS', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'H9_NSC')].columns
    ref.insert(0, 'H9-NSC', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'INSC fibroblast')].columns
    ref.insert(0, 'Fibroblast-derived NSC (Shahbazi et al.)', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'NSC008')].columns
    ref.insert(0, 'Fetal NSC (Li et al.)', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    cols = ref.loc[:, ref.columns.str.contains(r'fetal NSC rep')].columns
    ref.insert(0, 'Fetal NSC (Shahbazi et al.)', ref.loc[:, cols].mean(axis=1))
    ref.drop(cols, axis=1, inplace=True)

    ref_by_gene = ref.groupby(genes).sum()

    # now let's try clustering everything together
    abg = pd.concat((patient_data_by_gene, ref_by_gene), axis=1)

    # discard mitochondrial genes
    if remove_mt:
        idx = ~abg.index.isin(mt_ensg)
        abg = abg.loc[idx]
        # renorm
        if units == 'tpm':
            abg = abg.divide(abg.sum(), axis=1) * 1e6

    # discard genes expressed at low values
    idx = (abg > min_val).sum(axis=1) > min_n
    abg = abg.loc[idx]

    abg_log = np.log2(abg + eps)
    amad_log = transformations.median_absolute_deviation(abg_log).sort_values(ascending=False)

    if units == 'estimated_counts':
        # optionally could normalise here?
        pass

    row_colours_all = pd.DataFrame('gray', index=abg.columns, columns=[''])
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#96daff'
    row_colours_all.loc[row_colours_all.index.str.contains(r'H9-NSC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fetal')] = 'yellow'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fibroblast')] = '#fff89e'
    row_colours_all.loc[row_colours_all.index.str.contains('H1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('ES1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Hippocamp.* astro')] = '#e78ac3'
    row_colours_all.loc[row_colours_all.index.str.contains(r'ctx .*astro')] = '#b3b3b3'
    row_colours_all.loc[row_colours_all.index.str.contains(r'oligo')] = '#fccde5'

    for n_t in n_gene_try:
        fname = "all_samples_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t

        cm, mad_all = cluster_logdata_with_threshold(abg, n=n_t, eps=eps, col_colors=row_colours_all)
        cm.gs.update(bottom=0.3)
        cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "all_samples_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
        d = clustering.dendrogram_with_colours(
            abg_log.loc[amad_log.index[:n_t]],
            row_colours_all,
            fig_kws={'figsize': (10, 5.5)}
        )
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)



from load_data import rnaseq_data
from plotting import clustering
import references
from rnaseq import general, loader
from scripts.rnaseq import gtf_reader
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
import os
import re
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


if __name__ == '__main__':
    # units = 'estimated_counts'
    units = 'tpm'
    # remove_mt = True
    remove_mt = False
    outdir = unique_output_dir("salmon_insc_mouse", reuse_empty=True)
    n_gene_try = [1000, 2000, 3000, 5000][::-1]  # largest first, so we can reuse the MAD array

    if units == 'tpm':
        min_val = 1
        min_n = 4
        eps = .01
    elif units == 'estimated_counts':
        min_val = 10
        min_n = 4
        eps = .01

    if remove_mt:
        mt_ensg = set(gtf_reader.get_mitochondrial('GRCm38r88'))

    mouse_data = rnaseq_data.mouse_nsc_salmon(units=units)

    # aggregate to gene level
    mouse_data_by_gene = general.ensembl_transcript_quant_to_gene(mouse_data, tax_id=10090)

    # discard mitochondrial genes
    if remove_mt:
        idx = ~mouse_data_by_gene.index.isin(mt_ensg)
        mdbg = mouse_data_by_gene.loc[idx]
        # renorm
        if units == 'tpm':
            mdbg = mdbg.divide(mdbg.sum(), axis=1) * 1e6
    else:
        mdbg = mouse_data_by_gene

    # discard genes expressed at low values
    # idx = (mdbg > min_val).sum(axis=1) > min_n
    # mdbg = mdbg.loc[idx]

    if units == 'estimated_counts':
        # here we can normalise by library size if desired
        pass

    # ax = hist_logvalues(mouse_data_by_gene, thresholds=[min_val])
    # ax.figure.savefig(os.path.join(outdir, "log2_intensities_by_gene_with_min_tpm_threshold.png"), dpi=200)
    # ax.figure.savefig(os.path.join(outdir, "log2_intensities_by_gene_with_min_tpm_threshold.pdf"))

    mad = transformations.median_absolute_deviation(mdbg).sort_values(ascending=False)
    mdbg_log = np.log2(mdbg + eps)
    mad_log = transformations.median_absolute_deviation(mdbg_log).sort_values(ascending=False)
    row_colours = pd.DataFrame('gray', index=mdbg_log.columns, columns=[''])
    row_colours.loc[row_colours.index.str.contains(r'eNSC[0-9]med')] = '#66c2a5'
    row_colours.loc[row_colours.index.str.contains(r'eNSC[0-9]mouse')] = '#fc8d62'
    row_colours.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*mouse')] = '#8da0cb'
    row_colours.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*human')] = '#e78ac3'


    for n_t in n_gene_try:
        fname = "clustering_by_gene_corr_log_top%d_by_mad.{ext}" % n_t

        d = clustering.dendrogram_with_colours(mdbg_log.loc[mad_log.index[:n_t]], row_colours, fig_kws={'figsize': (10, 5.5)})
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        cm, _ = cluster_data_with_threshold(mdbg_log, n=n_t, mad=mad_log, col_colors=row_colours)
        cm.savefig(os.path.join(outdir, "clustermap_by_gene_corr_log_top%d_by_mad.png" % n_t), dpi=200)

    # bring in reference data
    ref_dats = [
        ('Pten/P53 study', rnaseq_data.mouse_gbm_pten_p53(source='salmon', units=units)),
        # ('Zhang et al., reprog.', rnaseq_data.gse78938_salmon(units=units)),
        ('Liu et al.', rnaseq_data.gse96950_salmon(units=units)),
        ('Wapinski et al.', rnaseq_data.gse43916_salmon(units=units)),
        # ('Friedmann-Morvinski et al.', rnaseq_data.gse73127_salmon(units=units)),
        # ('Friedmann-Morvinski et al.', rnaseq_data.gse64411_salmon(units=units)),
        # ('Zhang et al.', rnaseq_data.gse52564_salmon(units=units)),
        # ('Chen et al.', rnaseq_data.gse52125_salmon(units=units)),
        ('Yanez et al.', rnaseq_data.gse88982_salmon(units=units)),
        ('Lynch', rnaseq_data.gse78795_salmon(units=units)),
        ('Moyon et al.', rnaseq_data.gse66029_salmon(units=units)),
        ('Schmid et al.', rnaseq_data.gse75592_salmon(units=units)),
        ('Srinivasan et al.', rnaseq_data.gse75246_salmon(units=units)),
    ]

    # drop unneeded and rename
    for i, (auth, rd) in enumerate(ref_dats):
        # average over duplicates
        a = {}
        na = {}
        do_avg = False
        for col in rd.columns:
            if ' repl ' in col:
                do_avg = True
            t = re.sub(r' repl [0-9]*', '', col)
            if t in a:
                a[t] += rd.loc[:, col]
                na[t] += 1
            else:
                a[t] = rd.loc[:, col].copy()
                na[t] = 1
        for t in a:
            a[t] = a[t] / float(na[t])

        if do_avg:
            rd = pd.DataFrame(a)

        rd.columns = ["%s (%s)" % (col, auth) for col in rd.columns]
        ref_dats[i] = (auth, rd)

    ref = pd.concat([t[1] for t in ref_dats], axis=1)

    # just in case it's useful, stitch the data together and export before we do any removal
    tmp_ref_by_gene = general.ensembl_transcript_quant_to_gene(ref, tax_id=10090)
    tmp_abg = pd.concat((mouse_data_by_gene, tmp_ref_by_gene), axis=1)
    tmp_abg.to_excel(os.path.join(outdir, "tpm_values.xlsx"))

    # ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    ref = ref.loc[:, ~ref.columns.str.contains('Normal brain')]
    ref = ref.loc[:, ~ref.columns.str.contains('GBM')]
    ref = ref.loc[:, ~ref.columns.str.contains('TrNeuron')]
    ref = ref.loc[:, ~ref.columns.str.contains('TrAstrocyte')]
    ref = ref.loc[:, ~ref.columns.str.contains('Tumour')]
    ref = ref.loc[:, ~ref.columns.str.contains('MP and cMoP')]
    ref = ref.loc[:, ~ref.columns.str.contains('LPS')]
    ref = ref.loc[:, ~ref.columns.str.contains(r'day [148]')]

    # ref_by_gene = ref.groupby(genes).sum()
    ref_by_gene = general.ensembl_transcript_quant_to_gene(ref, tax_id=10090)

    # now let's try clustering everything together
    abg = pd.concat((mouse_data_by_gene, ref_by_gene), axis=1)

    # discard mitochondrial genes
    if remove_mt:
        idx = ~abg.index.isin(mt_ensg)
        abg = abg.loc[idx]
        # renorm
        if units == 'tpm':
            abg = abg.divide(abg.sum(), axis=1) * 1e6

    abg_log = np.log2(abg + eps)
    amad_log = transformations.median_absolute_deviation(abg_log).sort_values(ascending=False)

    if units == 'estimated_counts':
        # optionally could normalise here?
        pass

    # row_colours_all = pd.DataFrame('gray', index=abg.columns, columns=['Cell type', 'Study'])
    row_colours_all = pd.DataFrame('gray', index=abg.columns, columns=['Cell type',])

    # cell type
    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC'), 'Cell type'] = '#a5fff9' # pale turquoise(?)
    row_colours_all.loc[row_colours_all.index.str.contains(r'NSLC'), 'Cell type'] = '#a5fff9' # pale turquoise(?)
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC'), 'Cell type'] = '#1eb9d8' # light blue
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Nn]euron'), 'Cell type'] = '#6dada9'  # teal
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Aa]strocyte'), 'Cell type'] = '#ebff49' # yellow
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Oo]ligo'), 'Cell type'] = '#ffa8b3' # pale red
    row_colours_all.loc[row_colours_all.index.str.contains(r'OPC'), 'Cell type'] = '#b70017' # dark red
    row_colours_all.loc[row_colours_all.index.str.contains(r'MEF'), 'Cell type'] = '#f6ffaa' # pale yellow
    row_colours_all.loc[row_colours_all.index.str.contains(r'ESC'), 'Cell type'] = '#8b33dd' # pale purple
    row_colours_all.loc[row_colours_all.index.str.contains(r'[iI]PSC'), 'Cell type'] = '#8b33dd' # pale purple
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Mm]icroglia'), 'Cell type'] = '#ffd8af' # pale orange
    row_colours_all.loc[row_colours_all.index.str.contains(r'Yanez'), 'Cell type'] = '#ffa03d' # orange
    # these override previously-defined colours
    row_colours_all.loc[row_colours.index.str.contains(r'eNSC[0-9]med'), 'Cell type'] = '#96ff9d' # pale green
    row_colours_all.loc[row_colours.index.str.contains(r'eNSC[0-9]mouse'), 'Cell type'] = '#008408' # dark green
    row_colours_all.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*mouse'), 'Cell type'] = '#3543ff' # dark blue
    row_colours_all.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*human'), 'Cell type'] = '#c4c8ff' # pale blue

    for n_t in n_gene_try + [abg.shape[0]]:
        fname = "all_samples_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t

        # don't plot the clustermap with all genes
        if n_t < abg.shape[0]:
            cm, _= cluster_data_with_threshold(abg_log.loc[amad_log.index[:n_t]], col_colors=row_colours_all)
            cm.gs.update(bottom=0.3)
            cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "all_samples_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
        d = clustering.dendrogram_with_colours(
            abg_log.loc[amad_log.index[:n_t]],
            row_colours_all,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False
        )
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "all_samples_corrplot_log_top%d_by_mad.{ext}" % n_t
        cm = clustering.plot_correlation_clustermap(
            abg_log.loc[amad_log.index[:n_t]],
            row_colors=row_colours_all,
            n_gene=n_t,
        )
        plt.setp(cm.ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
        plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
        cm.gs.update(bottom=0.35, right=0.65)
        cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)


    # data subset: our data and Pten/P53 samples

    mouse_data = rnaseq_data.mouse_nsc_salmon(units=units)
    mouse_data_by_gene = general.ensembl_transcript_quant_to_gene(mouse_data, tax_id=10090)

    dat_pten = general.ensembl_transcript_quant_to_gene(
        rnaseq_data.mouse_gbm_pten_p53(source='salmon', units=units),
        tax_id=10090
    )
    abg = pd.concat((mouse_data_by_gene, dat_pten), axis=1)
    abg = abg.loc[:, ~abg.columns.str.contains(r'mDura[1-9NA]+mouse')]
    abg = abg.loc[:, ~abg.columns.str.contains(r'eNSC[356]mouse')]
    abg = abg.loc[:, ~abg.columns.str.contains(r'TG')]
    abg_log = np.log2(abg + eps)
    amad_log = transformations.median_absolute_deviation(abg_log).sort_values(ascending=False)

    if remove_mt:
        idx = ~abg.index.isin(mt_ensg)
        abg = abg.loc[idx]
        # renorm
        if units == 'tpm':
            abg = abg.divide(abg.sum(), axis=1) * 1e6

    row_colours_all = pd.DataFrame('gray', index=abg.columns, columns=['Cell type', ])

    # cell type
    row_colours_all.loc[row_colours_all.index.str.contains(r'WT'), 'Cell type'] = '#a5fff9'  # pale turquoise(?)
    row_colours_all.loc[row_colours_all.index.str.contains(r'TG'), 'Cell type'] = '#a5fff9'  # pale turquoise(?)
    row_colours_all.loc[row_colours_all.index.str.contains(r'GBM'), 'Cell type'] = '#ffa8b3'  # pale red
    # these override previously-defined colours
    row_colours_all.loc[row_colours.index.str.contains(r'eNSC[0-9]med'), 'Cell type'] = '#96ff9d'  # pale green
    # row_colours_all.loc[row_colours.index.str.contains(r'eNSC[0-9]mouse'), 'Cell type'] = '#008408'  # dark green
    row_colours_all.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*mouse'), 'Cell type'] = '#3543ff'  # dark blue
    # row_colours_all.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*human'), 'Cell type'] = '#c4c8ff'  # pale blue

    abg_log = np.log2(abg + eps)
    amad_log = transformations.median_absolute_deviation(abg_log).sort_values(ascending=False)

    if units == 'estimated_counts':
        # optionally could normalise here?
        pass

    for n_t in n_gene_try + [abg.shape[0]]:
        fname = "nsc_gbm_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t

        # don't plot the clustermap with all genes
        if n_t < abg.shape[0]:
            cm, _ = cluster_data_with_threshold(abg_log.loc[amad_log.index[:n_t]], col_colors=row_colours_all)
            cm.gs.update(bottom=0.3)
            cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "nsc_gbm_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
        d = clustering.dendrogram_with_colours(
            abg_log.loc[amad_log.index[:n_t]],
            row_colours_all,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False
        )
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "nsc_gbm_corrplot_log_top%d_by_mad.{ext}" % n_t
        cm = clustering.plot_correlation_clustermap(
            abg_log.loc[amad_log.index[:n_t]],
            row_colors=row_colours_all,
            n_gene=n_t,
        )
        plt.setp(cm.ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
        plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
        cm.gs.update(bottom=0.35, right=0.65)
        cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)


    # same samples but aggegate wherever possible
    abg.insert(0, 'mDura_mouse', abg.loc[:, abg.columns.str.contains('mDura')].mean(axis=1))
    abg.insert(0, 'eNSC_med', abg.loc[:, abg.columns.str.contains('eNSC')].mean(axis=1))
    abg.insert(0, 'GBM', abg.loc[:, abg.columns.str.contains('GBM')].mean(axis=1))
    abg.insert(0, 'WT', abg.loc[:, abg.columns.str.contains('WT')].mean(axis=1))

    abg_avg = abg.iloc[:, :4]

    # cell type
    row_colours_avg = pd.DataFrame('gray', index=abg_avg.columns, columns=['Cell type', ])
    row_colours_avg.loc[row_colours_avg.index.str.contains(r'WT'), 'Cell type'] = '#a5fff9'  # pale turquoise(?)
    row_colours_avg.loc[row_colours_avg.index.str.contains(r'GBM'), 'Cell type'] = '#ffa8b3'  # pale red
    row_colours_avg.loc[row_colours_avg.index.str.contains(r'eNSC'), 'Cell type'] = '#96ff9d'  # pale green
    row_colours_avg.loc[row_colours_avg.index.str.contains(r'mDura'), 'Cell type'] = '#3543ff'  # dark blue

    abg_avg_log = np.log2(abg_avg + eps)
    amad_avg_log = transformations.median_absolute_deviation(abg_avg_log).sort_values(ascending=False)

    for n_t in n_gene_try + [abg_avg_log.shape[0]]:
        fname = "nsc_gbm_avg_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t

        # don't plot the clustermap with all genes
        if n_t < abg_avg_log.shape[0]:
            cm, _ = cluster_data_with_threshold(abg_avg_log.loc[amad_avg_log.index[:n_t]], col_colors=row_colours_avg)
            cm.gs.update(bottom=0.3)
            cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "nsc_gbm_avg_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
        d = clustering.dendrogram_with_colours(
            abg_avg_log.loc[amad_avg_log.index[:n_t]],
            row_colours_avg,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False
        )
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        fname = "nsc_gbm_avg_corrplot_log_top%d_by_mad.{ext}" % n_t
        cm = clustering.plot_correlation_clustermap(
            abg_avg_log.loc[amad_avg_log.index[:n_t]],
            row_colors=row_colours_avg,
            n_gene=n_t,
        )
        plt.setp(cm.ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
        plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
        cm.gs.update(bottom=0.35, right=0.65)
        cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

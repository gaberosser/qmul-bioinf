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
from rnaseq import general, loader


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


def compare_mad_genes(data, samples1, samples2=None, max_ng=10000, step=50):
    """
    Test the overlap of two samples in terms of the top genes drawn (by descending MAD).
    :param data:
    :param samples1:
    :param samples2: Second sample set. If None, use all samples.
    :return:
    """
    d1 = data.loc[:, samples1]
    if samples2 is not None:
        d2 = data.loc[:, samples2]
    else:
        d2 = data

    mad1 = transformations.median_absolute_deviation(d1).sort_values(ascending=False)
    mad2 = transformations.median_absolute_deviation(d2).sort_values(ascending=False)

    ng = np.arange(step, max_ng + 1, step)

    iu = [(
        mad1.index[:i].intersection(mad2.index[:i]),
        mad2.index[:i].union(mad2.index[:i])
    ) for i in ng]
    iu_pct = np.array([
        t[0].size / float(t[1].size) * 100 for t in iu
    ])
    return ng, iu_pct



def hc_plot_dendrogram_vary_n_gene(
        data,
        row_colours,
        n_gene_arr=(1000, 2000, 3000, 5000, 10000),
        metric='correlation'
):
    """
    For each value in n_gene_arr, plot a dendrogram showing the result of hierarchical clustering of the data using
    that many genes (selected in descending MAD order)
    :param data: Cols are samples, rows are genes (or similar)
    :param row_colours: As passed to dendrogram routine
    :param n_gene_arr: The values to test
    :return:
    """
    mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
    fig_dict = {}
    for ng in n_gene_arr:
        the_dat = data.loc[mad.index[:ng]]
        d = clustering.dendrogram_with_colours(
            the_dat,
            row_colours,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False,
            metric=metric
        )
        fig_dict[ng] = d
    return fig_dict




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

    patient_obj = loader.load_by_patient(pids, source='salmon', units=units)
    patient_data = patient_obj.data

    # discard GBM and unused 016 iNSC
    patient_data = patient_data.loc[:, ~patient_data.columns.str.contains('GBM')]
    patient_data = patient_data.drop(['DURA061_NSC_N6_P4', 'DURA061_NSC_N1_P5'], axis=1)

    # discard mitochondrial genes
    if remove_mt:
        idx = ~patient_data.index.isin(mt_ensg)
        pdbg = patient_data.loc[idx]
        # renorm
        if units == 'tpm':
            pdbg = pdbg.divide(pdbg.sum(), axis=1) * 1e6
    else:
        pdbg = patient_data

    # discard genes expressed at low values
    idx = (pdbg > min_val).sum(axis=1) > min_n
    pdbg = pdbg.loc[idx]

    if units == 'estimated_counts':
        # here we can normalise by library size if desired
        pass

    ax = hist_logvalues(patient_data, thresholds=[min_val])
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

        d = clustering.dendrogram_with_colours(
            pdbg_log.loc[mad_log.index[:n_t]],
            row_colours,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False
        )
        d['fig'].savefig(os.path.join(outdir, fname_log.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # bring in reference data
    # IDs (if req), lab (appears in label), loader

    ref_dats = [
        (None, 'Barres et al.', loader.load_references('GSE73721', source='salmon', units=units),),
        (None, 'Caren et al.', loader.load_references('E-MTAB-3867', source='salmon', units=units),),
        (None, 'Yang et al.', loader.load_references('GSE80732', source='salmon', units=units),),
        (None, 'Shahbazi et al.', loader.load_references('GSE64882', source='salmon', units=units),),
        (None, 'Li et al.', loader.load_references('GSE84166', source='salmon', units=units),),
        (None, 'Duan et al.', loader.load_references('GSE61794', source='salmon', units=units),),
        (None, 'Bago et al.', loader.load_references('GSE92839', source='salmon', units=units),),
        (None, 'Kelley and Rinn', loader.load_references('GSE38993', source='salmon', units=units),),
        (None, 'ENCODE Wold', loader.load_references('encode_roadmap/ENCSR000EYP', source='salmon', units=units),),
        (['H1 PSC Costello_PSB'], 'ENCODE Costello', loader.load_references('encode_roadmap/ENCSR950PSB', source='salmon', units=units),),
        (['ENCSR000COU rep 1', 'ENCSR000COU rep 2'], 'ENCODE Gingeras', loader.load_references('encode_roadmap/ENCSR000COU', source='salmon', units=units),),
        (None, 'ENCODE Gingeras', loader.load_references('encode_roadmap/ENCSR490SQH', source='salmon', units=units),),
        (None, 'ENCODE Gingeras', loader.load_references('encode_roadmap/ENCSR244ISQ', source='salmon', units=units),),
        (['H1 PSC Ecker_WQY'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR670WQY', source='salmon', units=units),),
        (['H1 PSC Ecker_RSE'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR043RSE', source='salmon', units=units),),
        (['H1 NPC Ecker_XUX'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR977XUX', source='salmon', units=units),),
        (['H1 NPC Ecker_EET'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR572EET', source='salmon', units=units),),
    ]

    ref_arr = []
    ref_labels = []
    ref_cols = []
    batches = []

    for i, r in enumerate(ref_dats):
        the_dat = r[2].data
        if r[0] is not None:
            the_cols = r[0]
        else:
            the_cols = the_dat.columns
        ref_cols.extend(the_cols)
        ref_labels.extend(the_dat.columns)
        ref_arr.append(the_dat)
        batches.extend([r[1]] * the_dat.shape[1])

    ref = pd.concat(ref_arr, axis=1)
    ref.columns = ref_cols
    batches = pd.Series(batches, index=ref_cols)
    labels = pd.Series(ref_labels, ref_cols)
    # ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    # discard Barres irrelevant samples
    # discard immortalised cell line
    # discard fibroblasts (careful, need to be selective here)

    to_discard = [
        'Fetal ctx', 'tumor', 'ctx endo', 'whole cortex', 'myeloid', 'LM-NSC', r'^fibroblast', 'NHF1-hTERT fibroblasts',
        'foreskin fibroblast', 'lung fibroblast'
    ]
    for td in to_discard:
        the_idx = ~ref.columns.str.contains(td)
        ref = ref.loc[:, the_idx]
        batches = batches.loc[the_idx]
        labels = labels.loc[the_idx]

    # Take mean over various replicates (for simplicity)
    to_aggr = [
        ('Pollard NSC', 'Fetal NSC'),
        (r'H1-[12]', 'H1 PSC'),
        (r'H1-EPS', 'H1 EPSC'),
        (r'ES1-EPS', 'ES1 EPSC'),
        (r'H9_NSC', 'H9 NSC'),
        (r'INSC fibroblast', 'Fibroblast-derived iNSC'),
        (r'NSC008', 'Fetal NSC'),
        (r'fetal NSC rep', 'Fetal NSC'),
        (r'ReNcell CX', 'ReNcell CX NSC'),
        (r'NHF1-hTERT hiNSC', 'Fibroblast-derived iNSC'),
        (r'ENCSR000COU rep', 'H1 PSC'),
        (r'H1-hESC rep', 'H1 PSC'),
        (r'H7-hESC rep', 'H7 PSC'),
        (r'H9-hNPC rep', 'H9 NPC'),
        (r'H1 NPC Ecker', 'H1 NPC'),
        (r'H1 PSC Ecker', 'H1 PSC'),
    ]

    for srch, repl in to_aggr:
        idx = ref.columns.str.contains(srch)

        new_col = ref.loc[:, idx].mean(axis=1)
        new_batch = batches.loc[idx].iloc[0]
        new_label = repl

        ref = ref.loc[:, ~idx]
        ref.insert(ref.shape[1], repl, new_col, allow_duplicates=True)
        batches = batches.loc[~idx]
        batches = batches.append(pd.Series([new_batch], index=[repl]))
        labels = labels.loc[~idx]
        labels = labels.append(pd.Series([new_label], index=[repl]))

    # append study name to labels
    new_labels = ["%s (%s)" % (t) for t in zip(labels.index, batches.values)]

    ref.columns = new_labels
    batches.index = new_labels

    # now let's try clustering everything together
    abg = pd.concat((patient_data, ref), axis=1)

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
    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#96daff'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fibroblast')] = '#fff89e'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fetal')] = 'yellow'
    # row_colours_all.loc[row_colours_all.index.str.contains('H1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('ES1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('PSC')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Hippocamp.* astro')] = '#e78ac3'
    row_colours_all.loc[row_colours_all.index.str.contains(r'ctx .*astro')] = '#b3b3b3'
    row_colours_all.loc[row_colours_all.index.str.contains(r'oligo')] = '#fccde5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange

    plt_dict = hc_plot_dendrogram_vary_n_gene(abg_log, row_colours_all)
    for ng, x in plt_dict.items():
        fname = "all_samples_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # for n_t in n_gene_try:
    #     fname = "all_samples_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t
    #
    #     cm, mad_all = cluster_logdata_with_threshold(abg, n=n_t, eps=eps, col_colors=row_colours_all)
    #     cm.gs.update(bottom=0.3)
    #     cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
    #
    #     fname = "all_samples_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
    #     d = clustering.dendrogram_with_colours(
    #         abg_log.loc[amad_log.index[:n_t]],
    #         row_colours_all,
    #         fig_kws={'figsize': (5.5, 10)},
    #         vertical=False
    #     )
    #     d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    # remove Bares and replot
    # abg_nobarres = abg.loc[:, ~abg.columns.str.contains('Barres')]

    abg_nobarres_log = abg_log.loc[:, ~abg_log.columns.str.contains('Barres')]
    amad_nobarres_log = transformations.median_absolute_deviation(abg_nobarres_log).sort_values(ascending=False)
    row_colours_nobarres = row_colours_all.loc[abg_nobarres_log.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(abg_nobarres_log, row_colours_nobarres)
    for ng, x in plt_dict.items():
        fname = "no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')


    # kick out fetal samples
    x = abg_nobarres_log.loc[:, ~abg_nobarres_log.columns.str.contains('Fetal')]
    x = x.loc[:, ~x.columns.str.contains('Fibroblast-derived iNSC')]
    x = x.loc[:, ~x.columns.str.contains('ReNcell')]
    abg_nofetal_nobarres_log = x
    rc = row_colours_all.loc[x.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(x, rc)
    for ng, x in plt_dict.items():
        fname = "no_fetal_no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')


    # Encode samples, no Barres
    x = abg_nobarres_log.loc[:, ~abg_nobarres_log.columns.str.contains('Kelley')]
    x = x.loc[:, ~x.columns.str.contains('Duan')]
    x = x.loc[:, ~x.columns.str.contains('Bago')]
    abg_encode_nobarres_log = x

    rc = row_colours_all.loc[x.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(x, rc)
    for ng, x in plt_dict.items():
        fname = "encode_no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # Encode samples, no Barres, no fetal
    x = abg_encode_nobarres_log.loc[:, ~abg_encode_nobarres_log.columns.str.contains('Fetal')]
    x = x.loc[:, ~x.columns.str.contains('Fibroblast-derived iNSC')]
    x = x.loc[:, ~x.columns.str.contains('ReNcell')]
    abg_encode_nofetal_nobarres_log = x

    rc = row_colours_all.loc[x.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(x, rc)
    for ng, x in plt_dict.items():
        fname = "encode_no_fetal_no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # no Encode, no Barres
    x = abg_nobarres_log.loc[:, ~abg_nobarres_log.columns.str.contains('ENCODE')]
    abg_noencode_nobarres_log = x
    rc = row_colours_all.loc[x.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(x, rc)
    for ng, x in plt_dict.items():
        fname = "no_encode_no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # no Encode, no Barres, no fetal
    x = abg_noencode_nobarres_log.loc[:, ~abg_noencode_nobarres_log.columns.str.contains('Fetal')]
    x = x.loc[:, ~x.columns.str.contains('Fibroblast-derived iNSC')]
    x = x.loc[:, ~x.columns.str.contains('ReNcell')]
    abg_noencode_no_fetal_nobarres_log = x
    rc = row_colours_all.loc[x.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(x, rc)
    for ng, x in plt_dict.items():
        fname = "no_encode_no_fetal_no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

    plt.draw()
    plt.close('all')

    # for n_t in n_gene_try:
    #     fname = "no_barres_clustering_by_gene_log_corr_top%d_by_mad.{ext}" % n_t
    #
    #     cm, mad_all = cluster_logdata_with_threshold(abg_nobarres, n=n_t, eps=eps, col_colors=row_colours_nobarres)
    #     cm.gs.update(bottom=0.3)
    #     cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
    #
    #     fname = "no_barres_dendrogram_log_corr_top%d_by_mad.{ext}" % n_t
    #     d = clustering.dendrogram_with_colours(
    #         abg_nobarres_log.loc[amad_nobarres_log.index[:n_t]],
    #         row_colours_nobarres,
    #         fig_kws={'figsize': (5.5, 10)},
    #         vertical=False
    #     )
    #     d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)


    # more out of interest than anything else:
    # what proportion of the genes (ranked by MAD) are shared between the 'with Barres' and 'without Barres' sets?
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ng, iu_pct = compare_mad_genes(abg_log, samples1=abg_nobarres_log.columns)
    ax.plot(ng, iu_pct, label="With/without Barres")

    ng, iu_pct = compare_mad_genes(abg_log, samples1=abg_noencode_nobarres_log.columns, samples2=abg_nobarres_log.columns)
    ax.plot(ng, iu_pct, label="Without Barres, with/without Encode")

    ax.set_xlabel("Number of genes (descending MAD order)")
    ax.set_ylabel("Percentage common with and without Barres data")
    ax.legend(loc='lower right')
    fig.savefig(os.path.join(outdir, "gene_intersection_pct.png"), dpi=200)

    # now try with different sample choices


    # TODO: move this. PCA plot
    from plotting import pca
    from sklearn.decomposition import PCA

    p = PCA()
    y = p.fit_transform(abg_nobarres_log.transpose())

    subgroups = pd.Series('NA', index=abg_nobarres_log.columns)
    subgroups.loc[subgroups.index.str.contains('NPC')] = 'NPC'
    subgroups.loc[subgroups.index.str.contains('NSC')] = 'NPC'
    subgroups.loc[subgroups.index.str.contains('PSC')] = 'PSC'
    subgroups.loc[subgroups.index.str.contains('ESC')] = 'PSC'
    subgroups.loc[subgroups.index.str.contains('fetal')] = 'Fetal NSC'
    subgroups.loc[subgroups.index.str.contains('iNSC')] = 'iNSC'
    subgroups.loc[subgroups.index.str.contains(r'DURA[0-9]*_NSC')] = 'Dura iNSC'
    subgroups.loc[subgroups.index.str.contains('IPSC')] = 'Dura iPSC'
    sg_cmap = {
        'NPC': 'b',
        'PSC': 'g',
        'iNSC': 'y',
        'Dura iNSC': 'c',
        'IPSC': 'r',
        'Dura iPSC': 'm',
        'NA': 'k'
    }
    ax = pca.pca_plot_by_group_2d(y, subgroups, components=(0, 1), ellipses=False, colour_map=sg_cmap)
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "no_barres_pca_01.png"), dpi=200)
    ax = pca.pca_plot_by_group_2d(y, subgroups, components=(1, 2), ellipses=False, colour_map=sg_cmap)
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "no_barres_pca_12.png"), dpi=200)
    ax = pca.pca_plot_by_group_2d(y, subgroups, components=(0, 2), ellipses=False, colour_map=sg_cmap)
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "no_barres_pca_02.png"), dpi=200)
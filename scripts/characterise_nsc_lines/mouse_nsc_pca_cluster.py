from plotting import clustering, pca
from rnaseq import general, loader
from scripts.rnaseq import gtf_reader
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
from sklearn.decomposition import PCA
import os
import re


if __name__ == '__main__':
    source = 'salmon'
    # source = 'star'

    # units = 'estimated_counts'
    units = 'tpm'
    # units = 'cpm'
    # units = 'counts'

    # transform = 'vst'
    transform = 'log'

    # remove_mt = True
    remove_mt = False

    pca_add_sample_names = False

    outdir = unique_output_dir("mouse_nsc_pca_cluster", reuse_empty=True)
    n_gene_try = [1000, 2000, 3000, 5000][::-1]  # largest first, so we can reuse the MAD array

    if source == 'star':
        load_cls = loader.StarCountLoader
        load_kwargs = {}
    elif source == 'salmon':
        load_cls = loader.SalmonQuantLoader
        load_kwargs = {'units': units}
    else:
        raise ValueError("Unrecognised source %s" % source)

    if units == 'tpm':
        eps = .01
    elif units == 'estimated_counts':
        eps = .01
    elif units == 'cpm':
        eps = .01

    if source == 'star' and units == 'tpm':
        raise NotImplementedError("Unable to convert STAR counts into TPM.")

    if transform in {'vst', 'rlog'} and units not in {'counts', 'estimated_counts'}:
        raise AttributeError("vst and rlog transformations require counts as input")

    if remove_mt:
        mt_ensg = set(gtf_reader.get_mitochondrial('GRCm38r88'))

    loc = loader.RnaSeqFileLocations(
        root_dir=os.path.join(loader.RNASEQ_DIR, 'wtchg_p170390'),
        alignment_subdir='mouse',
        batch_id='wtchg_p170390'
    )

    kwargs = loc.loader_kwargs(source)
    kwargs.update(dict(
        tax_id=10090,
        samples=[u'mDura3N1human', u'mDura5N24Ahuman', u'mDura6N6human',
                 u'mDura3N1mouse', u'mDura5N24Amouse', u'mDura6N6mouse',
                 u'eNSC3mouse', u'eNSC5mouse', u'eNSC6mouse']
    ))
    kwargs.update(load_kwargs)
    obj1 = load_cls(**kwargs)

    loc = loader.RnaSeqFileLocations(
        root_dir=os.path.join(loader.RNASEQ_DIR, 'wtchg_p170506'),
        alignment_subdir='mouse',
        batch_id='wtchg_p170506'
    )

    kwargs = loc.loader_kwargs(source)
    kwargs.update(dict(
        tax_id=10090,
        samples=['eNSC3med', 'eNSC5med', 'eNSC6med',]
    ))
    kwargs.update(load_kwargs)
    obj2 = load_cls(**kwargs)

    our_obj = loader.loader.MultipleBatchLoader([obj1, obj2])

    our_dat = our_obj.data

    # discard mitochondrial genes
    if remove_mt:
        idx = ~our_dat.index.isin(mt_ensg)
        our_dat = our_dat.loc[idx]
        # renorm
        if units == 'tpm':
            our_dat = our_dat.divide(our_dat.sum(), axis=1) * 1e6

    if units == 'cpm':
        print "Converting to CPM. NB: this will not work if the input data are not counts!"
        # This line will have no effect if the input data are TPM
        our_dat = our_dat.divide(our_dat.sum(axis=0), axis=1) * 1e6

    if transform == 'vst':
        our_dat = transformations.vst_blind(our_dat)
    elif transform == 'rlog':
        our_dat = transformations.rlog_blind(our_dat)
    elif transform == 'log':
        our_dat = np.log2(our_dat + eps)

    mad = transformations.median_absolute_deviation(our_dat).sort_values(ascending=False)

    # PCA
    p = PCA()

    subgroups = pd.Series('NA', index=our_dat.columns)
    subgroups.loc[subgroups.index.str.contains(r'eNSC.med')] = 'eNSC regular media'
    subgroups.loc[subgroups.index.str.contains(r'eNSC.mouse')] = 'eNSC mouse induction media'
    subgroups.loc[subgroups.index.str.contains(r'mDura.*human')] = 'iNSC human media'
    subgroups.loc[subgroups.index.str.contains(r'mDura.*mouse')] = 'iNSC mouse media'

    sg_cmap = {
        'eNSC regular media': '#146dff',
        'eNSC mouse induction media': '#96beff',
        'iNSC human media': '#008702',
        'iNSC mouse media': '#7de07f',
        'NA': 'k'
    }

    for n_t in n_gene_try:
        y = p.fit_transform(our_dat.loc[mad.index[:n_t]].transpose())
        expl_var_ratio = p.explained_variance_ratio_
        ax = pca.pca_plot_by_group_2d(y, subgroups, components=(0, 1), ellipses=False, colour_map=sg_cmap)
        ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.4)
        ax.set_xlabel('PCA component 1 (%.1f%%)' % (expl_var_ratio[0] * 100.))
        ax.set_ylabel('PCA component 2 (%.1f%%)' % (expl_var_ratio[1] * 100.))
        ax.set_aspect('auto')
        ax.figure.tight_layout()
        ax.figure.savefig(os.path.join(outdir, "pca_top%d_by_mad.png" % n_t), dpi=200)
        for i in range(len(our_dat.columns)):
            ax.text(y[i, 0], y[i, 1], our_dat.columns[i], fontsize=10)
        ax.figure.tight_layout(pad=3)
        ax.figure.savefig(os.path.join(outdir, "pca_top%d_by_mad_with_names.png" % n_t), dpi=200)

    row_colours = pd.DataFrame('gray', index=our_dat.columns, columns=[''])
    row_colours.loc[row_colours.index.str.contains(r'eNSC[0-9]med')] = '#66c2a5'
    row_colours.loc[row_colours.index.str.contains(r'eNSC[0-9]mouse')] = '#fc8d62'
    row_colours.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*mouse')] = '#8da0cb'
    row_colours.loc[row_colours.index.str.contains(r'mDura.[AN0-9]*human')] = '#e78ac3'

    for n_t in n_gene_try:
        fname = "clustering_by_gene_corr_log_top%d_by_mad.{ext}" % n_t

        d = clustering.dendrogram_with_colours(our_dat.loc[mad.index[:n_t]], row_colours, fig_kws={'figsize': (10, 5.5)})
        d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)

        cm = clustering.plot_clustermap(
            our_dat.loc[mad.index[:n_t]],
            cmap='RdBu_r',
            metric='correlation',
            col_colors=row_colours
        )
        cm.gs.update(bottom=0.2)
        cm.savefig(os.path.join(outdir, "clustermap_by_gene_corr_log_top%d_by_mad.png" % n_t), dpi=200)
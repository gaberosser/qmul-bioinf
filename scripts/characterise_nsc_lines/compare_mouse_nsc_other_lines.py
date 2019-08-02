import os
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from plotting import clustering
from rnaseq import loader
from scripts.rnaseq import gtf_reader
from settings import RNASEQ_DIR
from stats import transformations
from utils.output import unique_output_dir


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
    # source = 'salmon'
    source = 'star'

    # units = 'estimated_counts'
    # units = 'tpm'
    units = 'counts'

    # transform = 'log'
    # transform = 'rlog'
    transform = 'vst'

    # remove_mt = True
    remove_mt = False

    outdir = unique_output_dir("salmon_insc_mouse", reuse_empty=True)
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
    if units == 'cpm':
        eps = .01
    elif units == 'estimated_counts':
        eps = 1
    elif units == 'counts':
        eps = 1

    if remove_mt:
        mt_ensg = set(gtf_reader.get_mitochondrial('GRCm38r88'))

    if source == 'star' and units == 'tpm':
        raise NotImplementedError("Unable to convert STAR counts into TPM.")

    if transform in {'vst', 'rlog'} and units not in {'counts', 'estimated_counts'}:
        raise AttributeError("vst and rlog transformations require counts as input")

    loc = loader.RnaSeqFileLocations(
        root_dir=os.path.join(loader.RNASEQ_DIR, 'wtchg_p170390'),
        alignment_subdir='mouse',
        batch_id='wtchg_p170390'
    )
    # kwargs = loc.loader_kwargs('salmon')
    kwargs = loc.loader_kwargs(source)
    kwargs.update(dict(
        tax_id=10090,
        samples=[u'mDura3N1human', u'mDura5N24Ahuman', u'mDura6N6human']
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

    ## TODO: move all of this to later

    # our_dat = our_obj.data
    #
    # # discard mitochondrial genes
    # if remove_mt:
    #     idx = ~our_dat.index.isin(mt_ensg)
    #     our_dat = our_dat.loc[idx]
    #     # renorm
    #     if units == 'tpm':
    #         our_dat = our_dat.divide(our_dat.sum(), axis=1) * 1e6
    #
    # if units == 'cpm':
    #     print "Converting to CPM. NB: this will not work if the input data are not counts!"
    #     # This line will have no effect if the input data are TPM
    #     our_dat = our_dat.divide(our_dat.sum(axis=0), axis=1) * 1e6
    #
    # if transform == 'vst':
    #     our_dat = transformations.vst_blind(our_dat)
    # elif transform == 'rlog':
    #     our_dat = transformations.rlog_blind(our_dat)
    # elif transform == 'log':
    #     our_dat = np.log2(our_dat + eps)

    # bring in reference data
    ref_names = [
        ('Brandner lab', 'brandner_mouse_gic'),
        ('Pten/P53 study', 'GBM_Pten_P53'),
        ('Zhang et al., reprog.', 'GSE78938'),
        ('Liu et al.', 'GSE96950'),
        ('Wapinski et al.', 'GSE43916'),
        ('Friedmann-Morvinski et al.', 'GSE73127'),
        ('Friedmann-Morvinski et al.', 'GSE64411/trimgalore'),
        ('Zhang et al.', 'GSE52564'),
        ('Chen et al.', 'GSE52125'),
        ('Yanez et al.', 'GSE88982'),
        ('Lynch', 'GSE78795'),
        ('Moyon et al.', 'GSE66029'),
        ('Schmid et al.', 'GSE75592'),
        ('Srinivasan et al.', 'GSE75246'),
    ]

    ref_objs = []
    for t in ref_names:
        loc = loader.RnaSeqFileLocations(
            root_dir=os.path.join(RNASEQ_DIR, t[1]),
            alignment_subdir='mouse',
            batch_id=t[0]
        )
        kwargs = loc.loader_kwargs(source)
        kwargs.update(dict(
            tax_id=10090,
        ))
        kwargs.update(load_kwargs)
        try:
            this_obj = load_cls(**kwargs)
            this_obj.meta.index = ["%s (%s)" % (x, t[0]) for x in this_obj.meta.index]
            this_obj.data.columns = ["%s (%s)" % (x, t[0]) for x in this_obj.data.columns]
            ref_objs.append(this_obj)

        except Exception as exc:
            print "Error: failed to load reference %s" % t[0]
            print repr(exc)

    ref_obj = loader.MultipleBatchLoader(ref_objs)

    # ref_obj = loader.load_references(
    #     [t[1] for t in ref_names],
    #     tax_id=10090,
    #     source=source,
    #     batch_names=[t[0] for t in ref_names],
    # )

    # remove unneeded samples

    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('Normal brain')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('GBM')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('TrNeuron')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('TrAstrocyte')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('Tumour')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('MP and cMoP')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains('LPS')]
    ref_obj.meta = ref_obj.meta.loc[~ref_obj.meta.index.str.contains(r'day [148]')]

    ref_obj.data = ref_obj.data.loc[:, ref_obj.meta.index]
    ref_obj.batch_id = ref_obj.meta.batch

    # aggregate
    ii, batch = ref_obj.batch_id.factorize()
    new_dat = {}
    new_meta = {}
    for i, b in enumerate(batch):
        the_meta = ref_obj.meta.loc[ii == i]
        the_dat = ref_obj.data.loc[:, ii == i]

        # average over duplicates
        to_keep = []
        a = {}
        m = {}
        na = {}
        do_avg = False
        for col in the_dat.columns:
            if ' repl ' in col:
                do_avg = True
            else:
                to_keep.append(col)
                continue
            t = re.sub(r' repl [0-9]*', '', col)
            if t in a:
                a[t] += the_dat.loc[:, col]
                na[t] += 1
            else:
                a[t] = the_dat.loc[:, col].copy()
                m[t] = the_meta.loc[col]
                na[t] = 1
        for t in a:
            a[t] = a[t] / float(na[t])

        if do_avg:
            # replace the affected columns
            the_meta = the_meta.loc[to_keep]
            the_dat = the_dat.loc[:, to_keep]
            for t in m:
                the_meta.loc[t] = m[t]
                the_dat.insert(0, t, a[t])
            the_dat = the_dat.loc[:, the_meta.index]

            new_dat[b] = the_dat
            new_meta[b] = the_meta

    for b in new_dat:
        idx = ref_obj.meta.batch != b
        ref_obj.meta = ref_obj.meta.loc[idx]
        ref_obj.data = ref_obj.data.loc[:, idx]

        ref_obj.meta = pd.concat((ref_obj.meta, new_meta[b]), axis=0)
        ref_obj.data = pd.concat((ref_obj.data, new_dat[b]), axis=1)
        ref_obj.batch_id = ref_obj.meta.batch

    all_obj = loader.MultipleBatchLoader([our_obj, ref_obj])

    dat = all_obj.data

    if transform == "vst":
        # cast to int as aggregation will have broken this
        dat = transformations.vst_blind(dat.astype(int))

    mad = transformations.median_absolute_deviation(dat)

    row_colours_all = pd.DataFrame('gray', index=dat.columns, columns=['Cell type',])

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
    row_colours_all.loc[row_colours_all.index.str.contains(r'eNSC[0-9]med'), 'Cell type'] = '#96ff9d' # pale green
    row_colours_all.loc[row_colours_all.index.str.contains(r'eNSC[0-9]mouse'), 'Cell type'] = '#008408' # dark green
    row_colours_all.loc[row_colours_all.index.str.contains(r'mDura.[AN0-9]*mouse'), 'Cell type'] = '#3543ff' # dark blue
    row_colours_all.loc[row_colours_all.index.str.contains(r'mDura.[AN0-9]*human'), 'Cell type'] = '#c4c8ff' # pale blue

    d = clustering.dendrogram_with_colours(
        dat.loc[mad.index[:3000]],
        row_colours_all,
        fig_kws={'figsize': (5.5, 10)},
        vertical=False
    )


    # cm = clustering.plot_correlation_clustermap(
    #     dat.loc[mad.index[:3000]],
    #     row_colors=row_colours_all,
    # )
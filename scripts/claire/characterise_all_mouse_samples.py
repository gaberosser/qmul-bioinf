from rnaseq import loader
from plotting import clustering, common, scatter
from stats import transformations
import pandas as pd
import numpy as np
import os, sys
from utils import output
import re
from copy import copy
from sklearn.decomposition import pca


def plot_pca(
        dat,
        colour_subgroups,
        p=None,
        components=(0, 1),
        marker_subgroups=None,
        ax=None,
        colour_map=None,
        marker_map=None,
        **kwargs
):
    if p is None:
        p = pca.PCA()
        pca_data = p.fit_transform(dat.transpose())
    else:
        pca_data = p.transform(dat.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    ax = scatter.scatter_with_colour_and_markers(
        pca_data[:, components],
        colour_subgroups=colour_subgroups,
        colour_map=colour_map,
        marker_subgroups=marker_subgroups,
        marker_map=marker_map,
        ax=ax,
        **kwargs
    )

    ax.set_xlabel("PCA component %s (%.1f%%)" % (components[0] + 1, variance_explained[components[0]]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (components[1] + 1, variance_explained[components[1]]))

    return p, ax


if __name__ == '__main__':
    min_val = 1
    n_above_min = 3
    n_gene_by_mad = 5000
    eps = 0.01  # offset to use when applying log transform
    source = 'salmon'
    quantile_norm = None

    outdir = output.unique_output_dir()

    obj1 = loader.load_references('wtchg_p170390', source=source, tax_id=10090)
    samples = ['eNSC%dmouse' % i for i in (3, 5, 6)] \
    + ['mDura%smouse' % i for i in ('3N1', '5N24A', '6N6')] \
    + ['mDura%shuman' % i for i in ('3N1', '5N24A', '6N6')]
    obj1.filter_samples(obj1.meta.index.isin(samples))

    obj2 = loader.load_references('wtchg_p170506', source=source, tax_id=10090)
    samples = ['eNSC%dmed' % i for i in (3, 5, 6)]
    obj2.filter_samples(obj2.meta.index.isin(samples))

    obj3 = loader.load_references('wtchg_p180443', source=source, tax_id=10090)

    obj = loader.loader.MultipleBatchLoader([obj1, obj2, obj3])

    the_dat = np.log2(obj.data + eps)
    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    studies = obj.meta.batch.copy()

    mouse_number = pd.Index([int(re.sub(r'^[a-zA-Z]*(?P<n>[356]).*', '\g<n>', t)) for t in obj.meta.index])

    p, ax = plot_pca(the_dat, mouse_number, marker_subgroups=studies)
    ax.figure.savefig(os.path.join(outdir, "all_mouse_nsc_by_individual.png"), dpi=200)

    prep_type = np.array(['med'] * len(obj.meta.index))
    ix = [bool(re.search(r'(mouse|mus)', t)) for t in obj.meta.index]
    prep_type[ix] = 'mouse'
    ix = [bool(re.search(r'hum', t)) for t in obj.meta.index]
    prep_type[ix] = 'human'
    prep_type = pd.Index(prep_type)

    p, ax = plot_pca(the_dat, prep_type, marker_subgroups=studies)
    ax.figure.savefig(os.path.join(outdir, "all_mouse_nsc_by_prep.png"), dpi=200)

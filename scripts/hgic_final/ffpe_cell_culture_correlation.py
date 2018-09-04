"""
Added on 4th Sept 2018

Modified version of scripts.publications.cruk_grant_jan_2018.ffpe_cell_culture_correlation.

Load RNA-Seq gene expression data and produce a pairwise correlation plot between FFPE and cell culture samples.

TODO (in original script first): Try similar approach with the methylation data?
"""

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import rnaseq.loader
from rnaseq import filter, NH_ID_TO_PATIENT_ID_MAP, general
from stats import transformations
from scipy import stats
import numpy as np
from utils import output
import os, sys
import consts


if __name__ == "__main__":
    source = 'salmon'
    pids = consts.PIDS
    apply_qn = True
    # apply_qn = False
    dist_metric = 'pearson'
    # dist_metric = 'spearman'
    remove_mt = True
    min_tpm = 1.
    eps = .1  # offset for log transform

    rna_ff_samples = [
        'NH15_1661DEF2C',
        'NH15_1877_SP1C',
        'NH15_2101_DEF1A',
        'NH16_270_DEF1Ereplacement',
        'NH16_616DEF1B',
        'NH16_677_SP1A',
        'NH16_2063_DEF1Areplacement',
        'NH16_2214DEF1A',
        'NH16_2255DEF1B2',
        'NH16_2806DEF3A1'
    ]

    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)

    if remove_mt:
        mt_ens = general.get_mitochondrial(9606)

    rna_cc_obj = rnaseq.loader.load_by_patient(pids, source=source, include_control=False)
    rna_ff_obj = rnaseq.loader.load_by_patient(pids, source=source, include_control=False, type='ffpe')

    # filter
    ix = rna_ff_obj.meta.index.isin(rna_ff_samples)
    rna_ff_obj.filter_samples(ix)

    ix = rna_cc_obj.meta.type == 'GBM'
    rna_cc_obj.filter_samples(ix)

    # add NH ID and patient ID to FFPE
    nh_id = rna_ff_obj.meta.index.str.replace(r'(_?)(DEF|SP).*', '')
    p_id = [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]
    rna_ff_obj.meta.insert(0, 'nh_id', nh_id)
    rna_ff_obj.meta.insert(0, 'patient_id', p_id)

    # extract data
    rna_ff_dat = rna_ff_obj.data.copy()
    rna_cc_dat = rna_cc_obj.data.copy()

    # update labels
    rna_ff_dat.columns = ["GBM%s" % t for t in p_id]

    cols = []
    for t in rna_cc_dat.columns:
        p, q = t.split('_')
        q = q.replace('n', ' & P')
        cols.append("%s (%s)" % (p, q))
    rna_cc_dat.columns = cols

    # filter low expression out
    rna_ff_dat = filter.filter_by_cpm(rna_ff_dat, min_cpm=min_tpm, min_n_samples=2)
    rna_cc_dat = filter.filter_by_cpm(rna_cc_dat, min_cpm=min_tpm, min_n_samples=2)

    # reduce to matching probes
    probes = rna_cc_dat.index.intersection(rna_ff_dat.index)

    if remove_mt:
        probes = probes[~probes.isin(mt_ens)]

    rna_ff_dat = np.log2(rna_ff_dat.loc[probes] + eps)
    rna_cc_dat = np.log2(rna_cc_dat.loc[probes] + eps)

    # QN
    if apply_qn:
        rna_ff_dat = transformations.quantile_normalisation(rna_ff_dat)
        rna_cc_dat = transformations.quantile_normalisation(rna_cc_dat)

    # correlation plot
    pdist = pd.DataFrame(index=rna_ff_dat.columns, columns=rna_cc_dat.columns, dtype=float)
    for ff in rna_ff_dat.columns:
        for cc in rna_cc_dat.columns:
            if dist_metric == 'pearson':
                pdist.loc[ff, cc] = stats.pearsonr(rna_ff_dat[ff], rna_cc_dat[cc])[0]
            elif dist_metric == 'spearman':
                pdist.loc[ff, cc] = stats.spearmanr(rna_ff_dat[ff], rna_cc_dat[cc]).correlation
            else:
                raise NotImplementedError("Unsupported distance metric %s." % dist_metric)
    pdist_z = pdist.subtract(pdist.mean(axis=1), axis=0).divide(pdist.std(axis=1), axis=0)

    # custom norming: divide by row sum AND col sum
    # preserve min and max for later
    pmin = pdist.values.flatten().min()
    pmax = pdist.values.flatten().max()

    pdist_s = pdist.divide(pdist.sum(axis=1), 0).divide(pdist.sum(axis=0), 1)
    # pdist_s = pdist_z

    pmin_s = pdist_s.values.flatten().min()
    pmax_s = pdist_s.values.flatten().max()
    pdist_s = (pdist_s - pmin_s) / (pmax_s - pmin_s)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax = sns.heatmap(pdist_z, vmin=-1.5, vmax=1.5, ax=ax)
    vmin = 0.1
    vmax = 0.9
    ax = sns.heatmap(
        pdist_s,
        ax=ax,
        cmap='bwr',
        vmin=vmin,
        vmax=vmax,
        cbar_kws={'ticks': [vmin, vmax], 'label': '%s correlation' % dist_metric.capitalize()}
    )

    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    ax.figure.tight_layout(rect=(0.05, 0.05, 1., 1.))

    # make arbitrary cax scaling
    cax = [t for t in ax.figure.get_axes() if t is not ax][0]
    cax.yaxis.set_ticklabels(['Low', 'High'])

    ax.figure.savefig(os.path.join(outdir, "ffpe_cc_correlation_heatmap.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "ffpe_cc_correlation_heatmap.tiff"), dpi=200)

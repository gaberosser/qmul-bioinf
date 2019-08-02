import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import rankdata

from load_data import rnaseq_data
from stats import transformations
from utils.output import unique_output_dir
from utils.reference_genomes import ensembl_to_gene_symbol, gene_symbol_to_ensembl

if __name__ == "__main__":
    outdir = unique_output_dir("tom_qpcr", reuse_empty=True)
    ref = 'GIBCO_NSC_P4'

    obj = rnaseq_data.all_hgic_loader(annotate_by="Ensembl Gene ID")
    dat = obj.data.loc[obj.data.index.str.contains('ENSG')]
    dat = dat.loc[:, ~obj.meta.index.str.contains('DURA')]
    # normalised version (by number of aligned reads)
    dat_n = dat.divide(dat.sum(axis=0), axis=1) * 1e6

    # remove any absent / mostly absent genes
    median_count = dat_n.median(axis=1).sort_values()
    keep_idx = median_count.loc[median_count != 0].index

    dat = dat.loc[keep_idx]
    dat_n = dat_n.loc[keep_idx]
    median_count = median_count.loc[keep_idx]

    # remove any genes that are (mostly) absent in NSC
    nsc_missing = dat.loc[:, ref] < 10.
    dat = dat.loc[~nsc_missing]
    dat_n = dat_n.loc[~nsc_missing]
    median_count= median_count.loc[~nsc_missing]

    # levels of HKG across GBM and GIBCO
    hkg = ['ATP5B', 'GAPDH', 'ACTB']
    hkg_ens = gene_symbol_to_ensembl(hkg)

    label_symbols = hkg + ['BMI1']
    label_ens = gene_symbol_to_ensembl(label_symbols)

    hkg_dat = dat_n.loc[hkg_ens, sorted(dat_n.columns)]
    hkg_dat.index = pd.Index(hkg, name='')
    hkg_dat_rel = hkg_dat.divide(hkg_dat.loc[:, ref], axis=0)
    ax = hkg_dat_rel.transpose().plot.bar()
    ax.set_ylim([0, 3.4])
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, 'housekeeping_levels.png'), dpi=200)

    # identifying stable HKG
    ranked_count = pd.Series(rankdata(median_count, method='ordinal'), index=median_count.index)
    ranked_perc = ranked_count / float(ranked_count.shape[0])
    mad = transformations.median_absolute_deviation(dat_n)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ranked_perc, np.log10(median_count))
    ax.scatter(ranked_perc.loc[label_ens], np.log10(median_count.loc[label_ens]), c='r')
    for g, e in zip(label_symbols, label_ens):
        ax.text(ranked_perc.loc[e], np.log10(median_count.loc[e]), g)
    ax.set_xlabel("Abundance percentile")
    ax.set_ylabel("Log10 normalised abundance")
    ax.figure.savefig(os.path.join(outdir, 'hkg_abundance.png'), dpi=200)

    # show the total variation using fill_between
    min_count = dat_n.min(axis=1)
    max_count = dat_n.max(axis=1)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill_between(ranked_perc, np.log10(min_count + 0.01), np.log10(max_count + 0.01))
    ax.scatter(ranked_perc.loc[label_ens], np.log10(median_count.loc[label_ens]), c='r')
    for g, e in zip(label_symbols, label_ens):
        ax.text(ranked_perc.loc[e], np.log10(median_count.loc[e]), g)
    ax.set_xlabel("Abundance percentile")
    ax.set_ylabel("Log10 normalised abundance")
    ax.figure.savefig(os.path.join(outdir, 'hkg_abundance_range.png'), dpi=200)

    # median absolute deviation *from GIBCO*
    ref_count_n = dat_n.loc[:, ref]
    med_dev_nsc = np.abs(dat_n.subtract(ref_count_n, axis=0).loc[:, ~dat_n.columns.str.contains('GIBCO')]).median(axis=1)
    med_dev_nsc_rel = med_dev_nsc / ref_count_n

    rel_dev_range = (max_count - min_count) / median_count

    # abundance vs variation
    # in advance, we define some thresholds
    mad_cutoff = 500.
    rel_cutoff = 0.3
    perc_cutoff = 0.99
    range_cutoff = 1.

    mad_candidates = ranked_perc.loc[
        (ranked_perc > perc_cutoff) & (mad < mad_cutoff)
    ].index

    rel_dev_candidates = med_dev_nsc_rel.loc[
        (ranked_perc > perc_cutoff) & (med_dev_nsc_rel < rel_cutoff)
    ].index

    range_candidates = rel_dev_range.loc[
        (ranked_perc > perc_cutoff) & (rel_dev_range < range_cutoff)
    ].index

    top1pct_idx = ranked_perc.loc[ranked_perc > perc_cutoff].sort_values().index

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ranked_perc.loc[top1pct_idx], mad.loc[top1pct_idx])
    ax.scatter(ranked_perc.loc[mad_candidates], mad.loc[mad_candidates], c='g')
    ax.scatter(ranked_perc.loc[hkg_ens], mad.loc[hkg_ens], c='r')
    for g, e in zip(hkg, hkg_ens):
        ax.text(ranked_perc.loc[e], mad.loc[e], g)
    ax.set_xlabel("Abundance percentile")
    ax.set_ylabel("Median absolute deviation")
    ax.set_title("%d genes meet MAD criteria" % mad_candidates.size)
    ax.figure.savefig(os.path.join(outdir, 'median_absolute_deviation.png'), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = ranked_perc.diff()[1]
    ax.bar(
        ranked_perc.loc[top1pct_idx],
        np.log10(max_count.loc[top1pct_idx] + 0.01) - np.log10(min_count.loc[top1pct_idx] + 0.01),
        width=width,
        bottom=np.log10(min_count.loc[top1pct_idx] + 0.01),
        color='b',
        zorder=1.,
        alpha=0.5
    )
    ax.scatter(ranked_perc.loc[range_candidates], np.log10(median_count.loc[range_candidates]), c='g', zorder=2)
    ax.scatter(ranked_perc.loc[hkg_ens], np.log10(median_count.loc[hkg_ens]), c='r', zorder=2)
    for g, e in zip(hkg, hkg_ens):
        ax.text(ranked_perc.loc[e], np.log10(median_count.loc[e]), g)
    ax.set_xlabel("Abundance percentile")
    ax.set_ylabel("Log10 normalised abundance")
    ax.figure.savefig(os.path.join(outdir, 'top_abundance_range.png'), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ranked_perc.loc[top1pct_idx], med_dev_nsc_rel.loc[top1pct_idx])
    ax.scatter(ranked_perc.loc[rel_dev_candidates], med_dev_nsc_rel.loc[rel_dev_candidates], c='g')
    ax.scatter(ranked_perc.loc[hkg_ens], med_dev_nsc_rel.loc[hkg_ens], c='r')
    for g, e in zip(hkg, hkg_ens):
        ax.text(ranked_perc.loc[e], med_dev_nsc_rel.loc[e], g)
    ax.set_ylim([0, 2])
    ax.set_xlabel("Abundance percentile")
    ax.set_ylabel("Relative median absolute difference from NSC")
    ax.set_title("%d genes meet relative NSC-MAD criteria" % rel_dev_candidates.size)
    ax.figure.savefig(os.path.join(outdir, 'relative_median_absolute_deviation_from_nsc.png'), dpi=200)

    final_candidates = ensembl_to_gene_symbol(
        mad_candidates.intersection(rel_dev_candidates).intersection(range_candidates)
    )
    print "Identified %d candidates" % final_candidates.size
    print '\n'.join(final_candidates)

    # now re-plot the first figure but with a subset of these
    new_hkg = ['GAPDH', 'ATP5B', 'ACTB', 'PPIA', 'H3F3B']
    new_hkg_ens = gene_symbol_to_ensembl(new_hkg)

    hkg_dat = dat_n.loc[new_hkg_ens, sorted(dat_n.columns)]
    hkg_dat.index = pd.Index(new_hkg, name='Housekeeping gene')

    hkg_dat_rel = hkg_dat.divide(hkg_dat.loc[:, ref], axis=0)
    cols = [ref] + sorted(hkg_dat_rel.columns[hkg_dat_rel.columns != ref])
    hkg_dat_rel = hkg_dat_rel.loc[:, cols]

    ax = hkg_dat_rel.transpose().plot.bar()
    ax.axhline(1.0, ls='--', c='k')
    ax.set_ylim([0, 3.4])
    ax.legend(loc='upper left')

    plt.tight_layout()
    leg = ax.get_legend()
    fr = leg.get_frame()
    leg.set_frame_on(True)
    fr.set_color('w')
    fr.set_alpha(0.4)
    ax.figure.savefig(os.path.join(outdir, 'hkg_levels.png'), dpi=200)

    # optionally add mean values
    ax.plot(hkg_dat_rel.mean(axis=0).values, 'ko')
    ax.figure.savefig(os.path.join(outdir, 'hkg_levels_with_mean.png'), dpi=200)

    # new set of HKG
    new_hkg = ['GAPDH', 'ATP5B', 'ACTB', 'PPIA', 'H3F3B']
    new_hkg_ens = gene_symbol_to_ensembl(new_hkg)

    hkg_dat = dat_n.loc[new_hkg_ens, sorted(dat_n.columns)]
    hkg_dat.index = pd.Index(new_hkg, name='Housekeeping gene')




    bmi_dat = dat_n.loc[gene_symbol_to_ensembl('BMI1')]

    bmi_vs_hkg = 1. / hkg_dat.divide(bmi_dat, axis=1)
    bmi_vs_hkg_rel = bmi_vs_hkg.divide(bmi_vs_hkg.loc[:, ref], axis=0)
    # column order
    cols = [ref] + sorted(bmi_vs_hkg_rel.columns[bmi_vs_hkg_rel.columns != ref])
    bmi_vs_hkg_rel = bmi_vs_hkg_rel.loc[:, cols]

    ax = bmi_vs_hkg_rel.transpose().plot.bar()
    ax.axhline(1.0, ls='--', c='k')
    ax.set_ylim([0, 2.5])
    ax.legend(loc='upper left')
    plt.tight_layout()
    leg = ax.get_legend()
    leg.set_title("BMI1 level, relative to")
    leg.set_frame_on(True)
    fr = leg.get_frame()
    fr.set_color('w')
    fr.set_alpha(0.4)
    ax.figure.savefig(os.path.join(outdir, 'BMI1_vs_HKG.png'), dpi=200)

    # optionally add mean values
    ax.plot(bmi_vs_hkg_rel.mean(axis=0).values, 'ko')
    ax.figure.savefig(os.path.join(outdir, 'BMI1_vs_HKG_with_mean.png'), dpi=200)

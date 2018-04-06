import os
import pickle
from utils import output
import settings
import numpy as np
import csv
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    indir = os.path.join(settings.DATA_DIR_NON_GIT, 'chipseq', 'wtchg_p170710')
    tss_dir = os.path.join(indir, 'human', 'bt2_alignment', 'transcript_tss')
    outdir = output.unique_output_dir("p170710_chipseq_analysis")

    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, index_col=0, header=0)

    traces = {}
    for pid in meta.index:
        the_name = meta.loc[pid, 'sample']
        the_fn = os.path.join(tss_dir, str(pid), "%s.trace" % str(pid))
        with open(the_fn, 'rb') as f:
            the_vals = np.array(f.readline().split(', ')).astype(float)
        traces[the_name] = the_vals

    traces = pd.DataFrame.from_dict(traces)
    pad = (traces.shape[0] - 1) / 2  # could check this is as expected
    coords = np.arange(-pad, pad + 1)
    traces.index = coords

    traces = traces.loc[:, meta.loc[:, 'sample']]

    # 1) One plot per input sample

    chip_target_colours = {
        'H3K36me3': '#1b9e77',
        'H3K4me3': '#d95f02',
        'H3K27ac': '#7570b3',
        'H3K27me3': '#e7298a',
        'BMI1': '#66a61e',
        'CHD7': '#e6ab02',
    }

    chip_sample_colours = {
        'GBM017': '#a6cee3',
        'GBM061': '#1f78b4',
        'GBM054': '#b2df8a',
        'iNSC017': '#33a02c',
        'iNSC054': '#fb9a99',
        'iNSC061': '#e31a1c',
        '3021Scr': '#fdbf6f',
        '3021shBMI1': '#ff7f00',
        '3021shCHD7': '#cab2d6',
        '3021shBMI1CHD7': '#6a3d9a',
    }

    # keep a copy of the enrichment
    ne = {}

    for the_name in sorted(chip_sample_colours.keys()):
        input_idx = (meta.input_sample == the_name) & (meta.loc[:, 'ChIP target'] == 'Input')
        if input_idx.sum() != 1:
            raise ValueError("Unable to find input data unambiguously for input sample %s" % the_name)
        else:
            input_idx = np.where(input_idx)[0][0]

        dat_input = traces.loc[:, meta.iloc[input_idx].loc['sample']]
        norm_input = meta.flagstat_reads_mapped.iloc[input_idx]

        chip_idx = (meta.input_sample == the_name) & (meta.loc[:, 'ChIP target'] != 'Input')
        dat_chip = traces.loc[:, meta.loc[chip_idx, 'sample']]
        norm_chip = meta.flagstat_reads_mapped[chip_idx]
        norm_chip.index = meta.loc[chip_idx, 'sample']

        chip_conditions = meta.loc[chip_idx, 'ChIP target']

        norm_factors = norm_input / norm_chip.astype(float)

        enrichment = dat_chip.divide(dat_input, axis=0).multiply(norm_factors, axis=1)
        enrichment.columns = chip_conditions

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for c in enrichment:
            ne.setdefault(c, {})[the_name] = enrichment.loc[:, c]
            ax.plot(enrichment.index, enrichment.loc[:, c].values, c=chip_target_colours[c], label=c)
        ax.legend(loc='upper right')
        ax.set_title(the_name)
        ax.set_xlabel('Position, relative to TSS')
        ax.set_ylabel('Normalised enrichment')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "%s.png" % the_name), dpi=200)

    plt.close('all')

    # 2) One plot per ChIP target

    for the_target in sorted(chip_target_colours.keys()):
        # get a list of samples to plot
        sample_names = meta.loc[meta.loc[:, 'ChIP target'] == the_target, 'input_sample'].values
        # retrieve relavent data
        to_plot = pd.DataFrame.from_dict(
            dict([(k, ne[the_target][k]) for k in sample_names])
        )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for c in to_plot:
            ax.plot(to_plot.index, to_plot.loc[:, c].values, c=chip_sample_colours[c], label=c)
        ax.legend(loc='upper right')
        ax.set_title(the_target)
        ax.set_xlabel('Position, relative to TSS')
        ax.set_ylabel('Normalised enrichment')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "%s.png" % the_target), dpi=200)
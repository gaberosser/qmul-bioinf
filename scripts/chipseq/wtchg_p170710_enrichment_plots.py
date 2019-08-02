import os
import pickle
from utils import output
import settings
import numpy as np
import csv
import copy
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    indir1 = os.path.join(settings.DATA_DIR, 'chipseq', 'wtchg_p170710')
    indir2 = os.path.join(settings.DATA_DIR, 'chipseq', 'wtchg_p170710_pilot')

    outdir = output.unique_output_dir("p170710_chipseq_analysis")

    tss_subdir = os.path.join('human', 'bt2_alignment', 'gene_tss_5000_Q10')

    meta_fn = os.path.join(indir1, 'sources.csv')
    meta = pd.read_csv(meta_fn, index_col=0, header=0)
    meta.insert(meta.shape[1], 'batch', 'P170710')
    meta_fn = os.path.join(indir2, 'sources.csv')
    meta = pd.concat((meta, pd.read_csv(meta_fn, index_col=0, header=0)), axis=0)
    meta.loc[meta.batch.isnull(), 'batch'] = 'pilot'
    idx = (meta.batch == 'pilot') & (meta.input_sample.str.contains('GBM054'))
    meta.loc[idx, 'input_sample'] = [
        "%s pilot" % t for t in meta.loc[idx, 'input_sample']
    ]

    traces = {}
    for batch, indir in zip(['pilot', 'P170710'], [indir2, indir1]):
        this_meta = meta.loc[meta.batch == batch]
        for pid in this_meta.index:
            the_name = this_meta.loc[pid, 'sample']
            the_fn = os.path.join(indir, tss_subdir, str(pid), "%s.trace" % str(pid))
            with open(the_fn, 'rb') as f:
                the_vals = np.array(f.readline().split(', ')).astype(float)
            traces[the_name] = the_vals

    traces = pd.DataFrame.from_dict(traces)
    pad = (traces.shape[0] - 1) / 2  # could check this is as expected
    coords = np.arange(-pad, pad + 1)
    traces.index = coords

    traces = traces.loc[:, meta.loc[:, 'sample']]

    # 1) One plot per input sample

    chip_target_styles = {
        'H3K36me3': {'color': '#1b9e77', 'lw': 2.},
        'H3K4me3': {'color': '#d95f02', 'lw': 2.},
        'H3K27ac': {'color': '#7570b3', 'lw': 2.},
        'H3K27me3': {'color': '#e7298a', 'lw': 2.},
        'BMI1': {'color': '#66a61e', 'lw': 2.},
        'CHD7': {'color': '#e6ab02', 'lw': 2.},
        'Input': {'color': 'gray', 'lw': 1., 'zorder': 999},
    }

    chip_sample_styles = {
        'GBM017': {'color': '#a6cee3'},
        'GBM061': {'color': '#1f78b4'},
        'GBM054': {'color': '#b2df8a'},
        'iNSC017': {'color': '#33a02c'},
        'iNSC054': {'color': '#fb9a99'},
        'iNSC061': {'color': '#e31a1c'},
        '3021Scr': {'color': '#fdbf6f'},
        '3021shBMI1': {'color': '#ff7f00'},
        '3021shCHD7': {'color': '#cab2d6'},
        '3021shBMI1CHD7': {'color': '#6a3d9a'},
        'GBM050': {'color': '#ffff99'},
        'GBM054 pilot': {'color': '#b15928'},
    }

    # keep a copy of the enrichment
    # we can compute this in a number of ways
    ne = {}
    ne_rl = {}
    ne_sub = {}

    ne_by_sample = {}
    ne_rl_by_sample = {}
    ne_sub_by_sample = {}

    for the_name in sorted(meta.input_sample.unique()):
        input_idx = (meta.input_sample == the_name) & (meta.loc[:, 'ChIP target'] == 'Input')
        if input_idx.sum() != 1:
            # raise ValueError("Unable to find input data unambiguously for input sample %s" % the_name)
            print "Unable to find input data unambiguously for input sample %s" % the_name
            continue
        else:
            input_idx = np.where(input_idx)[0][0]

        dat_input = traces.loc[:, meta.iloc[input_idx].loc['sample']]
        norm_input = meta.flagstat_reads_mapped_q10.iloc[input_idx]
        mean_input = dat_input.mean()

        chip_idx = (meta.input_sample == the_name) & (meta.loc[:, 'ChIP target'] != 'Input')
        dat_chip = traces.loc[:, meta.loc[chip_idx, 'sample']]
        norm_chip = meta.flagstat_reads_mapped_q10[chip_idx]
        norm_chip.index = meta.loc[chip_idx, 'sample']

        chip_conditions = meta.loc[chip_idx, 'ChIP target']

        norm_factors = norm_input / norm_chip.astype(float)

        enrichment = dat_chip.divide(dat_input, axis=0).multiply(norm_factors, axis=1)
        enrichment.columns = chip_conditions
        for c in enrichment:
            ne.setdefault(c, {})[the_name] = enrichment.loc[:, c]
            ne_by_sample.setdefault(the_name, {})[c] = ne[c][the_name]
        ne.setdefault('Input', {})[the_name] = pd.Series(1., index=enrichment.index)
        ne_by_sample[the_name]['Input'] = ne['Input'][the_name]

        enrichment_rl = dat_chip.multiply(norm_factors, axis=1) / mean_input
        enrichment_rl.columns = chip_conditions
        for c in enrichment_rl:
            ne_rl.setdefault(c, {})[the_name] = enrichment_rl.loc[:, c]
            ne_rl_by_sample.setdefault(the_name, {})[c] = ne_rl[c][the_name]
        ne_rl.setdefault('Input', {})[the_name] = dat_input / mean_input
        ne_rl_by_sample[the_name]['Input'] = ne_rl['Input'][the_name]

        enrichment_sub = dat_chip.multiply(norm_factors, axis=1).subtract(dat_input, axis=0)
        enrichment_sub.columns = chip_conditions
        for c in enrichment_sub:
            ne_sub.setdefault(c, {})[the_name] = enrichment_sub.loc[:, c]
            ne_sub_by_sample.setdefault(the_name, {})[c] = ne_sub[c][the_name]
        ne_sub.setdefault('Input', {})[the_name] = pd.Series(0., index=enrichment.index)
        ne_sub_by_sample[the_name]['Input'] = ne_sub['Input'][the_name]


    # 1) One plot per input sample

    the_ne_measures = {
        'rel': ne_by_sample,
        'rel_mean': ne_rl_by_sample,
        'sub': ne_sub_by_sample,
    }

    for the_method, this_ne in the_ne_measures.items():
        for the_name in this_ne:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for the_target, the_dat in sorted(this_ne[the_name].items(), key = lambda x: x[0]):
                ax.plot(the_dat.index, the_dat.values, label=the_target, **chip_target_styles[the_target])
            ax.legend(loc='upper right')
            ax.set_title(the_name)
            ax.set_xlabel('Position, relative to TSS')
            ax.set_ylabel('Normalised enrichment')
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "%s_%s.png" % (the_name, the_method)), dpi=200)

        plt.close('all')

    # 2) One plot per ChIP target
    the_ne_measures = {
        'rel': ne,
        'rel_mean': ne_rl,
        'sub': ne_sub,
    }

    for the_method, this_ne in the_ne_measures.items():
        for the_target in sorted(this_ne.keys()):
            if the_target == 'Input':
                continue
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for the_sample, the_dat in sorted(this_ne[the_target].items(), key=lambda x: x[0]):
                ax.plot(the_dat.index, the_dat.values, label=the_sample, **chip_sample_styles[the_sample])
                if the_method == 'rel_mean':
                    the_input = this_ne['Input'][the_sample]
                    the_input_style = copy.copy(chip_sample_styles[the_sample])
                    the_input_style['lw'] = 1.
                    the_input_style['ls'] = '--'
                    ax.plot(the_input.index, the_input.values, **the_input_style)
            if the_method == 'rel':
                ax.axhline(1.0, label='Input', **chip_target_styles['Input'])
            elif the_method == 'sub':
                ax.axhline(0.0, label='Input', **chip_target_styles['Input'])
            ax.legend(loc='upper right')
            ax.set_title(the_target)
            ax.set_xlabel('Position, relative to TSS')
            ax.set_ylabel('Normalised enrichment')
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "%s_%s.png" % (the_target, the_method)), dpi=200)
        plt.close('all')

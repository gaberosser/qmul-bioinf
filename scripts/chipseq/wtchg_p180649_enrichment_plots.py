import os
import pickle
from utils import output
import settings
import numpy as np
import csv
import copy
import itertools
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from plotting import common


if __name__ == '__main__':
    indir = os.path.join(settings.DATA_DIR_NON_GIT, 'chipseq', 'wtchg_p180649')

    outdir = output.unique_output_dir()

    tss_subdir = os.path.join('human', 'bt2_alignment', 'tss_enrichment_5000_Q10')

    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, index_col=0, header=0)

    traces = {}
    for pid in meta.index:
        the_name = meta.loc[pid, 'sample']
        the_fn = os.path.join(indir, tss_subdir, "%s.trace" % str(pid))
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

    cmap = plt.cm.Paired
    hgic_sample_order = [
        "%s%s" % t[::-1] for t in list(itertools.product(['018', '019', '026', '030', '031', '052'], ['GBM', 'iNSC']))
    ]

    chip_sample_colours= dict(zip(hgic_sample_order, cmap.colors))
    chip_sample_colours['1299'] = '#000000'
    chip_sample_colours['1299_shBMI1'] = '#4d4d4d'
    chip_sample_colours['1299_shCHD7'] = '#b3b3b3'
    chip_sample_colours['1299_shBMI1shCHD7'] = '#e6e6e6'

    # keep a copy of the enrichment
    # we can compute this in a number of ways
    ne = {}
    ne_rl = {}
    ne_sub = {}

    ne_by_sample = {}
    ne_rl_by_sample = {}
    ne_sub_by_sample = {}

    for the_name in sorted(meta.input_sample.unique()):
        input_idx = (meta.input_sample == the_name) & (meta.loc[:, 'chip_target'] == 'Input')
        if input_idx.sum() != 1:
            print "Unable to find input data unambiguously for input sample %s" % the_name
            continue
        else:
            input_idx = np.where(input_idx)[0][0]

        dat_input = traces.loc[:, meta.iloc[input_idx].loc['sample']]
        norm_input = meta.flagstat_reads_mapped_q10.iloc[input_idx]
        mean_input = dat_input.mean()

        chip_idx = (meta.input_sample == the_name) & (meta.loc[:, 'chip_target'] != 'Input')
        dat_chip = traces.loc[:, meta.loc[chip_idx, 'sample']]
        norm_chip = meta.flagstat_reads_mapped_q10[chip_idx]
        norm_chip.index = meta.loc[chip_idx, 'sample']

        chip_conditions = meta.loc[chip_idx, 'chip_target']

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
            fig = plt.figure(figsize=(10., 5.5))
            ax = fig.add_subplot(111)
            for the_sample, the_dat in sorted(this_ne[the_target].items(), key=lambda x: x[0]):
                ax.plot(the_dat.index, the_dat.values, label=the_sample, color=chip_sample_colours[the_sample])
                if the_method == 'rel_mean':
                    the_input = this_ne['Input'][the_sample]
                    the_input_style = {'color': chip_sample_colours[the_sample]}
                    the_input_style['lw'] = 1.
                    the_input_style['ls'] = '--'
                    ax.plot(the_input.index, the_input.values, **the_input_style)
            if the_method == 'rel':
                ax.axhline(1.0, label='Input', **chip_target_styles['Input'])
            elif the_method == 'sub':
                ax.axhline(0.0, label='Input', **chip_target_styles['Input'])
            # ax.legend(loc='upper right')
            common.legend_outside_axes(ax)
            ax.set_title(the_target)
            ax.set_xlabel('Position, relative to TSS')
            ax.set_ylabel('Normalised enrichment')
            fig.tight_layout()
            fig.subplots_adjust(right=0.8)
            fig.savefig(os.path.join(outdir, "%s_%s.png" % (the_target, the_method)), dpi=200)
        plt.close('all')

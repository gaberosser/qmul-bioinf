from methylation import dmr, loader, process
import pandas as pd
from utils import setops, excel, output
from stats import basic
from plotting import common
from sklearn.decomposition import pca
import collections
from matplotlib import pyplot as plt
import numpy as np
import os
from settings import GIT_LFS_DATA_DIR, OUTPUT_DIR
import multiprocessing as mp


if __name__ == "__main__":
    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 0.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu_permute',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    me_data_indir = os.path.join(OUTPUT_DIR, 'mb_methylation_data')
    de_results_indir = os.path.join(GIT_LFS_DATA_DIR, 'mb_de_bmi1_chd7')

    outdir = output.unique_output_dir("mb_de_dmr")
    norm_method = 'swan'

    obj = loader.IlluminaHumanMethylationLoader(
        base_dir=me_data_indir,
        meta_fn=os.path.join(me_data_indir, 'sources.csv'),
        norm_method=norm_method,
    )


    # obj = loader.load_by_patient(['3021', 'ICb1299'], norm_method=norm_method, include_control=False)

    # add condition and cell line column to meta
    meta = obj.meta
    # condition = pd.Series({
    #     '3021_1_Scr': 'scramble',
    #     '3021_1_shB': 'shBMI1',
    #     '3021_1_shC': 'shCHD7',
    #     '3021_1_shB+C': 'shBMI1shCHD7',
    #     'S': 'scramble',
    #     'B': 'shBMI1',
    #     'C': 'shCHD7',
    #     'B+C': 'shBMI1shCHD7',
    #     'ICb1299_Scr': 'scramble',
    #     'ICb1299_shBMI1': 'shBMI1',
    #     'ICb1299_shCHD7': 'shCHD7',
    #     'ICb1299_shBMI1CHD7': 'shBMI1shCHD7',
    #     'p62_3_shBmi1': 'shBMI1',
    #     'p62_3_shChd7': 'shCHD7',
    #     'p62_3_shB+C': 'shBMI1shCHD7',
    #     'p62_3_Scr': 'scramble',
    #
    # })
    # condition = condition.loc[meta.index]
    # meta.insert(0, 'condition', condition)
    #
    # cell_line = pd.Series('3021', index=meta.index)
    # cell_line[cell_line.index.str.contains('1299')] = 'ICb1299'
    # cell_line[cell_line.index.str.contains('p62')] = 'ICb1299'
    # meta.insert(0, 'cell_line', cell_line)

    anno = loader.load_illumina_methylationepic_annotation()

    me_data = obj.data.dropna()
    me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    dmr.add_merged_probe_classes(anno)

    # plot PCA

    p = pca.PCA()
    pca_dat = p.fit_transform(me_data.transpose())

    fig = plt.figure()
    ax = fig.add_subplot(111)

    marker_groups = meta.cell_line
    marker_dict = dict([
        ('3021', 'o'),
        ('ICb1299', 's'),
    ])

    colour_groups = meta.batch.copy()
    colour_dict = dict([
        ('2017-09-19', '#7fc97f'),
        ('2018-01-12', '#beaed4'),
        ('2018-03-19', '#fdc086'),
        ('2018-03-26', '#ffff99'),
        ('2018-04-09', '#386cb0'),
    ])

    for v1 in colour_dict:
        for v2 in marker_dict:
            this_ix = (colour_groups == v1) & (marker_groups == v2)
            ax.scatter(
                pca_dat[this_ix, 0],
                pca_dat[this_ix, 1],
                c=colour_dict[v1],
                marker=marker_dict[v2],
                edgecolor='k',
                s=40,
                label="%s / %s" % (v1, v2)
            )

    common.legend_outside_axes(ax, frameon=True, edgecolor='k', facecolor='w', framealpha=0.5)
    ax.set_xlabel("Principle component 1 (%.0f %%)" % (p.explained_variance_ratio_[0] * 100))
    ax.set_ylabel("Principle component 2 (%.0f %%)" % (p.explained_variance_ratio_[1] * 100))
    fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.75)
    fig.savefig(os.path.join(outdir, "methylation_pca.png"), dpi=200)

    # ECDF plot
    m_vals = np.linspace(-8, 8, 200)
    ecdf = basic.ecdf(me_data, m_vals)
    ecdf = pd.DataFrame(ecdf, index=m_vals, columns=me_data.columns)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for v1 in colour_dict:
        this_ix = (colour_groups == v1)
        ax.plot(m_vals, ecdf.loc[:, this_ix], c = colour_dict[v1], label=v1)
    hs, lbls = ax.get_legend_handles_labels()
    lbls, idx = np.unique(lbls, return_index=True)
    hs = [hs[i] for i in idx]
    ax.legend(hs, lbls, loc='lower right', frameon=True, edgecolor='k', facecolor='w', framealpha=0.5)

    ax.set_ylim([-0.01, 1.01])
    ax.set_xlabel("M value")
    ax.set_ylabel("ECDF")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "methylation_ecdf.png"), dpi=200)

    ## TODO: linear interp along y axis to identify regions that could be found DMR due to norming differences?
    ## TODO: apply this to the hGIC project

    # 1) DMR: All shBMI1 vs all scramble, etc... (aggregating cell lines)

    dmr_res_obj = dmr.DmrResults(anno=anno)
    dmr_res_obj.identify_clusters(**dmr_params)

    comparisons = [
        collections.OrderedDict([
            ('shCHD7', ['3021_1_shC', 'C', 'ICb1299_shCHD7', 'p62_3_shChd7']),
            ('scramble', ['3021_1_Scr', 'S', 'ICb1299_Scr', 'p62_3_Scr'])
        ]),
        collections.OrderedDict([
            ('shBMI1', ['3021_1_shB', 'B', 'ICb1299_shBMI1', 'p62_3_shBmi1']),
            ('scramble', ['3021_1_Scr', 'S', 'ICb1299_Scr', 'p62_3_Scr'])
        ]),
        collections.OrderedDict([
            ('shCHD7shBMI1', ['3021_1_shB+C', 'B+C', 'ICb1299_shBMI1CHD7', 'p62_3_shB+C']),
            ('scramble', ['3021_1_Scr', 'S', 'ICb1299_Scr', 'p62_3_Scr'])
        ]),
    ]

    dmr_aggregate = {}

    for c in comparisons:
        ttl = '-'.join(c.keys())
        the_obj = dmr_res_obj.copy()
        the_samples = c.values()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_aggregate[ttl] = the_obj
        print "%s: %d relevant DMRs, of which %d significant." % (
            ttl,
            len(the_obj.results_relevant),
            len(the_obj.results_significant)
        )

    # 2) DMR: shBMI1 vs scramble, etc for each cell line

    comparisons = [
        collections.OrderedDict([
            ('shCHD7_1299', ['ICb1299_shCHD7', 'p62_3_shChd7']),
            ('scramble_1299', ['ICb1299_Scr', 'p62_3_Scr'])
        ]),
        collections.OrderedDict([
            ('shCHD7_3021', ['3021_1_shC', 'C']),
            ('scramble_3021', ['3021_1_Scr', 'S'])
        ]),
        collections.OrderedDict([
            ('shBMI1_1299', ['ICb1299_shBMI1', 'p62_3_shBmi1']),
            ('scramble_1299', ['ICb1299_Scr', 'p62_3_Scr'])
        ]),
        collections.OrderedDict([
            ('shBMI1_3021', ['3021_1_shB', 'B']),
            ('scramble_3021', ['3021_1_Scr', 'S'])
        ]),
        collections.OrderedDict([
            ('shCHD7shBMI1_1299', ['ICb1299_shBMI1CHD7', 'p62_3_shB+C']),
            ('scramble_1299', ['ICb1299_Scr', 'p62_3_Scr'])
        ]),
        collections.OrderedDict([
            ('shCHD7shBMI1_3021', ['3021_1_shB+C', 'B+C']),
            ('scramble_3021', ['3021_1_Scr', 'S'])
        ]),
    ]

    dmr_separate = {}

    for c in comparisons:
        ttl = '-'.join(c.keys())
        the_obj = dmr_res_obj.copy()
        the_samples = c.values()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_separate[ttl] = the_obj
        print "%s: %d relevant DMRs, of which %d significant." % (
            ttl,
            len(the_obj.results_relevant),
            len(the_obj.results_significant)
        )

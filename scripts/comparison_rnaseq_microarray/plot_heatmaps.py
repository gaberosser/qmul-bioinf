from scripts.comparison_rnaseq_microarray import load_illumina_data, load_rnaseq_data, load_references, comparisons, consts
from plotting import heatmap
from microarray.process import aggregate_by_probe_set
import references
import pandas as pd
import numpy as np
from scipy import stats
import collections
from matplotlib import rc, pyplot as plt, gridspec as gridspec
import seaborn as sns
plt.interactive(True)
sns.set_style('white')


# standardised vmin
ZMAX = 3.
# addition for logging
eps = 1e-12
# aggregation method for microarray
AGGR_METHOD = 'median'

# load RNA-Seq healthy and MB
he_ct, he_tpm, he_meta = load_references.load_cerebellum_rnaseq_reference_data()
mb_tpm = load_rnaseq_data.load_rnaseq_cufflinks_gene_count_data(unit='tpm')

# change one gene name for ease of comparison
new_ix = np.array(he_tpm.index)
new_ix[new_ix == 'EYS'] = 'EGFL11'
he_tpm.index = new_ix
he_ct.index = new_ix

new_ix = np.array(mb_tpm.index)
new_ix[new_ix == 'EYS'] = 'EGFL11'
mb_tpm.index = new_ix

# also try log version
mb_tpm_log = np.log10(mb_tpm + eps)
he_tpm_log = np.log10(he_tpm + eps)

# concatenate
all_tpm = pd.concat((mb_tpm, he_tpm), axis=1, join='inner')
all_tpm_log = np.log10(all_tpm + eps)

# standardize using mean from pool of ALL data
# all_tpm_n = all_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)
# all_tpm_nlog = all_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)

# standardize using mean from pool of HEALTHY data
all_tpm_n = all_tpm.subtract(he_tpm.mean(axis=1), axis=0).divide(he_tpm.std(axis=1), axis=0)
all_tpm_nlog = all_tpm_log.subtract(he_tpm_log.mean(axis=1), axis=0).divide(he_tpm_log.std(axis=1), axis=0)

all_northcott = []
[all_northcott.extend(v) for _, v in consts.NORTHCOTT_GENES]

g = sns.clustermap(all_tpm_nlog.loc[all_northcott, :].dropna(), vmin=-3, vmax=3, row_cluster=False, col_cluster=True)
g.cax.set_visible(False)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

if True:
    fig, axs, cax = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        all_tpm_nlog,
        vmax=ZMAX,
        orientation='vertical',
        fig_kwargs={'figsize': [5, 8.5]}
    )

    fig, axs, cax = heatmap.grouped_expression_heatmap(
        consts.NORTHCOTT_GENES,
        all_tpm_nlog,
        vmax=ZMAX,
        cbar=False,
        orientation='vertical',
        fig_kwargs={'figsize': [5, 8.5]},
        heatmap_kwargs={'square': False}
    )


marray_data, pvals = load_illumina_data.load_normed_microarray_data(pval=None, return_pvals=True)
marray_data = marray_data.subtract(marray_data.min(axis=0))

HEALTHY_SAMPLES = dict(consts.SAMPLE_GROUPS_ZHANG)['Healthy cerebellum']

probe_set = load_illumina_data.load_illumina_array_library()
marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
marray_all = aggregate_by_probe_set(marray_ann, method=AGGR_METHOD)

# take mean over repeats
for sn in load_illumina_data.SAMPLE_NAMES:
    marray_all.loc[:, sn] = marray_all.loc[:, [sn, sn + '-R']].mean(axis=1)

keep_cols = list(HEALTHY_SAMPLES) + ['Pt1299', 'ICb1299-I', 'ICb1299-III', 'ICb1299-IV']

# marray_all = marray_all.loc[:, load_illumina_data.SAMPLE_NAMES]
marray_all = marray_all.loc[:, keep_cols]
marray_all = marray_all.dropna(axis=0, how='all')

marray_all_log = np.log10(marray_all + eps)

# standardise using pool of HEALTHY data
marray_he = marray_all.loc[:, HEALTHY_SAMPLES]
marray_he_log = marray_all_log.loc[:, HEALTHY_SAMPLES]

marray_all_n = marray_all.subtract(marray_he.mean(axis=1), axis=0).divide(marray_he.std(axis=1), axis=0)
marray_all_nlog = marray_all_log.subtract(marray_he_log.mean(axis=1), axis=0).divide(marray_he_log.std(axis=1), axis=0)

if False:

    fig, axs, cax = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        marray_all_nlog,
        vmax=ZMAX,
        cbar=False,
        orientation='vertical',
        fig_kwargs={'figsize': [5, 8.5]}
    )

    fig, axs, cax = heatmap.grouped_expression_heatmap(
        consts.NORTHCOTT_GENES,
        marray_all_nlog,
        vmax=ZMAX,
        cbar=False,
        orientation='vertical',
        fig_kwargs={'figsize': [5, 8.5]},
        heatmap_kwargs={'square': False}
    )
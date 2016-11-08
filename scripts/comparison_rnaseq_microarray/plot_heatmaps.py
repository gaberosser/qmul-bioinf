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

# standardize using mean from pool of ALL data
# all_tpm = mb_tpm.copy()
all_tpm = pd.concat((mb_tpm, he_tpm), axis=1, join='inner')
all_tpm_log = np.log10(all_tpm + eps)

all_tpm_n = all_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)
mb_tpm_n = mb_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)

all_tpm_nlog = all_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)
mb_tpm_nlog = mb_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)

fig, axs, cax = heatmap.grouped_expression_heatmap(
    consts.NANOSTRING_GENES,
    all_tpm_nlog,
    vmax=ZMAX,
    fig_kwargs={'figsize': [8.5, 5]}
)

fig, axs, cax = heatmap.grouped_expression_heatmap(
    consts.NANOSTRING_GENES,
    all_tpm_nlog,
    vmax=ZMAX,
    orientation='vertical',
    fig_kwargs={'figsize': [5, 8.5]}
)
from load_data import rnaseq_data
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from plotting import clustering, corr
import re

from scripts.paired_samples.astrocytes_comparison import plot_clustermap


data, meta = rnaseq_data.nih_gdc_gbm_preprocessed(units='counts')
data = data.loc[data.index.str.contains('ENSG')]

contains_arr = [
    'MES', re.compile(r'RTK_I$'), 'RTK_II'
]
col_colours = clustering.generate_colour_map_dict(
    meta,
    'subgroup',
    contains_arr,
    label='group',
    sample_names=data.columns,
    non_matching='gray'
)

# cg = plot_clustermap(data, z_score=0, col_colors=col_colours)
fig = plt.figure()
ax = fig.add_subplot(111)
z1 = hierarchy.linkage(data.transpose(), method='average', metric='euclidean')
r1 = hierarchy.dendrogram(z1, ax=ax)

fig = plt.figure()
ax = fig.add_subplot(111)
z2 = hierarchy.linkage(data.transpose(), method='single', metric='euclidean')
r2 = hierarchy.dendrogram(z2, ax=ax)


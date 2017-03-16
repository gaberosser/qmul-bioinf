from load_data import rnaseq_data
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy import stats
import references
from plotting import clustering, corr
import re
from scripts.rnaseq import gtf_reader

from scripts.paired_samples.astrocytes_comparison import plot_clustermap

rrna_ensg = set(gtf_reader.get_rrna())
mt_ensg = set(gtf_reader.get_mitochondrial())

data, meta = rnaseq_data.nih_gdc_gbm_preprocessed(units='counts')
data = data.loc[data.index.str.contains('ENSG')]
data = data.loc[~(
    data.index.isin(rrna_ensg) | data.index.isin(mt_ensg)
)]

# normalised version with counts / sum(counts)
datan = data.divide(data.sum(axis=0), axis=1)

contains_arr = [
    'MES', re.compile(r'RTK_I$'), 'RTK_II'
]
col_colours, legend_labels = clustering.generate_colour_map_dict(
    meta,
    'subgroup',
    contains_arr,
    label='group',
    sample_names=data.columns,
    non_matching='gray',
    group_names=['Mesenchymal',  'RTK I', 'RTK II']
)

clustering.dendrogram_with_colours(data, col_colours, legend_labels=legend_labels, metric='euclidean', method='average')
clustering.dendrogram_with_colours(data, col_colours, legend_labels=legend_labels, metric='euclidean', method='single')
clustering.dendrogram_with_colours(data, col_colours, legend_labels=legend_labels, metric='correlation', method='average')
clustering.dendrogram_with_colours(data, col_colours, legend_labels=legend_labels, metric='correlation', method='single')

# concatenate meta and data for plotting
# variables in columns

aa = pd.concat((datan.transpose(), meta), axis=1)
ax = sns.boxplot(y=references.gene_symbol_to_ensembl('PDGFRA'), x='subgroup', data=aa)
ax.set_xlabel('Methylation subgroup')
ax.set_ylabel('PDGFRa (proportion of reads)')
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
plt.tight_layout()

# test for significant differences
ens_pdgfra = references.gene_symbol_to_ensembl('PDGFRA')
gg = aa.loc[:, [ens_pdgfra]].groupby(aa.subgroup)
treatments = gg.apply(lambda x: x.values.flatten())

# ANOVA tests whether there are significant differences between the different treatments
print "Comparing PDGFRA in all methylation groups..."
print stats.f_oneway(*treatments)

# was it influenced by the smaller categories?
treatments2 = [datan.loc[ens_pdgfra, idx] for grp, idx in gg.groups.items() if 'GBM' in grp or 'IDH_HG' in grp]
print "Comparing PDGFRA in RTK I, RTK II, MES, IDH mut only..."
print stats.f_oneway(*treatments2) # apparently not!

# create group labels


# cg = plot_clustermap(data, z_score=0, col_colors=col_colours)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# z1 = hierarchy.linkage(data.transpose(), method='average', metric='euclidean')
# r1 = hierarchy.dendrogram(z1, ax=ax)
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# z2 = hierarchy.linkage(data.transpose(), method='single', metric='euclidean')
# r2 = hierarchy.dendrogram(z2, ax=ax)
#

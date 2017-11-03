from load_data import rnaseq_data
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy import stats
import references
from plotting import clustering, corr
import re
import os
from scripts.rnaseq import gtf_reader
from utils.output import unique_output_dir

from scripts.paired_samples.astrocytes_comparison import plot_clustermap

outdir = unique_output_dir("nih_gdc_tcga-gbm", reuse_empty=True)

rrna_ensg = set(gtf_reader.get_rrna())
mt_ensg = set(gtf_reader.get_mitochondrial())

data, meta = rnaseq_data.tcga_primary_gbm(units='counts')
data = data.loc[data.index.str.contains('ENSG')]

# check rRNA quantity then remove

rr = data.loc[rrna_ensg].sum(axis=0) / data.sum(axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax = rr.plot.bar(width=0.9, ax=ax)
ax.set_ylim([0, 0.1])
ax.set_ylabel("Proportion rRNA")
plt.tight_layout()
fig.savefig(os.path.join(outdir, 'proportion_rrna.pdf'))
fig.savefig(os.path.join(outdir, 'proportion_rrna.png'), dpi=200)

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

ens_pdgfra = references.gene_symbol_to_ensembl('PDGFRA')
ens_nf1 = references.gene_symbol_to_ensembl('NF1')
ens_cdkn2a = references.gene_symbol_to_ensembl('CDKN2A')

aa = pd.concat((datan.transpose(), meta), axis=1)

ax = sns.boxplot(y=ens_pdgfra, x='subgroup', data=aa)
ax.set_xlabel('Methylation subgroup')
ax.set_ylabel('PDGFRa (proportion of reads)')
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
plt.tight_layout()
ax.figure.savefig(os.path.join(outdir, "pdgfra_boxplot.png"), dpi=200)
ax.figure.savefig(os.path.join(outdir, "pdgfra_boxplot.pdf"))

ax = sns.boxplot(y=ens_nf1, x='subgroup', data=aa)
ax.set_xlabel('Methylation subgroup')
ax.set_ylabel('NF1 (proportion of reads)')
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
plt.tight_layout()
ax.figure.savefig(os.path.join(outdir, "nf1_boxplot.png"), dpi=200)
ax.figure.savefig(os.path.join(outdir, "nf1_boxplot.pdf"))

ax = sns.boxplot(y=ens_cdkn2a, x='subgroup', data=aa)
ax.set_xlabel('Methylation subgroup')
ax.set_ylabel('CDKN2A (proportion of reads)')
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
plt.tight_layout()


# test for significant differences
gg_pdgfra = aa.loc[:, [ens_pdgfra]].groupby(aa.subgroup)
treatments_pdgfra = gg_pdgfra.apply(lambda x: x.values.flatten())

# ANOVA tests whether there are significant differences between the different treatments
print "Comparing PDGFRA in all methylation groups..."
print stats.f_oneway(*treatments_pdgfra)

# was it influenced by the smaller categories?
print "Comparing PDGFRA in RTK I, RTK II, MES, IDH mut only..."
print stats.f_oneway(
    *[datan.loc[ens_pdgfra, idx] for grp, idx in gg_pdgfra.groups.items() if 'GBM' in grp or 'IDH_HG' in grp]
) # apparently not!

# nonparametric test between RTK I, RTK II
print "Wilcoxon rank sums test RTK I vs RTK II"
print stats.ranksums(gg_pdgfra.get_group('GBM_RTK_I'), gg_pdgfra.get_group('GBM_RTK_II'))

# direct t test between RTK I, RTK II
print "One way T test (normal assumption doesn't hold!)"
stats.ttest_ind(gg_pdgfra.get_group('GBM_RTK_I'), gg_pdgfra.get_group('GBM_RTK_II'))

gg_nf1 = aa.loc[:, [ens_nf1]].groupby(aa.subgroup)
treatments_nf1 = gg_nf1.apply(lambda x: x.values.flatten())

print "Comparing PDGFRA in all methylation groups..."
print stats.f_oneway(*treatments_nf1)

# was it influenced by the smaller categories?
print "Comparing PDGFRA in RTK I, RTK II, MES, IDH mut only..."
print stats.f_oneway(
    *[datan.loc[ens_nf1, idx] for grp, idx in gg_nf1.groups.items() if 'GBM' in grp or 'IDH_HG' in grp]
)

# nonparametric test between RTK I, RTK II
print "Wilcoxon rank sums test MES vs RTK I"
print stats.ranksums(gg_nf1.get_group('GBM_MES'), gg_nf1.get_group('GBM_RTK_I'))
print "Wilcoxon rank sums test MES vs RTK II"
print stats.ranksums(gg_nf1.get_group('GBM_MES'), gg_nf1.get_group('GBM_RTK_II'))


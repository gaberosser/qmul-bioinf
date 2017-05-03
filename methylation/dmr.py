from load_data import methylation_array
from methylation.process import m_from_beta, merge_illumina_probe_gene_classes
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

NORM_METHOD = 'swan'
anno = methylation_array.load_illumina_methylationepic_annotation()
b = methylation_array.gbm_nsc_methylationepic('swan')
m = m_from_beta(b)
anno.loc[:, 'merged_class'] = merge_illumina_probe_gene_classes(
    anno.loc[:, 'UCSC_RefGene_Group'], anno.loc[:, 'Relation_to_UCSC_CpG_Island']
)

# plot them
chr1 = anno.loc[anno.CHR == '1'].sort_values(by='MAPINFO', axis=0)

loc_from = 1000000
loc_to = loc_from + 500000
idx = (loc_from <= chr1.MAPINFO) & (chr1.MAPINFO < loc_to)

strand_y = chr1.Strand.where(chr1.Strand == 'F', 0)
strand_y.loc[strand_y == 'F'] = 1
strand_y = strand_y.astype(int)

col_mapper = {
    'tss': [1, 0, 0],
    'gene': [0, 1, 0],
    'island': [0, 0, 1]
}
colour_no_role = [0.5, 0.5, 0.5]

fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(111)

t = chr1.loc[idx]
for grp, c in col_mapper.items():
    x = t.MAPINFO[idx & t.merged_class.str.contains(grp)].values
    y = strand_y[idx & t.merged_class.str.contains(grp)].values
    ax.scatter(x, y, c=c, marker='|', label=grp, alpha=0.7, zorder=3)

x = t.MAPINFO[idx & (t.merged_class == '')].values
y = strand_y[idx & (t.merged_class == '')].values
ax.scatter(x, y, c=colour_no_role, marker='|', alpha=0.3, label='none', zorder=2)
ax.plot([loc_from, loc_to], [0, 0], 'k-', zorder=1, alpha=0.3)
ax.plot([loc_from, loc_to], [1, 1], 'k-', zorder=1, alpha=0.3)
ax.legend(loc='center right')
ax.yaxis.set_ticks([0, 1])
ax.yaxis.set_ticklabels(['R', 'F'])
ax.set_xlabel('Chromosomal coordinate')
ax.set_title('Chromosome 1: %d - %d' % (loc_from, loc_to))
fig.tight_layout()

# identify connected regions with same class
dm = 1.4
n_min = 4
d_max = 200



import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from load_data import microarray_data
from microarray import process
from plotting import bar
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq, load_xiaonan_microarray
from scripts.mb_subgroup_classifier.shrunken_centroids import NearestCentroidClassifier
from utils import log

logger = log.get_console_logger(__name__)

from scripts.comparison_rnaseq_microarray import consts

# define shrinkage delta values
n_core = 14
deltas = np.linspace(0., 12., 60)
# FLAT_PRIOR = False
FLAT_PRIOR = True
k = 10  # XV

# it's useful to maintain a list of known upregulated genes
nano_genes = []
for grp, arr in consts.NANOSTRING_GENES:
    nano_genes.extend(arr)
nano_genes.remove('EGFL11')
nano_genes.append('EYS')

# load Ncott data (285 non-WNT MB samples)

ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
    aggr_field='SYMBOL',
    aggr_method='max'
)
sort_idx = ncott_meta.subgroup.sort_values().index
ncott_meta = ncott_meta.loc[sort_idx]
ncott = process.yugene_transform(ncott.loc[:, sort_idx])

# X = ncott.copy()
# m = ncott_meta.copy()

# load Kool dataset
kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
    aggr_field='SYMBOL',
    aggr_method='max',
)
sort_idx = kool_meta.subgroup.sort_values().index
kool_meta = kool_meta.loc[sort_idx]
kool = process.yugene_transform(kool.loc[:, sort_idx])
kool_meta.loc[:, 'subgroup'] = (
    kool_meta.loc[:, 'subgroup'].str
        .replace('A', 'WNT')
        .replace('B', 'SHH')
        .replace('E', 'Group 3')
        .replace('C', 'Group 4')
        .replace('D', 'Group 4')
)

# X = kool.copy()
# m = kool_meta.copy()

# load Robinson dataset
robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
sort_idx = robi_meta.subgroup.sort_values().index
robi_meta = robi_meta.loc[sort_idx]
robi = process.yugene_transform(robi.loc[:, sort_idx])
robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')

X = robi.copy()
m = robi_meta.copy()

# if we lump groups C and D together, the classification becomes almost perfect (as you might expect):
# m.loc[:, 'subgroup'] = m.subgroup.replace('Group 4', 'Group C/D')
# m.loc[:, 'subgroup'] = m.subgroup.replace('Group 3', 'Group C/D')

# for statistical purposes, split the data into training/validation and test sets in the ratio 3:1

# shuffle columns first
n = X.columns.size
cols = X.columns[np.random.permutation(n)]
X = X.loc[:, cols]
m = m.loc[cols]

# now remove the top 25% for testing
cut_idx = int(n * 0.25)
X_test = X.iloc[:, :cut_idx]
m_test = m.iloc[:cut_idx]
X_train = X.iloc[:, cut_idx:]
m_train = m.iloc[cut_idx:]

# leave one out cross validation
n_err = []
for i in range(n):
    # leave out one
    this = []
    j = (np.arange(n) != i)
    tr = X.loc[:, j]
    mtr = m.loc[j]
    te = X.iloc[:, i:(i+1)]
    mte = m.iloc[i:(i+1)]
    obj = NearestCentroidClassifier(tr, mtr.subgroup, flat_prior=FLAT_PRIOR)
    for d in deltas:
        obj.shrink_centroids(d)
        c, nc, ni = obj.assess_test_data(te, true_labels=mte.subgroup)
        this.append(ni)
    n_err.append(this)

n_err = np.array(n_err)
mean_err = n_err.mean(axis=0)
min_delta = deltas[np.where(mean_err == mean_err.min())[0][-1]]

# load Xiao-Nan data and try to classify
p1299_sample_names = ('Pt1299', 'ICb1299-I', 'ICb1299-III', 'ICb1299-IV')
X_p1299, p1299_meta = load_xiaonan_microarray(yugene=True, gene_symbols=X.index, sample_names=p1299_sample_names)
X_cuff = load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=X.index)

clas_p1299 = pd.DataFrame(index=p1299_sample_names, columns=deltas)
clas_cuff = pd.DataFrame(index=X_cuff.columns, columns=deltas)

obj = NearestCentroidClassifier(X_train, m_train.subgroup, flat_prior=FLAT_PRIOR)
for d in deltas:
    obj.shrink_centroids(d)

    t, u = zip(*[(c, obj.class_probabilities(X_p1299.loc[:, c])) for c in X_p1299.columns])
    v = pd.Series(data=u, index=t)
    clas_p1299.loc[:, d] = v

    t, u = zip(*[(c, obj.class_probabilities(X_cuff.loc[:, c])) for c in X_cuff.columns])
    v = pd.Series(data=u, index=t)
    clas_cuff.loc[:, d] = v

# use the results to extract the relative probability of Grp3/4 in Xiao Nan 1299 samples
grps = ('Group 3', 'Group 4', 'SHH', 'WNT')
colours = ['r', 'b', 'c', 'k']
k_sample = 'ICb1299-IV'
y = np.array(
    list(
        clas_p1299.applymap(lambda x: [x.loc[g] for g in grps]).loc[k_sample].values
    )
).transpose()
fig = plt.figure(figsize=[12, 4])
ax = fig.add_subplot(111)
y = pd.DataFrame(y, columns=deltas, index=grps)
bar.stacked_bar_chart(deltas, y, colours=colours, width=0.8, ax=ax)
dt = deltas[1] - deltas[0]
ax.set_xlim([deltas[0] - dt * 0.5, deltas[-1] + dt * 0.5])
ax.set_ylim([0, 1])
ax.legend(loc='upper right', frameon=True)
ax.set_xlabel("Shrinkage parameter, $\Delta$")
ax.set_ylabel("Subgroup membership probability")
plt.tight_layout()



fig = plt.figure()
ax = fig.add_subplot(111)
bottom = np.zeros(len(deltas))
for i, g in enumerate(grps):
    p = clas_p1299.applymap(lambda x: x.loc[g]).loc[k_sample]
    ax.bar(deltas - 0.5 * dt, p, width=0.9 * dt, color=colours[i], label=g)
    bottom += p

X_xnan, xnan_meta = load_xiaonan_microarray(yugene=True, gene_symbols=X.index)
xnan_meta = xnan_meta.loc[~xnan_meta.loc[:, 'northcott classification'].isnull()]
X_xnan = X_xnan.loc[:, xnan_meta.index]
xnan_meta.loc[:, 'subgroup'] = xnan_meta.loc[:, 'northcott classification'].replace('C', 'Group 3').replace('D', 'Group 4')

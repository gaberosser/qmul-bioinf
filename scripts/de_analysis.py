import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import seaborn as sns
import numpy as np
from scripts.comparison_rnaseq_microarray.consts import NANOSTRING_GENES

NANOSTRING_GENES = list(NANOSTRING_GENES)
NANOSTRING_GENES[2] = ('Group C', ('IMPG2', 'GABRA5', 'EYS', 'NRL', 'MAB21L2', 'NPR3'),)
all_nano = []
[all_nano.extend(arr) for _, arr in NANOSTRING_GENES]

pc = (
    (1., 0.01, [0.40392157,  0.        ,  0.05098039]),
    (0.01, 0.001, [0.94666667,  0.2682353 ,  0.19607844]),
    (0.001, 0.0001, [0.98823529,  0.67154173,  0.56053827]),
    (0.0001, 0., [ 1.        ,  0.96078432,  0.94117647])
)

# fn = 'limma_de_1299all-healthy.csv'
# res = pd.read_csv(fn, index_col=0)

fn = 'deseq2_de_1299all-healthy.csv'
res = pd.read_csv(fn, index_col=0).dropna()

hkg = ['GAPDH', 'ACTB', 'H1FX', 'ATP2C1', 'LAMP1', 'CTNNA1', 'C14orf2']

# lfc = res.loc[hkg, 'log2FoldChange']
# pval =  res.loc[hkg, 'padj']
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# for p0, p1, c in pc:
#     idx = np.where((p1 < pval) & (pval <= p0))[0]
#     if len(idx):
#         rect = ax.bar(idx, lfc[idx], 1, color=c, label='%.4f < p < %.4f' % (p1, p0))
#
# ax.legend()
# ax.set_xticks(np.arange(len(hkg)) + 0.5)
# ax.set_xticklabels(hkg, rotation=45)
# plt.xlabel('HKG')
# plt.ylabel('Log2 FC')
# plt.show()
#
# fig.savefig('rnaseq_hkg.png', dpi=200)
# fig.savefig('rnaseq_hkg.pdf', dpi=200)

lfc = res.loc[all_nano, 'log2FoldChange']
pval =  res.loc[all_nano, 'padj']

fig = plt.figure()
ax = fig.add_subplot(111)

for p0, p1, c in pc:
    idx = np.where((p1 < pval) & (pval <= p0))[0]
    if len(idx):
        rect = ax.bar(idx, lfc[idx], 1, color=c, label='%.4f < p < %.4f' % (p1, p0))

ax.legend()
ax.set_xticks(np.arange(len(all_nano)) + 0.5)
ax.set_xticklabels(all_nano, rotation=45)
plt.xlabel('nanostring gene')
plt.ylabel('Log2 FC')
plt.show()

fig.savefig('rnaseq_nano.png', dpi=200)
fig.savefig('rnaseq_nano.pdf', dpi=200)


# lfc = res.loc[hkg, 'logFC']
# pval =  res.loc[hkg, 'adj.P.Val']
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# for p0, p1, c in pc:
#     idx = np.where((p1 < pval) & (pval <= p0))[0]
#     if len(idx):
#         rect = ax.bar(idx, lfc[idx], 1, color=c, label='%.4f < p < %.4f' % (p1, p0))
#
# ax.legend()
# ax.set_xticks(np.arange(len(hkg)) + 0.5)
# ax.set_xticklabels(hkg, rotation=45)
# plt.xlabel('HKG')
# plt.ylabel('Log2 FC')
# plt.show()
#
# fig.savefig('marr_hkg.png', dpi=200)
# fig.savefig('marr_hkg.pdf', dpi=200)

# lfc = res.loc[all_nano, 'logFC']
# pval = res.loc[all_nano, 'adj.P.Val']
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# for p0, p1, c in pc:
#     idx = np.where((p1 < pval) & (pval <= p0))[0]
#     if len(idx):
#         rect = ax.bar(idx, lfc[idx], 1, color=c, label='%.4f < p < %.4f' % (p1, p0))
#
# ax.legend()
# ax.set_xticks(np.arange(len(all_nano)) + 0.5)
# ax.set_xticklabels(all_nano, rotation=45)
# plt.xlabel('nanostring gene')
# plt.ylabel('Log2 FC')
# plt.show()
#
# fig.savefig('marr_nano.png', dpi=200)
# fig.savefig('marr_nano.pdf', dpi=200)
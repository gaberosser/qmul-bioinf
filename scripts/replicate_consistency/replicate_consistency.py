from load_data import rnaseq_data
import pandas as pd
from plotting import clustering
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns



def get_twodata(data, pattern, filter_idx=None):
    if filter_idx is not None:
        twodata = data.loc[filter_idx, obj.data.columns.str.contains(pattern)]
    else:
        twodata = data.loc[:, obj.data.columns.str.contains(pattern)]
    if twodata.shape[1] != 2:
        raise ValueError("Expected two samples for matching entry %s, got %d" % (m, twodata.shape[1]))
    return twodata


if __name__ == "__main__":
    obj = rnaseq_data.all_hgic_loader(annotate_by='Ensembl ID')
    # filters:
    # 1) real genes only
    # 2) only genes that have a count > 5 in 10 or more samples
    gene_idx = obj.data.index.str.contains('ENSG')
    idx = gene_idx & ((obj.data > 5).sum(axis=1) > 10)

    matches = (
        'GBM018',
        'GBM026',
        'GBM044',
        'DURA018',
        'DURA044',
    )

    # table of correlation coefficients for matching samples
    pair_corrcoeff = pd.Series(index=matches)
    for m in matches:
        twodata = get_twodata(obj.data, m, filter_idx=idx)
        pair_corrcoeff.loc[m] = twodata.corr().iloc[0, 1]

    # plots of log(absolute difference) vs proportion difference for each pair
    for i, m in enumerate(matches):
        twodata = get_twodata(obj.data, m, filter_idx=gene_idx)
        # filter out any low counts: at least one must be > 5
        twodata = twodata.loc[(twodata > 5).sum(axis=1) > 0]
        abs_diff = np.abs(twodata.iloc[:, 0] - twodata.iloc[:, 1])
        mean_val = twodata.sum(axis=1)
        rel_diff = abs_diff / mean_val
        out = sns.jointplot(np.log2(abs_diff + 1), rel_diff, alpha=0.3, stat_func=None)
        ax = out.ax_joint
        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.set_title(m)
        plt.tight_layout()

    cg = clustering.plot_correlation_clustermap(obj.data.loc[idx])
    cg.gs.update(bottom=0.25, right=0.7)
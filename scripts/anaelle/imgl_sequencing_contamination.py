from rnaseq import loader, filter
import numpy as np
from plotting import clustering
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from scipy.cluster import hierarchy as hc
from utils import output
import os


def pairwise_correlation(many, method='pearson'):
    if method == 'pearson':
        the_func = stats.pearsonr
    elif method == 'spearman':
        the_func = stats.spearmanr
    else:
        raise NotImplementedError("Unrecognised method %s" % method)

    cor = pd.DataFrame(index=many.index, columns=many.index)
    pval = pd.DataFrame(index=many.index, columns=many.index)

    for i, row in enumerate(many.index):
        for j in range(i, len(many.index)):
            col = many.index[j]
            c, p = the_func(many.loc[row], many.loc[col])
            cor.loc[row, col] = c
            cor.loc[col, row] = c
            pval.loc[row, col] = p
            pval.loc[col, row] = p

    return cor, pval



if __name__ == "__main__":
    min_tpm = 1.
    eps = 0.01
    outdir = output.unique_output_dir()
    # load just 1st lane
    obj = loader.load_references('wtchg_p180443/180911_K00150_0372_AHWV7TBBXX', tax_id=10090, source='salmon')
    ix = (obj.meta.species == 'mouse') | (obj.meta.index.str.contains(r'[iI]MGL'))
    dat = obj.data.loc[:, ix]
    dat_filt = filter.filter_by_cpm(dat, min_cpm=min_tpm, unless_cpm_gt=10.)
    log_dat_filt = np.log2(dat_filt + eps)

    # correlation clustermap
    row_colours = pd.DataFrame('g', index=log_dat_filt.columns, columns=['Sample type'])
    row_colours.loc[row_colours.index.str.contains('mDURA')] = 'k'
    row_colours.loc[row_colours.index.str.contains('mDURA5_NSCmus_N3BE50.2')] = 'y'
    row_colours.loc[row_colours.index.str.contains('mDURA6_NSCmus')] = 'y'

    # Spearman distance

    cor, pval = pairwise_correlation(log_dat_filt.transpose(), method='spearman')
    di = 1 - cor
    for i in range(di.shape[0]):
        di.iloc[i, i] = 0.
    di = hc.distance.squareform(di)

    cg = clustering.plot_correlation_clustermap(log_dat_filt, row_colors=row_colours, distance=di, method='average')
    cg.gs.update(bottom=0.3, right=0.7)
    cg.savefig(os.path.join(outdir, "log_tpm_spearman_distance_average_linkage.png"), dpi=200)

    # Euclidean distance

    # cg2 = clustering.plot_correlation_clustermap(
    #     log_dat_filt,
    #     row_colors=row_colours,
    #     method='complete',
    #     metric='euclidean',
    # )
    # cg2.gs.update(bottom=0.3, right=0.7)

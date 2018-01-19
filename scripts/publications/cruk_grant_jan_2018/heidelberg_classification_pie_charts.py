import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from load_data import rnaseq_data
import pandas as pd
import os
from utils import output
from settings import DATA_DIR_NON_GIT


if __name__ == "__main__":
    outdir = output.unique_output_dir("cruk_classification")
    ffpe_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-06-10_brandner', 'heidelberg_classifier', '2017_10', 'NH15-2101.calibrated_scores.csv')
    cc1_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'heidelberg_classifier', '2017_10', 'GBM019_P4_DNA_8-11-2016_CLEANED.calibrated_scores.csv')
    cc2_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-09-19', 'heidelberg_classifier', '2017_10', 'GBM019_NH16-2101_P6_FROM_8-11-2015_+_GBM019_NH16-2101_P3_FROM_26-10-2015.calibrated_scores.csv')

    # number of classes to include
    n = 4
    unioner = lambda x, y: set(x).union(y)

    ffpe = pd.read_csv(ffpe_fn, header=0, index_col=0).Score
    cc1 = pd.read_csv(cc1_fn, header=0, index_col=0).Score
    cc2 = pd.read_csv(cc2_fn, header=0, index_col=0).Score

    # get all classes in top N
    all_classes = sorted(reduce(unioner, [ffpe.index[:n], cc1.index[:n], cc2.index[:n]]))

    colours = {
        'GBM, MES': '#f4e842',
        'GBM, RTK I': '#0d680f',
        'GBM, RTK II': '#820505',
        'MB, G3': 'orange',
        'PLEX, PED B': '#4C72B0',
        'Other': 'gray'
    }

    top_n = pd.DataFrame(index=all_classes + ['Other'], columns=['FFPE', 'CC1', 'CC2'])
    top_n.loc[all_classes, 'FFPE'] = ffpe.loc[all_classes]
    top_n.loc[all_classes, 'CC1'] = cc1.loc[all_classes]
    top_n.loc[all_classes, 'CC2'] = cc2.loc[all_classes]
    top_n.loc['Other'] = 1 - top_n.sum()

    clist = [colours.get(k, 'gray') for k in all_classes + ['Other']]

    gs = GridSpec(1, 4, width_ratios=[2, 2, 2, 1])
    fig = plt.figure(figsize=(9, 2.7))
    axs = [fig.add_subplot(gs[i]) for i in range(4)]
    # fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
    axs[0].pie(top_n.loc[:, 'FFPE'].values, colors=clist)
    axs[1].pie(top_n.loc[:, 'CC1'].values, colors=clist)
    patches, texts = axs[2].pie(top_n.loc[:, 'CC2'].values, colors=clist)

    axs[0].set_aspect('equal')
    axs[0].text(0, 0.4, 'Tumour bulk (FFPE)', horizontalalignment='center')

    axs[1].set_aspect('equal')
    axs[1].text(0, 0.4, 'GIC culture 1', horizontalalignment='center')

    axs[2].set_aspect('equal')
    axs[2].text(0, 0.4, 'GIC culture 2', horizontalalignment='center')

    axs[3].axis('off')
    axs[3].legend(patches, top_n.index, loc='center right')

    gs.update(left=0, right=1, top=1.1, bottom=-.1, wspace=-.1)

    fig.savefig(os.path.join(outdir, "classification_pie_charts.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "classification_pie_charts.tiff"), dpi=200)
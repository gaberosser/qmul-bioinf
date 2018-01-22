import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from load_data import rnaseq_data
import pandas as pd
import os
from utils import output
from settings import DATA_DIR_NON_GIT


def get_top_n(data, n):
    all_classes = set()
    for c in data.columns:
        x = data.loc[:, c].sort_values(ascending=False)
        all_classes = all_classes.union(x.index[:n])
    top_n = pd.DataFrame(index=sorted(all_classes) + ['Other'], columns=data.columns)
    for c in data.columns:
        top_n.loc[all_classes, c] = data.loc[all_classes, c]
    top_n.loc['Other'] = 1 - top_n.sum()
    return top_n



if __name__ == "__main__":
    outdir = output.unique_output_dir("cruk_classification")
    ffpe_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-06-10_brandner', 'heidelberg_classifier', '2017_10', 'NH15-2101.calibrated_scores.csv')
    cc1_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'heidelberg_classifier', '2017_10', 'GBM019_P4_DNA_8-11-2016_CLEANED.calibrated_scores.csv')
    cc2_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-09-19', 'heidelberg_classifier', '2017_10', 'GBM019_NH16-2101_P6_FROM_8-11-2015_+_GBM019_NH16-2101_P3_FROM_26-10-2015.calibrated_scores.csv')

    # number of classes to include
    n = 3

    ffpe = pd.read_csv(ffpe_fn, header=0, index_col=0).Score
    ffpe.name = 'FFPE'
    cc1 = pd.read_csv(cc1_fn, header=0, index_col=0).Score
    cc1.name = 'CC1'
    cc2 = pd.read_csv(cc2_fn, header=0, index_col=0).Score
    cc2.name = 'CC2'
    data = pd.concat((ffpe, cc1, cc2), axis=1)


    # get all classes in top N
    colours = {
        'GBM, MES': '#f4e842',
        'GBM, RTK I': '#0d680f',
        'GBM, RTK II': '#820505',
        'MB, G3': 'orange',
        'PLEX, PED B': '#4C72B0',
        'Other': 'gray'
    }

    # three pie charts
    top_n = get_top_n(data, n)

    clist = [colours.get(k, 'gray') for k in top_n.index]

    gs = GridSpec(1, 4, width_ratios=[2, 2, 2, 1])
    fig = plt.figure(figsize=(9, 2.7))
    axs = [fig.add_subplot(gs[i]) for i in range(4)]

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

    fig.savefig(os.path.join(outdir, "classification_pie_charts_3.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "classification_pie_charts_3.tiff"), dpi=200)

    # three pie charts in a 2 x 2

    gs = GridSpec(2, 2)
    fig = plt.figure(figsize=(4.24, 4.11))
    top_axs = [fig.add_subplot(gs[0, i]) for i in range(2)]
    bottom_axs = [fig.add_subplot(gs[1, i]) for i in range(2)]
    top_axs[0].set_aspect('equal')
    [ax.set_aspect('equal') for ax in bottom_axs]

    top_axs[0].pie(top_n.loc[:, 'FFPE'].values, colors=clist)
    bottom_axs[0].pie(top_n.loc[:, 'CC1'].values, colors=clist)
    patches, texts = bottom_axs[1].pie(top_n.loc[:, 'CC2'].values, colors=clist)

    # axs[0].set_aspect('equal')
    top_axs[0].text(0, 0.4, 'Tumour bulk (FFPE)', horizontalalignment='center')

    # axs[1].set_aspect('equal')
    bottom_axs[0].text(0, 0.4, 'GIC culture 1', horizontalalignment='center')

    # axs[2].set_aspect('equal')
    bottom_axs[1].text(0, 0.4, 'GIC culture 2', horizontalalignment='center')

    top_axs[1].axis('off')
    top_axs[1].legend(patches, top_n.index, loc='center', fontsize=12)

    gs.update(left=0, right=1, top=1.1, bottom=-.1, wspace=-.15, hspace=-.35)

    fig.savefig(os.path.join(outdir, "classification_pie_charts_3_2x2.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "classification_pie_charts_3_2x2.tiff"), dpi=200)

    # two pie charts
    n = 3
    top_n = get_top_n(data.iloc[:, :2], n)
    clist = [colours.get(k, 'gray') for k in top_n.index]

    gs = GridSpec(2, 2, height_ratios=[3, 1])
    fig = plt.figure(figsize=(5.8, 3.4))
    pie_axs = [fig.add_subplot(gs[0, i]) for i in range(2)]
    leg_ax = fig.add_subplot(gs[1, :])
    # axs = [fig.add_subplot(gs[i]) for i in range(3)]
    pie_axs[0].pie(top_n.loc[:, 'FFPE'].values, colors=clist)
    patches, texts = pie_axs[1].pie(top_n.loc[:, 'CC1'].values, colors=clist)

    pie_axs[0].set_aspect('equal')
    pie_axs[0].text(0, 0.4, 'Tumour bulk (FFPE)', horizontalalignment='center')

    pie_axs[1].set_aspect('equal')
    pie_axs[1].text(0, 0.4, 'GIC culture 1', horizontalalignment='center')

    leg_ax.axis('off')
    leg_ax.legend(patches, top_n.index, loc='upper center')

    gs.update(left=0, right=1, top=1., bottom=0, wspace=-.2, hspace=-.3)

    fig.savefig(os.path.join(outdir, "classification_pie_charts_2.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "classification_pie_charts_2.tiff"), dpi=200)


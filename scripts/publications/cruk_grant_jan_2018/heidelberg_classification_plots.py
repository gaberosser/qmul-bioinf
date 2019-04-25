from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import ternary
import numpy as np
from load_data import rnaseq_data
import pandas as pd
import os
import re
import shutil
import zipfile
from utils import output
from plotting import common
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


def read_calibrated_scores_script():
    """
    Hold space for a copy-paste script used to read calibrated scores from their raw CSV files
    :return:
    """
    batches = [
        '2016-06-10_brandner',
        '2016-12-19_ucl_genomics',
        '2016-09-21_dutt',
        '2017-01-17_brandner',
        '2017-02-09_brandner',
        '2017-05-12',
        '2017-08-23',
        '2017-09-19',
        '2018-03-19'
    ]
    basedir = os.path.join(DATA_DIR_NON_GIT, 'methylation')
    for b in batches:
        print "*** %s ***" % b
        if os.path.isdir(os.path.join(basedir, b, 'heidelberg_classifier')):
            the_dir = os.path.join(basedir, b, 'heidelberg_classifier')
        else:
            the_dir = os.path.join(basedir, b, 'heidelberg')
        c_dir = sorted([os.path.join(the_dir, t) for t in os.listdir(the_dir) if os.path.isdir(os.path.join(the_dir, t))])[-1]
        flist = sorted([t for t in os.listdir(c_dir) if re.search(r'calibrated_', t)])
        if len(flist) == 0:
            print "WARNING: Batch %s has no matching files" % b
        for fn in flist:
            ff = os.path.join(c_dir, fn)
            df = pd.read_csv(ff, header=0, index_col=0)
            the_rows = sorted(df.index[df.index.str.contains(r'GBM') & df.index.str.contains(r'(RTK[ _]I$|RTK[ _]II$|MES)')])
            the_arr = df.iloc[:, -1].loc[the_rows].values
            print "%s, %s, %.3f, %.3f, %.3f" % ((fn, b) + tuple(the_arr))


def unzip_scores(the_dir):
    """
    Holding space for a copy-paste script used to extract scores from their zip file
    :return:
    """
    c_dir = sorted([os.path.join(the_dir, t) for t in os.listdir(the_dir) if os.path.isdir(os.path.join(the_dir, t))])[-1]
    flist = [t for t in os.listdir(c_dir) if re.search(r'\.zip', t)]
    for t in flist:
        z = zipfile.ZipFile(os.path.join(c_dir, t))
        fname = [u.filename for u in z.filelist if re.search(r'scores_cal.csv', u.filename)][0]
        z.extract(fname, path=c_dir)
        shutil.move(os.path.join(c_dir, fname), os.path.join(c_dir, t.replace('.zip', '.calibrated_scores.csv')))


if __name__ == "__main__":

    outdir = output.unique_output_dir("cruk_classification")
    in_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'classification.xlsx')
    ffpe_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-06-10_brandner', 'heidelberg_classifier', '2017_10', 'NH15-2101.calibrated_scores.csv')
    cc1_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'heidelberg_classifier', '2017_10', 'GBM019_P4_DNA_8-11-2016_CLEANED.calibrated_scores.csv')
    cc2_fn = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-09-19', 'heidelberg_classifier', '2017_10', 'GBM019_NH16-2101_P6_FROM_8-11-2015_+_GBM019_NH16-2101_P3_FROM_26-10-2015.calibrated_scores.csv')

    # number of classes to include
    n = 3

    # load summary sheet
    df = pd.read_excel(in_fn, sheet_name='GBM specific', header=0, index_col=None)
    # fix format of PIDs
    df.insert(0, 'pid', df.loc[:, 'Patient ID'].apply(lambda x: "%03.0f" % x))


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

    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']

    # pick out relevant results
    cc_ix = df.pid.isin(pids) & (df.Tissue == 'Primary culture')
    ff_ix = df.pid.isin(pids) & (df.Tissue == 'FFPE tumour')

    cc = df.loc[cc_ix].set_index('pid').sort_index()
    ff = df.loc[ff_ix].set_index('pid').sort_index()
    cc.insert(1, 'mean_passage', cc.Passage.astype(str).apply(lambda x: np.mean([float(t) for t in x.split(';')])))

    # let's make a TERNARY plot with the three components: RTK I, RTK II, MES
    fig, axs = plt.subplots(ncols=3, figsize=(10.6, 3.8))
    taxs = []

    for ax in axs:
        fig, tax = ternary.figure(scale=1., ax=ax)
        taxs.append(tax)
        tax.boundary(linewidth=1.0)
        tax.gridlines(multiple=0.2, color="gray", linewidth=1.0, alpha=0.2)
        tax.clear_matplotlib_ticks()

        ax = tax.get_axes()

        ax.text(-.02, -.02, 'Mesenchymal', horizontalalignment='left', verticalalignment='top')
        ax.text(1.02, -.02, 'RTK I', horizontalalignment='right', verticalalignment='top')
        ax.text(0.5, 0.5 * np.sqrt(3) + 0.02, 'RTK II', horizontalalignment='center', verticalalignment='bottom')

        ax.set_aspect('equal')
        ax.axis('off')

    # bundle them up
    tax_dict = {
        'RTK I': taxs[0],
        'RTK II': taxs[1],
        'MES': taxs[2]
    }

    ## TODO: the ternary plot has been updated and improved at scripts/hgic_final/heidelberg_classification_plots.py

    # each patient is a 'trajectory': FFPE -> early pass -> later passage
    bases = ['GBM_RTK_I', 'GBM_RTK_II', 'GBM_MES']
    # cmap_func = common.continuous_cmap(cc.mean_passage.max(), cmap='Blues')
    cmap = common.get_best_cmap(len(pids))
    ff_colour = '#ffffff'

    for p in pids:
        this_ff_score = ff.loc[p, bases]
        this_cc = cc.loc[p].sort_values(by='mean_passage')
        this_cc_score = this_cc.loc[:, bases]
        this_cc_pass = this_cc.loc[:, 'mean_passage']
        # this_colours = [ff_colour] + [cmap_func(x - 1) for x in this_cc_pass]
        this_colours = cmap[pids.index(p)]
        this_sizes = 20 + this_cc_pass.values ** 2.

        points = np.concatenate([[this_ff_score.values], this_cc_score.values], axis=0)

        tax = tax_dict[ff.loc[p, 'Result']]
        # here's how to annotate:
        # tax.annotate('GBM%s' % p, points[0], horizontalalignment='left')
        tax.scatter(points[:1], marker='^', c=this_colours, s=30, edgecolor='k', zorder=10, alpha=0.8, label=None)
        tax.scatter(points[1:], marker='o', c=this_colours, s=this_sizes, edgecolor='k', zorder=10, alpha=0.8, label=p)
        for i in range(2):
            tax.line(points[i], points[i+1], c='gray', lw=1., zorder=9)

    for tax in tax_dict.values():
        tax.legend()

    fig.subplots_adjust(left=0.0, right=1., wspace=0.)
    fig.savefig(os.path.join(outdir, 'ternary_plot_classification.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'ternary_plot_classification.tiff'), dpi=200)

    # TODO: do we need seaborn for the next few plots to improve their appearance?
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


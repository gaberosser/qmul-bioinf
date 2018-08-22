import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
import sys
import numpy as np
from plotting import common
from utils import output


if __name__ == '__main__':
    outdir = output.unique_output_dir("ipa_pathways_hgic_de")
    ## TODO: put this somewhere more portable!
    home = os.path.expanduser("~")

    # pathways by DE
    fn = os.path.join(home, 'Dropbox', 'research', 'qmul', 'data', 'hgic_project', 'ipa_de_all.xlsx')
    df = pd.read_excel(fn)

    subgroups = {
        'RTK I': ['018', '019', '030', '031'],
        'RTK II': ['017', '050', '054', '061'],
        'MES': ['026', '052']
    }

    pids_ordered = subgroups['RTK I'] + subgroups['RTK II'] + subgroups['MES']

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    # separate into patients
    pids = df.columns[::2].str.replace('GBM', '')
    res = {}
    for i, pid in enumerate(pids):
        sub_df = df.iloc[1:, (2 * i):(2 * i + 2)]
        sub_df.columns = ['pathway', '-log_pval']
        sub_df.set_index('pathway', inplace=True)
        res[pid] = sub_df.dropna()

    # top 30
    p30 = {}
    p30_all = set()
    for pid, df in res.items():
        p30[pid] = df.index[:30]
        p30_all.update(p30[pid])

    p30_all = sorted(p30_all)

    df_in30 = pd.DataFrame(False, index=p30_all, columns=pids)
    df_30 = pd.DataFrame(index=p30_all, columns=pids, dtype=float)

    for pid, d in res.items():
        df_30.loc[p30[pid], pid] = res[pid].loc[p30[pid]].values
        df_in30.loc[p30[pid], pid] = True

    # plot these with a sort of clustering output
    highlight_pathways = [
        'Aryl Hydrocarbon Receptor Signaling',
        'cAMP-mediated signaling',
        'Tec Kinase Signaling'
    ]

    colours = sns.color_palette("hls", 3)

    # matrix for plotting
    mat = df_in30.astype(int).copy()
    for i, h in enumerate(highlight_pathways):
        mat.loc[h, mat.loc[h] > 0] = i + 2

    cmap = [(1, 1, 1), (0, 0, 0)] + colours
    ix = df_in30.sum(axis=1).sort_values(ascending=False).index

    fig = plt.figure(figsize=(5.5, 6.7))
    ax = fig.add_subplot(111)

    ax = sns.heatmap(mat.loc[ix, pids_ordered], cmap=cmap, cbar=False, ax=ax)
    ticks = [np.where(ix == h)[0][0] for h in highlight_pathways]
    ax.yaxis.set_ticks(ticks)
    ax.yaxis.set_ticklabels(highlight_pathways)
    ax.tick_params(axis=u'both', which=u'both', length=0)

    for i, t in enumerate(ax.yaxis.get_ticklabels()):
        t.set_color(colours[highlight_pathways.index(t.get_text())])

    group_colours = {
        'RTK I': '#007319',
        'RTK II': '#890000',
        'MES': '#8900A9'
    }

    for i, t in enumerate(ax.xaxis.get_ticklabels()):
        grp = subgroups_lookup[t.get_text()]
        c = group_colours[grp]
        t.set_color(c)
        t.set_rotation(90)

    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "hgic_de_ipa_top30.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "hgic_de_ipa_top30.tiff"), dpi=200)

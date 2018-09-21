import pandas as pd
from matplotlib import pyplot as plt, colors
import seaborn as sns
import os
import sys
import numpy as np
from plotting import common
import collections
from utils import output, setops, excel
import consts
from settings import HGIC_LOCAL_DIR


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    top_n = 30

    # pathways to be highlighted
    highlight_pathways = [
        'Aryl Hydrocarbon Receptor Signaling',
        'cAMP-mediated signaling',
        'Tec Kinase Signaling'
    ]

    # pathways by DE
    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq',
        's0_individual_patients_direct_comparison',
        'ipa',
        'pathways',
        'full_de_all.xls'
    )
    res = pd.read_excel(fn, sheet_name=None, index_col=0)

    subgroups = consts.SUBGROUPS

    pids = consts.PIDS
    subgroup_ind = collections.OrderedDict([
        (k, pd.Index(pids).isin(subgroups[k])) for k in ['RTK I', 'RTK II', 'MES']
    ])

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    # top 30
    p_top = {}
    p_top_all = set()
    for pid, df in res.items():
        p_top[pid] = df.sort_values(by='-log_p', ascending=False).index[:top_n]
        p_top_all.update(p_top[pid])

    p_top_all = sorted(p_top_all)

    df_in_top = pd.DataFrame(False, index=p_top_all, columns=pids)
    df_top = pd.DataFrame(index=p_top_all, columns=pids, dtype=float)

    for pid, d in res.items():
        df_top.loc[p_top[pid], pid] = res[pid].loc[p_top[pid], '-log_p'].values
        df_in_top.loc[p_top[pid], pid] = True

    colours = sns.color_palette("hls", len(highlight_pathways))

    # matrix for plotting
    mat = df_in_top.astype(int).copy()
    for i, h in enumerate(highlight_pathways):
        mat.loc[h, mat.loc[h] > 0] = i + 2

    cmap = colors.LinearSegmentedColormap.from_list(
        'foo',
        colors=[(1, 1, 1), (0, 0, 0)] + colours,
        N=len(highlight_pathways) + 2
    )

    # index is ordered from shared pathways to individual pathways
    # this is actually the opposite order we require, but the heatmap will invert the y axis
    ix = df_in_top.sum(axis=1).sort_values(ascending=False).index

    fig = plt.figure(figsize=(5.5, 6.7))
    ax = fig.add_subplot(111)

    # NB reverse the index order
    ax = sns.heatmap(mat.loc[ix[::-1], pids], cmap=cmap, cbar=False, ax=ax)
    ticks = [np.where(ix == h)[0][0] for h in highlight_pathways]
    ax.yaxis.set_ticks(ticks)
    ax.yaxis.set_ticklabels(highlight_pathways, rotation=0)
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

    # add dashed lines dividing different set sizes
    path_cuts = [1, 2, 3, 5]
    n_pat = df_in_top.loc[ix].sum(axis=1)
    prev_cut = len(n_pat)
    new_ticks = []
    new_ticklabels = []

    for pc in path_cuts:
        this_cut = np.where(n_pat <= pc)[0]
        if len(this_cut) == 0:
            print "Warning: no pathways present with less than or equal to %d members" % pc
        else:
            ax.axhline(this_cut[0], c='gray', ls='--', lw=2.5)
            new_ticks.append((this_cut[0] + prev_cut) / 2.)
            new_ticklabels.append("<= %d" % pc)
            prev_cut = this_cut[0]
    # final addition for greater number of patients
    new_ticks.append(prev_cut / 2.)
    new_ticklabels.append("> %d" % path_cuts[-1])

    ax2 = ax.twinx()
    ax2.invert_yaxis()
    ax2.set_ylim(ax.get_ylim())
    ax2.yaxis.set_ticks(new_ticks)
    ax2.yaxis.set_ticklabels(new_ticklabels, rotation=90, color='gray')
    ax2.set_ylabel("Number of patients sharing pathway", color='gray')
    ax2.grid(False)

    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "hgic_de_ipa_top%d.png" % top_n), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "hgic_de_ipa_top%d.tiff" % top_n), dpi=200)

    # export a wideform dataframe containing all these pathways with log_p, etc.
    for_export = {}
    for pid in pids:
        for_export[pid] = res[pid].loc[p_top[pid]]

    vs, vc = setops.venn_from_arrays(*[p_top[pid] for pid in pids])
    out = setops.venn_set_to_wide_dataframe(
        for_export,
        vs,
        pids,
        full_data=res,
        cols_to_include=['-log_p', 'ratio', 'z'],
        static_cols_to_include=['genes']
    )
    # excel.pandas_to_excel(out, os.path.join(outdir, "ipa_de_top_%d_pathways.xlsx" % top_n))
    out.to_excel(os.path.join(outdir, "ipa_de_top_%d_pathways.xlsx" % top_n))


    """
    Note to myself:
    I did consider an UpSet plot here. However, with the full DE lists, the result isn't very edifying...
    With the exception of patient-specific pathways, all sets have 2 or fewer pathways. 
    """
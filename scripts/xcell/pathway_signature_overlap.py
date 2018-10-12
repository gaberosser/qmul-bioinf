import os
import pandas as pd
from settings import GIT_LFS_DATA_DIR, HGIC_LOCAL_DIR
from utils import setops, output
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import common
from scipy.ndimage import measurements
import numpy as np


if __name__ == '__main__':
    # minimum pval for pathways to be used
    alpha = 0.005
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05

    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    comparisons = ['syngeneic', 'gibco', 'h9']

    outdir = output.unique_output_dir()

    ipa_indir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways/'
    )
    xcell_sign_fn = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')

    xcell_s = pd.read_excel(xcell_sign_fn, header=0, index_row=0)
    xcell_signatures = {}
    for i, row in xcell_s.iterrows():
        xcell_signatures[row.Celltype_Source_ID] = set(row.iloc[2:].dropna().values)

    # load IPA pathway data
    # use a pre-assembled file of only highly significant pathways
    ipa_fn = os.path.join(ipa_indir, "full_de_ipa_results_significant.xlsx")
    ipa_res = pd.read_excel(ipa_fn, header=0, index_col=0)

    # load genes from raw reports
    file_patt = 'de_s2_{pid}_{cmp}.txt'
    ipa_pathway_signatures = {}
    for pid in pids:
        for c in comparisons:
            fn = os.path.join(ipa_indir, file_patt.format(pid=pid, cmp=c))
            this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
            this.columns = ['-logp', 'ratio', 'z', 'genes']
            this.index = [x.decode('utf-8') for x in this.index]
            for pw in this.index.intersection(ipa_res.index):
                this_list = set(this.loc[pw, 'genes'].split(','))
                if pw in ipa_pathway_signatures:
                    ipa_pathway_signatures[pw] = ipa_pathway_signatures[pw].union(this_list)
                else:
                    ipa_pathway_signatures[pw] = this_list

    # compare in a pairwise fashion
    so_ipa_not_cts = {}
    so_cts_not_ipa = {}
    so_both = {}
    all_ = [so_ipa_not_cts, so_cts_not_ipa, so_both]

    for pw, pw_arr in ipa_pathway_signatures.items():
        for a in all_:
            if pw not in a:
                a[pw] = {}
        for cts, cts_arr in xcell_signatures.items():
            m = len(pw_arr)
            n = len(cts_arr)
            its = len(pw_arr.intersection(cts_arr))
            so_ipa_not_cts[pw][cts] = m - its
            so_cts_not_ipa[pw][cts] = n - its
            so_both[pw][cts] = its

    so_ipa_not_cts = pd.DataFrame(so_ipa_not_cts)
    so_cts_not_ipa = pd.DataFrame(so_cts_not_ipa)
    so_both = pd.DataFrame(so_both)

    # now let's look at only the cell types of interest
    cell_types = sorted([
        'Basophils',
        'NKT',
        'CD4+ Tcm',
        'Mesangial cells',
        'mv Endothelial cells',
        'pDC',
        'iDC',
        'Memory B-cells',
        'B-cells',
        'Epithelial cells',
        'Melanocytes',
        'Class-switched memory B-cells',
        'Pericytes',
        'Osteoblast',
        'Th1 cells',
        'CD8+ naive T-cells',
        'CLP',
        'Smooth muscle',
        'Macrophages',
        'Macrophages M1',
        'Macrophages M2',
        'Neurons',
        'CD4+ T-cells',
        'HSC',
        'Mast cells',
        'Tregs',
    ])
    ix = []
    for ct in cell_types:
        ix.extend(so_both.index[so_both.index.str.contains(ct + '_')])

    pct_shared = (so_both.loc[ix] / (so_both.loc[ix] + so_ipa_not_cts.loc[ix]) * 100.).sort_index().transpose()

    col_colours = pd.Series(index=pct_shared.columns, name='Cell type')
    ix_lookup = pct_shared.columns.str.replace(r'(?P<ct>[^_]*)_.*', r'\g<ct>')
    a, b = ix_lookup.factorize()
    cc_cmap = common.get_best_cmap(a.max() + 1, cmap='jet')
    for i in range(a.max() + 1):
        col_colours[a == i] = cc_cmap[i]

    cg = sns.clustermap(
        pct_shared,
        cmap='Reds',
        row_cluster=False,
        col_cluster=False,
        # col_colors=col_colours,
        vmax=20.
    )
    # cg.gs.set_height_ratios([0.1, 0.005, 0.05, 0.9])
    cg.gs.set_height_ratios([0.1, 0.005, 0.9])
    cg.gs.set_width_ratios([0.04, 0.02, 0.8])
    cg.gs.update(top=0.99, left=0.02, right=0.67, bottom=0.15)

    # X labels: just one per group
    labels = a + 1
    com = measurements.center_of_mass(
        np.ones(len(ix_lookup)),
        labels,
        range(1, labels.max() + 1)
    )
    xticks = [t[0] for t in com]
    cg.ax_heatmap.xaxis.set_ticks(xticks)
    cg.ax_heatmap.xaxis.set_ticklabels(b)

    plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90, fontsize=6)
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), rotation=0, fontsize=6)

    # overlay outlines for values above threshold
    pct_max = 10.
    sp_masked = np.ma.masked_where(pct_shared < pct_max, np.zeros_like(pct_shared))
    cg.ax_heatmap.pcolor(
        sp_masked[::-1],
        edgecolors='k',
        facecolor='none',
        linewidths=1.,
        cmap='Greys_r'
    )
    cg.savefig(os.path.join(outdir, "genes_shared_pathways_cell_signatures.png"), dpi=200)
    cg.savefig(os.path.join(outdir, "genes_shared_pathways_cell_signatures.pdf"), dpi=200)
    cg.savefig(os.path.join(outdir, "genes_shared_pathways_cell_signatures.tiff"), dpi=200)

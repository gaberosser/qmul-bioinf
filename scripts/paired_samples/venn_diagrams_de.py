from matplotlib_venn import venn3, venn2
from matplotlib import pyplot as plt
import pandas as pd


if __name__ == "__main__":
    fn = '/home/gabriel/r_outputs/paired_analysis_de_rtkI.6/de_gene_lists.threeway.xls'
    dat = pd.read_excel(fn, sheetname=None)
    setnames = ('iNSC', 'eNSC_H9', 'eNSC_foetal')

    de_all = {}
    de_up = {}
    de_down = {}

    # all
    for i in range(1, 8):
        b = "{0:03b}".format(i)
        s = ' & '.join([setnames[j] for j in range(3) if b[j] == '1'])
        this = dat[s]
        de_all[b] = this.shape[0]
        d = this.loc[:, this.columns.str.contains('direction')]
        de_up[b] = (d == 'U').all(axis=1).sum()
        de_down[b] = (d == 'D').all(axis=1).sum()

    set_labels = ('hGIC:paired iPS NSC', 'hGIC:H9 NSC', 'hGIC:foetal NSC')

    fig, axs = plt.subplots(ncols=3)

    v_all = venn3(subsets=de_all, set_labels=set_labels, ax=axs[0])
    v_up = venn3(subsets=de_up, set_labels=None, ax=axs[1])
    v_down = venn3(subsets=de_down, set_labels=None, ax=axs[2])

    for v in [v_all, v_up, v_down]:
        plt.setp(v.patches, 'alpha', 0.4)

    v_all = venn3(subsets=de_all, set_labels=set_labels, ax=axs[0])
    v_up = venn3(subsets=de_up, set_labels=None, ax=axs[1])
    v_down = venn3(subsets=de_down, set_labels=None, ax=axs[2])

    for v in [v_all, v_up, v_down]:
        plt.setp(v.patches, 'facecolor', 'none')
        plt.setp(v.patches, 'edgecolor', 'k')


    # two way comparison
    de_all = {}
    de_up = {}
    de_down = {}

    # all
    for i in range(1, 4):
        b = "{0:02b}".format(i)
        s = ' & '.join([setnames[j] for j in range(2) if b[j] == '1'])
        this = dat[s]
        de_all[b] = this.shape[0]
        d = this.loc[:, this.columns.str.contains('direction')]
        de_up[b] = (d == 'U').all(axis=1).sum()
        de_down[b] = (d == 'D').all(axis=1).sum()

    set_labels = ('hGIC:paired iPS NSC', 'hGIC:H9 NSC')
    set_colors = ('b', 'r')
    fig, axs = plt.subplots(nrows=3)

    v_all = venn2(subsets=de_all, set_labels=set_labels, set_colors=set_colors, ax=axs[0])
    v_up = venn2(subsets=de_up, set_labels=None, set_colors=set_colors, ax=axs[1])
    v_down = venn2(subsets=de_down, set_labels=None, set_colors=set_colors, ax=axs[2])

    for v in [v_all, v_up, v_down]:
        plt.setp(v.patches, 'alpha', 0.4)

    v_all = venn2(subsets=de_all, set_labels=None, set_colors=set_colors, ax=axs[0])
    v_up = venn2(subsets=de_up, set_labels=None, set_colors=set_colors, ax=axs[1])
    v_down = venn2(subsets=de_down, set_labels=None, set_colors=set_colors, ax=axs[2])

    for v in [v_all, v_up, v_down]:
        plt.setp(v.patches, 'facecolor', 'none')
        plt.setp(v.patches, 'edgecolor', 'k')
        plt.setp(v.patches, 'zorder', '2')
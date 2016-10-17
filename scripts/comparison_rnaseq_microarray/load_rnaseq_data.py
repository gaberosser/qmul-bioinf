import pandas as pd
import numpy as np
import collections
import os

RNASEQ_GENE_COUNTS_DIR = 'data/rnaseq_GSE83696/'


def load_rnaseq_cufflinks_gene_count_data(unit='fpkm', return_ci=False):
    """
    Load in FPKM data from pre-existing Cufflinks run.
    Optionally convert units to TPM (values sum to 1). This is useful for comparison.
    NB: the confidence intervals are BROKEN, so don't use them.
    :param unit: Either 'fpkm' (do nothing) or 'tpm' (convert to TPM units)
    :param return_ci: If True return the confidence intervals. Not worth it: they are broken and meaningless.
    :return:
    """
    if unit not in ('fpkm', 'tpm'):
        raise ValueError("Unsupported unit %s. Supported options are tpm, fpkm.", unit)

    infiles = {
        'XZ1': 'XZ-1.genes.fpkm_tracking',
        'XZ2': 'XZ-2.genes.fpkm_tracking',
    }

    res = pd.DataFrame()
    res_ci_lo = pd.DataFrame()
    res_ci_hi = pd.DataFrame()

    for tag, fn in infiles.items():
        ff = os.path.join(RNASEQ_GENE_COUNTS_DIR, fn)
        t = pd.read_csv(ff, header=0, index_col=0, sep='\t')

        r = t.loc[t.FPKM_status == 'OK', 'FPKM']
        r = r.groupby(r.index).sum()
        res[tag] = r

        r_lo = t.loc[t.FPKM_status == 'OK', 'FPKM_conf_lo']
        r_lo = r_lo.groupby(r_lo.index).sum()
        res_ci_lo[tag] = r_lo

        r_hi = t.loc[t.FPKM_status == 'OK', 'FPKM_conf_hi']
        r_hi = r_hi.groupby(r_hi.index).sum()
        res_ci_hi[tag] = r_hi

        if unit == 'tpm':
            k = res[tag].sum()
            res[tag] /= k
            res_ci_lo[tag] /= k
            res_ci_hi[tag] /= k

    if return_ci:
        return res, res_ci_lo, res_ci_hi
    else:
        return res


def rnaseq_activity_in_groups(df):
    groups = (
        ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
        ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
        ('Group C', ('IMPG2', 'GABRA5', 'EYS', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
        ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
    )

    icb_cols = ('XZ1_FPKM', 'XZ2_FPKM')
    icb_cols_lo = tuple((t + '_conf_lo' for t in icb_cols))
    icb_cols_hi = tuple((t + '_conf_hi' for t in icb_cols))

    ci_mean = collections.OrderedDict()
    ci_lo = collections.OrderedDict()
    ci_hi = collections.OrderedDict()

    for g, arr in groups:
        ci_mean[g] = df.loc[arr, icb_cols].mean(axis=1)
        ci_lo[g] = df.loc[arr, icb_cols_lo].mean(axis=1)
        ci_hi[g] = df.loc[arr, icb_cols_hi].mean(axis=1)

    return ci_lo, ci_mean, ci_hi


def plot_rnaseq_activity_in_groups(df):

    from matplotlib import pyplot as plt
    plt.interactive(True)
    width = 0.8

    # df = load_rnaseq_cufflinks_gene_count_data()
    ci_lo, ci_centre, ci_hi = rnaseq_activity_in_groups(df)

    # identify global max/min for floor/ceil purposes
    ci_max = ci_min = 0
    for v in ci_lo.values():
        ci_min = min(ci_min, min(v[v > -np.inf]))
    for v in ci_hi.values():
        ci_max = max(ci_max, max(v[v < np.inf]))
    for v in ci_centre.values():
        ci_min = min(ci_min, min(v[v > -np.inf]))
        ci_max = max(ci_max, max(v[v < np.inf]))

    ci_min = np.floor(ci_min)
    ci_max = np.ceil(ci_max)

    fig, axes = plt.subplots(ncols=len(ci_lo), sharex=False, sharey=True)
    fig.subplots_adjust(wspace=0)

    first = True
    for ax, grp in zip(axes, ci_lo.keys()):
        x = np.arange(len(ci_lo[grp]))
        bottom = ci_lo[grp]
        height = ci_hi[grp] - ci_lo[grp]

        # two colour groups
        colours = []
        for b, h in zip(bottom, height):
            if b == -np.inf or h == np.inf:
                colours.append('red')
            else:
                colours.append('gray')

        bottom[bottom < ci_min] = ci_min
        height[height > ci_max] = ci_max

        ax.bar(x, height, width=width, bottom=bottom, color=colours, ec='k')
        ax.plot(x + width / 2., ci_centre[grp], 'ko', ms=5)
        # ax.plot([x[0] - 1., x[-1] + 2.], [0, 0], 'k--')
        ax.set(
            xticks=x + width / 2.,
            xticklabels=ci_centre[grp].index,
            xlabel=grp,
            xlim=[width - 1., len(ci_lo[grp])],
            ylim=[ci_min, ci_max]
        )
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=30, fontsize=12)
        ax.margins(0.05)  # Optional
        if first:
            ax.set(ylabel='FPKM')
            first = False

    plt.tight_layout(w_pad=0.)
    plt.show()


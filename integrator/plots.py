from matplotlib import pyplot as plt
import numpy as np
import os
from scipy import stats
from stats.basic import construct_contingency
from methylation import dmr
import seaborn as sns


def scatter_de_dmr_single(
        joint_de_dmr,
        de_fc_col='de_logfc',
        dmr_fc_col='dmr_median_delta',
        r_threshold=None,
        ax=None
):
    x = joint_de_dmr.loc[:, de_fc_col]
    y = joint_de_dmr.loc[:, dmr_fc_col]

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    idx_ul = ((x < 0) & (y > 0))
    ax.scatter(x.loc[idx_ul], y.loc[idx_ul], c='r')
    idx_br = ((x > 0) & (y < 0))
    ax.scatter(x.loc[idx_br], y.loc[idx_br], c='b')
    idx_od = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
    ax.scatter(x.loc[idx_od], y.loc[idx_od], c='gray')

    if r_threshold is not None:
        r = (x ** 2 + y ** 2) ** .5
        for i in np.where(idx_ul & (r > 8))[0]:
            ax.text(x.iloc[i], y.iloc[i], x.index[i], color='r')
        for i in np.where(idx_br & (r > 5))[0]:
            ax.text(x.iloc[i], y.iloc[i], x.index[i], color='b')

    ax.axhline(0, ls='--', c='k', alpha=0.4)
    ax.axvline(0, ls='--', c='k', alpha=0.4)
    ax.set_xlabel("DE logFC")
    ax.set_ylabel("DMR median delta")

    return ax


def scatter_plot_dmr_de(meth_de, fig_filestem, fig_titlestem='', outdir=None):



    for sid in meth_de:
        fig_title = ("%s %s" % (fig_titlestem, sid)).strip()
        fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, num=fig_title)
        for i, cls in enumerate(['all', 'tss', 'gene', 'island']):
            ax = axs.flat[i]

            # get values for ALL DMR clusters of this class
            x = meth_de[sid][cls].loc[:, 'logFC'].astype(float)
            y = meth_de[sid][cls].loc[:, 'me_mediandelta'].astype(float)

            # contingency table for Fisher's exact test
            conting = construct_contingency(x, y)
            logodds, fisherp = stats.fisher_exact(conting)

            # get values for DMR clusters that are ONLY in this class (no overlaps)
            if cls == 'all':
                print "%s - %s p = %.3e" % (
                    sid, cls, fisherp
                )
                ax.scatter(x, y, c='k')
                ax.axhline(0, c=0.4 * np.ones(3))
                ax.axvline(0, c=0.4 * np.ones(3))
                ttl = "%s; p={0}" % cls
                if fisherp < 0.001:
                    ttl = ttl.format('%.2e') % fisherp
                else:
                    ttl = ttl.format('%.3f') % fisherp

            else:
                cid_other = set(
                    np.unique(
                        np.concatenate([meth_de[sid][t].me_cid.values for t in dmr.CLASSES.difference({cls, })])
                    )
                )
                xu = x.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                yu = y.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                contingu = construct_contingency(xu, yu)
                logoddsu, fisherpu = stats.fisher_exact(contingu)

                print "%s - %s p = %.3e (incl overlaps), p = %.3e (excl overlaps)" % (
                    sid, cls, fisherp, fisherpu
                )

                ax.scatter(x, y, c='k')
                ax.scatter(xu, yu, c='b')
                ax.axhline(0, c=0.4 * np.ones(3))
                ax.axvline(0, c=0.4 * np.ones(3))

                ttl = "%s; p={0}; unique p={0}" % cls

                if fisherp < 0.001:
                    ttl = ttl.format('%.2e') % (fisherp, fisherpu)
                else:
                    ttl = ttl.format('%.3f') % (fisherp, fisherpu)

            ax.set_title(ttl)

        fig.text(0.5, 0.04, 'RNASeq DE logFC', ha='center')
        fig.text(0.04, 0.5, 'EPIC DMR median delta M', va='center', rotation='vertical')
        if outdir is not None:
            fn = os.path.join(outdir, "{0}_{1}.{{ext}}".format(fig_filestem, sid))
            fig.savefig(fn.format(ext='png'), dpi=200)
            fig.savefig(fn.format(ext='pdf'))
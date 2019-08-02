import os

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

from load_data import rnaseq_data
from rnaseq.filter import filter_by_cpm
from utils import output

if __name__ == '__main__':
    outdir = output.unique_output_dir("gmb061")
    pids = ['017', '050', '054', '061']
    # pids = [t for t in rnaseq_data.PATIENT_LOOKUP_STAR if t != 'GIBCO']
    obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    obj.data = obj.data.loc[obj.data.index.str.contains('ENSG')]
    dat_filt = filter_by_cpm(obj.data, min_n_samples=2)

    # look at the distribution of counts (ECDF)

    colours = plt.cm.jet(np.linspace(0, 1., obj.meta.shape[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    focus = 'DURA061_NSC_N6_P4'
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)
        cdf = the_data.cumsum() / float(the_data.sum())
        alpha = None if sname == focus else 0.3
        ax.plot(np.log2(the_data + 1), pct, c=colours[i], label=sname, alpha=alpha)
        # ax.plot(the_data, cdf, c=colours[i], label=sname)
    ax.legend(loc='lower right')
    ax.set_xlabel('Log2(count + 1)')
    ax.set_ylabel('Percentile')

    fig.savefig(os.path.join(outdir, 'percentile_vs_log_count.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_log_count.pdf'))

    # repeat and zoom in on the top 1% most expressed genes
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)

        ax.plot(the_data.loc[pct>0.99], pct[pct>0.99], c=colours[i], label=sname)
    ax.legend(loc='lower right')
    ax.set_xlabel('Count')
    ax.set_ylabel('Percentile')
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct.pdf'))

    # highlight DURA061 samples in the same plot
    colours = plt.cm.prism(np.linspace(0, 1., obj.meta.shape[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, sname in enumerate(obj.meta.index):
        the_data = dat_filt.loc[:, sname].sort_values()
        pct = stats.rankdata(the_data, method='max') / float(the_data.size)
        c = colours[i] if 'DURA061' in sname else 'gray'
        alpha = None if 'DURA061' in sname else 0.3
        label = sname if 'DURA061' in sname else None
        ax.plot(the_data.loc[pct > 0.99], pct[pct > 0.99], c=c, label=label, alpha=alpha)
    ax.legend(loc='lower right')
    ax.set_xlabel('Count')
    ax.set_ylabel('Percentile')
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct_focus061NSC.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'percentile_vs_count_top1pct_focus061NSC.pdf'))

    # scatterplots where replicates exist (all but NSC017)
    fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True)
    for i, pid in enumerate(pids):
        for j, typ in enumerate(('GBM', 'DURA')):
            the_idx = dat_filt.columns.str.contains("%s%s" % (typ, pid))
            if the_idx.sum() != 2:
                axs[j][i].set_title("%s%s (n=1)" % (typ, pid))
                continue
            the_data = np.log2(dat_filt.loc[:, the_idx] + 1.)
            axs[j][i].scatter(the_data.iloc[:, 0], the_data.iloc[:, 1])

            axs[j][i].set_title("%s%s r=%.3f" % (
                typ, pid, stats.linregress(the_data.iloc[:, 0], the_data.iloc[:, 1]).rvalue
            ))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'replicate_logcorr.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'replicate_logcorr.pdf'))

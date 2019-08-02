import pandas as pd
from rnaseq import loader, differential_expression, filter
from matplotlib import pyplot as plt
import seaborn as sns
from settings import GIT_LFS_DATA_DIR, DATA_DIR
from utils import output
import os


if __name__ == "__main__":
    outdir = output.unique_output_dir("small_rnaseq_overview")
    min_cpm = 1
    de_params = {
        'method': 'QLGLM',
        'lfc': 0,
        'fdr': 0.05
    }
    # load featurecounts summaries and generate a plot
    basedir = os.path.join(DATA_DIR, 'small_rnaseq')
    indirs = [
        'wtchg_p170389',
        'wtchg_p170645',
        'encode_roadmap/ENCSR511YYU',
        'GSE64744',
    ]
    study_names = ['P170389', 'P170645', 'ENCODE', 'GSE64744']
    metas = []
    for d, sn in zip(indirs, study_names):
        meta_fn = os.path.join(basedir, d, 'sources.csv')
        this_meta = pd.read_csv(meta_fn, header=0, index_col=0)
        if 'species' in this_meta:
            this_meta = this_meta.loc[this_meta.species != 'mouse']

        if 'tissue_type' in this_meta:
            # split into FFPE and cell culture
            this_meta_ff = this_meta.loc[this_meta.tissue_type == 'ffpe']
            this_meta = this_meta.loc[this_meta.tissue_type != 'ffpe']
            this_meta_ff.insert(this_meta_ff.shape[1], 'study', "%s (FFPE)" % sn)
            metas.append(this_meta_ff)

        this_meta.insert(this_meta.shape[1], 'study', sn)
        metas.append(this_meta)
    meta = pd.concat(metas, axis=0)

    meta.insert(meta.shape[1], 'pct_trimmed', 100 - meta.read_count / meta.read_count_pre_trimming.astype(float) * 100.)
    meta.insert(meta.shape[1], 'mapping_pct', meta.num_reads_mapped / meta.read_count.astype(float) * 100.)
    meta.insert(meta.shape[1], 'assignment_pct', meta.num_reads_assigned / meta.num_reads_mapped.astype(float) * 100.)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.boxplot(x='study', y='pct_trimmed', data=meta, ax=ax)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, 'pct_trimmed.png'), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.boxplot(x='study', y='mapping_pct', data=meta, ax=ax)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, 'mapping_pct.png'), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.boxplot(x='study', y='assignment_pct', data=meta, ax=ax)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, 'assignment_pct.png'), dpi=200)

    # now run differential expression on some of our samples
    fn = os.path.join(GIT_LFS_DATA_DIR, 'small_rnaseq_2018-04-02', 'p170645.quant')
    dat1 = pd.read_csv(fn, sep='\t', header=0, index_col=0, comment='#')
    dat1 = dat1.loc[:, dat1.columns.sort_values()]
    dat1.columns = dat1.columns.str.replace(
        r'.*/wtchg_p170645/human/bwa_alignment/',
        ''
    )
    dat1.columns = dat1.columns.str.replace('.sorted.bam', '')
    dat1 = dat1.loc[:, dat1.columns.str.contains(r'[0-9]')]
    m1 = pd.read_csv(
        os.path.join(basedir, 'wtchg_p170645', 'sources.csv'),
        header=0, index_col=0
    )
    m1.index = m1.index.astype(str)
    dat1.columns = m1.loc[dat1.columns, 'sample']

    fn = os.path.join(GIT_LFS_DATA_DIR, 'small_rnaseq_2018-04-02', 'p170389.quant')
    dat2 = pd.read_csv(fn, sep='\t', header=0, index_col=0, comment='#')
    dat2 = dat2.loc[:, dat2.columns.sort_values()]
    dat2.columns = dat2.columns.str.replace(
        r'.*/wtchg_p170389/human/bwa_alignment/',
        ''
    )
    dat2.columns = dat2.columns.str.replace('.sorted.bam', '')
    dat2 = dat2.loc[:, dat2.columns.str.contains(r'[0-9]')]
    m2 = pd.read_csv(
        os.path.join(basedir, 'wtchg_p170389', 'sources.csv'),
        header=0, index_col=0
    )
    m2.index = m2.index.astype(str)
    dat2.columns = m2.loc[dat2.columns, 'sample']

    dat = pd.concat((dat1, dat2), axis=1)
    m = pd.concat((m1, m2), axis=0).set_index('sample')

    m.loc[m.loc[:, 'type'] == 'GIC', 'type'] = 'GBM'

    # GBM vs iNSC DE
    # old style (one-by-one comparison)

    pids = ['017', '018', '019', '026', '030', '031', '044', '050', '052', '054']
    de_res_ind = {}

    for pid in pids:
        idx = m.index.str.contains(pid)
        the_dat = dat.loc[:, m.loc[idx].index.astype(str)]
        the_groups = m.loc[idx, 'type']
        the_groups.index = the_groups.index.astype(str)

        # filter
        the_dat = the_dat.loc[
            filter.filter_cpm_by_group(the_dat.divide(the_dat.sum(), axis=1) * 1e6, the_groups, min_cpm=min_cpm)
        ]

        the_comparison = ['GBM', 'iNSC']
        de_res_ind[pid] = differential_expression.run_one_de(
            the_dat,
            the_groups,
            the_comparison,
            add_gene_symbols=False,
            **de_params
        )

    print "Individual dispersion"
    for k, t in sorted(de_res_ind.items(), key=lambda x: x[0]):
        print "%s, %d DE miRNA features" % (k, t.shape[0])


    # new style
    de_res_ensembl = {}

    idx = m.index.str.contains('|'.join(pids))
    the_dat = dat.loc[:, m.loc[idx].index.astype(str)]
    the_groups = pd.Series(index=the_dat.columns)
    for pid in pids:
        the_groups[the_groups.index.str.contains('GBM') & the_groups.index.str.contains(pid)] = "GBM%s" % pid
        the_groups[the_groups.index.str.contains('NSC') & the_groups.index.str.contains(pid)] = "iNSC%s" % pid

    # initial filter
    the_cpm = the_dat.divide(the_dat.sum(), axis=1) * 1e6
    keep = (the_cpm > min_cpm).sum(axis=1) > 0
    the_dat = the_dat.loc[keep]

    for pid in pids:
        the_comparison = ('GBM%s' % pid, 'iNSC%s' % pid)
        the_de_res = differential_expression.run_one_de(
            the_dat,
            the_groups,
            the_comparison,
            add_gene_symbols=False,
            **de_params
        )
        # filter again
        this_cpm = the_cpm.loc[the_de_res.index, the_dat.columns.str.contains(pid)]
        this_groups = the_groups.loc[the_groups.index.str.contains(pid)]
        keep = filter.filter_cpm_by_group(this_cpm, this_groups, min_cpm=min_cpm)

        de_res_ensembl[pid] = the_de_res.loc[keep]

    print "Grouped dispersion"
    for k, t in sorted(de_res_ensembl.items(), key=lambda x: x[0]):
        print "%s, %d DE miRNA features" % (k, t.shape[0])

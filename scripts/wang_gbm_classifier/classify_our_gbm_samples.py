import pandas as pd
import numpy as np
from rnaseq import gsea, general, loader
from load_data import rnaseq_data
from utils.output import unique_output_dir
from settings import OUTPUT_DIR
import os
import references
from matplotlib import pyplot as plt
import seaborn as sns


SRC_MAP = {
    'star': 'counts',
    'salmon': 'tpm',
    'star/cufflinks': 'rpkm'
}


def prepare_gct_files_hgic(pids='all', outdir=None):
    """
    Prepare the GCT files required to perform classification of the hGIC samples:
    - hGIC FFPE
    - hGIC cell culture
    - Both combined
    In all cases, use FPKM units (cufflinks), TPM (salmon) and CPM (STAR).
    Use gene symbols as these are contained in the signatures.
    """
    ## FIXME: finish!
    if outdir is None:
        outdir = unique_output_dir("gct_files_for_wang")

    loaded = {}
    for typ in ('cell_culture', 'ffpe'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            this_dat = loader.load_by_patient(pids, type=typ, source=src, include_control=False)
            this_dat.meta = this_dat.meta.loc[~this_dat.meta.index.str.contains('DURA')]
            this_dat.data = this_dat.data.loc[:, this_dat.meta.index]
            loaded.setdefault(typ, {})[src] = this_dat

    for typ in ('cell_culture', 'ffpe'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            this_dat = loaded[typ][src]
            this_dat = references.translate_quantification_resolving_duplicates(
                this_dat.data,
                'Ensembl Gene ID',
                'Approved Symbol',
            )
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            gsea.data_to_gct(this_dat, fn)

    # combine
    for src in ('star', 'salmon', 'star/cufflinks'):
        this_dat1 = loaded['cell_culture'][src].data
        this_dat2 = loaded['ffpe'][src].data
        this_dat = pd.concat((this_dat1, this_dat2), axis=1)
        this_dat = references.translate_quantification_resolving_duplicates(
            this_dat,
            'Ensembl Gene ID',
            'Approved Symbol',
        )
        fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], 'all'))
        gsea.data_to_gct(this_dat, fn)

    infiles = []

    # FFPE
    ffpe_dat = rnaseq_data.load_salmon_by_patient_id('all', type='ffpe', include_control=False, units=units)
    # transcript -> gene
    ffpe_dat = general.ensembl_transcript_quant_to_gene(ffpe_dat)
    # ensembl -> gene symbol
    ffpe_dat = references.translate_quantification_resolving_duplicates(
        ffpe_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_ffpe_tpm.gct")
    gsea.data_to_gct(ffpe_dat, fn)
    infiles.append(fn)

    # GIC
    gic_dat = rnaseq_data.load_salmon_by_patient_id('all', type='cell_culture', include_control=False, units=units)
    # only keep GBM lines
    gic_dat = gic_dat.loc[:, gic_dat.columns.str.contains('GBM')]
    # transcript -> gene
    gic_dat = general.ensembl_transcript_quant_to_gene(gic_dat)
    # ensembl -> gene symbol
    gic_dat = references.translate_quantification_resolving_duplicates(
        gic_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_cc_tpm.gct")
    gsea.data_to_gct(gic_dat, fn)
    infiles.append(fn)

    # 3) Combined
    dat = gsea.combine_gct_files(*infiles)
    fn = os.path.join(outdir, "gbm_ffpe_and_cc_tpm.gct")
    gsea.data_to_gct(dat, fn)
    infiles.append(fn)

    return infiles


def prepare_gct_files_salmon(outdir=None, units='tpm'):
    """
    Prepare the GCT files required to perform classification:
    - Our GBM FFPE and cell culture samples
    - TCGA RNA-Seq cohort
    - Both combined
    In all cases, use FPKM units and gene symbols, as these are used by Wang
    """
    if outdir is None:
        outdir = unique_output_dir("gct_files_for_wang")

    infiles = []

    # FFPE
    ffpe_dat = rnaseq_data.load_salmon_by_patient_id('all', type='ffpe', include_control=False, units=units)
    # transcript -> gene
    ffpe_dat = general.ensembl_transcript_quant_to_gene(ffpe_dat)
    # ensembl -> gene symbol
    ffpe_dat = references.translate_quantification_resolving_duplicates(
        ffpe_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_ffpe_tpm.gct")
    gsea.data_to_gct(ffpe_dat, fn)
    infiles.append(fn)

    # GIC
    gic_dat = rnaseq_data.load_salmon_by_patient_id('all', type='cell_culture', include_control=False, units=units)
    # only keep GBM lines
    gic_dat = gic_dat.loc[:, gic_dat.columns.str.contains('GBM')]
    # transcript -> gene
    gic_dat = general.ensembl_transcript_quant_to_gene(gic_dat)
    # ensembl -> gene symbol
    gic_dat = references.translate_quantification_resolving_duplicates(
        gic_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_cc_tpm.gct")
    gsea.data_to_gct(gic_dat, fn)
    infiles.append(fn)

    # 3) Combined
    dat = gsea.combine_gct_files(*infiles)
    fn = os.path.join(outdir, "gbm_ffpe_and_cc_tpm.gct")
    gsea.data_to_gct(dat, fn)
    infiles.append(fn)

    return infiles


def prepare_gct_files_star(outdir=None):
    """
    Prepare the GCT files required to perform classification:
    - Our GBM FFPE and cell culture samples
    - TCGA RNA-Seq cohort
    - Both combined
    In all cases, use FPKM units and gene symbols, as these are used by Wang
    """
    if outdir is None:
        outdir = unique_output_dir("gct_files_for_wang")

    infiles = []

    # FFPE
    ffpe_dat = rnaseq_data.load_by_patient('all', type='ffpe', include_control=False)
    # ensembl -> gene symbol
    ffpe_dat = references.translate_quantification_resolving_duplicates(
        ffpe_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_ffpe_counts.gct")
    gsea.data_to_gct(ffpe_dat, fn)
    infiles.append(fn)

    # GIC
    gic_dat = rnaseq_data.load_by_patient('all', type='cell_culture', include_control=False)
    # only keep GBM lines
    gic_dat = gic_dat.loc[:, gic_dat.columns.str.contains('GBM')]
    # ensembl -> gene symbol
    gic_dat = references.translate_quantification_resolving_duplicates(
        gic_dat,
        'Ensembl Gene ID',
        'Approved Symbol',
    )
    fn = os.path.join(outdir, "gbm_cc_counts.gct")
    gsea.data_to_gct(gic_dat, fn)
    infiles.append(fn)

    # 3) Combined
    dat = gsea.combine_gct_files(*infiles)
    fn = os.path.join(outdir, "gbm_ffpe_and_cc_counts.gct")
    gsea.data_to_gct(dat, fn)
    infiles.append(fn)

    return infiles


def load_pvalue_results(fn):
    dat = pd.read_csv(fn, header=0, index_col=0, delimiter='\t')
    # only keep the p values
    ncol = dat.columns.size
    dat = dat.iloc[:, (ncol / 2):]
    dat.columns = dat.columns.str.replace('_pval', '')
    return dat


def simplicity_score(pvals):
    """
    For each sample (row), compute the simplicity score defined in Wang et al.
    :param pvals:
    :return:
    """
    # Rank the pvalues. This method chooses the first column it encounters in the event of a tie. This is fine as it
    # doesn't affect the outcome.
    n_cls = pvals.columns.size
    if n_cls < 2:
        raise AttributeError("Cannot compute a simplicity score with fewer than 2 classes")
    rnk = pvals.rank(axis=1, method='first')
    adds = pd.Series(index=pvals.index)
    adns = pd.Series(index=pvals.index)
    rng = pd.Series(index=pvals.index)
    for ix in pvals.index:
        p = pvals.loc[ix].values
        r = rnk.loc[ix].values
        p0 = p[r == 1]
        adds.loc[ix] = (p[r > 1] - p0).sum()
        this_adns = 0.
        for i in range(2, n_cls + 1):
            for j in range(i, n_cls + 1):
                this_adns += (p[r == j] - p[r == i])
        adns.loc[ix] = this_adns
        rng.loc[ix] = p[r == n_cls] - p0
    return (adds - adns) * rng / float(n_cls - 1)


def contingency_table(new, previous, vals=None, val_func=np.mean):
    """
    Previous values go on the INDEX, new values on the COLUMNS
    :param new:
    :param previous:
    :param vals:
    :return:
    """
    _, new_cls = new.factorize()
    new_cls = set(new_cls).difference({'None', 'Multi'})

    _, prev_cls = previous.factorize()
    prev_cls = set(prev_cls).difference({'None', 'Multi'})

    new_idx = list(new_cls) + ['Multi', 'None']
    prev_idx = list(prev_cls) + ['Multi', 'None']

    ctg = pd.DataFrame(
        index=prev_idx,
        columns=new_idx
    )

    for ix in ctg.index:

        if ix == "None":
            the_ids = previous.loc[previous.isnull()].index
        else:
            the_ids = previous.loc[previous == ix].index
        if len(the_ids) == 0:
            continue

        for col in ctg.columns:
            the_match = new.loc[the_ids]
            if col == "None":
                this_ix = the_match.isnull()
            else:
                this_ix = (the_match == col)
            if vals is None:
                # just count
                ctg.loc[ix, col] = this_ix.sum()
            else:
                # store values
                if val_func is None:
                    ctg.loc[ix, col] = vals.loc[this_ix.index[this_ix]].tolist()
                else:
                    ctg.loc[ix, col] = val_func(vals.loc[this_ix.index[this_ix]])

    return ctg


if __name__ == '__main__':
    alpha = 0.05
    outdir = unique_output_dir("wang_classification")
    n_perm = 1000

    ## TODO: tidy up and use the new function

    # prepare data
    # gct_files = prepare_gct_files_salmon(outdir=outdir)
    gct_files = prepare_gct_files_star(outdir=outdir)

    outdir = '/home/gabriel/python_outputs/gct_files_for_wang.0/'

    for typ in ('cell_culture', 'ffpe', 'all'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            gct_files.append(fn)
            gsea.wang_ssgsea_classification(fn)

    # run the algorithm on them
    p_res = []
    ss_res = []

    for fn in gct_files:
        # gsea.wang_ssgsea_classification(fn, n_perm)

        the_dir, the_stem = os.path.split(fn)
        outfn = os.path.join(the_dir, "p_result_%s.txt" % the_stem)
        if not os.path.exists(outfn):
            continue
        p_res.append(load_pvalue_results(outfn))
        ss_res.append(simplicity_score(p_res[-1]))

    p_res = {}
    ss_res = {}

    for typ in ('cell_culture', 'ffpe', 'all'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            the_dir, the_stem = os.path.split(fn)
            outfn = os.path.join(the_dir, "p_result_%s.txt" % the_stem)
            if not os.path.exists(outfn):
                continue
            this_pres = load_pvalue_results(outfn)
            p_res.setdefault(typ, {})[SRC_MAP[src]] = this_pres
            ss_res.setdefault(typ, {})[SRC_MAP[src]] = simplicity_score(this_pres)


    # indir = os.path.join(OUTPUT_DIR, 'wang_classification')
    # fn = os.path.join(indir, 'p_result_gbm_cc_and_ffpe_fpkm.gct.txt')
    # pvals_ours = load_pvalue_results(fn)
    # ss_ours = simplicity_score(pvals_ours)
    # nm_ours = (pvals_ours < alpha).sum(axis=1)
    # cls_ours = pd.Series(index=pvals_ours.index)
    # min_idx = np.argmin(pvals_ours.values, axis=1)
    # cls_ours.loc[nm_ours == 1] = pvals_ours.columns[min_idx[nm_ours == 1]]
    # cls_ours.loc[nm_ours > 1] = 'Multi'
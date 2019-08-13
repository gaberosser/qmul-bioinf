import os

import numpy as np
import pandas as pd

from hgic_consts import NH_ID_TO_PATIENT_ID_MAP
from rnaseq import gsea, loader
from scripts.hgic_final import consts
from utils import output, dictionary, excel, reference_genomes

SRC_MAP = {
    'star': 'counts',
    'salmon': 'tpm',
    'star/cufflinks': 'rpkm'
}


def nh_id_to_patient_id(arr):
    nh_id = pd.Index(arr).str.replace(r'(_?)(DEF|SP).*', '')
    return [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]


def prepare_gct_files_hgic(pids=consts.ALL_PIDS, outdir=None):
    """
    Prepare the GCT files required to perform classification of the hGIC samples:
    - hGIC FFPE
    - hGIC cell culture
    - Both combined
    In all cases, use FPKM units (cufflinks), TPM (salmon) and CPM (STAR).
    Use gene symbols as these are contained in the signatures.
    """
    if outdir is None:
        outdir = output.unique_output_dir()

    infiles = []

    loaded = {}
    for typ in ('cell_culture', 'ffpe'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            this_obj = loader.load_by_patient(pids, type=typ, source=src, include_control=False)
            this_obj.filter_samples(this_obj.meta.type == 'GBM')
            if typ == 'ffpe':
                # restrict to the 'best' versions (there are some duplicates where we tried twice)
                this_obj.filter_by_sample_name(consts.FFPE_RNASEQ_SAMPLES_ALL)
            this_dat = reference_genomes.translate_quantification_resolving_duplicates(
                this_obj.data,
                'Ensembl Gene ID',
                'Approved Symbol'
            )
            loaded.setdefault(typ, {})[src] = this_dat
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            gsea.data_to_gct(this_dat, fn)
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
    outdir = output.unique_output_dir("wang_classification")
    n_perm = 1000

    ## TODO: tidy up and use the new function

    # prepare data
    gct_files = prepare_gct_files_hgic(outdir=outdir)

    for typ in ('cell_culture', 'ffpe'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            gct_files.append(fn)
            gsea.wang_ssgsea_classification(fn, n_perm=n_perm)

    p_res = {}
    ss_res = {}

    for typ in ('cell_culture', 'ffpe'):
        for src in ('star', 'salmon', 'star/cufflinks'):
            fn = os.path.join(outdir, "%s_%s.gct" % (SRC_MAP[src], typ))
            the_dir, the_stem = os.path.split(fn)
            outfn = os.path.join(the_dir, "p_result_%s.txt" % the_stem)
            if not os.path.exists(outfn):
                continue
            this_pres = load_pvalue_results(outfn)
            p_res.setdefault(typ, {})[SRC_MAP[src]] = this_pres
            ss_res.setdefault(typ, {})[SRC_MAP[src]] = simplicity_score(this_pres)

    # export
    # easiest way is to flatten the dictionary, then combine
    export_p = dictionary.nested_dict_to_flat(p_res)
    export_ss = dictionary.nested_dict_to_flat(ss_res)
    to_export = {}

    for k in export_p:
        the_key = '_'.join(k)
        this = export_p[k].copy()
        this.insert(this.shape[1], 'Simplicity score', export_ss[k])
        if k[0] == 'ffpe':
            this.insert(this.shape[1], 'Patient ID', nh_id_to_patient_id(this.index))
        to_export[the_key] = this
    excel.pandas_to_excel(to_export, os.path.join(outdir, "wang_results.xlsx"))


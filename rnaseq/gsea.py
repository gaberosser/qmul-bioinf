import csv
import os
import re
import pandas as pd
from utils.log import get_console_logger
from utils import rinterface
import multiprocessing as mp
logger = get_console_logger(__name__)


def load_from_gct(infile):
    dat = pd.read_csv(infile, header=0, index_col=0, skiprows=2, delimiter='\t')
    dat = dat.drop('Description', axis=1)
    return dat


def data_to_gct(data, outfile):
    nrow, ncol = data.shape
    with open(outfile, 'wb') as f:
        c = csv.writer(f, delimiter='\t')
        # three header columns
        c.writerow(["#1.2"])
        c.writerow([nrow, ncol])
        c.writerow(['NAME', 'Description'] + data.columns.tolist())
        for name, vals in data.iterrows():
            c.writerow([name, "NA"] + vals.values.tolist())


def combine_gct_files(*infiles):
    all_dat = [load_from_gct(fn) for fn in infiles]
    idx = reduce(lambda x, y: x.intersection(y), [t.index for t in all_dat])
    n_lost_max = max([len(t.index.difference(idx)) for t in all_dat])
    if n_lost_max > 0:
        logger.warn("%d rows were discarded as they are not present in all samples.", n_lost_max)
    dat = pd.concat(all_dat, axis=1).dropna()
    return dat


def phenotypes_to_cls(groups, outfile):
    """

    :param groups: Iterable giving a group name (string) for each sample.
    :param outfile:
    :return:
    """
    nsample = len(groups)
    classes = sorted(set(groups))
    ncls = len(classes)
    # factorise the classes
    cls_map = dict([(t, i) for i, t in enumerate(classes)])
    with open(outfile, 'wb') as f:
        c = csv.writer(f, delimiter=' ')
        c.writerow([nsample, ncls, 1])
        c.writerow(['#'] + classes)
        c.writerow([cls_map[t] for t in groups])


def read_gmt_file(fn):
    """
    Parse a GMT file (from MSIGDB)
    :param fn:
    :return:
    """
    res = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            res[row[0]] = row[2:]
    return res


def create_gsea_params_file(outfn, rpt_label='foo', permute='gene_set', **kwargs):
    """
    Write a params file that can be used when running GSEA with the -param_file input arg.
    We apply a set of defaults - any supplied named / kw args will overwrite this.
    :param outfn:
    :param rpt_label: The label used as a prefix to the report files
    :param permute: The type of permutation:
    phenotype: permute across phenotypes
    gene_set: permute across genes [default]. This should be used for small numbers of samples.
    :param kwargs:
    :return:
    """
    defaults = {
        'collapse': 'false',
        'norm': 'meandiv',
        'nperm': 1000,
        'rnd_type': 'no_balance',
        'scoring_scheme': 'weighted',
        'rpt_label': 'RTK_II',
        'metric': 'Signal2Noise',
        'sort': 'real',
        'order': 'descending',
        'create_gcts': 'false',
        'create_svgs': 'false',
        'include_only_symbols': 'true',
        'make_sets': 'true',
        'median': 'false',
        'num': 1000,
        'plot_top_x': 20,
        'rnd_seed': 'timestamp',
        'save_rnd_lists': 'false',
        'set_max': 500,
        'set_min': 15,
        'zip_report': 'false',
        'gui': 'false'
    }
    defaults['rpt_label'] = rpt_label
    defaults['permute'] = permute
    defaults.update(kwargs)

    with open(outfn, 'wb') as f:
        for x in defaults.items():
            f.write("\t".join([str(t) for t in x]) + "\n")


def run_one_ssgsea(sample_data, gs, alpha=0.25, norm_by_gene_count=True, return_ecdf=False):
    # FIXME: CRUCIAL: which order do we rank in?
    # rank the sample in ascending order (1 corresponds to lowest expression)
    rs = sample_data.rank(method='average', ascending=True)  # ties resolved by averaging

    # sort in decreasing order
    # the most expressed genes come first in the list
    rs = rs.sort_values(ascending=False)

    # boolean vector for inclusion in gene set
    in_gene_set = rs.index.isin(gs).astype(float)
    out_gene_set = (~rs.index.isin(gs)).astype(float)

    # ECDF
    x_in = (rs * in_gene_set) ** alpha
    ecdf_in = (x_in.cumsum() / x_in.sum()).values

    # the ECDF for samples out is NOT weighted, which is strange
    ecdf_out = out_gene_set.cumsum() / out_gene_set.sum()

    # if we were to weight it, it would look like:
    # x_out = (rs * out_gene_set) ** alpha
    # ecdf_out = x_out.cumsum() / x_out.sum()

    # enrichment score is the difference in the integrals
    es = (ecdf_in - ecdf_out).sum()

    if norm_by_gene_count:
        es /= float(rs.shape[0])

    if return_ecdf:
        return (es, ecdf_in, ecdf_out)
    else:
        return es



def ssgsea(sample_data, gene_set, alpha=0.25, norm_by_gene_count=True, return_ecdf=False, threads=None):
    """
    Run single sample gene set enrichment analysis (ssGSEA) on the supplied data, following the details given in:

    Barbie, D.A., Tamayo, P., Boehm, J.S., Kim, S.Y., Moody, S.E., Dunn, I.F., Schinzel, A.C., Sandy, P., Meylan, E.,
    Scholl, C., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1.
    Nature 462, 108-112.

    See R/stats/ssgsea.R for a further example (found online)

    :param sample_data: Pandas Series
    :param gene_set: Dictionary. Each entry has a key giving the name of the gene set and value giving a list of genes.
    :param alpha: The weighting used to compute the ECDF. Default is 0.25 following Barbie et al. (2009)
    :param return_ecdf: If True, also return the two ECDFs being considered. Useful for plotting?
    :return:
    """
    b_one_sample = False
    if not isinstance(gene_set, dict):
        b_one_sample = True
        gene_set = {'result': gene_set}

    n = len(gene_set)

    pool = None
    if threads != 1 and n > 1:
        pool = mp.Pool(processes=threads)

    res = {}
    jobs = {}
    ssgsea_kwds = dict(
        alpha=alpha,
        norm_by_gene_count=norm_by_gene_count,
        return_ecdf=return_ecdf
    )

    for name, gs in gene_set.items():
        if pool is None:
            res[name] = run_one_ssgsea(
                sample_data,
                gs,
                **ssgsea_kwds
            )
        else:
            jobs[name] = pool.apply_async(
                run_one_ssgsea,
                args=(sample_data, gs),
                kwds=ssgsea_kwds
            )

    if pool is not None:
        pool.close()
        pool.join()
        for name, j in jobs.items():
            res[name] = j.get()

    if b_one_sample:
        return res['result']

    return res


    # # FIXME: CRUCIAL: which order do we rank in?
    # # rank the sample in ascending order (1 corresponds to lowest expression)
    # rs = sample_data.rank(method='average', ascending=True)  # ties resolved by averaging
    #
    # # sort in decreasing order
    # # the most expressed genes come first in the list
    # rs = rs.sort_values(ascending=False)
    #
    # # boolean vector for inclusion in gene set
    # in_gene_set = rs.index.isin(gene_set).astype(float)
    # out_gene_set = (~rs.index.isin(gene_set)).astype(float)
    #
    # # ECDF
    # x_in = (rs * in_gene_set) ** alpha
    # ecdf_in = (x_in.cumsum() / x_in.sum()).values
    #
    # # the ECDF for samples out is NOT weighted, which is strange
    # ecdf_out = out_gene_set.cumsum() / out_gene_set.sum()
    #
    # # if we were to weight it, it would look like:
    # # x_out = (rs * out_gene_set) ** alpha
    # # ecdf_out = x_out.cumsum() / x_out.sum()
    #
    # # enrichment score is the difference in the integrals
    # es = (ecdf_in - ecdf_out).sum()
    #
    # if norm_by_gene_count:
    #     es /= float(rs.shape[0])
    #
    # if return_ecdf:
    #     return (es, ecdf_in, ecdf_out)
    # else:
    #     return es


try:
    from rpy2.robjects import r
except ImportError:
    pass

def _wang_ssgsea_classification(gct_file, n_perm=1000):
    r("runSsGSEAwithPermutation")(gct_file, n_perm)


wang_ssgsea_classification = rinterface.RFunctionDeferred(_wang_ssgsea_classification, imports=['ssgsea.GBM.classification'])


def load_gsea_report_and_pathways(subdir, comparison='GBM', fdr=None, load_top_n_pathways=20):
    """
    Load a single GSEA report and (optionally) associated pathway details.
    :param subdir: Top level subdirectory. We search exhaustively within this to identify candidate GSEA reports.
    If >1 is found, a warning is issued and one is chosen arbitrarily.
    :param comparison: String giving the name of the class of interest.
    :param fdr: If supplied, we use this as a cutoff.
    :param load_top_n_pathways: The number of pathways to search for when loading details. To disable pathway loading,
    set this to None.
    :return:
    """
    # dicover the report file to load
    infiles = []
    patt = re.compile("gsea_report_for_%s.*\.xls" % comparison)
    for the_dir, _, the_files in os.walk(subdir):
        for f in the_files:
            if re.search(patt, f):
                infiles.append(os.path.join(the_dir, f))

    if len(infiles) == 0:
        raise AttributeError("No report files found in subdir %s" % subdir)
    elif len(infiles) > 1:
        logger.warning(
            "Found %d report files within subdirectory %s. Choosing one arbitrarily. "
            "This may NOT be the desired behaviour!",
            len(infiles),
            subdir
        )

    infile = infiles[0]
    indir = os.path.split(infile)[0]

    report = pd.read_csv(infile, sep='\t', header=0, index_col=None, usecols=[0, 3, 4, 5, 6, 7, 8])
    report.columns = [
        'pathway',
        'n_gene',
        'es',
        'nes',
        'nom_pval',
        'fdr',
        'fwer',
    ]

    # before filtering, we need to record the pathways
    pathways_to_load = report.pathway.copy()

    if fdr is not None:
        report = report.loc[report.fdr <= fdr]

    if load_top_n_pathways is not None:
        # load pathway details from separate files
        # in case we have filtered previously, run an intersection
        pathways_to_load = set(pathways_to_load.values[:load_top_n_pathways]).intersection(report.pathway.values)
        pathways = {}
        for p in pathways_to_load:
            infile = os.path.join(indir, "%s.xls" % p)
            if not os.path.isfile(infile):
                logger.error("Unable to find expected file %s. Skipping." % infile)
            else:
                the_path = pd.read_csv(infile, sep='\t', header=0, index_col=None, usecols=[1, 5, 6, 7, 8])
                the_path.columns = [
                    'name',
                    'rank_in_gene_list',
                    'rank_metric',
                    'running_es',
                    'core_enrichment'
                ]
                the_path['core_enrichment'] = (the_path['core_enrichment'] == 'Yes')
                pathways[p] = the_path
        return report, pathways
    else:
        return report
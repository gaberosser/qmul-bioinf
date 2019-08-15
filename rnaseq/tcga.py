import re

from utils import reference_genomes
import os
from urlparse import urljoin

import requests
import json
import pandas as pd

from utils.log import get_file_logger
from utils.output import unique_output_dir


API_ROOT = 'https://gdc-api.nci.nih.gov/'
LEGACY_API_ROOT = 'https://gdc-api.nci.nih.gov/legacy/'

FILES_ENDPOINT = urljoin(API_ROOT, 'files/')
LEGACY_FILES_ENDPOINT = urljoin(LEGACY_API_ROOT, 'files/')

DATA_ENDPOINT = urljoin(API_ROOT, 'data/')
LEGACY_DATA_ENDPOINT = urljoin(LEGACY_API_ROOT, 'data/')

CASE_ENDPOINT = urljoin(API_ROOT, 'cases/')

FILE_FIELDS = (
    'file_id',
    'file_name',
    'data_type',
    'cases.case_id',
    'cases.submitter_id'
)

def equal_query(field, value):
    return {
        "op": "=",
        "content": {
            "field": field,
            "value": value,
        }
    }


def in_query(field, arr):
    arr = list(arr)
    return {
        "op": "in",
        "content": {
            "field": field,
            "value": arr
        }
    }

def and_query(*args):
    qry = {
        "op": "and",
        "content": list(args)
    }
    return qry


def or_query(*args):
    qry = {
        "op": "or",
        "content": list(args)
    }
    return qry


qry_meth450 = equal_query("files.platform", "Illumina Human Methylation 450")
qry_trans = equal_query("files.data_type", "Gene Expression Quantification")
qry_gbm = equal_query("cases.project.project_id", "TCGA-GBM")
qry_primary = equal_query("cases.samples.sample_type", "Primary Tumor")


def download_data(file_id, outfile, legacy=False, create_dirs=True):
    endpoint = LEGACY_DATA_ENDPOINT if legacy else DATA_ENDPOINT
    response = requests.get(urljoin(endpoint, file_id))
    outdir = os.path.split(outfile)[0]

    if not os.path.isdir(outdir):
        if create_dirs:
            os.makedirs(outdir)
        else:
            raise AttributeError("Directory %s does not exist" % outdir)

    with open(outfile, 'wb') as fout:
        for blk in response.iter_content(1024):
            fout.write(blk)


def download_from_manifest(path_to_manifest, outdir=None, legacy=False):
    """
    Download all files from the provided manifest
    :param path_to_manifest:
    :param outdir: If None, create a unique output folder
    :return:
    """
    if outdir is None:
        outdir = unique_output_dir("nih_gdc_legacy")
    mani = pd.read_csv(path_to_manifest, sep='\t', header=0, index_col=0)
    for fid, row in mani.iterrows():
        outfile = os.path.join(outdir, row.filename)
        download_data(fid, outfile, legacy=legacy)


def get_meth450_case_ids(project='TCGA-GBM', sample_type='Primary Tumor'):
    """
    Get GBM case IDs associated with Illumina's Methylation 450k array.
    """
    queries = [equal_query("files.platform", "Illumina Human Methylation 450")]
    if project is not None:
        queries.append(equal_query("cases.project.project_id", project))
    if sample_type is not None:
        queries.append(equal_query("cases.samples.sample_type", sample_type))
    qry = {
        "filters": and_query(*queries),
        "format": "json",
        "fields": "case_id",
        "size": 10000,
    }

    response = requests.post(CASE_ENDPOINT, json=qry)
    return [t['case_id'] for t in response.json()['data']['hits']]


def get_rnaseq_case_ids(project='TCGA-GBM', sample_type='Primary Tumor'):
    queries = [equal_query("files.data_type", "Gene Expression Quantification")]
    if project is not None:
        queries.append(equal_query("cases.project.project_id", project))
    if sample_type is not None:
        queries.append(equal_query("cases.samples.sample_type", sample_type))
    qry = {
        "filters": and_query(*queries),
        "format": "json",
        "fields": "case_id",
        "size": 10000,
    }
    response = requests.post(CASE_ENDPOINT, json=qry)
    return [t['case_id'] for t in response.json()['data']['hits']]


def get_case_ids_with_paired_data():

    # Query 1: cases in TCGA-GBM with methylation data
    case_ids_meth = get_meth450_case_ids()

    # Query 2: cases in TCGA-GBM with transcriptome data
    case_ids_trans = get_rnaseq_case_ids()

    # find the intersecting case IDs

    return set(case_ids_meth).intersection(case_ids_trans)


def get_legacy_idat(case_ids):
    outdir = unique_output_dir("gdc-nih_methylation", reuse_empty=True)
    logger = get_file_logger("legacy_idats", os.path.join(outdir, "getter.log"))

    qry_case = in_query("cases.case_id", case_ids)

    qry_idat = equal_query("files.data_format", "idat")

    qry = {
        "filters": and_query(qry_primary, qry_case, qry_idat, qry_meth450),
        "format": "json",
        "fields": ','.join(FILE_FIELDS),
        "size": 10000
    }

    response = requests.post(LEGACY_FILES_ENDPOINT, json=qry)
    if response.status_code != 200:
        logger.error("Initial query failed: %s", response.content)
        raise ValueError("Query failed")

    res = response.json()['data']['hits']
    logger.info("Found %d idat files.", len(res))

    num_error = 0
    num_files = 0

    # we need to keep track of the files in order to write meta correctly
    meta = {}

    for r in res:
        if len(r['cases']) > 1:
            logger.error("File with ID %s has multiple case ID matches", r['file_id'])
        cid = r['cases'][0]['case_id']
        fid = r['file_id']
        fname = r['file_name']
        outfn = os.path.join(outdir, cid, fname)
        meta.setdefault(cid, [])
        logger.info("Case %s. File ID %s. Output path %s.", cid, fid, outfn)
        try:
            download_data(fid, outfn, legacy=True)
        except Exception:
            logger.exception("Failed to download %s for case id %s", fname, cid)
            num_error += 1
        else:
            logger.info("Downloaded case ID %s file ID %s to %s", cid, fid, outfn)
            meta[cid].append(r)
            num_files += 1

    logger.info("Downloaded %d files. Encountered %d errors.", num_files, num_error)

    num_meta = 0
    num_meta_errors = 0

    # write meta files
    for cid, arr in meta.iteritems():
        meta_fn = os.path.join(outdir, cid, 'meta.json')
        if os.path.exists(meta_fn):
            logger.error("Meta file already exists: %s", meta_fn)
            num_meta_errors += 1
        else:
            with open(meta_fn, 'wb') as f:
                json.dump(arr, f)
            num_meta += 1
    logger.info("Create %d meta files. Encountered %d errors.", num_meta, num_meta_errors)


def get_methylation_gene_expression_data(case_ids):
    outdir = unique_output_dir("gdc-nih_gene_expr", reuse_empty=True)
    logger = get_file_logger("nih_methylation_gene_counts", os.path.join(outdir, "getter.log"))

    # get relevant files for download

    qry_case = {
        "op": "in",
        "content": {"field": "cases.case_id", "value": list(case_ids)}
    }

    qry = {
        "filters": and_query(qry_primary, qry_case, or_query(qry_trans, qry_meth450)),
        "format": "json",
        "fields": ','.join(FILE_FIELDS),
        "size": 10000
    }
    response = requests.post(FILES_ENDPOINT, json=qry)
    if response.status_code != 200:
        logger.error("Initial query failed: %s", response.content)
        raise ValueError("Query failed")
    res = response.json()['data']['hits']

    meta = {}
    num_error = 0
    num_files = 0

    for r in res:
        if len(r['cases']) > 1:
            logger.error("File with ID %s has multiple case ID matches", r['file_id'])
        cid = r['cases'][0]['case_id']
        fid = r['file_id']
        fname = r['file_name']
        meta.setdefault(cid, {})
        if r['data_type'] == 'Gene Expression Quantification':
            if 'FPKM-UQ' in fname:
                continue
            elif 'FPKM' in fname:
                meta[cid]['fpkm'] = r
                outfn = os.path.join(outdir, cid, 'fpkm.gz')
            elif 'htseq.counts' in fname:
                meta[cid]['counts'] = r
                outfn = os.path.join(outdir, cid, 'counts.gz')
        elif r['data_type'] == 'Methylation Beta Value':
            meta[cid]['methylation'] = r
            outfn = os.path.join(outdir, cid, 'methylation.txt')
        try:
            download_data(fid, outfn)
        except Exception:
            logger.exception("Failed to download %s for case id %s", fname, cid)
            num_error += 1
        else:
            logger.info("Downloaded case ID %s file ID %s to %s", cid, fid, outfn)
            num_files += 1

    logger.info("Downloaded %d files. Encountered %d errors.", num_files, num_error)

    # run back through and write meta files

    num_meta = 0
    num_meta_errors = 0

    # write meta files
    for cid, d in meta.iteritems():
        meta_fn = os.path.join(outdir, cid, 'meta.json')
        if os.path.exists(meta_fn):
            logger.error("Meta file already exists: %s", meta_fn)
            num_meta_errors += 1
        else:
            with open(meta_fn, 'wb') as f:
                json.dump(d, f)
            num_meta += 1
    logger.info("Create %d meta files. Encountered %d errors.", num_meta, num_meta_errors)


def prepare_data_and_meta(indir, meta_fn=None):
    """
    :param meta_fn: Path to the Brennan Table S7 data, which gives expression subclasses. If None, assume it is
    one level *above* the indir.
    """
    if meta_fn is None:
        meta_fn = os.path.join(indir, '..', 'brennan_s7.csv')
    dirlist = os.listdir(indir)

    df = reference_genomes.conversion_table()
    df = df.loc[~df.loc[:, 'Ensembl Gene ID'].duplicated()].set_index('Ensembl Gene ID')

    dat = None  # this will hold all our data
    first = True

    for d in dirlist:
        if re.match(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}', d):
            the_dir = os.path.join(indir, d)
            the_file = os.listdir(the_dir)[0]  # assuming only one file per directory
            # load data and convert to a Series
            the_col = pd.read_csv(os.path.join(the_dir, the_file), sep='\t', index_col=0, header=None).iloc[:, 0]
            # only keep genes
            the_col = the_col.loc[the_col.index.str.contains('ENSG0')]
            # remove version numbers from the row names
            the_col.index = the_col.index.str.replace(r'\.[0-9]*', '')
            # rename
            col_name = the_file
            the_col.name = col_name

            # if this is the first file loaded, use it as a template
            if first:
                dat = pd.DataFrame(index=the_col.index)
                dat.loc[:, col_name] = the_col # the directory is the case ID
                first = False
            else:
                dat.loc[:, col_name] = the_col  # the directory is the case ID

    # now we can use the annotation provided by Brennan et al. and the manifest to add subgroup information
    # this is the old WHO subgrouping
    # we do this via the
    meta = pd.read_csv(meta_fn, index_col=0, header=0)

    fn = os.path.join(
        indir,
        'cases.json'
    )
    with open(fn, 'rb') as f:
        cases = json.load(f)
    cases_mani = pd.Series(
        [t['submitter_id'] for t in cases],
        index=[t['case_id'] for t in cases]
    )

    fn = os.path.join(
        indir,
        'files.json'
    )
    with open(fn, 'rb') as f:
        files = json.load(f)
    files_mani = pd.Series(
        [t['cases'][0]['case_id'] for t in files],
        index=[t['file_name'] for t in files],
    )

    # link filename -> case ID
    aa = files_mani.loc[dat.columns[:-1]]

    # I see duplicates here (2) - WHY? Probably the same tumour split into multiple samples?

    # report them - show all cases
    print "Duplicates detected: \n%s" % aa.loc[pd.Index(aa.values).duplicated(keep=False)].to_string()

    # pick one in each case
    # duplicated() doesn't mark the first entry by default
    aa = aa[~pd.Index(aa.values).duplicated()]

    dat = dat.loc[:, aa.index.intersection(dat.columns)]
    cols = aa.loc[dat.columns].values

    bb = cases_mani.loc[cols]
    # these should be unique...
    if pd.Index(bb.values).duplicated().any():
        raise AttributeError("Some case IDs are duplicated in the final dataset")
    dat.columns = bb.values

    # finally, only keep those cases that are also in the Brennan table
    not_in_meta = pd.Index(bb.values).difference(meta.index)

    if len(not_in_meta):
        this = bb.loc[pd.Index(bb.values).isin(not_in_meta)]
        print "Cases not in meta: \n%s" % this.to_string()

    meta = meta.loc[meta.index.intersection(bb.values)]
    dat = dat.loc[:, meta.index]

    return dat, meta

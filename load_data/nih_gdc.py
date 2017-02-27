import os
from urlparse import urljoin

import requests

from utils.log import get_file_logger
from utils.output import unique_output_dir

API_ROOT = 'https://gdc-api.nci.nih.gov/'
FILES_ENDPOINT = urljoin(API_ROOT, 'files/')
DATA_ENDPOINT = urljoin(API_ROOT, 'data/')
CASE_ENDPOINT = urljoin(API_ROOT, 'cases/')


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


def get_paired_methylation_gene_expression_data():
    outdir = unique_output_dir("gdc-nih_paired_methylation_gene_expr", reuse_empty=True)
    logger = get_file_logger("paired_methylation_gene_counts", os.path.join(outdir, "getter.log"))

    qry_meth450 = {
        "op": "=",
        "content": {
            "field": "files.platform",
            "value": "Illumina Human Methylation 450",
        }
    }

    qry_trans = {
        "op": "=",
        "content": {
            "field": "files.data_type",
            "value": "Gene Expression Quantification",
        }
    }

    qry_gbm = {
        "op": "=",
        "content": {
            "field": "cases.project.project_id",
            "value": "TCGA-GBM",
        }
    }

    qry_primary = {
        "op": "=",
        "content": {
            "field": "cases.samples.sample_type",
            "value": "Primary Tumor",
        }
    }

    # Query 1: cases in TCGA-GBM with methylation data

    qry = {
        "filters": and_query(qry_gbm, qry_meth450, qry_primary),
        "format": "json",
        "fields": "case_id",
        "size": 10000,
    }

    response = requests.post(CASE_ENDPOINT, json=qry)
    case_ids_meth = [t['case_id'] for t in response.json()['data']['hits']]

    # Query 2: cases in TCGA-GBM with transcriptome data

    qry = {
        "filters": and_query(qry_gbm, qry_trans, qry_primary),
        "format": "json",
        "fields": "case_id",
        "size": 10000,
    }

    response = requests.post(CASE_ENDPOINT, json=qry)
    case_ids_trans = [t['case_id'] for t in response.json()['data']['hits']]

    # find the intersecting case IDs

    case_ids = set(case_ids_meth).intersection(case_ids_trans)

    # get relevant files for download

    qry_case = {
        "op": "in",
        "content": {"field": "cases.case_id", "value": list(case_ids)}
    }

    qry = {
        "filters": and_query(qry_primary, qry_case, or_query(qry_trans, qry_meth450)),
        "format": "json",
        "fields": "file_id,file_name,data_type,cases.case_id",
        "size": 10000
    }
    response = requests.post(FILES_ENDPOINT, json=qry)

    # get paired data
    flist = {}
    for t in response.json()['data']['hits']:
        cid = t['cases'][0]['case_id']
        if t['data_type'] == 'Gene Expression Quantification':
            if 'FPKM-UQ' in t['file_name']:
                continue
            elif 'FPKM' in t['file_name']:
                flist.setdefault(cid, {})['fpkm'] = t['file_id']
            elif 'htseq.counts' in t['file_name']:
                flist.setdefault(cid, {})['counts'] = t['file_id']
        elif t['data_type'] == 'Methylation Beta Value':
            flist.setdefault(cid, {})['methylation'] = t['file_id']

    num_cases = 0
    num_files = 0
    num_error = 0
    for k, v in flist.items():
        if len(v) != 3:
            continue
        out_subdir = os.path.join(outdir, k)
        if not os.path.exists(out_subdir):
            os.makedirs(out_subdir)
        logger.info("Case %s. Output dir %s.", k, out_subdir)

        for fn, fid in v.items():
            ff = None
            if fn == 'counts':
                ff = os.path.join(out_subdir, fn) + '.gz'
            elif fn == 'fpkm':
                ff = os.path.join(out_subdir, fn) + '.gz'
            elif fn == 'methylation':
                ff = os.path.join(out_subdir, fn) + '.txt'
            else:
                logger.error("Unsupported data type %s", fn)
                raise ValueError("Unsupported data type %s" % fn)

            try:
                response = requests.get(urljoin(DATA_ENDPOINT, fid))
                with open(ff, 'wb') as fout:
                    for blk in response.iter_content(1024):
                        fout.write(blk)
            except Exception:
                logger.exception("Failed to download %s for case id %s", fn, k)
                num_error += 1
            else:
                logger.info("Downloaded %s data to %s", fn, ff)
                num_files += 1

    logger.info("Downloaded %d files. Encountered %d errors.", num_files, num_error)

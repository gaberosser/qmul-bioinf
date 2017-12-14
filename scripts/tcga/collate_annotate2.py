"""
Created on 13th Dec 2017
Aim: a generic approach to pull together the data that we download with the GDC client
"""
from glob import glob
import json
import os
import re
import pandas as pd
from utils.output import unique_output_dir
from settings import DATA_DIR_NON_GIT, DATA_DIR_NON_GIT2


file_types = [
    (r'\.htseq\.counts$', 'htseq'),
    (r'\.FPKM\.txt$', 'fpkm'),
    (r'\.FPKM-UQ\.txt$', 'fpkm-uq'),
    (r'HumanMethylation27', 'methylation_27k'),
    (r'HumanMethylation450', 'methylation_450k'),
]


def get_filetype(fname):
    for patt, ftype in file_types:
        if re.search(patt, fname):
            return ftype
    return 'unknown'


def get_meta(indir):
    meta_candidates = glob(os.path.join(indir, "metadata*.json"))
    if len(meta_candidates) == 1:
        meta_fn = meta_candidates[0]
        return meta_fn
    else:
        raise AttributeError("meta_fn not supplied and no single candidate found.")


def load_one_rnaseq(indir, fid, fname, col_name=None):
    """
    Load a single RNASeq result, formatting index correctly
    :param indir:
    :param fid:
    :param fname:
    :return: pd.Series
    """
    fn = os.path.join(indir, fid, fname)
    the_col = pd.read_csv(fn, sep='\t', index_col=0, header=None).iloc[:, 0]
    # only keep genes
    the_col = the_col.loc[the_col.index.str.contains('ENSG0')]
    # remove version numbers from the row names
    the_col.index = the_col.index.str.replace(r'\.[0-9]*', '')
    the_col.index.name = 'ensembl_gene_id'
    # rename
    if col_name is not None:
        the_col.name = col_name
    return the_col


def load_one_methylation(indir, fid, fname, col_name=None):
    fn = os.path.join(indir, fid, fname)
    the_col = pd.read_csv(fn, sep='\t', index_col=0, header=None).iloc[:, 0]
    # only keep genes
    the_col.index.name = 'probe_id'
    # rename
    if col_name is not None:
        the_col.name = col_name
    return the_col


def load_files(indir, loader, meta_fn=None, sample_types=('Primary Tumor',)):
    sample_types = set(sample_types)

    if meta_fn is None:
        # try to find an unambiguous meta file
        meta_fn = get_meta(indir)
        print "Using metadata in file %s" % meta_fn

    with open(meta_fn, 'rb') as f:
        meta = json.load(f)

    dat = {}
    info = {}

    for m in meta:
        fid = m['file_id']
        fname = m['file_name']

        ftype = get_filetype(fname)

        if len(m['cases']) > 1:
            raise AttributeError("Unexpectedly encountered >1 cases attributed to a single file: %s." % repr(m))
        c = m['cases'][0]
        cid = c['submitter_id']

        if len(c['samples']) > 1:
            print "Warning: %d samples attributed to case %s." % (len(c['samples']), cid)
            raise Exception

        for s in c['samples']:
            if s['sample_type'] in sample_types:
                dat.setdefault(ftype, None)
                if len(s['portions']) > 1:
                    print "Warning: %d samples attributed to case %s." % (len(s['portions']), cid)
                    raise Exception

                for p in s['portions']:
                    pid = p['submitter_id']
                    the_col = loader(indir, fid, fname, col_name=pid)

                    # if this is the first file loaded, use it as a template
                    if dat[ftype] is None:
                        dat[ftype] = pd.DataFrame(index=the_col.index)
                        dat[ftype].loc[:, pid] = the_col  # the directory is the case ID
                    else:
                        dat[ftype].loc[:, pid] = the_col  # the directory is the case ID

                    if pid not in info:
                        info[pid] = pd.Series({
                            'case_id': cid,
                            'file_id': m['file_id'],
                            'file_name': m['file_name'],
                            'sample_type': s['sample_type']
                        })

    info = pd.DataFrame(info).transpose()

    return info, dat


def save_to_disk(dat, filestem):
    # save to disk: CSV and Excel formats
    xl_writer = pd.ExcelWriter("%s.xlsx" % filestem)
    for ftype, t in dat.items():
        fn = os.path.join(outdir, "%s.%s.csv" % (ftype, filestem))
        t.to_csv(fn)
        t.to_excel(xl_writer, ftype)
    xl_writer.save()


if __name__ == "__main__":
    # RNA-Seq
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'raw')
    outdir = unique_output_dir("tcga_gbm_rnaseq", reuse_empty=True)
    rnaseq_info, rnaseq_dat = load_files(indir, loader=load_one_rnaseq)
    fout = os.path.join(outdir, "rnaseq")
    save_to_disk(rnaseq_dat, fout)

    # Methylation (27K)
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour', 'raw_27k')
    meth27_info, meth27_dat = load_files(indir, loader=load_one_methylation)
    meth27_dat = meth27_dat['methylation_27k'].dropna(how='all', axis=0)

    # Methylation (450K)
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour', 'raw_450k')
    meth450_info, meth450_dat = load_files(indir, loader=load_one_methylation)
    meth450_dat = meth450_dat['methylation_450k'].dropna(how='all', axis=0)


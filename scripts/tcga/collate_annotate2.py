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
    (r'\.htseq\.counts', 'htseq'),
    (r'\.FPKM\.txt', 'fpkm'),
    (r'\.FPKM-UQ\.txt', 'fpkm-uq'),
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


def load_one_microarray(indir, fid, fname, col_name=None):
    fn = os.path.join(indir, fid, fname)
    the_col = pd.read_csv(fn, sep='\t', index_col=0, header=0).iloc[1:, 0]
    # only keep genes
    the_col.index.name = 'gene_symbol'
    # rename
    if col_name is not None:
        the_col.name = col_name
    return the_col


def parse_meta(meta_fn):
    with open(meta_fn, 'rb') as f:
        meta = json.load(f)

    res = {}

    for m in meta:
        fid = m['file_id']
        fname = m['file_name']

        if len(m['cases']) > 1:
            raise AttributeError("Unexpectedly encountered >1 cases attributed to a single file: %s." % repr(m))
        c = m['cases'][0]
        cid = c['submitter_id']

        if len(c['samples']) > 1:
            raise AttributeError("%d samples attributed to case %s - expected 1." % (len(c['samples']), cid))

        s = c['samples'][0]

        if len(s['portions']) > 1:
            raise AttributeError("%d portions attributed to sample %s - expected 1." % (len(s['portions']), s))

        p = s['portions'][0]
        pid = p['submitter_id']

        this_a260_280 = []
        for a in p['analytes']:
            this_a260_280.append(a['a260_a280_ratio'])

        this_pct = {}
        for sl in p['slides']:
            for k in ('necrosis', 'normal_cells', 'stromal_cells', 'tumor_cells', 'tumor_nuclei'):
                the_key = "percent_%s" % k
                this_pct.setdefault(the_key, []).append(sl[the_key])

        if pid in res:
            if res[pid]['case_id'] != cid:
                print "case_id clash for %s: (%s, %s)" % (
                    pid, res[pid]['case_id'], pid
                )
            if res[pid]['sample_type'] != s['sample_type']:
                print "sample_type clash for %s: (%s, %s)" % (
                    pid, res[pid]['sample_type'], s['sample_type']
                )
            res[pid]['file_id'].append(fid)
            res[pid]['file_name'].append(fname)
        else:
            res[pid] = {
                'case_id': cid,
                'file_id': [fid],
                'file_name': [fname],
                'sample_type': s['sample_type'],
                'a260_280_ratio': this_a260_280,
            }
            res[pid].update(this_pct)

    return res


def load_files(
        indir,
        loader,
        meta_fn=None,
        sample_types=('Primary Tumor', 'Solid Tissue Normal'),
        filetype=None
):
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

        if filetype is None:
            ftype = get_filetype(fname)
        else:
            ftype = filetype

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


def save_to_disk(dat, filestem, compress=True):
    # save to disk: CSV and Excel formats
    xl_writer = pd.ExcelWriter("%s.xlsx" % filestem)
    for ftype, t in dat.items():
        fn = "%s.%s.csv" % (filestem, ftype)
        if compress:
            fn = "%s.%s.csv.gz" % (filestem, ftype)
            t.to_csv(fn, compression='gzip')
        else:
            fn = "%s.%s.csv" % (filestem, ftype)
            t.to_csv(fn)
        t.to_excel(xl_writer, ftype)
    xl_writer.save()


if __name__ == "__main__":
    outdir = unique_output_dir("tcga_gbm_data", reuse_empty=True)

    # RNA-Seq primary tumour
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'raw')
    rnaseq_info, rnaseq_dat = load_files(indir, loader=load_one_rnaseq)
    fout = os.path.join(outdir, "rnaseq")
    save_to_disk(rnaseq_dat, fout, compress=True)
    fout = os.path.join(outdir, "rnaseq.meta")
    rnaseq_info.to_csv("%s.csv" % fout)
    rnaseq_info.to_excel("%s.xlsx" % fout)

    # RNA-Seq solid tissue normal
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'solid_tissue_normal', 'raw')
    rnaseq_norm_info, rnaseq_norm_dat = load_files(indir, loader=load_one_rnaseq)
    fout = os.path.join(outdir, "rnaseq_normal")
    save_to_disk(rnaseq_norm_dat, fout, compress=True)
    fout = os.path.join(outdir, "rnaseq_normal.meta")
    rnaseq_norm_info.to_csv("%s.csv" % fout)
    rnaseq_norm_info.to_excel("%s.xlsx" % fout)

    # Methylation (27K)
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour', 'raw', '27k')
    meth27_info, meth27_dat = load_files(indir, loader=load_one_methylation)
    meth27_dat = meth27_dat['methylation_27k'].dropna(how='all', axis=0)
    meth27_dat.to_csv(os.path.join(outdir, "methylation.27k.csv.gz"), compression='gzip')
    meth27_dat.to_excel(os.path.join(outdir, "methylation.27k.xlsx"))
    meth27_info.to_csv(os.path.join(outdir, "methylation.27k.meta.csv"))
    meth27_info.to_excel(os.path.join(outdir, "methylation.27k.meta.xlsx"))

    # Methylation (450K)
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour', 'raw', '450k')
    meth450_info, meth450_dat = load_files(indir, loader=load_one_methylation)
    meth450_dat = meth450_dat['methylation_450k'].dropna(how='all', axis=0)
    meth450_dat.to_csv(os.path.join(outdir, "methylation.450k.csv.gz"), compression='gzip')
    # meth450_dat.to_excel(os.path.join(outdir, "methylation.450k.xlsx"))
    meth450_info.to_csv(os.path.join(outdir, "methylation.450k.meta.csv"))
    meth450_info.to_excel(os.path.join(outdir, "methylation.450k.meta.xlsx"))

    # Methylation (450K) solid tissue normal
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'solid_tissue_normal', 'raw', '450k')
    meth450_norm_info, meth450_norm_dat = load_files(indir, loader=load_one_methylation)
    meth450_norm_dat = meth450_norm_dat['methylation_450k'].dropna(how='all', axis=0)
    meth450_norm_dat.to_csv(os.path.join(outdir, "methylation_normal.450k.csv.gz"), compression='gzip')
    meth450_norm_info.to_csv(os.path.join(outdir, "methylation_normal.450k.meta.csv"))
    meth450_norm_info.to_excel(os.path.join(outdir, "methylation_normal.450k.meta.xlsx"))

    # arrays
    array_types = ['agilentg4502a_07_1', 'agilentg4502a_07_2', 'ht_hg_u133a']
    marr_dat = {}
    marr_info = {}
    root_dir = os.path.join(DATA_DIR_NON_GIT2, 'microarray', 'tcga_gbm', 'primary_tumour', 'raw')
    for at in array_types:
        marr_info[at], the_dat = load_files(os.path.join(root_dir, at), loader=load_one_microarray, filetype=at)
        marr_dat[at] = the_dat[at].dropna(how='all', axis=0)
    fout = os.path.join(outdir, 'microarray')
    save_to_disk(marr_dat, fout)
    fout = os.path.join(outdir, 'microarray.meta')
    save_to_disk(marr_info, fout, compress=False)

import os
from settings import DATA_DIR_NON_GIT, LOCAL_DATA_DIR
import csv
from gzip import GzipFile
import re


def get_ref_file(ref):
    if ref == 'GRCh38r90':
        infile = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90', 'gtf',
                              'Homo_sapiens.GRCh38.90.gtf.gz')
    elif ref == 'GRCh38r87':
        infile = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.release87', 'gtf',
                              'Homo_sapiens.GRCh38.87.gtf.gz')
    elif ref == 'GRCm38r88':
        infile = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'mouse', 'ensembl', 'GRCm38.p5.r88', 'gtf',
                              'Mus_musculus.GRCm38.88.gtf.gz')
    else:
        raise NotImplementedError("Unrecognised reference: %s" % ref)
    return infile


def get_rrna(ref='GRCh38r90'):
    infile = get_ref_file(ref)
    res = []
    with GzipFile(infile, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for r in c:
            if len(r) < 9:
                continue
            if 'rRNA' in r[8]:
                t = r[8].split('; ')[0]
                ensg = re.sub(r'gene_id "(?P<ensg>.*)"', r'\g<ensg>', t)
                res.append(ensg)

    # ensure it's unique
    res = list(set(res))
    return res


def get_mitochondrial(ref='GRCh38r90'):
    infile = get_ref_file(ref)
    res = []
    with GzipFile(infile, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for r in c:
            if len(r) < 9:
                continue
            if r[0] == 'MT' and r[2] == 'gene':
                t = r[8].split('; ')[0]
                ensg = re.sub(r'gene_id "(?P<ensg>.*)"', r'\g<ensg>', t)
                res.append(ensg)

    # ensure it's unique
    res = list(set(res))
    return res

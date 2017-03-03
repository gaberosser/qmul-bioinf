import os
from settings import DATA_DIR_NON_GIT
import csv
import re


def get_rrna():
    infile = os.path.join(DATA_DIR_NON_GIT, 'reference_genomes', 'ensembl', 'GRCh38', 'gtf', 'Homo_sapiens.GRCh38.87.gtf')
    res = []
    with open(infile, 'rb') as f:
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


if __name__ == '__main__':
    res = get_rrna()
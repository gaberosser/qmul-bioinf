import os
import pandas as pd
import numpy as np
from settings import DATA_DIR
import csv
import gzip

hu_anno_file = os.path.join(DATA_DIR, 'microarray_annotation', 'HuGene-1_1-st-v1.na36.hg19.transcript.csv.gz')
mo_anno_file = os.path.join(DATA_DIR, 'microarray_annotation', 'MoGene-1_1-st-v1.na36.mm10.transcript.csv.gz')
homolog_file = os.path.join(DATA_DIR, 'homologene', 'homologene.data')
homolog_header = [
    'hid',
    'taxonomy_id',
    'gene_id',
    'gene_symbol',
    'protein_gi',
    'protein_accessor'
]

mouse_tid = 10090
human_tid = 9606

cols_to_keep = [
    'probeset_id',
    'gene_assignment',
]

# Pandas struggles to read such a large CSV, so do it manually
with gzip.open(hu_anno_file, 'rb') as f:
    l = '#'
    while l[0] == '#':
        l = f.readline()
    c = csv.DictReader(f, fieldnames=l.split(','), delimiter=',')
    hu_anno = list(c)

## TODO: add index
hu_anno = pd.DataFrame(hu_anno)
mo_anno = pd.read_csv(mo_anno_file, comment='#', header=0)
homologs = pd.read_csv(homolog_file, header=None, names=homolog_header, sep='\t')

# map human <-> mouse
# get HID
hu_hid = homologs.loc[homologs.taxonomy_id == human_tid, 'hid']
mo_hid = homologs.loc[homologs.taxonomy_id == mouse_tid, 'hid']

# dedupe independently: no orthologs with the same hid
hu_hid = hu_hid[~hu_hid.duplicated()]
mo_hid = mo_hid[~mo_hid.duplicated()]

# reduce to common pool
joint_hid = np.intersect1d(hu_hid.values, mo_hid.values)

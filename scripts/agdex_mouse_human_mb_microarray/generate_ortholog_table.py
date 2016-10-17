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
    fieldnames = l.replace('"', '').replace('\n', '').split(',')
    c = csv.DictReader(f, fieldnames=fieldnames, delimiter=',')
    hu_anno = list(c)

hu_anno = pd.DataFrame(hu_anno)
hu_anno.set_index('probeset_id', inplace=True)
hu_anno.replace(to_replace='---', value=np.nan, inplace=True)

mo_anno = pd.read_csv(mo_anno_file, comment='#', header=0)
mo_anno.set_index('probeset_id', inplace=True)
mo_anno.replace(to_replace='---', value=np.nan, inplace=True)

homologs = pd.read_csv(homolog_file, header=None, names=homolog_header, sep='\t')

# map human <-> mouse
# get HID
hu_hid = homologs.loc[homologs.taxonomy_id == human_tid, 'hid']
mo_hid = homologs.loc[homologs.taxonomy_id == mouse_tid, 'hid']

# dedupe independently: no orthologs with the same hid
# hu_hid = hu_hid[~hu_hid.duplicated()]
# mo_hid = mo_hid[~mo_hid.duplicated()]

# reduce to common pool
joint_hid = np.intersect1d(hu_hid.values, mo_hid.values)

# only keep 1:1 mappings (but quantify this loss of m2m mappings)
hu_mappings = homologs.loc[(homologs.taxonomy_id == human_tid) & (homologs.hid.isin(joint_hid))]
hm_all = hu_mappings.shape[0]
hu_mappings = hu_mappings.loc[~hu_mappings.hid.duplicated()]
hm_uniq = hu_mappings.shape[0]

mo_mappings = homologs.loc[(homologs.taxonomy_id == mouse_tid) & (homologs.hid.isin(joint_hid))]
mo_all = mo_mappings.shape[0]
mo_mappings = mo_mappings.loc[~mo_mappings.hid.duplicated()]
mo_uniq = mo_mappings.shape[0]

# use gene symbols to map between homologene and annotations

def get_gene_symbol(x):
    return x.split(' // ')[1]

hu_anno.loc[~hu_anno.gene_assignment.isnull(), 'gene_symbol'] = \
    hu_anno.gene_assignment.loc[~hu_anno.gene_assignment.isnull()].apply(get_gene_symbol)

mo_anno.loc[~mo_anno.gene_assignment.isnull(), 'gene_symbol'] = \
    mo_anno.gene_assignment.loc[~mo_anno.gene_assignment.isnull()].apply(get_gene_symbol)

# lookup

hu_anno_match = hu_anno.loc[hu_anno.gene_symbol.isin(hu_mappings.gene_symbol.values)]
mo_anno_match = mo_anno.loc[mo_anno.gene_symbol.isin(mo_mappings.gene_symbol.values)]
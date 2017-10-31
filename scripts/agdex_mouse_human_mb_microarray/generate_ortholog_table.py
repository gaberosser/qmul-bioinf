import os
import pandas as pd
import numpy as np
from settings import GIT_LFS_DATA_DIR


homolog_file = os.path.join(GIT_LFS_DATA_DIR, 'homologene', 'homologeneb68.data')
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


def homologs(tid1, tid2, field='gene_symbol'):
    """
    Find matching properties between the two different taxonomies, specified by the IDs. The field tells us which
    attribute column to report.
    Homologene data is used to generate this.
    :param tid1: Taxonomy ID 1.
    :param tid2: Taxonomy ID 2
    :param field: Column name in homologene data.
    :return: DataFrame with two columns, named by taxonomy ID
    """
    homologs = pd.read_csv(homolog_file, header=None, names=homolog_header, sep='\t')
    # get HID
    hid1 = homologs.loc[homologs.taxonomy_id == tid1, 'hid']
    hid2 = homologs.loc[homologs.taxonomy_id == tid2, 'hid']

    joint_hid = np.intersect1d(hid1.values, hid2.values)

    # only keep 1:1 mappings (but quantify this loss of m2m mappings)
    map1 = homologs.loc[(homologs.taxonomy_id == tid1) & (homologs.hid.isin(joint_hid))]
    n1_all = map1.shape[0]
    map1 = map1.loc[~map1.hid.duplicated()]
    map1.set_index('hid', inplace=True)

    map2 = homologs.loc[(homologs.taxonomy_id == tid2) & (homologs.hid.isin(joint_hid))]
    n2_all = map2.shape[0]
    map2 = map2.loc[~map2.hid.duplicated()]
    map2.set_index('hid', inplace=True)

    # map on the field HID
    col_names = ["%s_%d" % (field, tid1), "%s_%d" % (field, tid2)]
    res = pd.DataFrame(columns=col_names, index=joint_hid)
    res.loc[:, col_names[0]] = map1.loc[:, field]
    res.loc[:, col_names[1]] = map2.loc[:, field]

    return res

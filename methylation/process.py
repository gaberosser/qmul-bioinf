import numpy as np
import pandas as pd
import re


def m_from_beta(beta):
    return np.log2(beta / (1 - beta))



def merge_illumina_probe_gene_classes(gene_classes, island_classes):
    """
    From DMRforPairs. Compute the merged probe classes of the Illumina Methylation array.
    :param gene_classes: Series containing the gene classes (UCSC_RefGene_Group in the default annotation file)
    :param island_classes: Series containing the island classes (Relation_to_UCSC_CpG_Island)
    :return: Series giving merged classes separated by semicolon where >1 applies
    """
    mapping = {
        'tss': 'TSS200|TSS1500',
        'gene': "Body|5'UTR|3'UTR|1stExon",
        'island': "Island|N_Shelf|N_Shore|S_Shelf|S_Shore"
    }

    # concatenate the strings in the two series
    c = gene_classes.str.cat(island_classes, sep=';', na_rep='')

    # instantiate new Series
    m = pd.Series(index=c.index, name='probe_group', dtype=str)

    for k, r in mapping.items():
        regex = re.compile(r)
        s = pd.Series(index=c.index, dtype=str)
        s.loc[c.str.contains(regex).fillna(False)] = k
        m = m.str.cat(s, sep=';', na_rep='')

    m = m.str.strip(';')
    m = m.str.replace(r';+', ';')

    return m
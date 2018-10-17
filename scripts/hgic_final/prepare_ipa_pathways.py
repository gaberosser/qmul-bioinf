"""
Here we load the genes involved in selected IPA pathways and export them (all together) to a single output file.
We do this for every pathway in the directory.
We remove any entry in the list does not have a valid human Entrez gene ID.
We convert IDs to Ensembl (where possible) for compatibility with downstream data.
"""

import os
import pandas as pd
from utils import output
import references
import csv
from settings import HGIC_LOCAL_DIR


# some of the pathways contain forbidden symbols; this list details which need renaming
TO_RENAME = {
    'hepatic_fibrosis_hepatic_stellate_cell_activation': 'Hepatic Fibrosis / Hepatic Stellate Cell Activation'
}


def smart_title(s):
    """
    Convert string s to title case, but leave existing uppercase words (e.g. Roman numerals)
    https://stackoverflow.com/questions/25512102/python-title-case-but-leave-pre-existing-uppercase
    :param s:
    :return:
    """
    return ' '.join(w if w.isupper() else w.capitalize() for w in s.split())


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    indir = os.path.join(HGIC_LOCAL_DIR, 'current/input_data/ipa_pathways/exported')
    flist = os.listdir(indir)

    processed = {}

    for fn in flist:
        name = fn.replace('.txt', '')
        name = TO_RENAME.get(name, smart_title(name.replace('_', ' ')))
        ff = os.path.join(indir, fn)
        this = pd.read_csv(ff, sep='\t', skiprows=1, header=0, index_col=0, usecols=range(9))

        lookup = this['Entrez Gene ID for Human'].dropna()
        if lookup.dtype != 'float':
            # this occurs when multiple Entrez IDs are given in one
            the_lookup = []
            for l in lookup.values:
                the_lookup.extend([float(t) for t in l.split('|')])
            lookup = the_lookup
        else:
            lookup = lookup.tolist()


        ens = references.entrez_to_ensembl(lookup)

        processed[name] = ens.dropna().values

    fout = os.path.join(outdir, 'ipa_exported_pathways.csv')
    with open(fout, 'wb') as f:
        c = csv.writer(f, lineterminator='\n')
        for k, v in processed.items():
            c.writerow([k] + v.tolist())

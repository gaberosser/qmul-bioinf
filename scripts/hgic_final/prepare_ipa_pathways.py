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
import unicodecsv as csv
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
    path_fn = os.path.join(indir, 'ipa_pathway_numbering.xlsx')

    output_gene_symbols = {}
    output_ens_id = {}

    pathways = pd.read_excel(path_fn, index_col=0, header=0)
    for path_id, row in pathways.iterrows():
        # get the export txt filename
        ff = os.path.join(indir, row.filename)
        this = pd.read_csv(ff, sep='\t', skiprows=1, header=0, index_col=None, usecols=range(9))

        # filter out non-genes
        this = this.loc[~this['Entrez Gene ID for Human'].isnull()]

        # The exports contain symbols and Entrez IDs

        # Entrez -> Ensembl
        lookup = this['Entrez Gene ID for Human']
        if lookup.dtype != 'float':
            # this occurs when multiple Entrez IDs are given in one
            the_lookup = []
            for l in lookup.values:
                the_lookup.extend([float(t) for t in l.split('|')])
            lookup = the_lookup
        else:
            lookup = lookup.tolist()

        ens = references.entrez_to_ensembl(lookup)
        output_ens_id[row['name']] = [path_id] + ens.dropna().values.tolist()

        # Entrez -> gene symbol
        gs = references.entrez_to_gene_symbol(lookup)
        output_gene_symbols[row['name']] = [path_id] + gs.dropna().values.tolist()

    # output to (UTF-8 encoded) CSV files

    fout = os.path.join(outdir, 'ipa_exported_pathways_symbols.csv')
    with open(fout, 'wb') as f:
        c = csv.writer(f, lineterminator='\n', encoding='utf-8')
        for k, v in output_gene_symbols.items():
            c.writerow([k] + v)

    fout = os.path.join(outdir, 'ipa_exported_pathways_ensembl_ids.csv')
    with open(fout, 'wb') as f:
        c = csv.writer(f, lineterminator='\n', encoding='utf-8')
        for k, v in output_ens_id.items():
            c.writerow([k] + v)

    # output to .rds file for import into R
    # this is for sharing with Erik Sulman's team
    from rpy2 import robjects
    to_r = [
        (v[0], robjects.StrVector(v[1:])) for v in output_gene_symbols.values()
    ]
    r_list = robjects.ListVector(to_r)
    robjects.r("saveRDS")(r_list, os.path.join(outdir, "signatures_for_ssgsea.rds"))

    # similar but maintain Ensembl IDS
    # this is for a similar comparison with our own data
    to_r = [
        (v[0], robjects.StrVector(v[1:])) for v in output_ens_id.values()
        ]
    r_list = robjects.ListVector(to_r)
    robjects.r("saveRDS")(r_list, os.path.join(outdir, "signatures_for_ssgsea_ensembl.rds"))

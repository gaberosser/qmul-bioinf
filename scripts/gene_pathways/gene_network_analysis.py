import networkx
import os
import pandas as pd
import settings
import references
import csv


pathway_files = {
    'all': 'msigdb.v6.1.symbols.gmt',
    'hallmark': 'h.all.v6.1.symbols.gmt',
    'kegg': 'c2.cp.kegg.v6.1.symbols.gmt',
    'go': 'c5.all.v6.1.symbols.gmt',
}


manual_gene_lookup = {
    'SDHD': 'ENSG00000204370',
    'CRIPAK': 'ENSG00000163945',
    'MAL2': 'ENSG00000147676'
}


if __name__ == "__main__":
    src = 'hallmark'

    indir = os.path.join(settings.GIT_LFS_DATA_DIR, 'msigdb')
    fn = os.path.join(indir, pathway_files[src])

    pathway_symbols = {}
    pathway_ens = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            pathway_symbols[row[0]] = row[2:]
            ee = references.gene_symbol_to_x_robust(row[2:], 'Ensembl Gene ID')
            if ee.isnull().any():
                still_missing = []
                for g in ee[ee.isnull()].index:
                    if g in manual_gene_lookup:
                        ee[g] = manual_gene_lookup[g]
                    else:
                        still_missing.append(g)
                if len(still_missing) > 0:
                    print "Pathway %s: Unable to find Ensembl ID for %d genes: %s" % (
                        row[0], len(still_missing), ', '.join(still_missing)
                    )
            pathway_ens[row[0]] = ee.values


    # now build a network
    # nodes are pathways, edges are genes
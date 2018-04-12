import networkx
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import settings
import references
import csv
import itertools


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

    # to get connectivity, we need to create the complementary dictionary (indexed by genes)
    gene_pathways = {}
    for pth, garr in pathway_symbols.iteritems():
        for g in garr:
            gene_pathways.setdefault(g, []).append(pth)

    # now build a network

    # 1) nodes are pathways, edges are genes
    graph_node_path = networkx.Graph(db=src)

    for pth in pathway_symbols:
        graph_node_path.add_node(pth)

    edges = {}
    for g in gene_pathways:
        conn_paths = gene_pathways[g]
        # graph_node_path.add_edges_from(itertools.combinations(conn_paths, 2), gene_symbol=g)
        for p1, p2 in itertools.combinations(conn_paths, 2):
            edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

    for (p1, p2), edge_attr in edges.iteritems():
        edge_attr['gene_count'] = len(edge_attr['genes'])
        graph_node_path.add_edge(p1, p2, **edge_attr)

    # plotting
    # we'll use the gene_count to define plotting parameters
    gc = []
    for (u, v, attr) in graph_node_path.edges(data=True):
        gc.append(attr['gene_count'])
    gc = np.array(gc)

    # need to reduce the range here somewhat
    ew = gc ** .5
    ew[gc <= 10] = 0.
    ec = gc.astype(float) / gc.max()
    pos = networkx.spring_layout(graph_node_path, iterations=200)
    networkx.draw_circular(
        graph_node_path,
        width=ew,
        edge_color=ec,
        alpha=0.8,
        edge_cmap=plt.cm.RdYlGn_r,
        edge_vmin=0, edge_vmax=1
    )
    # TODO: add node labels

    # 2) nodes are genes, edges are pathways
    graph_node_gene = networkx.Graph(db=src)

    for g in gene_pathways:
        graph_node_gene.add_node(g)

    for pth in pathway_symbols:
        conn_genes = pathway_symbols[pth]
        graph_node_gene.add_edges_from(itertools.combinations(conn_genes, 2), pathway=pth)

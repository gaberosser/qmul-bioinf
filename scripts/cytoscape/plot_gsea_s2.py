"""
Aim:

Visualise the GSEA pathway results with Cytoscape.

We're going to show the merged (union) pathway network for syngeneic and 2 x ref. comparisons, but use formatting
to distinguish nodes and edges that are reference-only or syngeneic-only.

We'll need one of these for each patient, I suppose?
"""

from apps.py2cytoscape.data.cyrest_client import CyRestClient
from apps.py2cytoscape.data import style

import networkx as nx
import os
import numpy as np
import itertools

import pandas as pd
from rnaseq import gsea
from settings import HGIC_LOCAL_DIR
import collections

from utils import setops, output


pathway_files = {
    'all': 'msigdb.v6.1.symbols.gmt',
    'hallmark': 'h.all.v6.2.symbols.gmt',
    'kegg': 'c2.cp.kegg.v6.2.symbols.gmt',
    'go': 'c5.all.v6.2.symbols.gmt',
}


if __name__ == '__main__':
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']

    # set a minimum pval for pathways to be used
    alpha = 0.05
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.1
    # small offset to avoid zero FDR
    # set this slightly lower than the smallest non-zero FDR value
    eps = 1e-6
    src = 'go'

    outdir = output.unique_output_dir()

    indir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/gsea/results/raw'
    )
    msigdb_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/gsea/msigdb/%s'
    ) % pathway_files[src]

    comparison_names = collections.OrderedDict([
        ('', 'syngeneic'),
        ('_h9_nsc', 'h9',),
        ('_gibco_nsc', 'gibco')
    ])

    # load C5 (all GO term pathways) for filtering and network analysis
    gmt = gsea.read_gmt_file(msigdb_fn)

    keep_pathways = gmt.keys()

    res = collections.OrderedDict()
    res_full = collections.OrderedDict()
    for pid in pids:
        for c in comparison_names:
            fn = os.path.join(indir, "%s%s.csv" % (pid, c))
            this = pd.read_csv(fn, sep='\t', header=0, index_col=0, usecols=[0, 3, 5, 7])
            this.columns = ['n_gene', 'nes', 'fdr']
            this = this.reindex(keep_pathways).dropna(how='all')
            res_full["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha_relevant]
            res["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha]

    # pathways_sign = sorted(setops.reduce_union(*[t.index for t in res.values()]))
    # pathways_rele = sorted(setops.reduce_union(*[t.index for t in res_full.values()]))

    # do everything with one PID for now
    pid = '018'

    # three networks to work with
    res_syn = res['%s_syngeneic' % pid]
    res_r1 = res['%s_h9' % pid]
    res_r2 = res['%s_gibco' % pid]

    all_pathways = setops.reduce_union(*[t.index for t in (res_syn, res_r1, res_r2)])

    p_to_g = dict([
        (p, gmt[p]) for p in all_pathways
    ])

    # to get connectivity, we need to create the complementary dictionary (indexed by genes)
    g_to_p = {}
    for p in all_pathways:
        for g in p_to_g[p]:
            g_to_p.setdefault(g, []).append(p)

    graph = nx.Graph(db=src)

    for pth in all_pathways:
        graph.add_node(pth, genes=p_to_g[pth], n_gene=len(p_to_g[pth]))

    edges = {}
    for g, p_arr in g_to_p.items():
        # graph_node_path.add_edges_from(itertools.combinations(conn_paths, 2), gene_symbol=g)
        for p1, p2 in itertools.combinations(p_arr, 2):
            edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

    for (p1, p2), edge_attr in edges.iteritems():
        edge_attr['gene_count'] = len(edge_attr['genes'])
        graph.add_edge(p1, p2, **edge_attr)




    p_to_p = {}
    for p, gs in p_to_g.items():
        pass


    # cy = CyRestClient()
    # # reset the session (in case something is already loaded)
    # cy.session.delete()
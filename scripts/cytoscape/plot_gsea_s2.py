"""
Aim:

Visualise the GSEA pathway results with Cytoscape.

We're going to show the merged (union) pathway network for syngeneic and 2 x ref. comparisons, but use formatting
to distinguish nodes and edges that are reference-only or syngeneic-only.

We'll need one of these for each patient, I suppose?
"""

from apps.py2cytoscape import cyrest
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
from plotting import common

from utils import setops, output


pathway_files = {
    'all': 'msigdb.v6.1.symbols.gmt',
    'hallmark': 'h.all.v6.2.symbols.gmt',
    'kegg': 'c2.cp.kegg.v6.2.symbols.gmt',
    'go': 'c5.all.v6.2.symbols.gmt',
}

if __name__ == '__main__':
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

        # minimum number of shared genes required to plot an edge
        min_n_gene_shared = 10

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

        res = collections.OrderedDict()
        res_full = collections.OrderedDict()
        for pid in pids:
            for c in comparison_names:
                fn = os.path.join(indir, "%s%s.csv" % (pid, c))
                this = pd.read_csv(fn, sep='\t', header=0, index_col=0, usecols=[0, 3, 5, 7])
                this.columns = ['n_gene', 'nes', 'fdr']
                this = this.reindex(gmt.keys()).dropna(how='all')
                res_full["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha_relevant]
                res["%s_%s" % (pid, comparison_names[c])] = this.loc[this.fdr < alpha]

        # functional API - the python bindings are incomplete here?
        cy = CyRestClient()
        # reset the session (in case something is already loaded)
        cy.session.delete()

        # command API - the python bindings are much better
        cy_cmd = cyrest.cyclient()

        # do everything with one PID for now
        pid = '018'
        for pid in pids:

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

            # we're going to use passthrough mapping to customise the node colour
            # we'll define 3 colourmaps, with -log10(p) assigning the shade:
            # greyscale for syn. and ref.
            # reds for ref. only
            # blues for syn. only
            # colours are defined by HEX values? Add these to the nodes
            logp_syn = -np.log10(res_syn.fdr + eps)
            logp_r1 = -np.log10(res_r1.fdr + eps)
            logp_r2 = -np.log10(res_r2.fdr + eps)
            vmax = max(
                logp_syn.max(),
                logp_r1.max(),
                logp_r2.max(),
            )

            cmap_both_func = common.continuous_cmap(0, vmax, cmap='Greys')
            cmap_syn_func = common.continuous_cmap(0, vmax, cmap='Blues')
            cmap_ref_func = common.continuous_cmap(0, vmax, cmap='Reds')
            vs, _ = setops.venn_from_arrays(*[t.index for t in (res_syn, res_r1, res_r2)])

            node_colours = {}
            for pth in vs['111'] + vs['101'] + vs['110']:
                node_colours[pth] = cmap_both_func(logp_syn[pth])
            for pth in vs['100']:
                node_colours[pth] = cmap_syn_func(logp_syn[pth])
            for pth in vs['011']:
                # mean P for refs
                node_colours[pth] = cmap_ref_func(0.5 * (logp_r1[pth] + logp_r2[pth]))
            for pth in vs['010']:
                node_colours[pth] = cmap_ref_func(logp_r1[pth])
            for pth in vs['001']:
                node_colours[pth] = cmap_ref_func(logp_r2[pth])

            graph = nx.Graph(db=src)

            for pth in all_pathways:
                graph.add_node(
                    pth,
                    genes=p_to_g[pth],
                    n_gene=len(p_to_g[pth]),
                    fill_colour=node_colours[pth],
                    comparisons=[v for u, v in zip([res_syn, res_r1, res_r2], ['Syngeneic', 'Ref. 1', 'Ref. 2']) if pth in u.index]
                )

            edges = {}
            for g, p_arr in g_to_p.items():
                # graph_node_path.add_edges_from(itertools.combinations(conn_paths, 2), gene_symbol=g)
                for p1, p2 in itertools.combinations(p_arr, 2):
                    edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

            for (p1, p2), edge_attr in edges.iteritems():
                edge_attr['gene_count'] = len(edge_attr['genes'])
                if edge_attr['gene_count'] >= min_n_gene_shared:
                    graph.add_edge(p1, p2, **edge_attr)

            # add network
            cy_net = cy.network.create_from_networkx(graph, collection=pid)

            cy_style = cy.style.create(pid)

            # passthrough node colour, label
            cy_style.create_passthrough_mapping(column='fill_colour', vp='NODE_FILL_COLOR', col_type='String')
            cy_style.create_passthrough_mapping(column='name', vp='NODE_LABEL', col_type='String')

            # node size based on number of genes
            vmax = max([len(v) for v in p_to_g.values()])
            node_size_map = style.StyleUtil.create_slope(min=0, max=vmax, values=(10, 50))
            cy_style.create_continuous_mapping(column='n_gene', vp='NODE_SIZE', col_type='Double', points=node_size_map)

            # edge width based on number of genes in common
            vmax = max([v['gene_count'] for v in edges.values()])
            edge_width_map = style.StyleUtil.create_slope(min=min_n_gene_shared, max=vmax, values=(0, 5))
            cy_style.create_continuous_mapping(column='gene_count', vp='EDGE_WIDTH', col_type='Double', points=edge_width_map)

            cy.style.apply(cy_style, network=cy_net)

            # layout
            cy_cmd.layout.force_directed(
                EdgeAttribute='gene_count',
                defaultSpringLength=70,
                defaultSpringCoefficient=50,
                maxWeightCutoff=vmax,
                network=pid
            )




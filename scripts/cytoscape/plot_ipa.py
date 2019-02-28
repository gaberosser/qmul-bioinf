"""
Aim:

Visualise the IPA core analysis results with Cytoscape.

We're going to show the merged (union) pathway network for syngeneic and 2 x ref. comparisons, but use formatting
to distinguish nodes and edges that are reference-only or syngeneic-only.

We'll need one of these for each patient, I suppose.

Then will investigate using a single combined network, highlighting only nodes that stand out as syngeneic only or
patient-specific.
"""

from apps.py2cytoscape import cyrest
from apps.py2cytoscape.data.cyrest_client import CyRestClient
from apps.py2cytoscape.data import style

import networkx as nx
import os
import numpy as np
import itertools

import pandas as pd
from cytoscape import cyto
from settings import HGIC_LOCAL_DIR
import collections
from plotting import common

from utils import setops, output, ipa
from scripts.hgic_final.analyse_xcell_results import IPA_PATHWAY_DIR, IPA_PATHWAY_FN
from scripts.hgic_final import consts


if __name__ == '__main__':
    pids = consts.PIDS
    # basic pvalue cutoff
    alpha = 0.05
    log_alpha = -np.log10(alpha)
    # more stringent cutoff, as used in the IPA S1/S2 analysis
    alpha_strict = 0.005
    log_alpha_strict = -np.log10(alpha_strict)
    # minimum number of shared genes required to plot an edge
    min_n_gene_shared = 10

    comparisons = ['syngeneic', 'h9', 'gibco']

    # load pathways, processed data and estimate signatures
    ipa_res = pd.read_excel(IPA_PATHWAY_FN)
    all_ipa = ipa.load_raw_reports(IPA_PATHWAY_DIR, "de_s2_{0}_{1}.txt", pids, comparisons)
    ipa_signatures = ipa.load_supported_signatures_from_raw(
        IPA_PATHWAY_DIR,
        "de_s2_{0}_{1}.txt",
        [pids, comparisons],
        pathways=ipa_res.index
    )

    cy_obj = cyto.CytoscapeSession()
    nx_graphs = {}

    # one network per patient:
    for pid in pids:
        this_ipa = [all_ipa[(pid, c)].loc[all_ipa[(pid, c)]['-logp'] >= log_alpha_strict] for c in comparisons]
        all_pathways = setops.reduce_union(*[t.index for t in this_ipa])

        p_to_g = {}
        for p in all_pathways:
            p_to_g[p] = setops.reduce_union(*[t.loc[p, 'genes'].split(',') if p in t.index else [] for t in this_ipa])

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
        logp_vals = [t['-logp'] for t in this_ipa]
        vmax = max([t.max() for t in logp_vals])

        # we need a lower offset for the non-grey colours, otherwise all the white shades look very similar
        vmin = -2
        cmap_both_func = common.continuous_cmap(0, vmax, cmap='Greys')
        cmap_syn_func = common.continuous_cmap(vmin, vmax, cmap='Blues')
        cmap_ref_func = common.continuous_cmap(vmin, vmax, cmap='Reds')
        vs, _ = setops.venn_from_arrays(*[t.index for t in this_ipa])

        node_significance = {}
        node_colours = {}
        node_attrs = {}
        for k, p_arr in vs.items():
            ix = [t == '1' for t in k]
            n = float(sum(ix))
            for pth in p_arr:
                m = 0
                for i, t in enumerate(ix):
                    if t:
                        m += logp_vals[i][pth]
                        node_attrs.setdefault(pth, {})["plogp_%s" % comparisons[i]] = logp_vals[i][pth]
                m /= n
                node_attrs.setdefault(pth, {})['mean_plogp'] = m
                if ix[0] and any(ix[1:]):
                    node_colours[pth] = cmap_both_func(m)
                elif ix[0]:
                    node_colours[pth] = cmap_syn_func(m)
                else:
                    node_colours[pth] = cmap_ref_func(m)

        graph = nx.Graph(db="%s_ipa" % pid)

        for pth in all_pathways:
            graph.add_node(
                pth,
                genes=sorted(p_to_g[pth]),
                n_gene=len(p_to_g[pth]),
                fill_colour=node_colours[pth],
                comparisons=[v for u, v in zip(this_ipa, comparisons) if pth in u.index],
                **node_attrs[pth]
            )

        edges = {}
        for g, p_arr in g_to_p.items():
            for p1, p2 in itertools.combinations(p_arr, 2):
                edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

        for (p1, p2), edge_attr in edges.iteritems():
            edge_attr['gene_count'] = len(edge_attr['genes'])
            if edge_attr['gene_count'] >= min_n_gene_shared:
                graph.add_edge(p1, p2, **edge_attr)

        nx_graphs[pid] = graph

    max_node_val = max([len(v) for v in p_to_g.values()])

    for name, graph in nx_graphs.items():
        this_net = cy_obj.add_networkx_graph(graph, name=name)
        this_net.passthrough_node_fill('fill_colour')
        this_net.passthrough_node_label('name')
        this_net.passthrough_node_size_linear('n_gene', xmax=max_node_val)
        this_net.passthrough_edge_width_linear('gene_count', xmin=min_n_gene_shared, ymin=0., ymax=5)

    layout_kwargs = dict(
        EdgeAttribute='gene_count',
        defaultSpringLength=70,
        defaultSpringCoefficient=50,
        maxWeightCutoff=max_node_val,
    )

    cy_obj.apply_layout(**layout_kwargs)

from cytoscape import cyto
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import networkx as nx
import itertools
import collections

from plotting import common
from utils import output, ipa, setops, dictionary
from settings import GIT_LFS_DATA_DIR
from scripts.hgic_final import ipa_results_s1_s2 as irss
from scripts.hgic_final import consts


def nx_graph_from_ipa_single(df, name=None, min_edge_count=0):
    p_to_g = {}
    for p, row in df.iterrows():
        p_to_g[p] = row.genes.split(',')

    # to get connectivity, we need to create the complementary dictionary (indexed by genes)
    g_to_p = {}
    for p, g_arr in p_to_g.items():
        for g in g_arr:
            g_to_p.setdefault(g, []).append(p)

    node_attrs = {}

    for p in p_to_g:
        node_attrs[p] = df.loc[p].drop('genes')

    graph = nx.Graph(name=name)

    for p in p_to_g.keys():
        graph.add_node(
            p,
            genes=sorted(p_to_g[p]),
            **node_attrs[p]
        )

    edges = {}
    for g, p_arr in g_to_p.items():
        for p1, p2 in itertools.combinations(p_arr, 2):
            edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

    for (p1, p2), edge_attr in edges.iteritems():
        edge_attr['gene_count'] = len(edge_attr['genes'])
        if edge_attr['gene_count'] >= min_edge_count:
            graph.add_edge(p1, p2, **edge_attr)

    return graph


def nx_graph_from_ipa_multiple(ipa_dict, name=None, min_edge_count=0, attr_cols=('genes', '-logp')):
    """
    Generate a single graph from multiple inputs.
    Each node is a pathway, which includes an attribute showing which patients are members.
    :param ipa_dict: Dictionary of IPA results dataframes
    :return:
    """
    p_to_g = collections.defaultdict(dict)
    p_to_g_union = {}

    # create a combined dataframe so we can add node attributes later
    this_comb = []

    for k, df in ipa_dict.items():
        this_df = pd.DataFrame(
            df[list(attr_cols)].values,
            index=df.index,
            columns=["%s_%s" % (k, t) for t in attr_cols]
        )
        this_comb.append(this_df)
        p_to_g[k] = {}
        for p, row in df.iterrows():
            this_genes = row.genes.split(',')
            p_to_g[k][p] = this_genes
            p_to_g_union.setdefault(p, set()).update(this_genes)

    this_comb = pd.concat(this_comb, axis=1, sort=True)
    node_attrs = {}

    for p in p_to_g_union:
        node_attrs[p] = this_comb.loc[p].dropna()
        node_attrs[p]['n_gene'] = len(p_to_g_union[p])

    g_to_p_union = dictionary.complement_dictionary_of_iterables(p_to_g_union)

    graph = nx.Graph(name=name)

    for p in p_to_g_union.keys():
        graph.add_node(
            p,
            all_genes=sorted(p_to_g_union[p]),
            **node_attrs[p]
        )

    edges = {}
    for g, p_arr in g_to_p_union.items():
        for p1, p2 in itertools.combinations(p_arr, 2):
            edges.setdefault((p1, p2), {}).setdefault('genes', []).append(g)

    for (p1, p2), edge_attr in edges.iteritems():
        edge_attr['gene_count'] = len(edge_attr['genes'])
        if edge_attr['gene_count'] >= min_edge_count:
            graph.add_edge(p1, p2, **edge_attr)

    return graph


if __name__ == '__main__':
    # set a minimum pval for pathways to be used
    alpha = 0.005
    plogalpha = -np.log10(alpha)
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    quantile = 0.99

    min_edge_count = 0

    pids = consts.PIDS

    outdir = output.unique_output_dir()

    indir = os.path.join(GIT_LFS_DATA_DIR, 'ipa_from_biplots')

    patient_colours = {
        '018': '#ccffcc',
        '019': '#4dff4d',
        '030': '#00cc00',
        '031': '#004d00',
        '017': '#ffcccc',
        '050': '#ff4d4d',
        '054': '#cc0000',
        '061': '#660000',
        '026': '#ff80ff',
        '052': '#800080'
    }

    cyto_nets = {}

    ## PC1, mean logFC (GIC vs iNSC)

    fn = os.path.join(indir, 'for_ipa_mean_logfc', 'pc1_mean.txt')
    this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
    this.columns = ['-logp', 'ratio', 'z', 'genes']
    # add ngenes column
    this.insert(3, 'n_gene', this.genes.str.split(',').apply(len))
    this.index = [x.decode('utf-8') for x in this.index]
    # NaN in the z attribute will cause a parsing error
    this.drop('z', axis=1, inplace=True)

    this = this.loc[this['-logp'] > plogalpha_relevant]

    cy_session = cyto.CytoscapeSession()
    gg = nx_graph_from_ipa_single(this, name='PC1 mean logFC', min_edge_count=min_edge_count)

    # add fill colours, deliberately saturating late to avoid deepest shade
    vmax = this['-logp'].max() * 1.2
    node_cmap = common.continuous_cmap(plogalpha_relevant, vmax, cmap='Reds')
    for k, n in gg.nodes.items():
        n['fill_colour'] = node_cmap(n['-logp'])

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)
    cyto_nets[gg.name] = this_net

    # formatting
    this_net.passthrough_node_label('name')
    this_net.passthrough_node_size_linear('n_gene')
    this_net.passthrough_edge_width_linear('gene_count', xmin=1., ymin=0.2, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')

    # this can be buggy if not applied last (?)
    this_net.passthrough_node_fill('fill_colour')

    ## PC2, separate logFCs

    all_ipa = ipa.load_raw_reports(
        os.path.join(indir, 'for_ipa_separate_logfc_pc2'),
        "{0}.txt",
        pids
    )

    # reduce to relevant results
    for k, df in all_ipa.items():
        all_ipa[k] = df.loc[df['-logp'] > plogalpha_relevant]

    gg = nx_graph_from_ipa_multiple(all_ipa, name='PC2 separate logFC')

    # add to node attributes (enabling pie chart)
    for k, n in gg.nodes.items():
        for p in pids:
            n[p] = int("%s_-logp" % p in n)

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    # formatting
    this_net.passthrough_node_label('name')
    this_net.passthrough_node_size_linear('n_gene')
    this_net.passthrough_edge_width_linear('gene_count', xmin=1., ymin=0.2, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])

    ## PC3, separate logFCs

    all_ipa = ipa.load_raw_reports(
        os.path.join(indir, 'for_ipa_separate_logfc_pc3'),
        "{0}.txt",
        pids
    )

    # reduce to relevant results
    for k, df in all_ipa.items():
        all_ipa[k] = df.loc[df['-logp'] > plogalpha_relevant]

    gg = nx_graph_from_ipa_multiple(all_ipa, name='PC3 separate logFC')

    # add to node attributes (enabling pie chart)
    for k, n in gg.nodes.items():
        for p in pids:
            n[p] = int("%s_-logp" % p in n)

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    # formatting
    this_net.passthrough_node_label('name')
    this_net.passthrough_node_size_linear('n_gene')
    this_net.passthrough_edge_width_linear('gene_count', xmin=1., ymin=0.2, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])

    layout_kwargs = dict(
        EdgeAttribute=None,
        defaultSpringLength=150,
        defaultSpringCoefficient=150,
        maxWeightCutoff=max([v['n_gene'] for v in gg.nodes.values()]),
    )

    cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))
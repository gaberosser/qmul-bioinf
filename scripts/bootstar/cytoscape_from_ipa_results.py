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

from settings import HGIC_LOCAL_DIR


if __name__ == '__main__':
    """
    Here we create a Cytoscape session for the visualisation of the IPA (and GO?) pathways identified as enriched
    in each patient. We'll be comparing between patients, hence it's a 'Strategy 1' approach.
    """
    # set a minimum pval for pathways to be used
    alpha = 0.01
    plogalpha = -np.log10(alpha)

    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    min_edge_count = 4

    pids = [
        '19', '31', '44', '49', '50', '52'
    ]

    patient_colours = dict(zip(pids, [
        '#1b9e77',
        '#d95f02',
        '#7570b3',
        '#e7298a',
        '#66a61e',
        '#e6ab02',
    ]))

    indir = '~/Documents/ipa_output_from_jb_20190313'

    outdir = output.unique_output_dir()

    # Cytoscape session
    cy_session = cyto.CytoscapeSession()
    cyto_nets = {}

    #######################################################
    # DE
    #######################################################
    ipa_de_res = ipa.load_raw_reports(indir, "{0}_pathways_BC.txt", pids)
    for k in ipa_de_res:
        ipa_de_res[k] = ipa_de_res[k].loc[ipa_de_res[k]['-logp'] >= plogalpha]

    gg = ipa.nx_graph_from_ipa_multiple(ipa_de_res, name='BC overlap pathways from DE genes', min_edge_count=min_edge_count)

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    # formatting
    this_net.passthrough_node_label('name')
    this_net.passthrough_node_size_linear('n_genes')
    this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])

    layout_kwargs = dict(
        EdgeAttribute='n_genes',
        defaultSpringLength=40,
        defaultSpringCoefficient=40,
        maxWeightCutoff=max([v['n_genes'] for v in gg.nodes.values()]),
    )

    cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))
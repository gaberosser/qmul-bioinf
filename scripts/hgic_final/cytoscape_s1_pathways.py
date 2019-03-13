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

    min_edge_count = 8

    pids = consts.PIDS

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

    de_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/ipa/pathways')
    dm_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/methylation/s0_individual_patients_direct_comparison/ipa/pathways')
    dedm_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq_methylation_combined/s0_individual_patients_direct_comparison/ipa/pathways')

    outdir = output.unique_output_dir()

    pids = consts.PIDS

    # Cytoscape session
    cy_session = cyto.CytoscapeSession()
    cyto_nets = {}

    #######################################################
    # DE
    #######################################################
    ipa_de_res = collections.OrderedDict()
    for pid in pids:
        fn = os.path.join(de_indir, "full_de_patient{pid}.xls".format(pid=pid))
        this_df = pd.read_excel(fn, skiprows=1, header=0, index_col=0)
        this_df.columns = ['-logp', 'ratio', 'z', 'genes']
        this_df.insert(3, 'n_gene', this_df.genes.str.split(',').apply(len))
        # filter to include only relevant pathways
        ipa_de_res[pid] = this_df.loc[this_df['-logp'] >= plogalpha]

    gg = ipa.nx_graph_from_ipa_multiple(ipa_de_res, name='IPA from DE genes', min_edge_count=min_edge_count)

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

    #######################################################
    # DM
    #######################################################
    ipa_dm_res = ipa.load_raw_reports(dm_indir, "{0}.txt", pids)
    # filter to include only significant pathways
    for k in ipa_dm_res:
        ipa_dm_res[k] = ipa_dm_res[k].loc[ipa_dm_res[k]['-logp'] >= plogalpha]

    gg = ipa.nx_graph_from_ipa_multiple(ipa_dm_res, name='IPA from DM genes', min_edge_count=min_edge_count)

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
        defaultSpringLength=150,
        defaultSpringCoefficient=150,
        maxWeightCutoff=max([v['n_genes'] for v in gg.nodes.values()]),
    )

    cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))
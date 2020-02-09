from cytoscape import cyto
import os
import pandas as pd
import numpy as np
import csv
import time
from matplotlib import pyplot as plt
import seaborn as sns
import networkx as nx
import itertools
import collections
from scipy.spatial import distance

from plotting import common
from utils import output, ipa, setops, dictionary
from settings import GIT_LFS_DATA_DIR
from scripts.hgic_final import ipa_results_s1_s2 as irss
from scripts.hgic_final import consts

from settings import HGIC_LOCAL_DIR


def plot_network(g, x, y, ax=None):
    if ax is None: ax = plt.gca()
    # nodes
    ax.scatter(x, y, marker='o', edgecolor='k', facecolor='none', s=40)
    # edges
    for (k1, k2), e in g.edges.items():
        i = g.nodes.keys().index(k1)
        j = g.nodes.keys().index(k2)
        ax.plot(
            [x[i], x[j]],
            [y[i], y[j]],
            color='0.5',
            linestyle='-',
            lw=0.1 * e['n_genes'],
            alpha=0.5
        )


if __name__ == '__main__':
    """
    Here we create a Cytoscape session for the visualisation of the IPA (and GO?) pathways identified as enriched
    in each patient. We'll create a separate network for each of the patients, hence it's a 'Strategy 0' approach.
    """
    # set a minimum pval for pathways to be used
    alpha = 0.01
    plogalpha = -np.log10(alpha)

    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    min_edge_count = 8

    label_top_n_nodes = 15

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

    outdir = output.unique_output_dir()

    pids = consts.PIDS

    # generate standalone legend
    legend_dict = collections.OrderedDict([
        (pid, {'class': 'patch', 'edgecolor': 'k', 'facecolor': patient_colours[pid], 'linewidth': 1.})
        for pid in pids
    ])
    legend_dict = {'': legend_dict}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('off')
    fig.set_facecolor('w')
    common.add_custom_legend(ax, legend_dict, loc='center')
    fig.set_size_inches((1.4, 2.4))
    fig.savefig(os.path.join(outdir, "cytoscape_legend.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "cytoscape_legend.tiff"), dpi=200)

    # load exported pathways
    ipa_pathway_export_fn = os.path.join(HGIC_LOCAL_DIR, 'current/input_data/ipa_pathways/ipa_exported_pathways_symbols.csv')
    ipa_pathways = {}
    with open(ipa_pathway_export_fn, 'r') as f:
        c = csv.reader(f)
        for row in c:
            ipa_pathways[row[0].decode('utf-8')] = row[2:]  # the second entry is the 'anonymisation number' for sending to externals


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

        # drop 'z' column before generating graph as this is often undefined
        gg = ipa.nx_graph_from_ipa_single(ipa_de_res[pid].drop('z', axis=1), name=pid, min_edge_count=min_edge_count)

        # order by plogp
        plogp = dict([(k, v['-logp']) for k, v in gg.nodes.items()])

        to_label = sorted(plogp, key=lambda x: -plogp[x])[:label_top_n_nodes]

        for k, v in gg.nodes.iteritems():
            if k in to_label:
                v['name_vis_top%d' % label_top_n_nodes] = k
                v['z_top%s' % label_top_n_nodes] = 10.
            else:
                v['name_vis_top%d' % label_top_n_nodes] = ''
                v['z_top%d' % label_top_n_nodes] = 5.

        this_net = cy_session.add_networkx_graph(gg, name=gg.name)
        cyto_nets[gg.name] = this_net

        # this_net.passthrough_node_size_linear('n_genes')
        this_net.passthrough_node_size_linear('-logp')
        this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
        this_net.set_node_border_width(0.)
        this_net.set_edge_colour('#b7b7b7')
        this_net.set_node_fill_colour(patient_colours[pid])
        this_net.set_node_transparency(255)
        this_net.passthrough_node_label('name')

        layout_kwargs = dict(
            EdgeAttribute='n_genes',
            defaultSpringLength=150,
            defaultSpringCoefficient=150,
            maxWeightCutoff=max([v['n_genes'] for v in gg.nodes.values()]),
        )

        cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))

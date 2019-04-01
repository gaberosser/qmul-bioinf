from cytoscape import cyto
import os
import pandas as pd
import numpy as np
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

    de_pathways_manually_labelled = [
        'Chondroitin Sulfate Biosynthesis (Late Stages)',
        'Chondroitin Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis (Late Stages)',
        'Dermatan Sulfate Biosynthesis',
        u'T Helper Cell Differentiation',
        u'Nur77 Signaling in T Lymphocytes',
        u'Cdc42 Signaling',
        u'B Cell Development',
    ]
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
    # add mean plogp
    mean_plogp = {}
    for k, v in gg.nodes.iteritems():
        all_logp = [v['%s_-logp' % pid] for pid in pids if '%s_-logp' % pid in v]
        if len(all_logp) > 0:
            v['plogp_mean'] = np.mean(all_logp)
            mean_plogp[k] = v['plogp_mean']

    # add conditional node labels
    # these should have a higher z value than the unlabelled nodes
    to_label = sorted(mean_plogp, key=lambda x:-mean_plogp[x])[:label_top_n_nodes]
    to_label.extend(de_pathways_manually_labelled)
    for k, v in gg.nodes.iteritems():
        if k in to_label:
            v['name_vis'] = k
            v['z'] = 10.
        else:
            v['z'] = 5.

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    # formatting
    this_net.passthrough_node_label('name_vis')
    # this_net.passthrough_node_size_linear('n_genes')
    this_net.passthrough_node_size_linear('plogp_mean')
    this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)
    this_net._create_passthrough_mapping('z', 'NODE_Z_LOCATION', col_type='Double')

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])


    # try to write my own node layout algo
    n = len(gg.nodes)
    default_length = 50.
    min_length = 10.
    spring_const = 50.
    group_const = 50.
    n_iter_max = 100
    delta_t = 0.1

    idx_to_node = pd.Series(gg.nodes)
    node_to_idx = pd.Series(range(n), index=gg.nodes)

    # pre-compute pairwise edge weights and group similarity values
    group_similarity = np.zeros([n] * 2)

    for i, k1 in enumerate(gg.nodes):
        for j, k2 in enumerate(gg.nodes):
            if j <= i:
                continue
            members_i = [pid for pid in pids if gg.nodes[k1][pid] == 1]
            members_j = [pid for pid in pids if gg.nodes[k2][pid] == 1]
            n_common = len(set(members_i).intersection(members_j))
            if n_common == 0:
                continue
            n_tot = float(len(members_i) + len(members_j))
            # n_union = float(len(set(members_i).union(members_j)))
            group_similarity[i, j] = 2 * n_common / n_tot
    group_similarity = group_const * group_similarity

    edge_similarity = np.zeros([n] * 2)

    for k, e in gg.edges.iteritems():
        i = node_to_idx[k[0]]
        j = node_to_idx[k[1]]
        edge_similarity[i, j] = e['n_genes']
    eq_lengths = default_length - edge_similarity
    eq_lengths[eq_lengths < min_length] = min_length

    # initial locations: random
    d = (1 - .5 ** (1 / float(n))) ** .5
    l = default_length / d
    xy = np.random.rand(n, 2) * l

    iter_count = 0
    m = 1e6  # max movement - used as convergence indicator
    tol = 2.
    while m > tol:
        if iter_count > n_iter_max:
            break

        # pairwise distances
        # dr = distance.squareform(distance.pdist(xy))

        dx = reduce(
            lambda x, y: x - y,
            np.meshgrid(xy[:, 0], xy[:, 0])
        )
        dy = reduce(
            lambda x, y: x - y,
            np.meshgrid(xy[:, 1], xy[:, 1])
        )
        dr = (dx ** 2 + dy ** 2) ** .5

        # unit vectors
        dv = np.array([dx, dy])

        # forces (as vectors)


        iter_count += 1


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
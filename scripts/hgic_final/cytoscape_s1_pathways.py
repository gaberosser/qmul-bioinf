from cytoscape import cyto
import os
import pandas as pd
import numpy as np
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


def compute_group_similarity(nx_g, pids=consts.PIDS):
    idx_to_node = pd.Series(nx_g.nodes)
    n = len(nx_g.nodes)
    # precompute indices of the pairwise comparisons
    triu_ind = np.triu_indices(n, k=1)

    group_similarity = []

    for i, j in zip(*triu_ind):
        k_i = idx_to_node[i]
        k_j = idx_to_node[j]

        # group similarity
        members_i = [pid for pid in pids if nx_g.nodes[k_i][pid] == 1]
        members_j = [pid for pid in pids if nx_g.nodes[k_j][pid] == 1]
        n_common = len(set(members_i).intersection(members_j))
        n_tot = float(len(members_i) + len(members_j))
        group_similarity.append((2 * n_common / n_tot) ** 2)

    return np.array(group_similarity)


def network_layout_force_directed_algorithm(
        cyto_net,
        group_similarity,
        edge_weight_col='n_genes',
        default_edge_length=50.,
        min_edge_length=10.,
        spring_const=50.,
        group_const=300.,
        node_repulsion_const=4e5,
        contraction_const=1e5,
        n_iter_max=1000,
        delta_t=0.1,
        animated_gif_out=None,
        update_each_frame=True,
):
    generate_animation = animated_gif_out is not None
    if generate_animation:
        update_each_frame = True

    nx_g = cyto_net.nx_graph
    n = len(nx_g.nodes)

    # it will be useful to precompute indices of the pairwise comparisons
    triu_ind = np.triu_indices(n, k=1)

    idx_to_node = pd.Series(nx_g.nodes)

    # pre-compute pairwise edge weights and group similarity values
    edge_similarity = []

    for i, j in zip(*triu_ind):
        k_i = idx_to_node[i]
        k_j = idx_to_node[j]

        # edge_similarity
        if (k_i, k_j) in gg.edges:
            edge_similarity.append(nx_g.edges[(k_i, k_j)][edge_weight_col])
        else:
            edge_similarity.append(0)

    edge_similarity = np.array(edge_similarity)

    eq_lengths = default_edge_length - edge_similarity
    eq_lengths[eq_lengths < min_edge_length] = min_edge_length

    # we need a boolean multiplier to encode connectivity
    k_e = np.ones_like(eq_lengths) * edge_similarity.astype(bool) * spring_const

    # initial locations: random
    d = (1 - .5 ** (1 / float(n))) ** .5
    l = default_edge_length / d
    xy = np.random.rand(2, n) * l

    # apply layout
    x = dict(zip(gg.nodes.keys(), 10 * xy[0]))
    y = dict(zip(gg.nodes.keys(), 10 * xy[1]))
    this_net.add_node_column('x', x)
    this_net.add_node_column('y', y)
    this_net.passthrough_node_position('x', 'y')
    this_net.view_fit_content()

    if generate_animation:
        images = []
        import PIL
        from StringIO import StringIO

    iter_count = 0
    m = 1e6  # max movement - used as convergence indicator
    tol = 2.
    while m > tol:
        if iter_count > n_iter_max:
            break

        # centre of mass
        com = xy.mean(axis=1)

        dx = xy[0, triu_ind[1]] - xy[0, triu_ind[0]]
        dy = xy[1, triu_ind[1]] - xy[1, triu_ind[0]]
        dr = (dx ** 2 + dy ** 2) ** .5

        # unit pairwise vectors
        dv = np.array([dx, dy]) / dr

        # unit vector from com
        ux = xy[0] - com[0]
        uy = xy[1] - com[1]
        dist_from_com = (ux ** 2 + uy ** 2) ** .5
        ux /= dist_from_com
        uy /= dist_from_com
        u = np.array([ux, uy])

        # forces (as vectors)
        # regular node interaction via edges
        f_e = k_e * (eq_lengths - dr) * -dv

        # node repulsion
        # apply a small offset in the denominator to avoid very large jumps
        f_n = node_repulsion_const / (dr ** 2 + 1.) * -dv

        # group attraction
        f_g = group_similarity * group_const * dr * dv

        # pairwise forces in a squareform
        f_tot = delta_t ** 2 / 2. * (f_e + f_n + f_g)

        # 'open out' pairwise format, remembering to change the sign on the lower diagonal
        # if we have any non-symmetric terms, this won't work
        f_tot_x = distance.squareform(f_tot[0])
        f_tot_x[np.tril_indices_from(f_tot_x, k=1)] *= -1

        f_tot_y = distance.squareform(f_tot[1])
        f_tot_y[np.tril_indices_from(f_tot_y, k=1)] *= -1

        # sum for each node
        force = np.array([
            f_tot_x.sum(axis=1),
            f_tot_y.sum(axis=1),
        ])

        # contraction force
        f_c = contraction_const * np.exp(-(dist_from_com ** 2)) * -u

        force += delta_t ** 2 / 2. * f_c

        d_xy = force * delta_t ** 2 * 0.5

        # summary metric for extent of movement
        m = (d_xy ** 2).sum(axis=0).max()

        xy = xy + d_xy

        print "Iteration %d, m = %.2f" % (iter_count, m)

        # optionally update the plot
        if update_each_frame:
            x = dict(zip(gg.nodes.keys(), 10 * xy[0]))
            y = dict(zip(gg.nodes.keys(), 10 * xy[1]))
            this_net.add_node_column('x', x)
            this_net.add_node_column('y', y)
            this_net.view_fit_content()

        if generate_animation:
            images.append(PIL.Image.open(StringIO(this_net.obj.get_png(height=2000))))

        iter_count += 1

    # apply final layout
    x = dict(zip(gg.nodes.keys(), 10 * xy[0]))
    y = dict(zip(gg.nodes.keys(), 10 * xy[1]))
    this_net.add_node_column('x', x)
    this_net.add_node_column('y', y)
    this_net.passthrough_node_position('x', 'y')

    if generate_animation:
        images[0].save(
            animated_gif_out,
            save_all=True,
            append_images=images[1:],
            duration=100,  # units: ms
            loop=0
        )


if __name__ == '__main__':
    """
    Here we create a Cytoscape session for the visualisation of the IPA (and GO?) pathways identified as enriched
    in each patient. We'll be comparing between patients, hence it's a 'Strategy 1' approach.
    """
    generate_animation = True

    # set a minimum pval for pathways to be used
    alpha = 0.01
    plogalpha = -np.log10(alpha)

    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    min_edge_count = 8

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

    # add a few different versions of conditional node labels
    # in all cases these should have a higher z value than the unlabelled nodes

    # 1. Manually listed pathways plus top N
    to_label = sorted(mean_plogp, key=lambda x:-mean_plogp[x])[:label_top_n_nodes]
    to_label.extend(de_pathways_manually_labelled)
    for k, v in gg.nodes.iteritems():
        if k in to_label:
            v['name_vis_1'] = k
            v['z_1'] = 10.
        else:
            v['name_vis_1'] = ''
            v['z_1'] = 5.

    # 2. Manually listed pathways only
    for k, v in gg.nodes.iteritems():
        if k in de_pathways_manually_labelled:
            v['name_vis_2'] = k
            v['z_2'] = 10.
        else:
            v['name_vis_2'] = ''
            v['z_2'] = 5.

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    # this_net.passthrough_node_size_linear('n_genes')
    this_net.passthrough_node_size_linear('plogp_mean')
    this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])

    group_similarity = compute_group_similarity(gg)

    animated_gif_fn = None
    if generate_animation:
        animated_gif_fn = os.path.join(outdir, "net_layout_animation_de.gif")

    network_layout_force_directed_algorithm(
        this_net,
        group_similarity,
        animated_gif_out=animated_gif_fn
    )

    # cycle through node labelling types and export images

    # no labels
    this_net.export_png(os.path.join(outdir, "cytoscape_de_network_no_node_labels.png"), height=2000)

    # label version 1
    this_net.passthrough_node_label('name_vis_1')
    this_net._create_passthrough_mapping('z_1', 'NODE_Z_LOCATION', col_type='Double')
    this_net.export_png(os.path.join(outdir, "cytoscape_de_network_selected_and_topn_node_labels.png"), height=2000)

    # label version 2
    this_net.passthrough_node_label('name_vis_2')
    this_net._create_passthrough_mapping('z_2', 'NODE_Z_LOCATION', col_type='Double')
    this_net.export_png(os.path.join(outdir, "cytoscape_de_network_selected_node_labels.png"), height=2000)

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

    group_similarity = compute_group_similarity(gg)

    animated_gif_fn = None
    if generate_animation:
        animated_gif_fn = os.path.join(outdir, "net_layout_animation_dm.gif")

    network_layout_force_directed_algorithm(
        this_net,
        group_similarity,
        animated_gif_out=animated_gif_fn
    )

    # layout_kwargs = dict(
    #     EdgeAttribute='n_genes',
    #     defaultSpringLength=150,
    #     defaultSpringCoefficient=150,
    #     maxWeightCutoff=max([v['n_genes'] for v in gg.nodes.values()]),
    # )
    #
    # cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))

    # NB can convert the animated gifs generated here into MP4 (more efficient and versatile) on the command line:
    # ffmpeg -i net_layout_animation_dm.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" net_layout_animation_dm.mp4
    # Source: https://unix.stackexchange.com/a/294892/253111
    # TODO: consider getting there directly with openCV
    # Source: https://blog.extramaster.net/2015/07/python-pil-to-mp4.html

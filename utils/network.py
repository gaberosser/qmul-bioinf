import networkx as nx
import itertools
from utils import dictionary, setops
import collections


def edges_from_node_members(node_dict, member_key='members'):
    # get the complementary dict
    member_node_dict = dictionary.complement_dictionary_of_iterables(node_dict)
    edges = {}
    for m, k_arr in member_node_dict.items():
        # connect every pair of nodes
        for p1, p2 in itertools.combinations(k_arr, 2):
            edges.setdefault((p1, p2), {}).setdefault(member_key, []).append(m)
    return edges


def add_edges_to_graph(graph, edges, member_key='members', min_edge_count=0):
    """

    :param graph: networkx graph.
    :param edges: Computed by `edges_from_node_members`
    :param member_key:
    :param min_edge_count: Only edges with this many shared members or more are added
    :return:
    """
    member_count_key = 'n_%s' % member_key
    for (p1, p2), edge_attr in edges.iteritems():
        edge_attr[member_count_key] = len(edge_attr[member_key])
        if edge_attr[member_count_key] >= min_edge_count:
            graph.add_edge(p1, p2, **edge_attr)


def nx_graph_from_overlapping_members(node_member_dict, member_key='members', name=None, min_edge_count=0, node_attrs=None):
    """
    Generate a networkx.Graph from the input dictionary, which is keyed by node names. Values are 'members' of that
    node, which are used to form the edges in the network.
    :param node_member_dict:
    :param member_key: This is the name used for the node and edge attributes which list the members. It might be more
    context-specific to rename it (e.g. 'genes').
    :param name:
    :param min_edge_count: Optionally specify a minimum overlap between nodes, below which edges are not created.
    :param node_attrs: If supplied, this is a dict of dicts (or similar), keyed by node names. The attributes contained
    in it are associated with nodes in the resulting graph.
    :return:
    """
    if node_attrs is None:
        node_attrs = {}

    member_count_key = 'n_%s' % member_key
    graph = nx.Graph(name=name)

    for k in node_member_dict.keys():
        node_attr = node_attrs.get(k, {})
        if not isinstance(node_attr, dict):
            node_attr = dict(node_attr)
        node_attr.update({member_key: node_member_dict[k], member_count_key: len(node_member_dict[k])})
        graph.add_node(
            k,
            **node_attr
        )

    edges = edges_from_node_members(node_member_dict, member_key=member_key)
    add_edges_to_graph(graph, edges, member_key=member_key, min_edge_count=min_edge_count)

    return graph


def nx_graph_from_multiple_member_dicts(
        dict_of_node_member_dict,
        member_key='members',
        participants_key='participants',
        name=None,
        min_edge_count=0,
        dict_of_node_attrs=None,
        reduction_method='union',
        add_boolean_indicator=True,
        use_namespace_for_individuals=True
):
    """
    Generate a single networkx graph to summarise multiple 'member dictionary' inputs.
    We do this by reducing the information so that an edge reflects a summary of all instances with that edge.
    :param dict_of_node_member_dict: Dict of dicts, each of which is keyed by node name and has values giving
    members of that node. This is the same format required by `nx_graph_from_overlapping_members`.
    :param id_key: The key name used to identify the network ID.
    :param member_key: This is the key used for the node and edge attributes which list the members. It might be more
    context-specific to rename it (e.g. 'genes').
    :param participants_key: The key used to list the participants (top level - e.g. patients)
    :param name:
    :param min_edge_count: Optionally specify a minimum overlap between nodes, below which edges are not created.
    :param dict_of_node_attrs: Dictionary of dictionary of dictionaries, levels as follows:
    - network
        - node
            - attribute dictionary
    :param reduction_method: The method by which we summarise multiple lists of members in order to connect edges
    :param add_boolean_indicator: If True (default), we attach an attribute to each node for each possible member that
    is a Boolean indicator of that member's participation
    :return:
    """
    member_count_key = 'n_%s' % member_key

    if reduction_method == 'union':
        reduce_fun = setops.reduce_union
    elif reduction_method == 'intersection':
        reduce_fun = setops.reduce_intersection
    else:
        raise AttributeError("Unsupported reduction_method %s." % reduction_method)

    sub_graphs = dict([
        (
            k,
            nx_graph_from_overlapping_members(
                v,
                member_key=member_key,
                name=k,
                min_edge_count=min_edge_count,
                node_attrs=dict_of_node_attrs.get(k) if dict_of_node_attrs is not None else None,
            )
        ) for k, v in dict_of_node_member_dict.iteritems()
    ])

    graph = nx.Graph(name=name)
    # nodes = collections.defaultdict(dict)
    node_members = collections.defaultdict(dict)

    # combine subgraphs
    for k, g in sub_graphs.items():
        for node_name, node_attrs in g.nodes.items():
            node_members[node_name][k] = node_attrs[member_key]

    # combine subgraphs
    # prepare node attrs
    node_attrs = collections.defaultdict(dict)

    for k, v in dict_of_node_attrs.items():
        for node_name, attrs in v.items():
            # rename attrs
            new_attrs = {}
            for x, y in attrs.items():
                if use_namespace_for_individuals:
                    the_key = "%s::%s_%s" % (k, k, x)
                else:
                    the_key = "%s_%s" % (k, x)
                new_attrs[the_key] = y
            # new_attrs = dict([
            #     (the_key, y) for x, y in attrs.items()
            # ])
            node_attrs[node_name].update(new_attrs)
            node_attrs[node_name].setdefault(participants_key, []).append(k)

    # add boolean indicator if required
    if add_boolean_indicator:
        for nn in node_attrs:
            for k in sub_graphs:
                if use_namespace_for_individuals:
                    the_key = "%s::%s" % (k, k)
                else:
                    the_key = k
                node_attrs[nn][the_key] = int(k in node_attrs[nn][participants_key])

    # reduce node members
    node_members_reduced = {}
    for node_name, coll in  node_members.items():
        node_members_reduced[node_name] = sorted(reduce_fun(*coll.values()))

    edges = edges_from_node_members(node_members_reduced, member_key=member_key)

    # add nodes and attributes to graph
    for node_name, md in node_members.items():
        node_attrs[node_name][member_key] = node_members_reduced[node_name]
        node_attrs[node_name][member_count_key] = len(node_members_reduced[node_name])

        graph.add_node(
            node_name,
            **node_attrs[node_name]
        )

    add_edges_to_graph(graph, edges, member_key=member_key, min_edge_count=min_edge_count)

    return graph

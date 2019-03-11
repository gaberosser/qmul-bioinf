import datetime
import os
import csv
import itertools
import pandas as pd
import collections
import networkx as nx
from utils import dictionary, network


def results_to_ipa_format(
        blocks,
        outdir,
        incl_cols=('logFC', 'FDR'),
        identifier='Ensembl',
):
    """
    Create IPA-compatible text files for each of the items in blocks.
    These are suitable for batch upload.
    :param blocks: Dict containing the results to export. Each element is treated as a separate file.
    """
    incl_cols = list(incl_cols)

    for k, bl in blocks.iteritems():
        fn = os.path.join(outdir, "%s.txt" % k)
        header = [
            ['Key', 'Value'],
            ['observation_name', k],
            ['date_created', datetime.datetime.now().isoformat()]
        ]
        if identifier is not None:
            header += [['identifier_types', identifier]]
        with open(fn, 'wb') as f:
            c = csv.writer(f, delimiter='\t')
            # meta header
            c.writerows(header)
            c.writerow(['Data_begins_here'])
            # data column header
            c.writerow(['ID'] + incl_cols)
            # reduced block
            reduced_block = bl.loc[:, incl_cols]
            c.writerows(reduced_block.itertuples())


def load_raw_reports(indir, file_pattern, *args):
    """
    Load all raw reports matching file_pattern in the specified directory.
    Reports are processed to clean up the column naming, then packaged into a dictionary.
    :param indir: Input directory.
    :param file_pattern: A string supporting .format() operations (i.e. {0}, {1}, etc).
    The arrays in *args will be passed in, so these should match.
    :param args: One or more arrays containing variables that will be substituted into filename_format using .format().
    :return: Flat dictionary, keys given as tuples matching the input *args.
    """
    res = collections.OrderedDict()
    for tup in itertools.product(*args):
        fn = os.path.join(indir, file_pattern.format(*tup))
        this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
        this.columns = ['-logp', 'ratio', 'z', 'genes']
        # add ngenes column
        this.insert(3, 'n_gene', this.genes.str.split(',').apply(len))
        this.index = [x.decode('utf-8') for x in this.index]
        k = tup
        if len(k) == 1:
            k = k[0]
        res[k] = this
    return res


def load_supported_signatures_from_raw(indir, file_pattern, args, pathways=None):
    """
    Based on all the raw reports, compute the union of the features involved in each of the pathways.
    For example, this can be used to estimate a number of the DE genes in each pathway.
    :param indir:
    :param file_pattern:
    :param args: Iterable of iterables, which will be passed into the input *args of `load_raw_reports`
    :param pathways: If supplied, only these pathways are included
    :return: Dictionary, keyed by pathway name. Values are lists of features.
    """
    res = load_raw_reports(indir, file_pattern, *args)
    ipa_pathway_signatures = {}
    for this in res.values():
        if pathways is None:
            this_pw_list = this.index
        else:
            this_pw_list = this.index.intersection(pathways)
        for pw in this_pw_list:
            this_list = set(this.loc[pw, 'genes'].split(','))
            if pw in ipa_pathway_signatures:
                ipa_pathway_signatures[pw] = ipa_pathway_signatures[pw].union(this_list)
            else:
                ipa_pathway_signatures[pw] = this_list

    return ipa_pathway_signatures


def nx_graph_from_ipa_single(df, name=None, min_edge_count=0):
    """
    Generate a networkx graph from IPA results. Nodes are pathways, edges represent similarity between pathways due to
    shared genes.
    :param df: A single result from load_raw_results. Index is pathways names. Must have a column named 'genes' giving
    the deregulated genes in a string, separated by commas.
    :param name: Name to give the graph, optional.
    :param min_edge_count: If supplied, this is applied as a cutoff: edges are only created when the number of shared
    genes is equal to or greater than this value.
    :return: networkx.Graph object
    """
    p_to_g = {}
    node_attrs = {}

    for p, row in df.iterrows():
        p_to_g[p] = row.genes.split(',')
        node_attrs[p] = dict(df.loc[p].drop('genes'))

    return network.nx_graph_from_overlapping_members(
        p_to_g,
        member_key='genes',
        name=name,
        min_edge_count=min_edge_count,
        node_attrs=node_attrs
    )


def nx_graph_from_ipa_multiple(
        ipa_dict,
        name=None,
        min_edge_count=0,
        attr_cols=('genes', '-logp', 'ratio', 'n_gene'),
        reduction_method='union'
):
    """
    Generate a single networkx graph from multiple inputs.
    Each node is a pathway, which includes an attribute showing which patients are members.
    :param ipa_dict: Dictionary of IPA results dataframes, as returned by `load_raw_reports`. Keys are used to name
    attribute columns.
    :param name: Optionally proivde a name for the graph.
    :param min_edge_count: See `nx_graph_from_ipa_single`.
    :param attr_cols: Tuple giving the columns that will be added to the node attributes.
    :return: networkx.Graph object. Nodes will have attributes named `{key}_{attr_name}`, where {key} is given by the
    input dictionary.
    """
    p_to_g = {}
    node_attrs = {}

    for k, df in ipa_dict.items():
        p_to_g[k] = {}
        node_attrs[k] = {}
        for p, row in df.iterrows():
            this_genes = row.genes.split(',')
            p_to_g[k][p] = this_genes
            node_attrs[k][p] = dict(df.loc[p][list(attr_cols)])

    return network.nx_graph_from_multiple_member_dicts(
        p_to_g,
        member_key='genes',
        participants_key='patients',
        name=name,
        min_edge_count=min_edge_count,
        dict_of_node_attrs=node_attrs,
        reduction_method=reduction_method,
        use_namespace_for_individuals=False
    )

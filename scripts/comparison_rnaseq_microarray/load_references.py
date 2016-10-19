from settings import DATA_DIR
import os
import dill
import pandas as pd
import logging
import operator

# (re)set logger
logger = logging.getLogger("load_references")
logger.handlers = []
sh = logging.StreamHandler()
fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(fmt)
logger.addHandler(sh)
logger.setLevel(logging.INFO)


MB_GENE_GROUPS = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1')),
)
MB_GENES = reduce(operator.add, [list(t[1]) for t in MB_GENE_GROUPS])

MB_ENTREZ_GROUPS = (
    ('WNT', (11197, 3371, 2571, 27123, 2018)),
    ('SHH', (27295, 2138, 64399, 474, 6422)),
    ('Group C', (50939, 2558, 346007, 4901, 10586, 4883)),
    ('Group D', (3736, 8320, 202559, 221662, 137970, 4938)),
)
MB_ENTREZ = reduce(operator.add, [list(t[1]) for t in MB_ENTREZ_GROUPS])


def get_structure_ids_by_parent(ontol, parent_id):
    """
    Recurse through the DataFrame ontol, retrieving structure IDs that are directly beneath parent_id in the
    nested hierarchy
    :param ontol: DataFrame indexed by structure ID and with a column named parent_structure_id
    :param parent_id: The parent ID
    :return: set of structure IDs
    """
    struct_ids = {parent_id}
    this_set = {parent_id}
    while True:
        filt = ontol[ontol.parent_structure_id.isin(this_set)]
        if filt.empty:
            break
        else:
            struct_ids.update(filt.index)
            this_set = set(filt.index)
    return struct_ids


def prepare_cerebellum_microarray_reference_data(outfile=None):

    INDIR = os.path.join(DATA_DIR, 'allen_human_brain_atlas/microarray')
    DONOR_NUMBERS = [
        9861,
        10021,
        12876,
        14380,
        15496,
        15697
    ]
    PARENT_STRUCT_ID = 4696

    if outfile is None:
        outfile = os.path.join(INDIR, 'cerebellum_expression_dataframes.pkl')

    # load probe library
    probe_fn = os.path.join(INDIR, 'Probes.csv')
    probes = pd.read_csv(probe_fn, index_col=0)
    # keep only those probes with an Entrez ID
    probes = probes.dropna(axis=0, subset=['entrez_id'])

    # load ontology library
    ontol_fn = os.path.join(INDIR, 'Ontology.csv')
    ontol = pd.read_csv(ontol_fn, index_col=0)

    struct_ids = get_structure_ids_by_parent(ontol, PARENT_STRUCT_ID)

    expression = pd.DataFrame()
    sample_meta = pd.DataFrame()

    for dn in DONOR_NUMBERS:
        logger.info("Processing donor %d", dn)
        this_dir = os.path.join(INDIR, "donor%d" % dn)
        expre_fn = os.path.join(this_dir, 'MicroarrayExpression.csv')
        sampl_fn = os.path.join(this_dir, 'SampleAnnot.csv')
        pacal_fn = os.path.join(this_dir, 'PACall.csv')

        sampl = pd.read_csv(sampl_fn)
        pacal = pd.read_csv(pacal_fn, header=None, index_col=0).astype(bool)
        expre = pd.read_csv(expre_fn, header=None, index_col=0)

        # set sample IDs
        sampl_idx = pd.Index(['%d_%d' % (dn, i) for i in range(sampl.shape[0])])
        sampl.index = sampl_idx
        expre.columns = sampl_idx
        pacal.columns = sampl_idx

        # filter sample annotation and expression for recognized probes
        pacal = pacal.loc[probes.index]
        expre = expre.loc[probes.index]

        # filter expression by PA call
        # this replaces all non-significant results with NaN
        expre = expre[pacal]

        # filter sample annotation by ontology
        sampl_idx2 = sampl[sampl.structure_id.isin(struct_ids)].index

        # filter expression by sample
        expre = expre[sampl_idx2]

        # concatenate along axis 1
        expression = pd.concat([expression, expre], axis=1)

        # add sample metadata to the list
        this_meta = sampl.loc[sampl_idx2]
        this_meta['donor_id'] = dn
        sample_meta = sample_meta.append(this_meta)

        logger.info("Completed donor %d", dn)

    # prepend gene symbol and entrez ID to the total expression dataframe
    expression = pd.concat([probes.loc[expression.index, ['gene_symbol', 'entrez_id']], expression], axis=1)

    logger.info("Saving results to %s", outfile)
    with open(outfile, 'wb') as f:
        dill.dump(
            {'expression': expression, 'meta': sample_meta},
            f)
    logger.info('Complete')

    return expression, sample_meta


def load_cerebellum_microarray_reference_data():

    INDIR = os.path.join(DATA_DIR, 'allen_human_brain_atlas/microarray')
    infile = os.path.join(INDIR, 'cerebellum_expression_dataframes.pkl')
    if not os.path.exists(infile):
        logger.info("Unable to find pre-prepared pickled file %s. Recomputing.", infile)
        return prepare_cerebellum_microarray_reference_data()
    else:
        logger.info("Loading from pre-prepared pickled file %s.", infile)
        with open(infile, 'rb') as f:
            res = dill.load(f)
        return res['expression'], res['meta']


def _microarray_probe_to_gene(marray_data, lookup=None, method='max', groupby='gene_symbol'):
    grp_data = marray_data.groupby(groupby)
    if method == 'max':
        data = grp_data.max()
    elif method == 'mean':
        data = grp_data.mean()
    elif method == 'sum':
        data = grp_data.sum()
    else:
        raise ValueError("Method %s not supported", method)

    if lookup is not None:
        return data.loc[lookup]
    else:
        return data


def microarray_gene_markers(marray_data, genes=None, method='max'):
    """
    Load the raw microarray Allen reference data then convert to a gene-centric dataframe.
    Many genes are represented on multiple probes.
    Following Zhao et al (PLoS One 2014), we take the most active gene whenever we encounter ambiguity (method='max').
    Also can take the sum, or mean, etc... specified by the method param.
    Optionally permit providing a list of genes of interest, otherwise all are computed.
    :param marray_data: pandas dataframe containing normalised intensity data. Must have a column 'gene_symbol'
    :param genes: Optional list of genes to use, otherwise use all.
    :param method: 'max', 'mean', 'sum'
    :return:
    """
    data = marray_data.drop('entrez_id', axis=1)
    return _microarray_probe_to_gene(data, lookup=genes, method=method, groupby='gene_symbol')


def microarray_entrez_markers(marray_data, entrez_ids=None, method='max'):
    """
    Analogous function to microarray_gene_markers, only entrez IDs are used instead of symbols
    :param marray_data:
    :param entrez_ids:
    :param method:
    :return:
    """
    data = marray_data.drop('gene_symbol', axis=1)
    return _microarray_probe_to_gene(data, lookup=entrez_ids, method=method, groupby='entrez_id')


def prepare_cerebellum_rnaseq_reference_data(outfile=None):
    DONOR_NUMBERS = [
        9861,
        10021
    ]
    INDIR = os.path.join(DATA_DIR, 'allen_human_brain_atlas/rnaseq')
    PARENT_STRUCT_ID = 4696

    if outfile is None:
        outfile = os.path.join(INDIR, 'cerebellum_expression_dataframes.pkl')

    # load gene library
    # unnecessary unless we want Entrez IDs
    # genes_fn = os.path.join(INDIR, 'Genes.csv')
    # genes = pd.read_csv(genes_fn, header=0, index_col=0)

    # load ontology library
    ontol_fn = os.path.join(INDIR, 'Ontology.csv')
    ontol = pd.read_csv(ontol_fn, index_col=0)

    struct_ids = get_structure_ids_by_parent(ontol, PARENT_STRUCT_ID)

    expression_cts = pd.DataFrame()
    expression_tpm = pd.DataFrame()
    sample_meta = pd.DataFrame()

    for dn in DONOR_NUMBERS:
        logger.info("Processing donor %d", dn)
        this_dir = os.path.join(INDIR, "donor%d" % dn)
        tpm_fn = os.path.join(this_dir, 'RNAseqTPM.csv.gz')
        cts_fn = os.path.join(this_dir, 'RNAseqCounts.csv.gz')
        sampl_fn = os.path.join(this_dir, 'SampleAnnot.csv.gz')

        sampl = pd.read_csv(sampl_fn)
        tpm = pd.read_csv(tpm_fn, header=None, index_col=0)
        cts = pd.read_csv(cts_fn, header=None, index_col=0)

        # set sample IDs
        sampl_idx = pd.Index(['%d_%d' % (dn, i) for i in range(sampl.shape[0])])
        sampl.index = sampl_idx
        tpm.columns = sampl_idx
        cts.columns = sampl_idx

        # filter sample annotation by ontology
        sampl_idx2 = sampl[sampl['ontology_structure_id'].isin(struct_ids)].index

        # filter expression by sample
        tpm = tpm[sampl_idx2]
        cts = cts[sampl_idx2]

        # concatenate along axis 1
        expression_tpm = pd.concat([expression_tpm, tpm], axis=1)
        expression_cts = pd.concat([expression_cts, cts], axis=1)

        # add sample metadata to the list
        this_meta = sampl.loc[sampl_idx2]
        this_meta['donor_id'] = dn
        sample_meta = sample_meta.append(this_meta)

        logger.info("Completed donor %d", dn)

    # give expresion index a meaningful name
    expression_tpm.index.name = 'gene_symbol'
    expression_cts.index.name = 'gene_symbol'

    logger.info("Saving results to %s", outfile)
    with open(outfile, 'wb') as f:
        dill.dump(
            {
                'tpm': expression_tpm,
                'count': expression_cts,
                'meta': sample_meta
            },
            f)
    logger.info('Complete')

    return expression_cts, expression_tpm, sample_meta


def load_cerebellum_rnaseq_reference_data():

    INDIR = os.path.join(DATA_DIR, 'allen_human_brain_atlas/rnaseq')
    infile = os.path.join(INDIR, 'cerebellum_expression_dataframes.pkl')
    if not os.path.exists(infile):
        logger.info("Unable to find pre-prepared pickled file %s. Recomputing.", infile)
        return prepare_cerebellum_rnaseq_reference_data()
    else:
        logger.info("Loading from pre-prepared pickled file %s.", infile)
        with open(infile, 'rb') as f:
            res = dill.load(f)
        return res['count'], res['tpm'], res['meta']


def plot_microarray_gene_markers():
    from matplotlib import pyplot as plt
    import seaborn as sns
    plt.interactive(True)

    expr, meta = load_cerebellum_microarray_reference_data()

    markers = microarray_gene_markers(expr, genes=MB_GENES).transpose()
    fig, axs = plt.subplots(ncols=len(MB_GENE_GROUPS), sharey=True, figsize=(10, 6))
    for i in range(len(MB_GENE_GROUPS)):
        grp, arr = MB_GENE_GROUPS[i]
        ax = axs[i]
        subdata = markers.loc[:, markers.columns.isin(arr)].fillna(0)
        subdata.boxplot(ax=ax, rot=90)
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(0.5, -.21)
        if i == 0:
            ax.set_ylabel('Normalised intensity')
    ylim = list(ax.get_ylim())
    ylim[0] = -0.4
    ax.set_ylim(ylim)
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.99, wspace=0.05, hspace=0.)


def plot_microarray_entrez_markers():
    from matplotlib import pyplot as plt
    import seaborn as sns
    plt.interactive(True)

    expr, meta = load_cerebellum_microarray_reference_data()

    markers = microarray_entrez_markers(expr, entrez_ids=MB_ENTREZ).transpose()
    fig, axs = plt.subplots(ncols=len(MB_ENTREZ_GROUPS), sharey=True, figsize=(10, 6))
    for i in range(len(MB_ENTREZ_GROUPS)):
        grp, arr = MB_ENTREZ_GROUPS[i]
        ax = axs[i]
        subdata = markers.loc[:, markers.columns.isin(arr)].fillna(0)
        subdata.boxplot(ax=ax, rot=90)
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(0.5, -.21)
        if i == 0:
            ax.set_ylabel('Normalised intensity')
    ylim = list(ax.get_ylim())
    ylim[0] = -0.4
    ax.set_ylim(ylim)
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.2, top=0.99, wspace=0.05, hspace=0.)


if __name__ == '__main__':
    mb_marr_ref, mb_marr_meta = load_cerebellum_microarray_reference_data()
    mb_marr_ref_by_entrez = microarray_entrez_markers(mb_marr_ref, entrez_ids=MB_ENTREZ)
    # plot_microarray_entrez_markers()
    # plot_microarray_gene_markers()
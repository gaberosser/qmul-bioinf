import datetime
import os
import csv
import itertools
import pandas as pd
import collections


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
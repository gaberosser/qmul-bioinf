import datetime
import os
import csv


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
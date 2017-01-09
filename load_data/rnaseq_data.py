import pandas as pd
import os
import references
from settings import DATA_DIR, DATA_DIR_NON_GIT

INDEX_FIELDS = (
    'Approved Symbol',
    'Entrez Gene ID',
    'RefSeq IDs',
    'Ensembl Gene ID'
)


def gse83696(index_by='Ensembl Gene ID'):
    """
    Data are initially indexed by Ensembl ID. Coversion is carried out using HGNC data, if required.
    Index field options are: Approved Symbol, Entrez Gene ID, RefSeq IDs
    :param index_by:
    :return:
    """
    # TODO: convert this into a generic loader for htseq-count outputs
    indir = os.path.join(DATA_DIR, 'rnaseq_GSE83696', 'htseq-count')
    samples = [
        ('XZ1', 'XZ-1.count'),
        ('XZ2', 'XZ-2.count'),
        ('XZ3', 'XZ-3.count'),
        ('XZ4', 'XZ-4.count'),
        ('XZ5', 'XZ-5.count'),
        ('XZ6', 'XZ-6.count'),
        ('XZ7', 'XZ-7.count'),
        ('XZ8', 'XZ-8.count'),
    ]
    df = pd.DataFrame()
    for sn, fn in samples:
        t = pd.read_csv(os.path.join(indir, fn), sep='\t', index_col=0, header=None).iloc[:, 0]
        df.loc[:, sn] = t

    if index_by is not None and index_by != 'Ensembl Gene ID':
        new_idx = references.translate(df.index, to_field=index_by, from_field='Ensembl Gene ID')
        new_idx.dropna(inplace=True)
        df = df.loc[new_idx.index]
        df.index = new_idx.values
    return df


def paired_samples_gene_counts(units='counts', annotate_by=None, annotation_type='protein_coding'):
    """
    :param units: One of 'counts', 'fpkm', 'tpm'
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation. Otherwise
      add all supported annotations
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    :return:
    """
    supported_annot = ('Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs')
    if annotate_by is not None and annotate_by not in supported_annot:
        raise ValueError("Unrecognised annotation requested. Supported options are %s" % ', '.join(supported_annot))

    supported_units = ('counts', 'fpkm', 'tpm')
    if units not in supported_units:
        raise ValueError("Unrecognised units requested. Supported options are %s" % ', '.join(supported_units))

    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    infiles = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'GBM018',
        'GBM019',
        'GBM024',
        'GBM026',
        'GBM031',
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA024N28_NSC',
        'DURA026N31D_NSC',
        'DURA031N44B_NSC',
    )

    first_run = True
    res = None
    for fn, mfn in zip(infiles, metafiles):
        meta = pd.read_csv(mfn, header=0, index_col=0)
        dat = pd.read_csv(fn, comment='#', header=0, index_col=0, sep='\t')

        # keep only counts for each requested sample
        codes_to_keep = meta.loc[meta.loc[:, 'sample'].isin(samples)].index
        files_to_keep = ['%s.bam' % t for t in codes_to_keep]
        sample_names = meta.loc[codes_to_keep, 'sample'].values

        if units in ('fpkm', 'tpm'):
            nreads = meta.loc[codes_to_keep, 'read_count']
            nreads.index = sample_names
            lengths = dat.Length

        dat = dat.loc[:, files_to_keep]
        dat.columns = sample_names

        if units in ('fpkm', 'tpm'):
            dat = dat.divide(nreads, axis=1).divide(lengths, axis=0) * 1e9
            if units == 'tpm':
                dat = dat.divide(dat.sum(axis=0), axis=1)

        if first_run:
            res = dat.copy()
            first_run = False
        else:
            res += dat

    # load genenames data for annotation
    df = references.conversion_table(type=annotation_type)
    df.set_index('Ensembl Gene ID', inplace=True)
    annot = df.loc[res.index.intersection(df.index), ['Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs']]

    # add annotations (where available)
    # take a record of the columns first
    cols = res.columns
    res = pd.concat((res, annot), axis=1)

    if annotate_by is not None:
        # drop any rows that do not have an annotation
        res.dropna(axis=0, subset=[annotate_by], inplace=True)
        # set index then drop other annotations
        res.set_index(annotate_by, inplace=True)
        res = res.loc[:, cols]

    return res
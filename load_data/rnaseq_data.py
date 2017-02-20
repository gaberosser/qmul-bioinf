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


def featurecounts(
        count_files,
        metafiles,
        samples=None,
        units='counts',
        annotate_by='all',
        annotation_type='protein_coding'):
    """
    Iterate over one or more lanes of gene count data computed using featureCounts, summing data over the lanes.

    :param count_files: Iterable. Each entry is the path to a featureCounts output.
    :param metafiles: Iterable, same length as count_files. Each entry is the path to a CSV file. Each row details a
    sample. The column 'sample' gives the identifier. The first (index) column is the file stem in the counts file.
    :param samples: If supplied, these are the samples to keep. Everything else is dropped.
    :param units: One of 'counts', 'fpkm', 'tpm'
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
    If 'all' add all supported annotations.
    If None, add no extra annotations.
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    :return:
    """
    supported_annot = ('Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs', 'all')
    if annotate_by is not None and annotate_by not in supported_annot:
        raise ValueError("Unrecognised annotation requested. Supported options are %s" % ', '.join(supported_annot))

    supported_units = ('counts', 'fpkm', 'tpm')
    if units not in supported_units:
        raise ValueError("Unrecognised units requested. Supported options are %s" % ', '.join(supported_units))

    first_run = True
    res = None
    lengths = None
    nreads = None
    meta = None
    for fn, mfn in zip(count_files, metafiles):
        meta = pd.read_csv(mfn, header=0, index_col=0)
        dat = pd.read_csv(fn, comment='#', header=0, index_col=0, sep='\t')

        # keep only counts for each requested sample
        if samples is not None:
            codes_to_keep = meta.loc[meta.loc[:, 'sample'].isin(samples)].index
        else:
            # keep all samples
            codes_to_keep = meta.index

        files_to_keep = ['%s.bam' % t for t in codes_to_keep]
        sample_names = meta.loc[codes_to_keep, 'sample'].values
        lengths = dat.Length

        dat = dat.loc[:, files_to_keep]
        dat.columns = sample_names

        if first_run:
            res = dat.copy()
            # get total reads and length of each transcript
            nreads = meta.loc[codes_to_keep, 'read_count']
            nreads.index = sample_names
            first_run = False
        else:
            # sum counts and total count tally
            res += dat
            nr = meta.loc[codes_to_keep, 'read_count']
            nr.index = sample_names
            nreads += nr

    # normalise (if required) AFTER combining all count files
    # we are assuming that the transcripts do not vary between files
    if units == 'fpkm':
        res = res.divide(nreads, axis=1).divide(lengths, axis=0) * 1e9
    elif units == 'tpm':
        rpk = res.divide(lengths, axis=0)
        res = rpk.divide(rpk.sum(axis=0), axis=1) * 1e6

    if annotate_by is None:
        return res, meta
    else:
        # load genenames data for annotation
        df = references.conversion_table(type=annotation_type)
        df.set_index('Ensembl Gene ID', inplace=True)

        if annotate_by == 'all':
            annot = df.loc[res.index.intersection(df.index), ['Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs']]
        else:
            annot = df.loc[res.index.intersection(df.index), annotate_by]
        # take a record of the columns first
        # cols = res.columns

        # add annotation columns
        res = pd.concat((res, annot), axis=1)

        # if one annotation was requested, set that as the index
        if annotate_by != 'all':
            # drop any rows that do not have an annotation
            res.dropna(axis=0, subset=[annotate_by], inplace=True)
            # set index
            res.set_index(annotate_by, inplace=True)
            # res = res.loc[:, cols]

        return res, meta


def gbm_paired_samples(units='counts', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
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
    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


def gbm_astrocyte_nsc_samples(units='counts', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA018N2_ASTRO_DAY12',
        'DURA019N8C_ASTRO_DAY12',
    )
    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


def mb_zhao_cultures(units='counts', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        '1078',
        '1595',
        '1487'
    )
    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


def brainrnaseq_preprocessed():
    """
    Load the pre-processed results from www.brainrnaseq.org
    These are in units of FPKM, annotated by gene. Annoyingly, the gene symbols are for mouse.
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE73721')
    infile = os.path.join(indir, 'fpkm', 'fpkm.csv')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    data = pd.read_csv(infile, header=0, index_col=0)

    return data, meta

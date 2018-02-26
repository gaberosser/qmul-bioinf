import pandas as pd
import os
import warnings
import glob
import random
import collections
import numpy as np
import gzip
import csv
from settings import LOCAL_DATA_DIR


GTF_SOURCES = ('ensembl', 'havana', 'ensembl_havana')


def get_reference_genome_directory(tax_id, version):
    if tax_id == 9606:
        if version is None:
            return os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90')
        else:
            raise NotImplementedError("Not yet supporting multiple versions.")
    elif tax_id == 10090:
        if version is None:
            return os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'mouse', 'ensembl', 'GRCm38.p5.r88')
        else:
            raise NotImplementedError("Not yet supporting multiple versions.")
    else:
        raise ValueError("Unsupported tax_id: %d" % tax_id)


def reference_genome_chrom_lengths(tax_id=9606, version=None, fn=None):
    """
    Get the lengths of each chromosome in the reference genome.
    :param tax_id: Taxonomy ID
    :param version: Not currently supported, but will allow the specification of more than one version.
    :param fn: If supplied, this directly overrides any other inputs and is used.
    """
    if fn is None:
        indir = os.path.join(get_reference_genome_directory(tax_id, version), 'fa')
        flist = glob.glob(os.path.join(indir, '*.fai'))
        if len(flist) == 0:
            raise AttributeError("No .fai files found in directory %s." % indir)
        elif len(flist) > 1:
            warnings.warn("Found >1 .fai files. Using an arbitrary choice.")
        fn = flist[0]

    fai = pd.read_csv(fn, sep='\t', header=None, index_col=0)
    return fai.iloc[:, 0]


def random_genomic_interval(chrom_lengths, n_bp):
    """
    Generate a random genomic interval based on the input chromosome lengths
    """
    the_chr = chrom_lengths.index[random.randint(0, chrom_lengths.size - 1)]
    l = chrom_lengths.loc[the_chr]
    the_start = random.randint(0, l - n_bp)
    the_end = the_start + n_bp - 1
    return (the_chr, the_start, the_end)


def get_promoter_regions(
        tax_id=9606,
        version=None,
        fn=None,
        upstream_dist=1000,
        sources=GTF_SOURCES,
        zero_based=True
):
    if fn is None:
        indir = os.path.join(get_reference_genome_directory(tax_id, version), 'gtf')
        flist = glob.glob(os.path.join(indir, '*.gtf'))
        is_gz = False
        if len(flist) == 0:
            flist = glob.glob(os.path.join(indir, '*.gtf.gz'))
            if len(flist) == 0:
                raise AttributeError("Directory %s contains no .gtf or .gtf.gz files" % indir)
            is_gz = True
        if len(flist) > 1:
            warnings.warn("Found >1 .gtf or .gtf.gz files. Using an arbitrary choice: %s." % flist[0])
        fn = flist[0]
    else:
        is_gz = fn[-3:].lower() == '.gz'

    res = []

    if is_gz:
        fstream_func = gzip.open
    else:
        fstream_func = open

    with fstream_func(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for i, row in enumerate(c):
            if len(row) == 1:
                continue
            chr = row[0]
            src = row[1]

            if sources is not None and src not in sources:
                continue

            typ = row[2]
            blk = row[8].split('; ')

            attr = dict([t.split(' ', 1) for t in blk])
            if (
                typ == 'exon'
                and attr.get('exon_number') == '"1"'
                and attr.get('gene_biotype') == '"protein_coding"'
            ):
                # strip speech marks
                attr = dict([(k, v.strip('"')) for k, v in attr.items()])
                strand = row[6]

                if strand == '+':
                    start = int(row[3])
                    stop = int(row[4])
                    tss_minus = start - upstream_dist
                elif strand == '-':
                    start = int(row[4])
                    stop = int(row[3])
                    tss_minus = start + upstream_dist
                else:
                    raise ValueError("We expect strand to be + or - but it isn't here: %s" % str(row))

                if zero_based:
                    start -= 1
                    stop -= 1

                attr['chr'] = chr
                attr['source'] = src
                attr['type'] = typ
                attr['start'] = start
                attr['stop'] = stop
                attr['strand'] = strand
                attr['promoter_region_start'] = tss_minus
                attr['promoter_region_end'] = start

                res.append(attr)

    return res

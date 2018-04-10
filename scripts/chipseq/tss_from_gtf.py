import gzip
import os
import re
import csv
import collections
from settings import LOCAL_DATA_DIR
from utils import output


LIMIT = int(1e12)
SOURCES = {'ensembl', 'havana', 'ensembl_havana'}


def get_feature_tss_from_gtf(fn, filter_func, name_field, distance=2500, sources=None, limit=None):
    """
    :param filter_func: Function that takes two inputs: type (GTF field 3) and attributes (final GTF field, pre-split)
    and returns True if the entry should be included
    :param name_field: The attributes field to use for the name
    """
    if sources is None:
        sources = SOURCES

    if os.path.splitext(fn)[1].lower() == '.gz':
        opener = gzip.open
    else:
        opener = open

    seen = collections.defaultdict(list)
    regions = []

    with opener(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        line_no = 1
        for row in c:
            if len(row) == 1:
                continue
            if limit is not None and line_no > limit:
                break
            chr = row[0]
            src = row[1]
            typ = row[2]
            strand = row[6]
            blk = row[8].split('; ')

            if strand == '+':
                start = int(row[3])
                # stop = int(row[4])
            elif strand == '-':
                start = int(row[4])
                # stop = int(row[3])
            else:
                raise ValueError("Unexpected strand value %s" % strand)

            attr = dict([t.split(' ') for t in blk])

            if src in sources and filter_func(typ, attr):
                name = attr[name_field].strip('"')

                # BED format is 0-based
                # GTF format is 1-based
                # NB samtools depth does not include the final base (?) so only subtract 1 from the start
                tss_minus = max(start - distance - 1, 0)
                tss_plus = start + distance

                # only include this TSS if it inhabits a unique region
                if (chr, tss_minus) not in seen:
                    regions.append([chr, tss_minus, tss_plus, strand])

                # add the name to a list
                seen[(chr, tss_minus)].append(name)

            line_no += 1

    return regions, seen


def get_transcript_tss_from_gtf(fn, distance=2500, sources=None, limit=None):
    filter_func = lambda typ, attr: \
        typ == 'transcript' and attr.get('transcript_biotype') == '"protein_coding"'
    name_field = 'transcript_name'
    return get_feature_tss_from_gtf(fn, filter_func, name_field, distance=distance, sources=sources, limit=limit)


def get_gene_tss_from_gtf(fn, distance=2500, sources=None, limit=None):
    filter_func = lambda typ, attr: \
        typ == 'gene' and attr.get('gene_biotype') == '"protein_coding"'
    name_field = 'gene_name'
    return get_feature_tss_from_gtf(fn, filter_func, name_field, distance=distance, sources=sources, limit=limit)


def write_bed_file(region_data, names, fn):
    """
    :param region_data: Iterable of regions, each defined by an iterable of (chrom., start coord, end coord, strand)
    :param names: Dictionary indexed by (chr, start_coord). Values either give the name as a string or a list, in
    which case names are joined with a semicolon
    """
    with open(fn, 'wb') as f:
        c = csv.writer(f, delimiter='\t')
        for row in region_data:
            ix = (row[0], row[1])
            the_names = names[ix]
            if not isinstance(the_names, str):
                the_names = ';'.join(the_names)
            c.writerow(
                [
                    row[0],
                    row[1],
                    row[2],
                    the_names,
                    '.',  # empty score field
                    row[3]
                ]
            )


if __name__ == "__main__":

    distance = 5000 # distance from TSS to include

    fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.release87',
        'gtf',
        'Homo_sapiens.GRCh38.87.gtf.gz'
    )
    outdir = output.unique_output_dir("chipseq_analysis")

    reg, names = get_gene_tss_from_gtf(fn, distance=distance, sources=SOURCES)
    fn_out = os.path.join(outdir, 'gene_tss_pad_%d.bed' % distance)
    write_bed_file(reg, names, fn_out)

    reg, names = get_transcript_tss_from_gtf(fn, distance=distance, sources=SOURCES)
    fn_out = os.path.join(outdir, 'transcript_tss_pad_%d.bed' % distance)
    write_bed_file(reg, names, fn_out)

import gzip
import os
import csv
import collections
import numpy as np
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


def opener(fn, *args, **kwargs):
    if os.path.splitext(fn)[1].lower() == '.gz':
        return gzip.open(fn, *args, **kwargs)
    else:
        return open(fn, *args, **kwargs)


def coverage_reader(cov_fn, exclude=None, include=None):
    if isinstance(exclude, str):
        exclude = {exclude}
    elif exclude is not None:
        exclude = set(list(exclude))

    if isinstance(include, str):
        include = {include}
    elif include is not None:
        include = set(list(include))

    with opener(cov_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        current_ch = None
        current_res = []
        current_count = 0
        for line in c:
            if exclude is not None and line[0] in exclude:
                continue
            if include is not None and line[0] not in include:
                continue
            if line[0] != current_ch:
                if current_count > 0:
                    print "Complete"
                    yield (current_ch, np.array(current_res).astype(np.uint32))
                # new chromosome
                current_ch = line[0]
                current_count = 0
                print "Start reading depth for chromosome %s" % current_ch
                current_res = []
            current_res.append(line[1:])
            current_count += 1
            if (current_count % 1000000) == 0:
                print "Line %d" % current_count
        # yield final result
        print "Reached end of file"
        yield (current_ch, np.array(current_res).astype(np.uint32))


def depth_to_trace(depth, bed_regions):
    """
    :param depth: Numpy array with N + 1 columns, where N is the number of samples.
    The first col corresponds to genomic coordinate. The remaining cols give the coverage of each sample.
    :param bed_regions: List of regions, used to call samtools depth and generate the depth data. Each region is
    specified by 4 parameters: start coord, end coord, name, strand
    """
    bp = depth[:, 0]
    n_trace = len(bed_regions)
    n_sample = depth.shape[1] - 1
    n_depth = depth.shape[0]

    # use these indices to carry out a fast check whether to skip a region
    lower = depth[:, 0].min()
    upper = depth[:, 0].max()

    # infer the trace size from the first traces
    trace_length = max([
        reg[1] - reg[0] for reg in bed_regions
    ])

    this_traces = np.zeros((trace_length, n_sample, n_trace), dtype=np.uint32)

    # some slots may not get filled, so we need to keep a running tally
    n = 0
    skipped = 0
    features = []

    # get the coverage for each sample in each tss region
    for i, (start, end, name, strand) in enumerate(bed_regions):
        if i % 500 == 0:
            print i
        # account for 0-indexing of BED file
        start += 1
        # only need to correct the start coord - the GTF file deliberately includes one extra base at the end

        # if any part of the region is not within the coordinates of the depth array, skip
        if (end > upper) or (start < lower):
            # region is definitely not within the depth array
            skipped += 1
            continue

        # NB: searchsorted will return the array length if the element is not found
        ix0 = np.searchsorted(bp, start)

        ix1 = ix0 + trace_length
        if ix1 > n_depth:
            # this means that we are looking beyond the end of the depth array
            print "Bed region (%d, %d) extends beyond the depth array" % (start, end)
            skipped += 1
            continue
        else:
            the_trace = depth[ix0:ix1, 1:]
            if strand == '-':
                the_trace = the_trace[::-1]
            this_traces[:, :, n] = the_trace
            n += 1
            features.append(name)

    print "Computed %d traces. Skipped %d." % (n, skipped)
    return this_traces[..., :n], features


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

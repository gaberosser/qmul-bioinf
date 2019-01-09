import pandas as pd
import numpy as np
import os
import warnings
import glob
import collections
import pysam
import gffutils
import gzip
import csv
import random
import subprocess
import re
import multiprocessing as mp
from settings import LOCAL_DATA_DIR
from utils import log
_logger = log.get_console_logger()

GTF_SOURCES = ('ensembl', 'havana', 'ensembl_havana')


def _single_lookup(region, gtf_fn):
    db = GtfAnnotation(gtf_fn=gtf_fn, logger=None)
    return list(db.region(region=region))


def multiple_region_lookup(regions, db, njobs=None):
    res = {}
    if njobs != 1:
        pool = mp.Pool(processes=njobs)
        jobs = {}
        for reg in regions:
            jobs[reg] = pool.apply_async(_single_lookup, args=(reg, db.gtf_fn))
        pool.close()
        pool.join()
        for k, v in jobs.items():
            res[k] = v.get()
    else:
        for reg in regions:
            res[reg] = _single_lookup(reg, db.gtf_fn)
    return res


class GtfAnnotation(gffutils.FeatureDB):
    def __init__(self, gtf_fn, db_fn=None, logger=_logger, **kwargs):
        if not os.path.isfile(gtf_fn):
            raise ValueError("Unable to find the specified GTF file %s." % gtf_fn)
        if db_fn is None:
            # place DB file in same location as GTF
            fstem, ext = os.path.splitext(gtf_fn)
            if ext.lower() == 'gz':
                fstem, ext = os.path.splitext(fstem)
            db_fn = "%s.gffutils_db" % fstem

        self.gtf_fn = gtf_fn
        self.db_fn = db_fn
        self.logger = logger

        if os.path.isfile(db_fn):
            self.log_info("Using existing DB file %s", db_fn)
            # logger.info("Using existing DB file %s", db_fn)
        else:
            # logger.info("Creating new DB file %s", db_fn)
            self.log_info("Creating new DB file %s", db_fn)
            gffutils.create_db(gtf_fn, db_fn, disable_infer_genes=True, disable_infer_transcripts=True)

        super(GtfAnnotation, self).__init__(self.db_fn, **kwargs)

    def log_info(self, *args):
        if self.logger is not None:
            self.logger.info(*args)


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


def get_overlapping_paired_reads(bam_fn, chrom):
    """
    Iterate over the supplied BAM file (which must be sorted for this to work), finding paired reads within the
    specified chromosome that overlap.
    This is a bit niche, but useful for debugging double-counting issues, etc.
    """
    s = pysam.AlignmentFile(bam_fn, 'rb')
    rd_seen = set()
    for rd in s.fetch(chrom):
        if rd.is_unmapped:
            continue
        if rd in rd_seen:
            continue
        if rd.is_read1:
            rd1 = rd
            rd2 = s.mate(rd)
            rd_seen.add(rd2)
        else:
            rd2 = rd
            rd1 = s.mate(rd)
            rd_seen.add(rd1)

        if rd1.is_reverse == rd2.is_reverse:
            # this is strange, skip
            continue

        if rd1.is_reverse:
            rd_f = rd2
            rd_r = rd1
        else:
            rd_f = rd1
            rd_r = rd2

        # check overlap
        if rd_r.pos <= rd_f.pos + rd_f.alen:
            yield (rd_f, rd_r)


def random_sampling_alignments(
        bam_fn,
        p=0.1,
        proper_pair=False,
        min_qual=None,
        include_unmapped=False
):
    """
    Iterate over reads at random intervals from the supplied file
    :param bam_fn: BAM, SAM or similar alignment file
    :param p: Probability of returning a read
    :return:
    """
    s = pysam.AlignmentFile(bam_fn, 'rb')
    for rd in s:
        if rd.is_unmapped and not include_unmapped:
            continue
        if proper_pair and not rd.is_proper_pair:
            continue
        if min_qual is not None and rd.mapq < min_qual:
            continue
        if random.random() < p:
            yield rd


def samtools_random_sampling_bam(
    bam_fn,
    frac_reads=None,
    est_n_reads=None,
    max_n_reads=None,
    seed=1,
    samtools_args=tuple()
):
    """
    Iterate over reads at random intervals from the supplied file
    :param bam_fn: BAM, SAM or similar alignment file
    :param frac_reads: The fraction of reads to sample
    :param est_n_reads: The (approx) number of reads to sample. Returning the exact number is *not* guaranteed.
    :param max_n_reads: If supplied, stop the output as soon as this number is reached.
    :param seed: The random seed to use
    :param samtools_args: Passed directly to samtools, for example:
    ['-q', 10] to include only alignments with quality >= 10
    :return: Read generator
    """
    if frac_reads is not None and est_n_reads is not None:
        raise ValueError("Must specify one of frac_reads OR est_n_reads, but not both.")

    if frac_reads is None and est_n_reads is None:
        raise ValueError("Must specify one of frac_reads OR est_n_reads, but not both.")

    if est_n_reads is not None:
        est_total = estimate_number_of_bam_reads(bam_fn)
        _logger.info("Estimated line count in bam file is %d" % est_total)
        frac_reads = est_n_reads / float(est_total)

    frac_as_pct = frac_reads * 100.
    if frac_as_pct < 1e-6:
        frac_as_pct = 1e-6
    pct_as_str = '{:09.6f}'.format(frac_as_pct).replace('.', '')

    cmd = (
        "samtools", "view",
        "-s", "%d.%s" % (seed, pct_as_str)
      ) + tuple([str(t) for t in samtools_args]) + (bam_fn,)
    _logger.info(' '.join(cmd))
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
    )

    for i, line in enumerate(p.stdout):
        if max_n_reads is not None and (i + 1) > max_n_reads:
            raise StopIteration
        yield line


def get_mean_insert_size_pe_subsample(
        bam_fn,
        p=None,
        n=100000,
        proper_pair=True,
        min_qual=10,
):
    """
    Estimate the mean insert length (length of sequence between adapters) by subsampling the supplied file.
    Optionally apply filters by setting them in kwargs (passed to random_sampling_alignments)
    NB we must limit ourselves to alignments on one strand to avoid double counting and avoid negative lengths.
    :param bam_fn:
    :param n: Stop reading file after this many reads have been observed
    :param p: Probability of outputting a single read. Automatically selected if absent.
    :param proper_pair, min_qual: Sensible values set.
    :return:
    """
    if p is None:
        # estimate number of forward reads
        estimated_rd_count = estimate_number_of_bam_reads(bam_fn) / 2.
        # we don't know what proportion of reads will be filtered, so let's be conservative and ensure we don't
        # run to the end of the file too soon.
        p = n / float(estimated_rd_count) * 10.
    it = random_sampling_alignments(bam_fn=bam_fn, p=p, proper_pair=proper_pair, min_qual=min_qual)
    curr = 0.
    count = 0
    for rd in it:
        if not rd.is_reverse:
            curr += rd.tlen
            count += 1
            if count == n:
                break
    return curr / float(n)


def estimate_number_of_bam_reads(fn, bytes=10000000):
    """
    Estimate the number of reads in a BAM file based on the first N bytes (or exactly, if an idex is available).
    This is an underestimation, as the header is included (but should be a minor factor).
    Adapted from https://www.biostars.org/p/1890/
    :param fn:
    :param bytes: Number of bytes to use to estimate [100Mb]
    :return:
    """
    DEVNULL = open(os.devnull, 'wb')
    bai_fn = "%s.bai" % fn
    # if a BAM index is available, this is both faster and perfectly correct.
    if os.path.isfile(bai_fn):
        cmd = "samtools idxstats {fn} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}'".format(
            fn=fn
        )
        return int(
            subprocess.check_output(
                cmd,
                stderr=DEVNULL,
                shell=True
            )
        )

    total_size = int(
        subprocess.check_output("ls -ln '" + fn + "' | awk '{print $5}'", stderr=subprocess.STDOUT, shell=True)
    )
    cmd = "head -c {bytes_used} {fn} | samtools view - | wc -l".format(
        bytes_used = bytes, fn = fn
    )
    reads_in_subsample = int(subprocess.check_output(
        cmd,
        stderr=DEVNULL,
        shell=True
    ))
    return int((reads_in_subsample / float(bytes)) * total_size)


def get_feature_sequence_lengths(fa_file, features=None):
    from Bio import SeqIO
    feat_lens = collections.OrderedDict()

    with open(fa_file, 'rb') as f:
        fa_reader = SeqIO.parse(f, 'fasta')
        for feat in fa_reader:
            the_id = feat.id
            if features is not None and the_id not in features:
                continue
            feat_lens[the_id] = len(feat.seq)

    return feat_lens


def write_bed_file(region_data, fn):
    """
    :param region_data: Dict with values defining regions as an iterable of (chrom., start coord, end coord, strand)
    Keys are the unique names for each region
    :param fn: Filename to write BED file to
    """
    with open(fn, 'wb') as f:
        c = csv.writer(f, delimiter='\t')
        for name, row in region_data.items():
            if not isinstance(name, str):
                name = ';'.join(name)
            c.writerow(
                [
                    row[0],
                    row[1],
                    row[2],
                    name,
                    '.',  # empty score field
                    row[3]
                ]
            )


def motif_locations(fa_file, motif='CG', features=None):
    """
    Find the locations of all occurrences of the supplied motif (default: CpG sites).
    :param fa_file:
    :param motif: This is used to search for the motif. Can contain regex characters. NB case sensitive!
    :param features: If supplied, this is an iterable containing the names of features to include
    :return:
    """
    from Bio import SeqIO
    res = collections.OrderedDict()
    feat_lens = collections.OrderedDict()

    with open(fa_file, 'rb') as f:
        fa_reader = SeqIO.parse(f, 'fasta')
        for feat in fa_reader:
            the_id = feat.id
            if features is not None and the_id not in features:
                continue
            the_seq = str(feat.seq)
            it = re.finditer(motif, the_seq)
            start_coords = np.array([t.start() for t in it])
            res[the_id] = start_coords
            feat_lens[the_id] = len(the_seq)
    return res, feat_lens


def cg_content_windowed(fa_file, motif='CG', window_size=20000, features=None):
    """
    Summarise the CpG (or some other motif) content in the supplied fasta sequence
    :param fa_file:
    :param motif: This is used to search for the motif. Can contain regex characters. NB case sensitive!
    :param window_size: The window used to summarise
    :param features: If supplied, this is an iterable containing the names of features to include
    :return: Dictionary, keyed by feature ID (e.g. chromosome). Values are pd.Series, index is start coord, value is
    CpG count.
    """
    start_coords, feat_lens = motif_locations(fa_file, motif=motif, features=features)
    binned_counts = collections.OrderedDict()
    for feat in start_coords.keys():
        edges = range(1, feat_lens[feat] + 1, window_size)
        if edges[-1] != feat_lens[feat]:
            edges.append(feat_lens[feat])
        counts, _ = np.histogram(start_coords[feat], edges)
        binned_counts[feat] = pd.Series(counts, index=edges[:-1])

    return binned_counts, feat_lens


def gc_fraction_from_sequence(seq):
    return len(re.findall(r'[GC]', seq)) / float(len(seq))


def estimate_gc_content_from_bam(
    bam_fn,
    frac_reads=0.001,
    est_n_reads=None,
    proper_pair=True,
    min_qual=None,
    n_reads_max=None,
    include_unmapped=False,
    paired_end=True,
    extra_args=tuple()
):
    """
    Estimate the GC content in a BAM file by sampling randomly
    :param bam_fn:
    :param frac_reads: The fraction of reads to use when estimating
    :param proper_pair:
    :param min_qual:
    :return: Generator that outputs %GC for reads as they are encountered
    """
    # setup flags where required
    gc_frac_args = []

    if paired_end:
        if include_unmapped:
            gc_frac_args.extend(['-f', 1]) # paired read
        elif proper_pair:
            gc_frac_args.extend(['-f', 3])  # paired read, mapped in proper pair
        else:
            gc_frac_args.extend(['-f', 3, '-F', 4])  # paired read, NOT unmapped

    else:
        if not include_unmapped:
            gc_frac_args.extend(['-F', 4])  # this read is NOT unmapped

    if min_qual is not None:
        gc_frac_args.extend(['-q', min_qual])

    it = samtools_random_sampling_bam(
        bam_fn,
        frac_reads=frac_reads,
        est_n_reads=est_n_reads,
        max_n_reads=n_reads_max,
        samtools_args=tuple(gc_frac_args) + extra_args
    )
    for rd in it:
        seq = rd.split('\t')[9]
        yield gc_fraction_from_sequence(seq)

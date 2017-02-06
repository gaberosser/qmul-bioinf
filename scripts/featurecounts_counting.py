#!/usr/bin/env python
import subprocess
import os
import sys
import re
import csv
from gzip import GzipFile
from log import get_file_logger, get_console_logger
logger = get_console_logger(__name__)


FEATURECOUNTS_CMD = 'featureCounts'



if __name__ == "__main__":
    """
    Usage: featurecounts_counting.py path_to_bams path_to_GTF path_to_output passed_to_featureCounts
    If no parameters are aupplied to pass to featureCounts, a default set are used that work for
    Illumina paired end RNA-Seq.
    """
    read_dir = sys.argv[1]
    ref_fn = sys.argv[2]
    out_dir = sys.argv[3]
    fc_args = sys.argv[4:]

    if len(fc_args) == 0:
        fc_args = [
            '-p',  # paired end
            '-B',  # only count read pairs that are both successfully aligned
            # '-P',  # Check validity of paired-end distance
            # '-d', '50',  # Minimum template length
            # '-D', '600',  # Maximum template length
            '-s', '2',  # Reverse-stranded
            '-T', '12',  # Parallel threads
        ]


    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'featureCounts.log')

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    outfile = os.path.join(out_dir, 'counts.txt')
    logger = get_file_logger(__name__, log_fn)

    # BAM file discovery
    rr = re.compile(r'\.bam$', flags=re.IGNORECASE)
    # flist = [os.path.join(read_dir, t) for t in os.listdir(read_dir) if re.search(rr, t)]
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]

    # the column names in the output file are the FULL PATH to the BAMs
    # this looks much nicer if we chdir now and use the filename only
    os.chdir(read_dir)
    cmd = [
        FEATURECOUNTS_CMD,
        '-a', ref_fn,
        '-o', outfile,
    ] + fc_args + flist
    stderr_fn = os.path.join(out_dir, 'featureCounts.stderr')
    logger.info("Running featureCounts, piping stderr > %s", stderr_fn)
    logger.info("%s", ' '.join(cmd))
    with open(stderr_fn, 'wb') as stderr:
        subprocess.call(cmd, stderr=stderr)

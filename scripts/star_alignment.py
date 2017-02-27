#!/usr/bin/env python
import os
import re
import subprocess
import sys

sys.path.append(os.path.dirname(__file__) + '/..')

from gzip import GzipFile

import numpy as np
from Bio import SeqIO

from utils.log import get_file_logger, get_console_logger

logger = get_console_logger(__name__)

STAR_CMD = 'STAR'


if __name__ == "__main__":
    """
    Usage: star_alignment.py path_to_reads path_to_genomeindex path_to_output passed_to_STAR
    """
    # sys.path.append(os.get)

    read_dir = sys.argv[1]
    ref_fn = sys.argv[2]
    out_dir = sys.argv[3]
    star_args = sys.argv[4:]

    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'star_alignment')

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = get_file_logger(__name__, log_fn)
    
    # fastq.gz file discovery
    rr = re.compile(r'\.fastq(\.gz)?$', flags=re.IGNORECASE)
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
    # check for existing output and identify pairs of files
    fl = {}
    for t in flist:
        is_gz = t[-2:].lower() == 'gz'
        base = re.sub(r'_[12]\.fastq(\.gz)?', '', t)
        read_num = int(re.sub(r'\.fastq(\.gz)?', '', t)[-1])
        out_prefix = os.path.join(out_dir, base)
        # if SAM output file exists, log warning and skip
        if os.path.isfile("%sLog.final.out" % out_prefix):
            logger.warn("File already exists: %s. Skipping.", base)
            continue

        if base not in fl:
            fl[base] = {
                read_num: os.path.join(read_dir, t),
                'out_prefix': out_prefix,
                'gzipped': is_gz,
            }
        else:
            fl[base][read_num] = os.path.join(read_dir, t)

    # log the filelist
    logger.info("Found %d fastq pairs: %s.", len(fl), ', '.join(fl.keys()))

    # start the jobs
    for base, d in fl.iteritems():
        # create a temporary stderr file (StringIO doesn't work)
        stderr_fn = os.path.join(out_dir, '%s.stderr' % base)
        cmd = [STAR_CMD] + star_args + [
            '--genomeDir',
            ref_fn,
            '--readFilesIn %s %s' % (d[1], d[2]),
            '--outFileNamePrefix',
            d['out_prefix'],
            '--outSAMtype BAM Unsorted SortedByCoordinate',  # needs to be included like this to avoid quoting spaces
            '--outSAMstrandField',
            'intronMotif',
            '--quantMode',
            'GeneCounts',
        ]
        if d['gzipped']:
            cmd += ['--readFilesCommand', 'zcat']
        logger.info("Running alignment on sample %s: %s", base, ' '.join(cmd))
        logger.info("Sending stderr to temporary file %s", stderr_fn)
        with open(stderr_fn, 'wb') as stderr:
            try:
                subprocess.call(cmd, stderr=stderr)
            except Exception as exc:
                logger.exception("Alignment failed.")
        logger.info("Alignment complete")


#!/usr/bin/env python
import os
import re
import subprocess
import sys

sys.path.append(os.path.dirname(__file__) + '/../../')

from gzip import GzipFile

import numpy as np
from Bio import SeqIO

from utils.log import get_file_logger, get_console_logger

logger = get_console_logger(__name__)

QMAP_CMD = 'qualimap'

if __name__ == "__main__":
    """
    Usage: qualimap_rnaseq_bam_qc.py path_to_bams path_to_gtf path_to_output passed_to_function
    """
    # sys.path.append(os.get)

    read_dir = sys.argv[1]
    ref_fn = sys.argv[2]
    out_dir = sys.argv[3]
    qmap_args = sys.argv[4:]

    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'qualimap')

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = get_file_logger(__name__, log_fn)

    # fastq.gz file discovery
    rr = re.compile(r'\.bam$', flags=re.IGNORECASE)
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
    # check for existing output and identify files
    fl = {}
    for t in flist:
        base = re.sub(r'\.bam$', '', t)
        out_subdir = os.path.join(out_dir, base)
        # if output folder exists, log warning and skip
        if os.path.isdir(out_subdir):
            logger.warn("Folder already exists: %s. Skipping.", out_subdir)
            continue
        else:
            os.makedirs(out_subdir)

        fl[base] = {
            'outdir': out_subdir,
            'bam': os.path.join(read_dir, t),
        }

    # log the filelist
    logger.info("Found %d BAM files: %s.", len(fl), ', '.join(fl.keys()))

    # start the jobs
    for base, d in fl.iteritems():
        # create a temporary stderr file (StringIO doesn't work)
        stderr_fn = os.path.join(out_dir, '%s.stderr' % base)
        cmd = [QMAP_CMD] + ['rnaseq'] + qmap_args + [
            '-bam', d['bam'],
            '-gtf', ref_fn,
            '-outdir', d['outdir'],
        ]

        logger.info("Running qualimap on sample %s: %s", base, ' '.join(cmd))
        logger.info("Sending stderr to temporary file %s", stderr_fn)
        with open(stderr_fn, 'wb') as stderr:
            try:
                subprocess.call(cmd, stderr=stderr)
            except Exception as exc:
                logger.exception("Qualimap failed.")
        logger.info("Qualimap complete")

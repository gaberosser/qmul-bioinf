#!/usr/bin/env python
import os
import re
import subprocess
import sys

sys.path.append(os.path.dirname(__file__) + '/..')
from utils.log import get_file_logger, get_console_logger


logger = get_console_logger(__name__)



if __name__ == "__main__":
    """
    Usage: apply_fastqc.py path_to_reads output_dir passed_to_fastqc
    """
    # sys.path.append(os.get)

    read_dir = sys.argv[1]
    out_dir = sys.argv[2]

    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'fastqc')

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
        base = re.sub(r'_[12]\.fastq(\.gz)?', '', t)
        read_num = int(re.sub(r'\.fastq(\.gz)?', '', t)[-1])
        if base not in fl:
            fl[base] = {
                read_num: os.path.join(read_dir, t),
            }
        else:
            fl[base][read_num] = os.path.join(read_dir, t)

    # log the filelist
    logger.info("Found %d fastq pairs: %s.", len(fl), ', '.join(fl.keys()))

    # start the jobs
    for base, d in fl.iteritems():
        CMD = ["fastqc", "-o", out_dir, '--noextract', d[1], d[2]]
        logger.info("Running fastqc on sample %s: %s", base, ' '.join(CMD))
        try:
            subprocess.call(CMD)
        except Exception as exc:
            logger.exception("Failed.")
        logger.info("Fastqc complete")
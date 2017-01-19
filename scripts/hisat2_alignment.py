#!/usr/bin/env python
import subprocess
import os
import sys
import re
import csv
from gzip import GzipFile
from Bio import SeqIO
from log import get_file_logger, get_console_logger
import numpy as np
from StringIO import StringIO
logger = get_console_logger(__name__)

HISAT_CMD = 'hisat2'


# TODO: add proper argument parsing, which is especially useful for --threads, which is used in more than one place.


def main(read_dir, ref_dir, prealign_dir, *args):
    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(ref_dir):
        raise ValueError("Could not find specified reference directory %s" % ref_dir)
    if not os.path.isfile(prealign_dir):
        raise ValueError("Could not find specified prealignment directory %s" % prealign_dir)


def get_read_length(read_fn, n_read=10000):
    """
    Get the read length from the .fastq.gz file
    :param read_fn: .fastq.gz file
    :return:
    """
    with GzipFile(read_fn, mode='rb') as f:
        h = SeqIO.QualityIO.FastqGeneralIterator(f)
        i = 0
        l = []
        while i < n_read:
            try:
                t = h.next()
                l.append(len(t[1]))
                i += 1
            except StopIteration:
                logger.warning("Requested %d reads but reached the end of the file after %d", n_read, i)
    return int(np.round(np.mean(l)))


if __name__ == "__main__":
    """
    Usage: hisat2_alignment.py path_to_reads path_to_bt2 path_to_output passed_to_hisat2
    """
    read_dir = sys.argv[1]
    ref_fn = sys.argv[2]
    out_dir = sys.argv[3]
    hisat_args = sys.argv[4:]

    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'hisat2_alignment')

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = get_file_logger(__name__, log_fn)
    
    # fastq.gz file discovery
    rr = re.compile(r'\.fastq\.gz$', flags=re.IGNORECASE)
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
    # check for existing output and identify pairs of files
    fl = {}
    for t in flist:
        base = re.sub(r'_[12]\.fastq\.gz', '', t)
        read_num = int(re.sub(r'\.fastq\.gz', '', t)[-1])
        sam_out = os.path.join(out_dir, "%s.sam" % base)
        # if SAM output file exists, log warning and skip
        if os.path.isfile(sam_out):
            logger.warn("File already exists: %s. Skipping.", sam_out)
            continue

        if base not in fl:
            fl[base] = {
                read_num: os.path.join(read_dir, t),
                'sam_out': sam_out,
            }
        else:
            fl[base][read_num] = os.path.join(read_dir, t)

    # log the filelist
    logger.info("Found %d fastq pairs: %s.", len(fl), ', '.join(fl.keys()))

    # start the jobs
    for base, d in fl.iteritems():
        # create a temporary stderr file (StringIO doesn't work)
        stderr_fn = os.path.join(out_dir, '%s.stderr' % base)
        cmd = [HISAT_CMD] + hisat_args + [
            '-x',
            ref_fn,
            '-1',
            d[1],
            '-2',
            d[2],
            '-S',
            d['sam_out']
        ]
        logger.info("Running alignment on sample %s: %s", base, ' '.join(cmd))
        logger.info("Sending stderr to temporary file %s", stderr_fn)
        with open(stderr_fn, 'wb') as stderr:
            try:
                subprocess.call(cmd, stderr=stderr)
            except Exception as exc:
                logger.exception("Alignment failed.")
        logger.info("Alignment complete")

        # convert to BAM and delete the original SAM
        cmd = ['samtools', 'view', '--threads', '8', '-b', d['sam_out']]
        bam_out = re.sub(r'\.sam$', '\.bam', d['sam_out'])
        logger.info("Converting %s -> %s", d['sam_out'], bam_out)
        with open(bam_out, 'wb') as stdout:
            logger.info("%s", ' '.join(cmd))
            try:
                subprocess.call(cmd, stdout=stdout)
            except Exception as exc:
                logger.exception("Conversion failed.")
        logger.info("Conversion complete.")

        # delete SAM file
        cmd = ['rm', d['sam_out']]
        logger.info("Deleting original SAM file %s", d['sam_out'])
        try:
            subprocess.call(cmd)
        except Exception as exc:
            logger.exception("Deletion failed.")

        logger.info("All complete")

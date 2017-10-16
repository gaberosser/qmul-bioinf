#!/usr/bin/env python
import os
import re
import subprocess
import sys
import argparse

sys.path.append(os.path.dirname(__file__) + '/..')

from utils.log import get_file_logger, get_console_logger
from utils.output import unique_output_dir

logger = get_console_logger(__name__)

"""
If we want to use cufflinks after this, need to use --dta-cufflinks to add splicing tags.
Not sure how this affects the alignment, though?
"""

BT2_CMD = 'bowtie2'
SAMTOOLS_CMD = 'samtools'


if __name__ == "__main__":
    """
    Usage: bowtie2_alignment.py path_to_reads path_to_bt2 path_to_output passed_to_hisat2
    """

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.set_defaults(sort=True)

    required.add_argument("-x", "--reference", help="Bowtie2-compiled reference (stem, no extensions).", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        test_out = os.path.join(args.read_dir, 'bt2')
        if not os.path.isdir(test_out):
            os.makedirs(test_out)
            args.out_dir = test_out
        else:
            args.out_dir = unique_output_dir('bt2', root_output_dir=args.read_dir)
        sys.stderr.write("Created output directory %s\n" % args.out_dir)

    read_dir = args.read_dir
    ref_fn = args.reference
    out_dir = args.out_dir

    # put global log output into the output directory
    log_fn = os.path.join(out_dir, 'bowtie2_alignment')

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    logger = get_file_logger(__name__, log_fn)

    # fastq.gz file discovery
    ## TODO: support SE reads, too?
    rr = re.compile(r'\.fastq(\.gz)?$', flags=re.IGNORECASE)
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
    # check for existing output and identify pairs of files
    ignore = set()
    fl = {}
    for t in flist:
        base = re.sub(r'_[12]\.fastq(\.gz)?', '', t)
        if base in ignore:
            continue
        read_num = int(re.sub(r'\.fastq(\.gz)?', '', t)[-1])
        bam_out = os.path.join(out_dir, "%s.bam" % base)
        # if SAM output file exists, log warning and skip
        if os.path.isfile(bam_out):
            logger.warn("File already exists: %s. Skipping.", bam_out)
            # add to ignore list so we'll skip its pair
            ignore.add(base)
            continue

        if base not in fl:
            fl[base] = {
                read_num: os.path.join(read_dir, t),
                'bam_out': bam_out,
            }
        else:
            fl[base][read_num] = os.path.join(read_dir, t)

    # log the filelist
    logger.info("Found %d fastq pairs: %s.", len(fl), ', '.join(fl.keys()))

    # start the jobs
    for base, d in fl.iteritems():
        # create a temporary stderr file (StringIO doesn't work)
        stderr_fn = os.path.join(out_dir, '%s.stderr' % base)
        cmd_tuple = (
            BT2_CMD, '-x', ref_fn, '-1', d[1], '-2', d[2]
        ) + tuple(extra)
        bam_tuple = (
            SAMTOOLS_CMD, 'view', '-bS'
        )

        logger.info("Running alignment on sample %s: %s", base, ' '.join(cmd_tuple))
        logger.info("Sending stderr to temporary file %s", stderr_fn)
        with open(stderr_fn, 'wb') as stderr, open(d['bam_out'], 'wb') as stdout:
            try:
                # subprocess.call(cmd, stderr=stderr)
                bt2_proc = subprocess.Popen(cmd_tuple, stderr=stderr, stdout=subprocess.PIPE)
                st_proc = subprocess.Popen(bam_tuple, stdin=bt2_proc.stdout, stdout=stdout)
                st_proc.wait()
            except Exception as exc:
                logger.exception("Alignment failed.")

        logger.info("Alignment complete")



    logger.info("All complete")

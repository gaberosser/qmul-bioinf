#!/usr/bin/env python
import subprocess
import os
import sys
import csv
from log import get_file_logger

INSERT_SH = 'scripts/bash/insert_size_from_bam.sh'


def main(read_dir, ref_dir, prealign_dir, *args):
    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(ref_dir):
        raise ValueError("Could not find specified reference directory %s" % ref_dir)
    if not os.path.isfile(prealign_dir):
        raise ValueError("Could not find specified prealignment directory %s" % prealign_dir)


def get_read_length(read_fn):
    """
    Get the read length from the .fastgz file
    :param read_fn:
    :return:
    """
    pass



if __name__ == "__main__":
    """
    Usage: hisat2_alignment.py path_to_reads path_to_bt2 path_to_insert_size.csv passed_to_hisat2
    """
    read_dir = sys.argv[1]
    ref_dir = sys.argv[2]
    prealign_dir = sys.argv[3]
    hisat_args = sys.argv[4:]

    # put log output into the output directory
    log_fn = os.path.join(prealign_dir, 'hisat2_alignment')

    logger = get_file_logger(__name__, log_fn)

    if not os.path.isdir(read_dir):
        raise ValueError("Could not find specified read directory %s" % read_dir)
    if not os.path.isdir(ref_dir):
        raise ValueError("Could not find specified reference directory %s" % ref_dir)
    if not os.path.isfile(prealign_dir):
        raise ValueError("Could not find specified prealignment directory %s" % prealign_dir)

    ins_fn = os.path.join(prealign_dir, 'insert_size.csv')
    if os.path.exists(ins_fn):
        logger.info("Using existing insert sizes CSV %s", ins_fn)
    else:
        logger.info("Unable to find insert sizes CSV %s so will create now", ins_fn)
        logger.error("The process of calling this from python is SLOW and I don't know why!")
        with open(ins_fn, 'rb') as f:
            subprocess.call([INSERT_SH, prealign_dir], stdout=)


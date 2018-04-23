from chipseq import loader, annotation, feature_enrichment
from settings import LOCAL_DATA_DIR, RNASEQ_DIR
from utils import output
import os
import collections
import pandas as pd
import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    run_type = 'default'
    N_peak = 100
    window_size = 10001

    sname = '3021'

    peaks = loader.load_macs2_by_patient(sname, run_type=run_type)
    base_dirs = [k.base_dir for k in loader.INPUTS_LOOKUP[sname]]

    # extract N most confident peaks from each, excluding those that are too large
    for c in peaks.meta.index:
        # find the associated bam file - we need this to extract the trace(s)

        the_dat = peaks.data[c]
        the_dat = the_dat[(the_dat.end - the_dat.start) < ((window_size - 1) / 2.)].sort_values(
            '-log10q', ascending=False
        )[:N_peak]

        the_centres = the_dat.start + the_dat.rel_peak_pos
        the_chroms = the_dat.chrom




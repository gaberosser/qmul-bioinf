import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
from settings import OUTPUT_DIR


if __name__ == '__main__':
    base_dir = os.path.join(OUTPUT_DIR, 'dmr_without_classes.3')
    outdir = base_dir
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    required_comparisons = {
        '017': ['hypo', 'hyper'],
        '019': ['hyper'],
        '030': ['hypo', 'hyper'],
        '031': ['hyper'],
        '050': ['hyper'],
        '054': ['hyper'],
        '061': ['hyper'],
        '026': ['hypo', 'hyper'],
        '052': ['hyper'],
    }

    res = {}

    for p in pids:
        for typ in required_comparisons[p]:
            # load patient-specific result; use this for results tbl
            # load full result for comparison; add this in showing enrichment
            # create plot (?) with motif (?)
            pass
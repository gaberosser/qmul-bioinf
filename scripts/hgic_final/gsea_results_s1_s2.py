import os
from settings import HGIC_LOCAL_DIR
import pandas as pd
from utils import output, setops, excel
from plotting import common
import consts
from matplotlib import pyplot as plt, patches, collections, gridspec
import seaborn as sns
import numpy as np


if __name__ == "__main__":
    # load all GSEA results - use both the NSC and GBM excel files. Keep the FDR and the NES (for 'direction')
    # (can use the CSV summarised files for convenience here)
    # Carry out a very similar analysis to that run in ipa_results_s1_s2, including the same plot.
    indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/s0_individual_patients_direct_comparison/gsea/results/raw')
    pass

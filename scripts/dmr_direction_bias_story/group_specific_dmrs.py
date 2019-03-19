from plotting import bar, common, pie, polar
from methylation import loader, dmr, process
import pandas as pd
from stats import nht
from utils import output, setops, genomics, log, dictionary
import multiprocessing as mp
import os
import collections
import pickle
import numpy as np
from scipy import stats, cluster
import matplotlib
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from matplotlib import cm
from sklearn.neighbors import KernelDensity

import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts
from scripts.dmr_direction_bias_story import analyse_dmr_direction_and_distribution as addd

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we seek to identify DMRs that distinguish the hypo- and hypermethylated groups.

This is initially carried out directly using the GIC-iNSC results.

We then query these DMRs to determine whether they have different location distributions.

Possible directions:
- Run a direct DM comparison to identify DMRs without the need for a comparator
- Look for shared but discordant between the two groups
- Link to DE and look for a change in the concordance relative to 'all DMRs' (try this with TSS only)
"""
from load_data import rnaseq_data
import os
from settings import DATA_DIR_NON_GIT, DATA_DIR_NON_GIT2
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    # load our iNSC, iPSC
    pids = rnaseq_data.PATIENT_LOOKUP_CC.keys()

    obj = rnaseq_data.load_by_patient(pids, include_control=False)
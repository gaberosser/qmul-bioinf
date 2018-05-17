from methylation import loader
import pandas as pd
import os
# import glob
# import re
# from settings import GIT_LFS_DATA_DIR, DATA_DIR_NON_GIT
#
# from utils import rinterface, log
# from rpy2 import robjects
# from rpy2.robjects import pandas2ri, r
# pandas2ri.activate()
# logger = log.get_console_logger(__name__)


idat_dir = '/home/gabriel/data/methylation/2016-06-10_brandner/idat/'
meta_fn = '/home/gabriel/data/methylation/2016-06-10_brandner/sources.csv'
outdir = '/home/gabriel/data/methylation/2016-06-10_brandner/beta'

if not os.path.exists(outdir):
    os.makedirs(outdir)

func = loader.idat_files_to_beta
beta = func(
    idat_dir,
    meta_fn,
    outdir=outdir,
    samples=None,
    array_type='EPIC',
    name_col='Sample_Name',
    annotation=None
)
import os
from settings import GIT_LFS_DATA_DIR
import pandas as pd
import numpy as np

ANNOTATION_DIR = os.path.join(GIT_LFS_DATA_DIR, 'microarray_annotation')


def load_from_r_format(library):
    """
    Load annotation data from a TSV outputted by R.
    :param library: String giving the name of the library required.
    :return: Pandas dataframe
    """
    ## TODO: could add the R function here directly, allowing computing of annotations ON THE FLY, rather than a two-
    ## step process where we generate in R first. Just a thought.
    if library == "mogene10sttranscriptcluster.db":
        infile = os.path.join(ANNOTATION_DIR, 'mogene10sttranscriptcluster.tsv')
    else:
        raise ValueError("Unsupported microarray library %s", library)

    res = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    res = res.replace("NA", np.nan)

    return res
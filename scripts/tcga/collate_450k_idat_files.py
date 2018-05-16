import os
import json
import pandas as pd
from settings import DATA_DIR_NON_GIT

"""
Scenario:

We have downloaded a number of IDAT files from TCGA (legacy portal) using the `gdc-client` (or similar).
We have access to the manifesto (used to download files) and the 'cart' metadata JSON
IDAT files stored in one directory per file, so we need to reunite pairs.
We'll also associate each with the case UUID and submitter ID and save this as new metadata.
"""

indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '')
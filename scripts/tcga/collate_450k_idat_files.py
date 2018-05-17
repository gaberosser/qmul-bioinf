import os
import shutil
import json
import pandas as pd
from settings import DATA_DIR_NON_GIT
from collate_annotate2 import parse_meta, get_meta
from utils import output

"""
Scenario:

We have downloaded a number of IDAT files from TCGA (legacy portal) using the `gdc-client` (or similar).
We have access to the manifesto (used to download files) and the 'cart' metadata JSON
IDAT files stored in one directory per file, so we need to reunite pairs.
We'll also associate each with the case UUID and submitter ID and save this as new metadata.
"""

indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour', 'idat')
outdir = output.unique_output_dir('tcga_450k')

meta_fn = get_meta(indir)
res = parse_meta(meta_fn)
meta_export = {}

for k, v in res.items():
    if len(v['file_name']) != 2:
        raise AttributeError("Require 2 idat files, found %d." % len(v['file_name']))
    if len(v['file_id']) != 2:
        raise AttributeError("Require 2 input folders, found %d." % len(v['file_id']))
    f1 = os.path.join(indir, v['file_id'][0], v['file_name'][0])
    f2 = os.path.join(indir, v['file_id'][1], v['file_name'][1])
    out_subdir = os.path.join(outdir, v['case_id'])
    if not os.path.exists(out_subdir):
        os.makedirs(out_subdir)
    for d, fn in zip(v['file_id'], v['file_name']):
        if '_Red.idat' in fn:
            fn_out = "%s_Red.idat" % k
        elif '_Grn.idat' in fn:
            fn_out = "%s_Grn.idat" % k
        else:
            raise AttributeError("Unexpected input filename: %s (case %s)" % (fn, k))
        shutil.copyfile(os.path.join(indir, d, fn), os.path.join(out_subdir, fn_out))
    meta_export[k] = {
        'a260_280_ratio': ';'.join([str(t) for t in v['a260_280_ratio']]),
        'case_id': v['case_id'],
        'filestem': k,
        'sample_type': v['sample_type'],
    }
    for t in ('necrosis', 'normal_cells', 'stromal_cells', 'tumor_cells', 'tumor_nuclei'):
        meta_export[k]["percent_%s" % t] = ';'.join([str(x) for x in v["percent_%s" % t]])
meta_export = pd.DataFrame.from_dict(meta_export).transpose()
meta_export.to_csv(os.path.join(outdir, "sources.csv"))

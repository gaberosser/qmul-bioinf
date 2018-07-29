from settings import GIT_LFS_DATA_DIR
import pandas as pd
from methylation import loader
from utils import output, excel, dictionary, setops
from glob import glob
import os


if __name__ == '__main__':
    anno = loader.load_illumina_methylationepic_annotation(split_genes=False)

    # 1. Annotate DMPs and re-export to Excel

    dmp_fns = glob(os.path.join(GIT_LFS_DATA_DIR, 'mb_dmp', '*.xlsx'))
    print "Found %d relevant input (DMP) files: %s" % (len(dmp_fns), ', '.join(dmp_fns))
    outdir = output.unique_output_dir("mb_dmps")
    res = {}
    anno_cols = ['chr', 'strand', 'loc', 'gene', 'relation']

    for fn in dmp_fns:
        base = os.path.splitext(os.path.basename(fn))[0]
        res[base] = {}
        dat = pd.read_excel(fn, sheet_name=None)
        for cmp, df in dat.items():
            this_anno = anno.reindex(df.index)
            this_res = []
            # gene - relation unique values
            for ix, row in this_anno.iterrows():
                df_row = df.loc[ix].tolist()
                if pd.isna(row.UCSC_RefGene_Name):
                    gr = [('', '')]
                else:
                    g = row.UCSC_RefGene_Name.split(';')
                    r = row.UCSC_RefGene_Group.split(';')
                    gr = sorted(set(zip(g, r)))
                for t in gr:
                    this_res.append(
                        [ix] + df_row + [row.CHR, row.Strand, row.MAPINFO, t[0], t[1]]
                    )
            res[base][cmp] = pd.DataFrame(this_res, columns=['probe_id'] + df.columns.tolist() + anno_cols)

        # save to Excel
        out_fn = os.path.join(outdir, os.path.basename(fn))
        excel.pandas_to_excel(res[base], out_fn, write_index=False)

    # 2.1 Look for common DMPs

    # 'Refold' the previous results dictionary
    res_flat = dictionary.nested_dict_to_flat(res)
    res_by_cmp = dict([(k[::-1], v) for k, v in res_flat.items()])
    res_by_cmp = dictionary.flat_dict_to_nested(res_by_cmp)
    common_dmps = {}

    for cmp, d in res_by_cmp.items():
        common_dmps[cmp] = setops.reduce_intersection(*[t.probe_id for t in d.values()])

    # 2.2 Look for common genes

    # 'Refold' the previous results dictionary
    common_genes = {}

    for cmp, d in res_by_cmp.items():
        common_genes[cmp] = setops.reduce_intersection(*[t.gene for t in d.values()])

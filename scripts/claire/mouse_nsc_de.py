from load_data import rnaseq_data
from rnaseq import differential_expression
from utils.output import unique_output_dir
import references
import os
import pandas as pd


if __name__ == '__main__':
    outdir = unique_output_dir("mouse_NSC_DE", reuse_empty=True)

    lfc = 1
    fdr = 0.01

    dat = rnaseq_data.mouse_nsc_star()
    the_groups = pd.Series('eNSC', index=dat.columns)
    the_groups[dat.columns.str.contains('mDura')] = 'iNSC'
    the_contrast = 'iNSC - eNSC'

    res = differential_expression.edger_glmqlfit(dat, the_groups, the_contrast, lfc=lfc, fdr=fdr)

    gs = references.ensembl_to_gene_symbol(res.index, tax_id=10090)
    res.insert(0, "Gene symbol", gs.values)

    out_fn = os.path.join(outdir, "de_insc-ensc.xlsx")

    xl_writer = pd.ExcelWriter(out_fn)
    res.to_excel(xl_writer, "DE results", index=True)
    xl_writer.save()
from load_data import rnaseq_data
from rnaseq import general, gsea
from utils.output import unique_output_dir
import references
import os
import pandas as pd


if __name__ == '__main__':
    outdir = unique_output_dir("mouse_gsea_files", reuse_empty=True)
    dat = rnaseq_data.mouse_nsc_salmon()
    dat = general.ensembl_transcript_quant_to_gene(dat, tax_id=10090)
    idx = dat.columns.str.contains(r'eNSC[0-9]med') | dat.columns.str.contains(r'mDura[0-9AN]*human')
    dat = dat.loc[:, idx]
    the_groups = pd.Series('eNSC', index=dat.columns)
    the_groups[dat.columns.str.contains('mDura')] = 'iNSC'

    # now switch from Ensembl to gene symbol and capitalize (why?)
    gs = references.ensembl_to_gene_symbol(dat.index, tax_id=10090)
    gs = gs.str.upper()
    gs = gs.loc[~gs.index.duplicated()]

    gs.dropna(inplace=True)
    dat = dat.loc[gs.index]
    dat.index = gs

    # this leaves some duplicate values
    # we'll take the average
    dupe_idx = dat.index[dat.index.duplicated()]
    dupe_map = dat.index.isin(dupe_idx)
    dupes = dat.loc[dupe_map]
    dat = dat.loc[~dupe_map]
    dupes_mean = dupes.groupby(dupes.index).mean()
    dat = dat.append(dupes_mean)

    out_fn = os.path.join(outdir, "cv_mouse_nsc_salmon_tpm.{ext}")
    gsea.data_to_gct(dat, out_fn.format(ext='gct'))
    gsea.phenotypes_to_cls(the_groups, out_fn.format(ext='cls'))
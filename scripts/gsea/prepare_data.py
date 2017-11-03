from load_data import rnaseq_data
from utils.output import unique_output_dir
from rnaseq import gsea
import os
import references


if __name__ == '__main__':
    outdir = unique_output_dir("gsea_data")
    # load all data
    pids = ['017', '050', '054', '061', '018', '052', '044']
    obj = rnaseq_data.load_by_patient(pids, include_control=False)
    dat = obj.get_fpkm()
    # convert to gene symbols
    idx = references.ensembl_to_gene_symbol(dat.index).dropna()
    idx = idx.loc[~idx.index.duplicated()]
    dat = dat.loc[idx.index]
    dat.index = idx

    # write individual patient data
    for pid in pids:
        the_idx = obj.meta.index[obj.meta.index.str.contains(pid)]
        the_dat = dat.loc[:, the_idx]
        the_classes = obj.meta.loc[the_idx, 'type'].values
        out_fn = os.path.join(outdir, "%s_fpkm.{ext}" % pid)
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

    # write grouped RTK II data
    rtkii_pids = ['017', '050', '054', '061']
    the_idx = obj.meta.index.str.contains(r'|'.join(rtkii_pids))
    the_dat = dat.loc[:, the_idx]
    the_classes = obj.meta.loc[the_idx, 'type'].values
    out_fn = os.path.join(outdir, "rtkii_fpkm.{ext}")
    gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
    gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

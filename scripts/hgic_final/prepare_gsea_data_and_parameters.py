from rnaseq import loader, general
from rnaseq import gsea
import os
import pandas as pd
from utils import output
import consts


def ens_index_to_gene_symbol(df):
    general.add_gene_symbols_to_ensembl_data(df)
    tmp = df['Gene Symbol'].dropna()
    df = df.loc[tmp.index]
    df.set_index('Gene Symbol', inplace=True)
    return df


if __name__ == '__main__':

    # load all data
    pids = consts.PIDS
    units = 'tpm'

    outdir = output.unique_output_dir()

    out_subdir = os.path.join(outdir, units)
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
        print "Created output subdirectory %s" % out_subdir

    source_by_units = {
        'tpm': 'salmon',
        'counts': 'star',
        'fpkm': 'star/cufflinks'
    }

    obj = loader.load_by_patient(pids, source=source_by_units[units], include_control=True)

    # set gibco aside
    dat_gibco = obj.data.loc[:, obj.data.columns.str.contains('GIBCO')]
    dat_gibco = ens_index_to_gene_symbol(dat_gibco)

    # drop any cell types other than GBM and iNSC
    ix = obj.meta['type'].isin(['GBM', 'iNSC'])
    # drop unneeded GBM061 samples
    ix = ix & (~obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    obj.filter_samples(ix)

    # convert to gene symbols
    dat = ens_index_to_gene_symbol(obj.data)

    # load reference dataset(s)
    ref_obj = loader.load_references('GSE61794', strandedness='u', source=source_by_units[units])
    ix = ref_obj.meta.index.str.contains('NSC')
    ref_obj.filter_samples(ix)

    # convert to gene symbols
    dat_h9 = ens_index_to_gene_symbol(ref_obj.data)

    # write single patient syngeneic comparison data
    for pid in pids:
        the_idx = dat.columns.str.contains(pid)
        the_dat = dat.loc[:, the_idx]
        the_classes = pd.Series('GBM', index=the_dat.columns)
        the_classes.loc[the_classes.index.str.contains('NSC')] = 'iNSC'
        out_fn = os.path.join(out_subdir, "%s.{ext}" % (pid))
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

    # params
    for pid in pids:
        out_fn = os.path.join(out_subdir, "%s.{ext}" % (pid))
        gsea.create_gsea_params_file(out_fn.format(ext='params'), rpt_label=pid)

    # write single patient reference comparison data
    refs = {
        'gibco_nsc': dat_gibco,
        'h9_nsc': dat_h9
    }
    for pid in pids:
        ix = dat.columns.str.contains(pid) & dat.columns.str.contains('GBM')
        n_ix = ix.sum()
        for rnm, rd in refs.items():
            the_dat = pd.concat(
                (dat.loc[:, ix], rd),
                axis=1
            )
            the_classes = pd.Series(['GBM'] * n_ix + ['NSC'] * rd.shape[1], index=the_dat.columns)
            out_fn = os.path.join(out_subdir, "%s_%s.{ext}" % (pid, rnm))
            gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
            gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

    # params
    for pid in pids:
        for rnm, rd in refs.items():
            out_fn = os.path.join(out_subdir, "%s_%s.{ext}" % (pid, rnm))
            gsea.create_gsea_params_file(out_fn.format(ext='params'), rpt_label="%s_%s" % (pid, rnm))

    """
    At this point, we need to run GSEA. The script to do this is in gsea_run_all_comparisons.sh
    """
    
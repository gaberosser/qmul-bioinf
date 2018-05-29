from rnaseq import loader, general
from utils.output import unique_output_dir
from rnaseq import gsea
import os
import references
import pandas as pd
from settings import LOCAL_DATA_DIR


""" 
Running GSEA from the command line:
java -Xmx16192m -cp /opt/gsea/gsea-3.0.jar xtools.gsea.Gsea -res /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/tpm_salmon/rtkii.gct -cls /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/tpm_salmon/rtkii.cls#GBM_versus_iNSC -gmx /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/msigdb.v6.1.symbols.gmt -collapse false -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label RTK_II -metric Signal2Noise -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 1000 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /home/gabriel/gsea_home/output/rtk_II -gui false

"""


if __name__ == '__main__':
    outdir = unique_output_dir("gsea_data", reuse_empty=True)
    # load all data
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    units = 'tpm'

    source_by_units = {
        'tpm': 'salmon',
        'counts': 'star',
        'fpkm': 'star/cufflinks'
    }

    obj = loader.load_by_patient(pids, source=source_by_units[units], include_control=False)

    # drop any cell types other than GBM and iNSC
    ix = obj.meta['type'].isin(['GBM', 'iNSC'])
    ix = ix & (~obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    dat = obj.data.loc[:, ix.index[ix]]

    # convert to gene symbols
    general.add_gene_symbols_to_ensembl_data(dat)
    tmp = dat['Gene Symbol'].dropna()
    dat = dat.loc[tmp.index]
    dat.set_index('Gene Symbol', inplace=True)

    # write individual patient data
    for pid in pids:
        # the_idx = obj.meta.index[obj.meta.index.str.contains(pid)]
        the_idx = dat.columns.str.contains(pid)
        the_dat = dat.loc[:, the_idx]
        # the_classes = obj.meta.loc[the_idx, 'type'].values
        the_classes = pd.Series('GBM', index=the_dat.columns)
        the_classes.loc[the_classes.index.str.contains('NSC')] = 'iNSC'
        out_fn = os.path.join(outdir, "%s_%s.{ext}" % (pid, units))
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

    # load reference dataset(s)
    ext_ref_name = 'GSE61794'
    ext_ref_strandedness = 'u'

    ref_obj = loader.load_references(ext_ref_name, strandedness=ext_ref_strandedness)
    ref_obj.meta = ref_obj.meta.loc[ref_obj.meta.index.str.contains('NSC')]
    ref_obj.data = ref_obj.data.loc[:, ref_obj.meta.index]

    # TODO: update from here - export S2 comparisons too

    if False:
        # write grouped RTK II data
        rtkii_pids = ['017', '050', '054']
        the_idx = dat.columns.str.contains(r'|'.join(rtkii_pids))
        the_dat = dat.loc[:, the_idx]
        # the_classes = obj.meta.loc[the_idx, 'type'].values
        the_classes = pd.Series('GBM', index=the_dat.columns)
        the_classes.loc[the_classes.index.str.contains('NSC')] = 'iNSC'

        out_fn = os.path.join(outdir, "rtkii_fpkm.{ext}")
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

        # write grouped RTK I data
        rtki_pids = ['019', '030', '031']
        the_idx = dat.columns.str.contains(r'|'.join(rtki_pids))
        the_dat = dat.loc[:, the_idx]
        # the_classes = obj.meta.loc[the_idx, 'type'].values
        the_classes = pd.Series('GBM', index=the_dat.columns)
        the_classes.loc[the_classes.index.str.contains('NSC')] = 'iNSC'

        out_fn = os.path.join(outdir, "rtki_fpkm.{ext}")
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

from load_data import rnaseq_data
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
    pids = ['017', '050', '054', '019', '030', '031']
    # obj = rnaseq_data.load_salmon_by_patient_id(pids, include_control=False)
    # dat = obj.get_fpkm()

    dat = rnaseq_data.load_salmon_by_patient_id(pids, include_control=False)

    # drop any iPSC
    dat = dat.loc[:, ~dat.columns.str.contains('IPSC')]

    idx = dat.index.str.replace(r'.[0-9]+$', '')
    dat.index = idx

    fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90',
                      'gene_to_transcript.txt')
    gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

    # aggregate to gene level
    genes = gene_transcript.loc[dat.index, 'Gene stable ID']
    dat = dat.groupby(genes).sum()

    # convert to gene symbols
    idx = references.ensembl_to_gene_symbol(dat.index).dropna()
    idx = idx.loc[~idx.index.duplicated()]
    dat = dat.loc[idx.index]
    dat.index = idx

    # write individual patient data
    for pid in pids:
        # the_idx = obj.meta.index[obj.meta.index.str.contains(pid)]
        the_idx = dat.columns.str.contains(pid)
        the_dat = dat.loc[:, the_idx]
        # the_classes = obj.meta.loc[the_idx, 'type'].values
        the_classes = pd.Series('GBM', index=the_dat.columns)
        the_classes.loc[the_classes.index.str.contains('NSC')] = 'iNSC'
        out_fn = os.path.join(outdir, "%s_fpkm.{ext}" % pid)
        gsea.data_to_gct(the_dat, out_fn.format(ext='gct'))
        gsea.phenotypes_to_cls(the_classes, out_fn.format(ext='cls'))

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

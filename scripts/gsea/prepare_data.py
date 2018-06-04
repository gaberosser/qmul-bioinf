from rnaseq import loader, general
from utils.output import unique_output_dir
from rnaseq import gsea
import os
import references
import pandas as pd
from settings import OUTPUT_DIR


def ens_index_to_gene_symbol(df):
    general.add_gene_symbols_to_ensembl_data(df)
    tmp = df['Gene Symbol'].dropna()
    df = df.loc[tmp.index]
    df.set_index('Gene Symbol', inplace=True)
    return df


if __name__ == '__main__':
    outdir = os.path.join(OUTPUT_DIR, "gsea_data")

    # load all data
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    units = 'tpm'

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
    At this point, we need to run GSEA.

    ##############################################
    # 1. Running directly from the command line: #
    ##############################################

    java -Xmx16192m -cp /opt/gsea/gsea-3.0.jar xtools.gsea.Gsea -res /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/tpm_salmon/rtkii.gct -cls /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/tpm_salmon/rtkii.cls#GBM_versus_iNSC -gmx /home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/msigdb.v6.1.symbols.gmt -collapse false -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label RTK_II -metric Signal2Noise -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 1000 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /home/gabriel/gsea_home/output/rtk_II -gui false

    #############################
    # 2. Use a parameters file: #
    #############################

    # See rnaseq.gsea.create_gsea_params_file.
    # Specify with -params_file

    for i in *.params;
        do b=$(basename $i .params);
        if [[ $b = *"nsc"* ]]; then
            c="#GBM_versus_NSC";
        else
            c="#GBM_versus_iNSC";
        fi;
        java -Xmx16192m -cp /opt/gsea/gsea-3.0.jar xtools.gsea.Gsea \
        -res ${b}.gct -cls ${b}.cls${c} -gmx $gmx \
        -out ./${b} \
        -param_file ${b}.params

        STATUS=$?

        if [ $STATUS == 0 ]; then
            echo "SUCCESS: ${b}"
            mv ${b}.gct ${b}.cls ${b}.params ${b}/
        else
            echo "FAILED: ${b}"
        fi
    done

    #######################################
    # 3. parallel method for many files:  #
    #######################################

    runGsea() {
        # how do I make this globally accessible within parallel? export it?
        gmx="/home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/msigdb.v6.1.symbols.gmt"
        b=$(basename $1 .params)

        if [[ $b = *"nsc"* ]]; then
            c="#GBM_versus_NSC";
        else
            c="#GBM_versus_iNSC";
        fi;

        java -Xmx16192m -cp /opt/gsea/gsea-3.0.jar xtools.gsea.Gsea \
        -res ${b}.gct -cls ${b}.cls${c} -gmx $gmx \
        -out ./${b} \
        -param_file ${b}.params

        STATUS=$?

        if [ $STATUS == 0 ]; then
            echo "SUCCESS: ${b}"
            mv ${b}.gct ${b}.cls ${b}.params ${b}/
        else
            echo "FAILED: ${b}"
        fi
    }
    # export so that the parallel env has the function in scope
    export -f runGsea

    parallel -j 3 runGsea ::: *.params

    """
    
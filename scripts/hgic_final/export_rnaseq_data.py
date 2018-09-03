from rnaseq import loader, filter, general
from utils import output
import os
import sys
import consts


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)

    star_cc_obj = loader.load_by_patient(
        consts.PIDS,
        type='cell_culture',
        source='star',
        include_control=True,
    )
    salmon_cc_obj = loader.load_by_patient(
        consts.PIDS,
        type='cell_culture',
        source='salmon',
        include_control=True,
    )
    star_ff_obj = loader.load_by_patient(
        consts.PIDS,
        type='ffpe',
        source='star',
        include_control=True,
    )
    salmon_ff_obj = loader.load_by_patient(
        consts.PIDS,
        type='ffpe',
        source='salmon',
        include_control=True,
    )

    # combine cell culture and FFPE
    star_all_obj = loader.loader.MultipleBatchLoader([
        star_cc_obj,
        star_ff_obj
    ])

    salmon_all_obj = loader.loader.MultipleBatchLoader([
        salmon_cc_obj,
        salmon_ff_obj
    ])

    # STAR gene counts
    dat = star_all_obj.data.copy()
    general.add_gene_symbols_to_ensembl_data(dat)
    dat.to_excel(os.path.join(outdir, "star_counts.xlsx"))

    # convert to cpm
    dat = star_all_obj.data.copy()
    dat = dat.divide(dat.sum(axis=0), axis=1) * 1e6
    general.add_gene_symbols_to_ensembl_data(dat)
    dat.to_excel(os.path.join(outdir, "star_cpm.xlsx"))

    # Salmon TPM
    dat = salmon_all_obj.data.copy()
    general.add_gene_symbols_to_ensembl_data(dat)
    dat.to_excel(os.path.join(outdir, "salmon_tpm.xlsx"))

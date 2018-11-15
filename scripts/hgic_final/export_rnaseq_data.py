from rnaseq import loader, filter, general
from utils import output
import os
import sys
import consts


def add_genes_and_save(dat, fn, type='excel'):
    general.add_gene_symbols_to_ensembl_data(dat)
    if type == 'excel':
        save_func = lambda x: x.to_excel
    elif type == 'csv':
        save_func = lambda x: x.to_csv
    else:
        raise NotImplementedError("Unsupported type %s." % type)

    save_func(dat)(fn)


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    keep_samples = consts.S1_RNASEQ_SAMPLES_FB + consts.S1_RNASEQ_SAMPLES_INSC + consts.S1_RNASEQ_SAMPLES_GIC + \
        consts.S1_RNASEQ_SAMPLES_IPSC + consts.S1_RNASEQ_SAMPLES_IAPC

    star_cc_obj = loader.load_by_patient(
        consts.PIDS,
        type='cell_culture',
        source='star',
        include_control=True,
    )
    ix = star_cc_obj.meta.index.isin(keep_samples)
    star_cc_obj.filter_samples(ix)

    salmon_cc_obj = loader.load_by_patient(
        consts.PIDS,
        type='cell_culture',
        source='salmon',
        include_control=True,
    )
    ix = salmon_cc_obj.meta.index.isin(keep_samples)
    salmon_cc_obj.filter_samples(ix)

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

    add_genes_and_save(
        star_cc_obj.data.copy(),
        os.path.join(outdir, "cell_culture_star_counts.xlsx")
    )
    star_cc_obj.meta.to_excel(os.path.join(outdir, "cell_culture_star_meta.xlsx"))

    dat = star_cc_obj.data.copy()
    dat = dat.divide(dat.sum(axis=0), axis=1) * 1e6
    add_genes_and_save(
        dat,
        os.path.join(outdir, "cell_culture_star_cpm.xlsx")
    )

    add_genes_and_save(
        salmon_cc_obj.data.copy(),
        os.path.join(outdir, "cell_culture_salmon_tpm.xlsx")
    )
    salmon_cc_obj.meta.to_excel(os.path.join(outdir, "cell_culture_salmon_meta.xlsx"))

    add_genes_and_save(
        star_ff_obj.data.copy(),
        os.path.join(outdir, "ffpe_star_counts.xlsx")
    )
    star_ff_obj.meta.to_excel(os.path.join(outdir, "ffpe_star_meta.xlsx"))

    dat = star_ff_obj.data.copy()
    dat = dat.divide(dat.sum(axis=0), axis=1) * 1e6
    add_genes_and_save(
        dat,
        os.path.join(outdir, "ffpe_star_cpm.xlsx")
    )

    add_genes_and_save(
        salmon_ff_obj.data.copy(),
        os.path.join(outdir, "ffpe_salmon_tpm.xlsx")
    )
    salmon_ff_obj.meta.to_excel(os.path.join(outdir, "ffpe_salmon_meta.xlsx"))

    # since it may be useful, also export a full set of FFPE data
    # (i.e. including patients not in the current cohort)
    ## FIXME: we need an 'exclusion list' here (or a full inclusion list), otherwise broken samples will be included
    if False:
        star_ff_obj_all = loader.load_by_patient(
            'all',
            type='ffpe',
            source='star',
            include_control=False,
        )

        add_genes_and_save(
            star_ff_obj_all.data.copy(),
            os.path.join(outdir, "ffpe_star_counts_all_samples.xlsx")
        )
        star_ff_obj_all.meta.to_excel(os.path.join(outdir, "ffpe_star_meta_all_samples.xlsx"))

        dat = star_ff_obj_all.data.copy()
        dat = dat.divide(dat.sum(axis=0), axis=1) * 1e6
        add_genes_and_save(
            dat,
            os.path.join(outdir, "ffpe_star_cpm_all_samples.xlsx")
        )

        salmon_ff_obj_all = loader.load_by_patient(
            'all',
            type='ffpe',
            source='salmon',
            include_control=False,
        )

        add_genes_and_save(
            salmon_ff_obj_all.data.copy(),
            os.path.join(outdir, "ffpe_salmon_tpm_all_samples.xlsx")
        )
        salmon_ff_obj_all.meta.to_excel(os.path.join(outdir, "ffpe_salmon_meta_all_samples.xlsx"))
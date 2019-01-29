from methylation import loader, process
import consts
import os
import collections
import pandas as pd
import numpy as np
from utils import output
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP


def expand_annotation_by_gene_and_cluster(
    df,
    drop_geneless=False,
    gene_name_col='UCSC_RefGene_Name',
    gene_rel_col='UCSC_RefGene_Group',
):
    if drop_geneless:
        df = df.loc[~df[gene_name_col].isnull()]

    # make (gene, relation) lists
    gr = []
    for k, row in df[[gene_name_col, gene_rel_col]].iterrows():
        gr.append(set(zip(*row.values)))

    gr = pd.DataFrame(gr, index=df.index)

    if not drop_geneless:
        idx = gr[0].isnull()
        gr.loc[idx, 0] = [('!', '!')] * idx.sum()

    gr = gr.stack().reset_index(level=-1, drop=True)
    gr = pd.concat([gr.apply(lambda x: x[0]), gr.apply(lambda x: x[1])], axis=1)
    gr.columns = ['gene', 'relation']
    gr.insert(0, 'probe_id', gr.index)


    # pivot table around the relation
    tbl = gr.pivot_table(index=['probe_id', 'gene'], columns='relation', fill_value=0, aggfunc='size').reset_index()
    tbl.columns.name = None

    # extract other attributes from original dataframe
    attrs = df.drop([gene_name_col, gene_rel_col], axis=1).loc[tbl.probe_id]

    # combine them
    res = pd.concat((tbl.set_index('probe_id'), attrs), axis=1).reset_index()

    # if maintaining geneless results, deprotect these now
    if not drop_geneless:
        new_gene_col = res.loc[:, 'gene'].replace('!', np.nan)  # don't use None as the replacement variable!
        res.loc[:, 'gene'] = new_gene_col
        # drop the pivot column, but only if it actually exists
        if '!' in res.columns:
            res.drop('!', axis=1, inplace=True)
        # quick sanity check
        cols = sorted(set(gr.relation.unique()).difference(['!']))
        ix = res.gene.isnull()
        if not (res.loc[ix, cols].sum(axis=1) == 0).all():
            raise ValueError("Some 'geneless' clusters still have gene relations defined.")

    res.set_index('index', inplace=True)
    return res


if __name__ == "__main__":
    """
    Here, we simply load the methylation data and export it in an efficient manner (restricting the floating point
    bit depth to save space).

    We also export the annotation and metadata separately.
    """
    norm_method = 'swan'
    # the float format is used when exporting to Excel - it reduces the file size by restricting the precision
    float_format = '%.2f'
    outdir = output.unique_output_dir()
    anno = loader.load_illumina_methylationepic_annotation()
    obj_cc = loader.load_by_patient(
        consts.PIDS,
        type='cell_culture',
        norm_method=norm_method,
        reduce_to_common_probes=False
    )
    obj_ff = loader.load_by_patient(consts.PIDS, type='ffpe', norm_method=norm_method, reduce_to_common_probes=False)
    # add useful patient ID column to metadata
    obj_ff.meta.insert(0, 'patient_id', [NH_ID_TO_PATIENT_ID_MAP[k] for k in obj_ff.meta.index])

    # export methylation data
    obj_cc.data.to_excel(
        os.path.join(outdir, "methylation_beta_cell_culture.xlsx"),
        float_format=float_format
    )
    obj_ff.data.to_excel(
        os.path.join(outdir, "methylation_beta_ffpe.xlsx"),
        float_format=float_format
    )

    # export metadata
    obj_cc.meta.to_excel(
        os.path.join(outdir, "metadata_cell_culture.xlsx")
    )
    obj_ff.meta.to_excel(
        os.path.join(outdir, "metadata_ffpe.xlsx")
    )

    # prepare annotation data
    df = expand_annotation_by_gene_and_cluster(anno.copy(), drop_geneless=False)
    # export
    df.to_excel(os.path.join(outdir, 'epic_array_probe_annotation.xlsx'))



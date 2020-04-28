import pandas as pd
import os
import pickle
import collections
import numpy as np

from scripts.hgic_final import two_strategies_combine_de_dmr as tscdd, consts
from utils import setops, output, log, druggable_genome, excel, reference_genomes
from settings import INTERMEDIATE_DIR
from rnaseq import differential_expression, loader as rnaseq_loader
from methylation import dmr, annotation_gene_to_ensembl
from integrator import rnaseq_methylationarray
from plotting import common
logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir(reuse_empty=True)

    min_cpm = 1.
    # Set this to True to include reference methylation data
    # This will limit the number of available probes (to 450K)
    include_external_dm_refs = False

    de_params = consts.DE_PARAMS
    dmr_params = consts.DMR_PARAMS

    norm_method_s1 = 'swan'

    pids = consts.PIDS

    subgroups = consts.SUBGROUPS

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

    if not os.path.isdir(DE_LOAD_DIR):
        raise AttributeError("To run this script, we require pre-computed DE results residing in %s" % DE_LOAD_DIR)

    if not os.path.isdir(DMR_LOAD_DIR):
        raise AttributeError("To run this script, we require pre-computed DMR results residing in %s" % DMR_LOAD_DIR)

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # some quantities relating to set membership

    # load methylation
    # For S1, we just need the paired comparison. This is important - any other samples lead to a change in the final
    # probe list (any NA rows are dropped from ALL samples).
    me_obj, anno = tscdd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tscdd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    rnaseq_obj = obj = rnaseq_loader.load_by_patient(pids)

    rnaseq_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    # only keep the required syngeneic samples for this analysis
    dat_s1 = rnaseq_obj.data
    meta_s1 = rnaseq_obj.meta

    the_hash = tscdd.de_results_hash(meta_s1.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise AttributeError("Unable to load pre-computed DE results, expected at %s" % fn)

    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

    # get the joint table
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    # run the dgidb lookup against all genes
    # have to chunk this operation to avoid error
    all_genes = sorted(setops.reduce_union(*[t.gene.values for t in joint_de_dmr_s1.values()]))
    dgi_all = druggable_genome.dgidb_lookup_drug_gene_interactions(all_genes)

    # manually resolve a few known ambiguities
    ambig = {
        'ELTD1': 'ADGRL4',
        'ODZ3': 'TENM3'
    }
    for k, v in ambig.items():
        x = [t for t in dgi_all['ambiguous'][k] if t['geneName'] == v][0]
        dgi_all['interactions'][k] = x['interactions']

    de_dmr_by_member = [joint_de_dmr_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member)

    # define short and long list

    # long list
    ss = setops.specific_sets(pids)
    ps_de_dm_long = collections.OrderedDict([
        (pid, venn_set[ss[pid]]) for pid in pids
    ])

    ps_de_dm_long_list = setops.reduce_union(
        *ps_de_dm_long.values()
    )

    # short list
    vs_dm, vc_dm = setops.venn_from_arrays(*[dmr_res_s1[pid].results_significant.keys() for pid in pids])
    vs_de, vc_de = setops.venn_from_arrays(*[de_res_s1[pid]['Gene Symbol'].dropna() for pid in pids])

    ps_dm = {
        pid: vs_dm[v] for pid, v in ss.items()
    }
    ps_dm_genes = {}
    for pid, cids in ps_dm.items():
        ps_dm_genes[pid] = sorted(setops.reduce_union(*[[t[0] for t in dmr_res_s1.clusters[c].genes] for c in cids]))

    ps_de = {
        pid: vs_de[v] for pid, v in ss.items()
    }

    ps_de_dm_short = {}
    for pid in pids:
        this_genes = set(ps_de[pid]).intersection(ps_dm_genes[pid])
        this_clusters = []
        for c in ps_dm[pid]:
            this = [t[0] for t in dmr_res_s1.clusters[c].genes]
            if len(this_genes.intersection(this)):
                this_clusters.extend([(c, g) for g in this_genes.intersection(this)])
        # ps_de_dm_short[pid] = joint_de_dmr_s1[pid].reindex(this_clusters).dropna(axis=0, how='all')
        ps_de_dm_short[pid] = this_clusters

    ps_de_dm_short_list = setops.reduce_union(
        *ps_de_dm_short.values()
    )

    # sanity check: everything on the shortlist is also on the longlist
    assert len(ps_de_dm_short_list.difference(ps_de_dm_long_list)) == 0, "Expected: shortlist is a subset of the longlist"

    # promiscuous list (DE and DM in all)
    all_dm = vs_dm[''.join(['1'] * len(pids))]
    all_de = vs_de[''.join(['1'] * len(pids))]

    all_dm_genes = sorted(setops.reduce_union(*[[t[0] for t in dmr_res_s1.clusters[c].genes] for c in all_dm]))

    this_genes = set(all_de).intersection(all_dm_genes)
    all_de_dm = []
    for c in all_dm:
        this = [t[0] for t in dmr_res_s1.clusters[c].genes]
        all_de_dm.extend([(c, g) for g in this_genes.intersection(this)])

    # table export
    all_dm_genes = sorted(setops.reduce_union(
        *[[t[0] for t in dmr_res_s1.clusters[c].genes] for c in dmr_res_s1.clusters]
    ))
    all_dm_ens = annotation_gene_to_ensembl.gene_to_ens(all_dm_genes).dropna()

    gene_to_dm_cluster = {}
    for c, cl in dmr_res_s1.clusters.items():
        for g in cl.genes:
            gene_to_dm_cluster.setdefault(g[0], set()).add(c)

    ens_to_dm_cluster = {}
    for c, cl in dmr_res_s1.clusters.items():
        for g in cl.genes:
            if g[0] in all_dm_ens:
                e = all_dm_ens[g[0]]
                ens_to_dm_cluster.setdefault(e, set()).add(c)


    all_de_ens_with_dm = sorted(setops.reduce_union(
        *[de_res_full_s1[pid].reindex(all_dm_ens).dropna().index for pid in pids]
    ))
    ens_to_gs = reference_genomes.ensembl_to_gene_symbol(all_de_ens_with_dm)

    # single mega table
    single_de_dm_df = []
    de_fdr = {}
    dm_fdr = {}

    # plus separate DGIdb dump
    dgi_db_df = []

    # use to determine which relations are relevant
    all_relations = sorted(
        setops.reduce_union(*[[t[1] for t in cl.genes] for cl in dmr_res_s1.clusters.values()])
    )


    for e in all_de_ens_with_dm:
        g = ens_to_gs[e]

        related_dm_clusters = ens_to_dm_cluster[e]

        related_dm_genes = all_dm_ens.index[all_dm_ens == e]

        # DGIdb results
        this_n_interactions = len(dgi_all['interactions'].get(g, []))

        for d in dgi_all['interactions'].get(g, []):
            dgi_db_df.append([e, g, d['drugName'], d['drugChemblId'], ';'.join(d['interactionTypes'])])

        this_block = []
        # DE results are the same across this whole block
        this_de_fdr = []
        this_de_lfc = []
        for pid in pids:
            this_de_lfc.append(de_res_full_s1[pid].reindex([e])['logFC'].squeeze())
            this_de_fdr.append(de_res_full_s1[pid].reindex([e])['FDR'].squeeze())
        de_fdr[e] = this_de_fdr

        for cid in related_dm_clusters:
            # get relevant meth - gene relations
            this_rels = set([t[1] for t in dmr_res_s1.clusters[cid].genes if all_dm_ens.get(t[0], '') == e])
            rel_chunk = [int(t in this_rels) for t in all_relations]

            this_dm_delta_m = []
            this_dm_fdr = []
            for pid in pids:
                this_dm_delta_m.append(dmr_res_full_s1[pid][cid]['median_change'])
                this_dm_fdr.append(dmr_res_full_s1[pid][cid].get('padj'))
            dm_fdr[cid] = this_dm_fdr

            # what 'list' does it belong to?
            aa = set([(cid, dm_g) for dm_g in related_dm_genes])
            ns = len(aa.intersection(ps_de_dm_short_list))
            nl = len(aa.intersection(ps_de_dm_long_list))
            na = len(aa.intersection(all_de_dm))

            mmbr = ''
            if ns > 0:
                logger.info("Shortlisted genes: %s", ', '.join([str(t) for t in aa.intersection(ps_de_dm_short_list)]))
                mmbr = 'shortlist'
            elif nl > 0:
                logger.info("Longlisted genes: %s",
                            ', '.join([str(t) for t in aa.intersection(ps_de_dm_long_list)]))
                mmbr = 'longlist'
            elif na > 0:
                logger.info("Promiscuous genes: %s",
                            ', '.join([str(t) for t in aa.intersection(all_de_dm)]))
                mmbr = 'common'

            # create block entry
            this_block.append(
                [e, g, cid] + rel_chunk + this_de_lfc + this_dm_delta_m + [this_n_interactions, mmbr]
            )

        single_de_dm_df.extend(this_block)

    # create a df
    cols = ['ensembl_gene_id', 'gene_symbol', 'dmr_cluster'] \
           + all_relations \
           + ['%s_de_logfc' % pid for pid in pids] \
           + ['%s_dm_delta_m' % pid for pid in pids] \
           + ['num_DGIdb_interactions', 'list_member']
    single_de_dm_df = pd.DataFrame(single_de_dm_df, columns=cols)

    de_fdr = pd.DataFrame(de_fdr, index=pids).transpose().fillna(1.)
    dm_fdr = pd.DataFrame(dm_fdr, index=pids).transpose().fillna(1.)

    dgi_cols = [
        'ensembl_gene_id', 'gene_symbol', 'name', 'chembl_id', 'interaction_type'
    ]
    dgi_db_df = pd.DataFrame(dgi_db_df, columns=dgi_cols)

    # export to Excel and format it
    out_fn = os.path.join(outdir, "de_dmr_dgi_results_wide.xlsx")
    sheet_name = "DE DMR results"

    writer = pd.ExcelWriter(out_fn, engine='xlsxwriter')
    single_de_dm_df.to_excel(writer, sheet_name=sheet_name)
    dgi_db_df.to_excel(writer, sheet_name='DGIdb')

    workbook = writer.book
    worksheet = writer.sheets[sheet_name]
    formatter = excel.ExcelFormatter(writer)

    # merge gene-related data in some columns
    # this has to happen first, as it overwrites data and formatting
    vcentred_fmt_dict = {'valign': 'vcenter'}
    formatter.add_format('vcentre', vcentred_fmt_dict)
    vcentre_fmt = formatter._formats['vcentre']

    col_to_merge = ['ensembl_gene_id', 'gene_symbol'] + ['%s_de_logfc' % pid for pid in pids] + ['num_DGIdb_interactions']
    col_to_merge_letters = [excel.get_column_letter(single_de_dm_df, c) for c in col_to_merge]
    dupe_ix = single_de_dm_df['ensembl_gene_id'].duplicated()
    dupes = single_de_dm_df['ensembl_gene_id'][dupe_ix].unique()
    for e in dupes:
        row_ix = single_de_dm_df.index[single_de_dm_df['ensembl_gene_id'] == e]
        first_row = excel.get_row_number(single_de_dm_df, row_ix[0])
        last_row = excel.get_row_number(single_de_dm_df, row_ix[-1])
        for c, col in zip(col_to_merge_letters, col_to_merge):
            val = single_de_dm_df.loc[row_ix, col].fillna('').unique()
            if len(val) != 1:
                raise ValueError(
                    "Found %d unique values in column %s for gene %s." % (len(val), col, e)
                )
            val = val[0]
            merge_range = '%s%s:%s%s' % (c, first_row, c, last_row)
            worksheet.merge_range(merge_range, val, cell_format=vcentre_fmt)

    # comments (DGIdb)
    comments = dgi_db_df.groupby('ensembl_gene_id').apply(
        lambda x: '\n'.join([
            t['name'] + ("" if t.interaction_type == "" else " (%s)" % t.interaction_type)
            for _, t in x.iterrows()
        ])
    )
    col_letter = excel.get_column_letter(single_de_dm_df, 'num_DGIdb_interactions')
    for ix, row in single_de_dm_df[single_de_dm_df['num_DGIdb_interactions'] > 0].iterrows():
        row_letter = excel.get_row_number(single_de_dm_df, ix)
        worksheet.write_comment('%s%s' % (col_letter, row_letter), comments.loc[row['ensembl_gene_id']])

    # shading (logFC)
    green_cmap = common.continuous_cmap(-8, 0, 'Greens_r')
    red_cmap = common.continuous_cmap(0, 8, 'Reds')

    breaks_de = list(zip(range(-5, 5), range(-4, 6)))
    breaks_de[0] = (None, breaks_de[0][1])
    breaks_de[-1] = (breaks_de[-1][0], None)
    cmap_de_formats = []
    for lo, hi in breaks_de:
        if lo is None or lo < 0:
            this_fmt_dict = {'bg_color': green_cmap(hi - 0.5)}
            this_fmt_name = '{}_{}'.format(lo or 'LOW', hi)
            formatter.add_format(this_fmt_name, this_fmt_dict)
            cmap_de_formats.append(this_fmt_name)
        else:
            this_fmt_dict = {'bg_color': red_cmap(lo + 0.5)}
            this_fmt_name = '{}_{}'.format(lo, hi or 'HIGH')
            formatter.add_format(this_fmt_name, this_fmt_dict)
            cmap_de_formats.append(this_fmt_name)

    # range to apply
    cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_logfc')]
    c0 = excel.get_column_letter(single_de_dm_df, cols[0])
    c1 = excel.get_column_letter(single_de_dm_df, cols[-1])
    r0 = excel.get_row_number(single_de_dm_df, single_de_dm_df.index[0])
    r1 = excel.get_row_number(single_de_dm_df, single_de_dm_df.index[-1])
    range_str = "%s%d:%s%d" % (c0, r0, c1, r1)

    formatter.add_quantitative_colourmap('de', breaks_de, cmap_de_formats)
    formatter.conditional_format_colourmap(range_str, 'de', sheet_name=sheet_name)

    # shading (delta M)
    # same colourmap as for DE
    # range to apply
    cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_delta_m')]
    c0 = excel.get_column_letter(single_de_dm_df, cols[0])
    c1 = excel.get_column_letter(single_de_dm_df, cols[-1])
    r0 = excel.get_row_number(single_de_dm_df, single_de_dm_df.index[0])
    r1 = excel.get_row_number(single_de_dm_df, single_de_dm_df.index[-1])
    range_str = "%s%d:%s%d" % (c0, r0, c1, r1)

    formatter.add_quantitative_colourmap('dm', breaks_de, cmap_de_formats)
    formatter.conditional_format_colourmap(range_str, 'dm', sheet_name=sheet_name)

    # borders around significant results
    outlined_fmt_dict = {
        'border': 1
    }
    formatter.add_format('border', outlined_fmt_dict)
    fmt = formatter._formats['border']

    # DE
    cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_logfc')]
    col_nums = [single_de_dm_df.columns.tolist().index(t) for t in cols]
    for row_num, (ix, row) in enumerate(single_de_dm_df.iterrows()):
        for j, (col_num, col) in enumerate(zip(col_nums, cols)):
            e = row['ensembl_gene_id']
            pid = col.replace('_de_logfc', '')
            this_fdr = de_fdr.loc[e, pid]
            if this_fdr < de_params['fdr']:
                formatter.apply_format_one_cell(row_num + 1, col_num + 1, 'border', sheet_name=sheet_name)

    # DM
    cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_delta_m')]
    col_nums = [single_de_dm_df.columns.tolist().index(t) for t in cols]
    for row_num, (ix, row) in enumerate(single_de_dm_df.iterrows()):
        for j, (col_num, col) in enumerate(zip(col_nums, cols)):
            cid = row['dmr_cluster']
            pid = col.replace('_dm_delta_m', '')
            this_fdr = dm_fdr.loc[cid, pid]
            if this_fdr < dmr_params['alpha']:
                formatter.apply_format_one_cell(row_num + 1, col_num + 1, 'border', sheet_name=sheet_name)

    # set text colour depending on concordancy
    conc_fmt_dict = {
        'bold': True
    }
    disc_fmt_dict = {
        'italic': True,
        'font_color': '#666666'
    }
    formatter.add_format('concordant', conc_fmt_dict)
    formatter.add_format('discordant', disc_fmt_dict)

    cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_de_')]
    matching_cols = single_de_dm_df.columns[single_de_dm_df.columns.str.contains('_dm_')]
    col_nums = [single_de_dm_df.columns.tolist().index(t) for t in cols]
    matching_col_nums = [single_de_dm_df.columns.tolist().index(t) for t in matching_cols]

    for row_num, (ix, row) in enumerate(single_de_dm_df.iterrows()):
        for j, (col_num, col) in enumerate(zip(col_nums, cols)):
            logfc = row[col]
            other_col = matching_cols[j]
            other_col_num = matching_col_nums[j]

            delta_m = row[other_col]
            if logfc is None or delta_m is None:
                # cannot define concordancy
                continue

            if np.sign(logfc) == np.sign(delta_m):
                # discordant
                formatter.apply_format_one_cell(row_num + 1, col_num + 1, 'discordant', sheet_name=sheet_name)
                formatter.apply_format_one_cell(row_num + 1, other_col_num + 1, 'discordant', sheet_name=sheet_name)
            else:
                formatter.apply_format_one_cell(row_num + 1, col_num + 1, 'concordant', sheet_name=sheet_name)
                formatter.apply_format_one_cell(row_num + 1, other_col_num + 1, 'concordant', sheet_name=sheet_name)

    # place a dividing line to separate DM and DE
    first_dm_col = col_nums[0] + 1
    left_border_fmt_dict = {
        'left': 6
    }
    formatter.add_format('divider', left_border_fmt_dict)
    for col_num in [col_nums[0], matching_col_nums[0]]:
        for row_num in range(single_de_dm_df.shape[0]):
            formatter.apply_format_one_cell(row_num + 1, col_num + 1, 'divider', sheet_name=sheet_name)

    writer.save()





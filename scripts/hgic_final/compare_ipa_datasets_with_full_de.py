import pandas as pd
import os
import numpy as np


def compare(curr_df, ipa_df):
    if ipa_df.index.sort_values().equals(curr_df.index.sort_values()):
        print "Indexes match"
    else:
        print "Indexes DO NOT match"
        ix = ipa_df.index.difference(curr_df.index)
        print "In IPA dataset and not in current (%d): %s" % (
            len(ix),
            ', '.join(ix)
        )
        ix = curr_df.index.difference(ipa_df.index)
        print "In current and not in IPA dataset (%d): %s" % (
            len(ix),
            ', '.join(ix)
        )

    # check remaining values agree
    df2 = curr_df.reindex(ipa_df.index)
    b = np.isclose(df2["%s_logFC" % p], df["%s_logFC" % p])
    if b.all():
        print "The logFC values are close in all cases."
    else:
        print "The logFC values are not close in %d cases." % (~b).sum()
    b = np.isclose(np.log10(df2["%s_FDR" % p]), np.log10(df["%s_FDR" % p]))
    if b.all():
        print "The FDR values are close in all cases."
    else:
        print "The FDR values are not close in %d cases." % (~b).sum()


if __name__ == '__main__':
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    # # full DE:
    # ipa_indir = os.path.join(
    #     os.path.expanduser('~'),
    #     'Dropbox',
    #     'research',
    #     'qmul',
    #     'hGIC_project',
    #     'current',
    #     'core_pipeline',
    #     'rnaseq',
    #     'ipa',
    #     'data_uploaded'
    # )
    # # path to most current full DE results:
    # curr_fn = os.path.join(
    #     os.path.expanduser('~'),
    #     'Dropbox',
    #     'research',
    #     'qmul',
    #     'hGIC_project',
    #     'current',
    #     'core_pipeline',
    #     'rnaseq',
    #     'full_de.xlsx'
    # )
    # filestem = 'full_de_patient{pid}.txt'

    # patient specific (S1):
    ipa_indir = os.path.join(
        os.path.expanduser('~'),
        'Dropbox',
        'research',
        'qmul',
        'hGIC_project',
        'current',
        'core_pipeline',
        'rnaseq',
        's1_inter_patient',
        'ipa',
        'data_uploaded'
    )
    # path to most current S1 patient specific results:
    curr_fn = os.path.join(
        os.path.expanduser('~'),
        'Dropbox',
        'research',
        'qmul',
        'hGIC_project',
        'current',
        'core_pipeline',
        'rnaseq',
        's1_inter_patient',
        'patient_specific_de.xlsx'
    )
    filestem = 'patient_specific_de_{pid}.txt'


    full_de = pd.read_excel(curr_fn)
    for p in pids:
        this_fn = os.path.join(ipa_indir, filestem.format(pid=p))
        df = pd.read_csv(this_fn, sep='\t', header=0, index_col=0).sort_values(by='%s_FDR' % p)
        this_full_de = full_de.loc[
            full_de[p] == 'Y',
            full_de.columns.str.contains(p) | full_de.columns.str.contains('Gene')
        ].sort_values(by='%s_FDR' % p)

        print p
        compare(this_full_de, df)


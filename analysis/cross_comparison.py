from utils import setops
import pandas as pd


# TODO: this is named inaccurately - it does a lot more!
def compute_cross_comparison_correction(res, samples, external_refs, set_type='pair_only'):
    """
    Compute the _correction_ list of features for the supplied results. These are the features that are
    EITHER present in every reference comparison but no cross-comparisons (set_type='ref_only')
    OR present in no reference comparison but all cross-comparisons (set_type='pair_only')
    :param res: Dictionary containing comparison results. Each comparison is keyed by the tuple (i, j), where i and j
    are the IDs of the two groups being compared. Values are iterables of unique feature identifiers (e.g. gene IDs,
    DMR cluster IDs).
    :param samples: The core sample list, without including external references.
    :param external_refs: A list of external reference sample names.
    :param set_type: See description.
    :return: Iterable of feature IDs
    """

    members_rows = samples
    members_cols = members_rows + external_refs

    the_venn_set = pd.DataFrame(index=members_rows, columns=members_cols)
    for i in members_rows:
        p = res[(i, i)]
        for j in members_cols:
            r = res[(i, j)]
            x, _ = setops.venn_from_arrays(p, r)
            if set_type == 'pair_only':
                kset = '10'
            elif set_type == 'ref_only':
                kset = '01'
            else:
                raise AttributeError("set_type must be 'pair_only' or 'ref_only'.")
            the_venn_set.loc[i, j] = x[kset]

    # For each reference, get the features that are pair only in that reference and not in any of the iNSC
    vs_diff = pd.DataFrame(index=members_rows, columns=external_refs)
    for i in members_rows:
        for j in external_refs:
            the_ref = the_venn_set.loc[i, j]
            all_else = the_venn_set.loc[i, members_rows]
            union_all_else = setops.reduce_union(*all_else)
            vs_diff.loc[i, j] = sorted(set(the_ref).difference(union_all_else))

    # Intersection down the columns gives us a correction list for each reference
    vs_specific_to_ref = vs_diff.apply(lambda x: setops.reduce_intersection(*x))

    # Intersection across the references gives us a final list that need correcting
    vs_specific_to_all_refs = setops.reduce_intersection(*vs_specific_to_ref)

    return {
        'specific_to_each_ref': vs_specific_to_ref,
        'specific_to_all_refs': vs_specific_to_all_refs,
        'venn_set': the_venn_set,
        'ref_diff_set': vs_diff
    }
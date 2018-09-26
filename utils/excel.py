import pandas as pd


def pandas_to_excel(sheets, fn, write_index=True, **kwargs):
    """

    :param blocks: Dictionary containing the different comparisons to save. Values are pandas dataframes.
    :param fn: Output file
    :param write_index: If True (default) then the index is written as the first column. May wish to disable this if it
    has no meaning.
    :return:
    """
    xl_writer = pd.ExcelWriter(fn)
    # sort the keys for a more sensible order
    keys = sorted(sheets.keys())
    for k in keys:
        bl = sheets[k]
        bl.to_excel(xl_writer, k, index=write_index, **kwargs)
    xl_writer.save()
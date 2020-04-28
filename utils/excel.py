import pandas as pd
import string
from utils import log
import xlsxwriter

logger = log.get_console_logger()


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


def colnum_string(n):
    """
    https://stackoverflow.com/a/23862195/1474448
    :param n:
    :return:
    """
    s = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        s = chr(65 + remainder) + s
    return s


def get_column_letter(df, colname, index=True):
    """
    Get the letter of the column in the supplied dataframe with the specified name.
    :param write_index: If True (default), index/row names will be in column A, so we need to offset by one.
    """
    if index:
        try:
            n_idx = len(df.index.levels)
        except AttributeError:
            n_idx = 1
    else:
        n_idx = 0

    i = df.columns.tolist().index(colname) + n_idx + 1
    return colnum_string(i)


def get_row_number(df, index_name, header=True):
    """
    Get the row number of the DataFrame with the supplied index
    :param df:
    :param index_name:
    :param header: Has the header been written?
    :return:
    """
    if header:
        try:
            n_head = len(df.columns.levels)
        except AttributeError:
            n_head = 1
    else:
        n_head = 0

    return df.index.tolist().index(index_name) + 1 + n_head


class ExcelFormatter(object):
    def __init__(self, writer):
        self.writer = writer
        self.book = writer.book
        self._formats = {}
        self._format_dicts = {}
        self._custom_colourmaps = {}

    @property
    def formats_inverse(self):
        # FIXME: issue warning if this dictionary is poorly formed (same format appears under >1 name)
        return {v: k for k, v in self._formats.items()}

    def get_sheet(self, name=None):
        if name is None:
            sheet = self.book.worksheets()[0]
            if len(self.writer.sheets) > 1:
                logger.warning(
                    "No sheet name specified but the workbook contains %d sheets. Arbitrarily choosing the first (%s).",
                    len(self.writer.sheets),
                    sheet.name
                )
        else:
            sheet = self.writer.sheets[name]
        return sheet

    def add_format(self, name, format_dict):
        self._formats[name] = self.book.add_format(format_dict)
        self._format_dicts[name] = format_dict

    def add_quantitative_colourmap(self, name, bins, formats):
        """
        Create a custom quantitative colourmap that links existing static formats to bins
        :param name:
        :param bins: Iterable of (min, max) defining bin edges. These are on the half-open interval [min, max).
        A value of None in first or last position means an open interval.
        These don't have to be in order and we don't check that the request makes sense.
        :param formats:
        :return:
        """
        # small offset for ranges
        eps = 1e-12
        assert len(bins) == len(formats), "Length of bins does not match length of formats"
        the_cmap = []

        for b, f in zip(bins, formats):
            assert len(b) == 2, "Invalid bin interval {}".format(str(b))
            assert not (b[0] is None and b[1] is None), "Undefined interval (both limits are None)."
            assert f in self._formats, "Requested format {} has not been defined: run add_format().".format(f)

            if b[0] is None:
                fmt = {
                    'type': 'cell',
                    'criteria': '<',
                    'value': b[1],
                    'format': self._formats[f]
                }
            elif b[1] is None:
                fmt = {
                    'type': 'cell',
                    'criteria': '>=',
                    'value': b[0],
                    'format': self._formats[f]
                }
            else:
                fmt = {
                    'type': 'cell',
                    'criteria': 'between',
                    'minimum': b[0],
                    'maximum': float(b[1]) - eps,
                    'format': self._formats[f]
                }
            the_cmap.append(fmt)

        self._custom_colourmaps[name] = the_cmap

    def add_qualitative_colourmap(self, name, values, formats):
        assert len(values) == len(formats), "Length of bins does not match length of formats"
        the_cmap = []
        for val, f in zip(values, formats):
            fmt = {
                'type': 'cell',
                'criteria': 'equal to',
                'value': '"{}"'.format(val),
                'format': self._formats[f]
            }
            the_cmap.append(fmt)

        self._custom_colourmaps[name] = the_cmap

    def conditional_format_single(self, condition_dict, range_str, format_name, sheet_name=None):
        assert format_name in self._formats, \
            "Format name {} has not been defined: run add_format() first.".format(format_name)
        sheet = self.get_sheet(sheet_name)
        fmt = dict(condition_dict)
        fmt['format'] = self._formats[format_name]
        sheet.conditional_format(range_str, fmt)

    def conditional_format_colourmap(self, range_str, format_name, sheet_name=None):
        assert format_name in self._custom_colourmaps, \
            "Format name {} has not been defined: run add_{{quantitative,qualitative}}_colourmap() first.".format(
                format_name,

            )
        sheet = self.get_sheet(sheet_name)
        the_cmap = self._custom_colourmaps[format_name]
        for fmt in the_cmap:
            sheet.conditional_format(range_str, fmt)

    def apply_format_one_cell(self, row, col, format_name, sheet_name=None, keep_existing_format=True):
        # FIXME: see issue with string values
        sheet = self.get_sheet(sheet_name)
        # get
        try:
            x = sheet.table[row][col]
            fmt = x.format
            if hasattr(x, 'number'):
                val = x.number
            elif hasattr(x, 'string'):
                ## FIXME: the shared string table seems only to be populated at write time, so this doesn't work
                val = sheet.str_table.string_array[x.string]
            elif x.__class__.__name__ == 'Blank':
                val = None
            else:
                raise NotImplementedError("Unrecognised value type: %s", type(x))

        except KeyError:
            val = None
            fmt = None

        if fmt is None:
            fmt = self._formats[format_name]
        else:
            if keep_existing_format:
                # attempt to combine formats
                existing_fmt_name = self.formats_inverse[fmt]
                existing_fmt_dict = self._format_dicts[existing_fmt_name]
                new_fmt_dict = self._format_dicts[format_name]
                # new name
                combined_fmt_name = '_'.join([existing_fmt_name, format_name])

                if combined_fmt_name not in self._formats:
                    # make new combined format
                    combined_fmt_dict = dict(existing_fmt_dict, **new_fmt_dict)
                    self.add_format(combined_fmt_name, combined_fmt_dict)
                    logger.info(
                        "Combined previous format %s with new format %s under new name %s.",
                        existing_fmt_name,
                        format_name,
                        combined_fmt_name
                    )

                fmt = self._formats[combined_fmt_name]
            else:
                logger.warn("Existing format will be overwritten")
                fmt = self._formats[format_name]

        sheet.write(row, col, val, fmt)

    def merge_range(self, first_row, first_col, last_row, last_col, sheet_name=None, overwrite=False):
        sheet = self.get_sheet(sheet_name)
        # FIXME: this would be nice, but the shared_string table issue gets in the way
        raise NotImplementedError




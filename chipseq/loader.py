import os
from load_data import loader
import pandas as pd
from settings import CHIPSEQ_DIR



class ChIPSeqFileLocations(object):
    def __init__(self, root_dir, alignment_subdir=None, batch_id=None, tax_id=9606):
        self.root_dir = root_dir
        self.alignment_subdir = alignment_subdir
        if batch_id is None:
            self.batch_id = os.path.split(self.root_dir)[-1]
        else:
            self.batch_id = batch_id
        self.tax_id = tax_id

        self.base_dir = self.root_dir if alignment_subdir is None else os.path.join(self.root_dir, self.alignment_subdir)
        self.meta_file = os.path.join(self.root_dir, 'sources.csv')

        common_kwds = {
            'tax_id': self.tax_id,
            'batch_id': self.batch_id,
            'meta_fn': self.meta_file
        }

        self.params = {}

        self.params['default'] = {
            'base_dir': os.path.join(self.base_dir, 'macs2', 'default'),
            'run_type': 'default'
        }

        self.params['broad'] = {
            'base_dir': os.path.join(self.base_dir, 'macs2', 'broad'),
            'run_type': 'broad'
        }

        for k, v in self.params.items():
            v.update(common_kwds)

    def loader_kwargs(self, typ='default'):
        try:
            return self.params[typ]
        except KeyError:
            raise NotImplementedError("Unrecognised type: %s" % typ)


wtchg_p170710 = ChIPSeqFileLocations(
    root_dir=os.path.join(CHIPSEQ_DIR, 'wtchg_p170710'),
    alignment_subdir='human/bt2_alignment',
)

wtchg_p170710_pilot = ChIPSeqFileLocations(
    root_dir=os.path.join(CHIPSEQ_DIR, 'wtchg_p170710_pilot'),
    alignment_subdir='human/bt2_alignment',
)


MACS2_SAMPLE_LOOKUP = {
    '050': {wtchg_p170710_pilot: [
        'GBM050H3K4me3',
        'GBM050H3K36me3',
        'GBM050H3K27me3',
        'GBM050H3K27ac'
    ]},
    '017': {wtchg_p170710: [
        'GBM17H3K27ac',
        'GBM17H3K36me3',
        'GBM17H3K4me3',
        'Dura17H3K27ac',
        'Dura17H3K36me3',
        'Dura17H3K4me3',
    ]},
    '061': {wtchg_p170710: [
        'GBM61H3K27ac',
        'GBM61H3K36me3',
        'GBM61H3K4me3',
        'Dura61H3K27ac',
        'Dura61H3K36me3',
        'Dura61H3K4me3'
    ]},
    '3021': {wtchg_p170710: [
        '3021_1ScraBMI1',
        '3021_1ScraCHD7',
        '3021_1ScraK27',
        '3021_1ScraK4',
    ]},
    '3021_shBMI1': {wtchg_p170710: [
        '3021_1shBaCHD7',
        '3021_1shBaK27',
        '3021_1shBaK4',
    ]},
    '3021_shCHD7': {wtchg_p170710: [
        '3021_1shCaBMI1',
        '3021_1shCaK27',
        '3021_1shCaK4',
    ]},
    '3021_shBMI1shCHD7': {wtchg_p170710: [
        '3021_1shB+CaK27',
        '3021_1shB+CaK4',
    ]},
    '054': {
        wtchg_p170710_pilot: [
            'GBM054H3K4me3',
            'GBM54H3K36me3',
            'GBM054H3K27me3',
            'GBM054H3K27ac',
        ],
        wtchg_p170710: [
            'GBM54H3K27ac_nonSS',
            'GBM54H3K36me3_nonSS',
            'GBM54H3K4me3_nonSS',
            'Dura54H3K27ac',
            'Dura54H3K36me3',
            'Dura54H3K4me3',
        ]},
}


class MACS2PeaksLoader(loader.MultipleFileLoader):
    meta_col_filename = 'sample'
    row_indexed = False

    def __init__(self, run_type='default', *args, **kwargs):
        self.run_type = run_type
        if run_type == 'default':
            self.file_pattern = "*_peaks.narrowPeak"
            self.data_columns = ['chrom', 'start', 'end', 'peak_name', 'disp_q', 'null', 'fc', '-log10p', '-log10q',
                                 'rel_peak_pos']
        elif run_type == 'broad':
            self.file_pattern = "*_peaks.broadPeak"
            self.data_columns = ['chrom', 'start', 'end', 'peak_name', 'disp_q', 'null', 'fc', '-log10p', '-log10q']
        else:
            raise NotImplementedError("Unsupported run type %s." % run_type)
        super(MACS2PeaksLoader, self).__init__(*args, **kwargs)

    def load_one_file(self, fn):
        dat = pd.read_csv(fn, sep='\t', index_col=None, header=None)
        dat.columns = self.data_columns
        return dat

    def generate_input_path(self, fname):
        """
        Given the filename from the meta file, generate the path to the actual data (e.g. constructing subdir structure)
        """
        if self.file_pattern[0] == '*':
            to_append = self.file_pattern[1:]
        else:
            to_append = self.file_pattern

        return os.path.join(self.base_dir, "".join([fname, to_append]))


def load_macs2_by_patient(
        patient_ids,
        run_type='default',
        **kwargs
):
    """
    Load all MACS2 peak data associated with the patient ID(s) supplied
    :param patient_ids: Iterable or single int or char
    :param run_type: Specify the type of MACS2 run (default, broad)
    :param kwargs: Passed to the loader
    :return:
    """
    # ensure patient IDs are in correct form
    if patient_ids == 'all':
        patient_ids = MACS2_SAMPLE_LOOKUP.keys()
    elif not hasattr(patient_ids, '__iter__'):
        patient_ids = [patient_ids]

    # precompute the loaders required to avoid reloading multiple times
    by_loader = {}
    for pid in patient_ids:
        d = MACS2_SAMPLE_LOOKUP[pid]
        for ldr, slist in d.items():
            by_loader.setdefault(ldr, []).extend(slist)

    objs = []
    for ldr, samples in by_loader.items():
        the_kwargs = dict(ldr.loader_kwargs(run_type))
        the_kwargs.update(kwargs)
        objs.append(
            MACS2PeaksLoader(
                samples=samples,
                **the_kwargs
            )
        )

    if len(objs) > 1:
        res = loader.MultipleBatchLoader(objs)
    else:
        res = objs[0]
        # make samples column the meta index
        # res.meta.set_index('sample', inplace=True)

    # TODO: apply original ordering?

    return res
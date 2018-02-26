import math
import re
import os
from pyscript import jobs, sge

class TrimGaloreSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 120


class TrimGaloreBase(jobs.ArrayJob):
    title = 'trim_galore'
    required_args = ['read_dir']

    ext = r'fastq(\.gz)?'
    cleanup_regex = [
        (r'_[12]$', ''),
    ]

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$NAME', None),
        ('$FASTQ', '! -z'),
    ]
    param_delim = ':'
    core_cmd = "trim_galore -o {out_dir} $FASTQ"

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class TrimGalorePEBase(TrimGaloreBase, jobs.PEFastqBartsMultiLaneMixin):
    file_sep = ' '  # the character used to separate files of the same read number in different lanes

    def setup_params(self, read_dir, *args, **kwargs):
        super(TrimGalorePEBase, self).setup_params(read_dir, *args, **kwargs)
        # the format here is interleaved reads, e.g. LANE1_1 LANE1_2 LANE2_1 LANE2_2 etc...
        new_params = []
        for p in self.params:
            sep1 = p[1].split(' ')
            sep2 = p[2].split(' ')
            new_fq = [self.file_sep.join(t) for t in zip(sep1, sep2)]
            new_params.append(
                [p[0], self.file_sep.join(new_fq)] + p[3:]
            )
        self.params = new_params

    def set_default_extras(self):
        if '--paired' not in self.extra_args:
            self.extra_args.extend(
                ['--paired']
            )


class TrimGalorePEApocrita(TrimGaloreSgeRequirements, TrimGalorePEBase):
    pass


class TrimGalorePEBash(jobs.BashArrayJobMixin, TrimGalorePEBase):
    pass



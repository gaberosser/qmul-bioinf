import math
import re
import os
from pyscript import jobs, sge


class VCFRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 30


class VCFStatsBase(jobs.ArrayJob):
    title = 'vcf_stats'
    required_args = ['filename', 'n_lines_per_shard']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$FILE', '! -z'),
        ('$CONTIG', '! -z'),
        ('$STARTPOS', None),
        ('$ENDPOS', None),
    ]
    param_delim = ','

    core_cmd = """

    """

    # def set_default_extras(self):
    #     if '--gzip' not in self.extra_args:
    #         self.extra_args.extend(
    #             ['--gzip']
    #         )
    #     if '--split-files' not in self.extra_args:
    #         self.extra_args.extend(
    #             ['--split-files']
    #         )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['filename'], self.args['n_lines_per_shard'])
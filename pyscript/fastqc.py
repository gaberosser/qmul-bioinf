import math
import re
import os
from pyscript import jobs, sge

class FastQCSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 30


class FastQCBase(jobs.ArrayJob):
    title = 'fastqc'
    required_args = ['read_dir', 'threads']

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
    core_cmd = "fastqc -o {out_dir} $FASTQ"

    def set_default_extras(self):
        if '--noextract' not in self.extra_args:
            self.extra_args.extend(
                ['--noextract']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])

    def setup_params(self, read_dir, *args, **kwargs):
        """
        Generate the parameters array, run names and check output subdirs (if necessary)
        :return:
        """
        self.params = []
        self.run_names = []
        rr = re.compile(r'\.{ext}$'.format(ext=self.ext), flags=re.IGNORECASE)
        flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
        # check for existing output and identify pairs of files
        for t in flist:
            base = jobs.filename_to_name(t, self.ext, cleanup_regex_arr=self.cleanup_regex)

            if self.include is not None:
                if base not in self.include:
                    self.logger.info("Skipping file %s as it is not in the list of included files", base)
                    continue
            if self.exclude is not None:
                if base in self.exclude:
                    self.logger.info("Skipping file %s as it is in the list of excluded files", base)
                    continue

            self.params.append([base, os.path.abspath(os.path.join(read_dir, t))])
            self.run_names.append(base)


class FastQCApocrita(FastQCSgeRequirements, FastQCBase):
    pass


class FastQCBash(jobs.BashArrayJobMixin, FastQCBase):
    pass
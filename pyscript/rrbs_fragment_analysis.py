import math
import os
from pyscript import jobs, sge


class SgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        ram_per_core = 1.
        return "%dG" % ram_per_core

    @property
    def runtime_mins(self):
        return 60


class RRBSFragmentBase(jobs.ArrayJob):
    title = 'rrbs_fragment_analysis'
    required_args = ['bam_dir', 'threads', 'bed_fn', 'saf_fn']
    fc_args = ['-p']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$NAME', None),
        ('$BAM', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'

    @property
    def core_cmd(self):
        cmd = "samtools depth -a -b {cg_bed_fn} {bam_fn} | gzip > {outfn}\n"
        cmd += "featureCounts -a {saf_fn} -T 12 -F SAF {bam_fn} -o {outfn} {extra}"
        core_cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -1 $READS1 -2 $READS2 -o $SUBDIR {extra}"

    def set_default_extras(self):
        # default operation is PE reads, but can override this with --se
        if '--se' in self.extra_args:
            self.extra_args.pop(self.extra_args.index('--se'))
        else:
            self.extra_args.extend(
                ['-p']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['bam_dir'])

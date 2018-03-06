import math
from pyscript import jobs, sge


class BwaSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "2G"

    @property
    def runtime_mins(self):
        return 60


class BwaSEBase(jobs.ArrayJob):
    title = 'bwa_se'
    required_args = ['read_dir', 'threads', 'index']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$NAME', None),
        ('$READS', '! -z'),
        ('$OUTFILE', None),
    ]
    param_delim = ':'
    core_cmd = 'bwa aln {extra} -t {threads} {index} $READS > "$OUTFILE.sai"\n'
    core_cmd += 'bwa samse {index} "$OUTFILE.sai" $READS | samtools view -b > "$OUTFILE.bam"\n'
    core_cmd += 'samtools sort -@ {threads} "$OUTFILE.bam" > "$OUTFILE.sorted.bam"\n'
    core_cmd += 'rm "$OUTFILE.sai" "$OUTFILE.bam"\n'


    def set_default_extras(self):
        if '--srna' in self.extra_args:
            self.extra_args.remove('--srna')
            self.extra_args.extend(
                ['-n', '1', '-o', '0', '-e', '0', '-k', '1']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class BwaSEApocrita(BwaSgeRequirements, BwaSEBase, jobs.SEFastqFileIteratorMixin):
    pass


class BwaSEBash(jobs.BashArrayJobMixin, BwaSEBase, jobs.SEFastqFileIteratorMixin):
    pass


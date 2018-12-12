import math
from pyscript import jobs, sge


class Bowtie2SgeRequirements(sge.ApocritaArrayJobMixin):
    min_ram_gb = 4.
    @property
    def ram_per_core(self):
        nthread = int(self.args['threads'])
        # aim for minimum
        gb_per_core = int(math.ceil(self.min_ram_gb / float(nthread)))
        print "Num thread: %d" % nthread
        print "Gb per core: %d" % gb_per_core
        return "%dG" % gb_per_core

    @property
    def runtime_mins(self):
        return 120


class Bowtie2SEBase(jobs.ArrayJob):
    title = 'bt2_se'
    required_args = ['read_dir', 'threads', 'index']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$READS', '! -z'),
        ('$OUTFILE', None),
    ]
    param_delim = ':'
    core_cmd = 'bowtie2 -p {threads} -x {index} -U $READS {extra} | samtools view -b > ${{OUTFILE}}.bam 1> ${{OUTFILE}}.log\n'
    core_cmd += 'samtools sort -@ {threads} "${{OUTFILE}}.bam" > "${{OUTFILE}}.sorted.bam"\n'
    core_cmd += 'rm "${{OUTFILE}}.bam"\n'

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class Bowtie2SEApocrita(Bowtie2SgeRequirements, Bowtie2SEBase, jobs.SEFastqFileIteratorMixin):
    pass


class Bowtie2SEBash(jobs.BashArrayJobMixin, Bowtie2SEBase, jobs.SEFastqFileIteratorMixin):
    pass


class Bowtie2PEBase(jobs.ArrayJob):
    title = 'bt2_pe'
    required_args = ['read_dir', 'threads', 'index']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$OUTFILE', None),
    ]
    param_delim = ':'
    core_cmd = 'bowtie2 -p {threads} -x {index} -1 $READS1 -2 $READS2 {extra} | samtools view -b > ${{OUTFILE}}.bam\n'
    core_cmd += 'samtools sort -@ {threads} "${{OUTFILE}}.bam" > "${{OUTFILE}}.sorted.bam"\n'
    core_cmd += 'rm "${{OUTFILE}}.bam"\n'


    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class Bowtie2PEApocrita(Bowtie2SgeRequirements, Bowtie2PEBase, jobs.SEFastqFileIteratorMixin):
    pass


class Bowtie2PEBash(jobs.BashArrayJobMixin, Bowtie2PEBase, jobs.SEFastqFileIteratorMixin):
    pass


class Bowtie2MultilaneIlluminaPEApocrita(Bowtie2SgeRequirements, Bowtie2PEBase, jobs.PEFastqIlluminaMultiLaneMixin):
    file_sep = ','  # the character used to separate files of the same read number in different lanes


class Bowtie2MultilaneIlluminaPEBash(jobs.BashArrayJobMixin, Bowtie2PEBase, jobs.PEFastqIlluminaMultiLaneMixin):
    file_sep = ','  # the character used to separate files of the same read number in different lanes


def multilane_run(run_type):
    import argparse
    import sys
    import os

    if run_type == 'bash':
        cls = Bowtie2MultilaneIlluminaPEBash
    elif run_type == 'apocrita':
        cls = Bowtie2MultilaneIlluminaPEApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')

    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    required.add_argument("-x", "--index", help="Directory of pre-computed BT2 index", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'bt2_alignment')
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    if args.include is not None:
        args.include = args.include.split(',')
    if args.exclude is not None:
        args.exclude = args.exclude.split(',')

    obj = cls(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()
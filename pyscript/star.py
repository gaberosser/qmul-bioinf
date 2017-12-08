import math
from pyscript import jobs, sge


class StarSgeRequirements(sge.ApocritaArrayJobMixin):
    min_ram_gb = 32.
    estimated_runtime_per_core = 720.
    @property
    def ram_per_core(self):
        nthread = int(self.args['threads'])
        # aim for minimum 32Gb
        gb_per_core = int(math.ceil(self.min_ram_gb / float(nthread)))
        print "Num thread: %d" % nthread
        print "Gb per core: %d" % gb_per_core
        return "%dG" % gb_per_core

    @property
    def runtime_mins(self):
        nthread = int(self.args['threads'])
        # roughly 60 mins with 12 cores
        return self.estimated_runtime_per_core / float(nthread)


class StarBase(jobs.ArrayJob):
    title = 'star_multilane_pe'
    required_args = ['read_dir', 'threads', 'genomeDir']
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
        ('$GZ', None)
    ]
    param_delim = ':'

    core_cmd = "STAR {extra} --genomeDir {genomeDir} " + \
    "--readFilesIn $READS1 $READS2 --outFileNamePrefix $SUBDIR --runThreadN {threads}"

    def set_default_extras(self):
        if '--outSAMstrandField' not in self.extra_args:
            self.extra_args.extend(
                ['--outSAMstrandField', 'intronMotif']
            )
        if '--quantMode' not in self.extra_args:
            self.extra_args.extend(
                ['--quantMode', 'GeneCounts']
            )
        if '--outSAMtype' not in self.extra_args:
            self.extra_args.extend(
                ['--outSAMtype', 'BAM SortedByCoordinate']
            )
        if '--readFilesCommand' not in self.extra_args:
            self.extra_args.extend(
                ['$([ "$GZ" = "gz" ] && echo "--readFilesCommand zcat" || echo "")']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class StarMultilanePEApocrita(StarSgeRequirements, StarBase, jobs.PEFastqIlluminaMultiLaneMixin):
    pass


class StarMultilanePEBash(jobs.BashArrayJobMixin, StarBase, jobs.PEFastqIlluminaMultiLaneMixin):
    pass


class StarPEApocrita(StarBase, jobs.PEFastqFileIteratorMixin, StarSgeRequirements):
    pass


class StarPEBash(StarBase, jobs.PEFastqFileIteratorMixin, jobs.BashArrayJobMixin):
    pass


class StarSEApocrita(StarBase, jobs.SEFastqFileIteratorMixin, StarSgeRequirements):
    pass


class StarSEBash(StarBase, jobs.SEFastqFileIteratorMixin, jobs.BashArrayJobMixin):
    pass
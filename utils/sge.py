import os
import re
import subprocess
import math
import datetime
import csv

from utils.log import get_file_logger
from utils.output import unique_output_dir


def sge_submission_header(work_dir=None, threads=1, ram_per_core="1G", runtime="1:0:0", arr_size=1):
    header = [
        '#!/bin/sh',
        '#$ -j y # Join stdout and stderr'
    ]
    if work_dir is None:
        header.append('#$ -cwd # Use current directory as working directory')
    else:
        header.append("#$ -wd {work_dir} # Set the working directory for the job to the current directory".
                      format(work_dir=work_dir))

    header.append("#$ -pe smp {threads} # Request {threads} CPU cores".format(threads=threads))
    header.append("#$ -l h_rt={runtime}  # Runtime".format(runtime=runtime))
    header.append("#$ -l h_vmem={ram} # Request RAM / core".format(ram=ram_per_core))

    if arr_size is not None and arr_size > 1:
        header.append("#$ -t 1-{nfile}".format(nfile=arr_size))

    header.append(
        "# source bashrc to get correct path\n"
        ". $HOME/.bashrc"
    )

    return '\n'.join(header)


def sge_array_params_boilerplate(params_fn, argnames):
    sh = ['INPUTS=$(sed -n "${{SGE_TASK_ID}}p" {params_fn})'.format(params_fn=params_fn)]
    for i, nm in enumerate(argnames):
        sh.append('%s=$(echo $INPUTS | cut -d , -f %d)' % (nm, i + 1))
    return '\n'.join(sh)


def sge_tracking_files_boilerplate(submitted_fn, completed_fn):
    """
    Generate two pieces of code:
    - One to output a task number to the 'submitted tasks' file
    - One to check the exit code ($STATUS) and mark as complete if it is OK
    :param submitted_fn:
    :param completed_fn:
    :return:
    """
    sh_submit = 'printf "${SGE_TASK_ID}\\n" >> %s' % submitted_fn
    sh_complete = 'if [ $STATUS == 0 ]; then printf "${SGE_TASK_ID}\\n" >> %s; fi' % completed_fn
    return sh_submit, sh_complete


class SgeJob(object):
    title = None
    core_cmd = None
    log_dir = os.path.join(os.environ['HOME'], 'log')
    required_args = []
    require_empty_outdir = False
    create_outdir = False

    def __init__(
            self,
            out_dir=os.path.abspath('.'),
            extra_args=tuple(),
            **kwargs
    ):
        self.out_dir = os.path.abspath(out_dir)
        # initiate logger
        self.now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

        self.logger = None
        self.set_logger()

        self.args = {}
        self.extra_args = []
        self.set_args(*extra_args, **kwargs)

        if not self.check_inputs():
            raise Exception("Failed input checking. See logs for more details.")

        self.script_fn = os.path.join(self.out_dir, "%s.%s.sh" % (self.title, self.now_str))
        self.conf_fn = os.path.join(self.out_dir, "%s.%s.conf" % (self.title, self.now_str))
        self.submitted_fn = os.path.join(self.out_dir, "submitted.%s.%s" % (self.title, self.now_str))
        self.completed_fn = os.path.join(self.out_dir, "completed.%s.%s" % (self.title, self.now_str))

        self.write_conf_file()

        # write a basic opening log entry
        self.logger.info(
            "Job %s. Initialised logger from directory %s. The output and working directory is %s.",
            self.title,
            os.path.abspath(os.path.curdir),
            self.out_dir,
        )

        self.n_tasks = 1
        self.prepare_submission()

    def set_logger(self):
        self.logger = get_file_logger(
            name=self.title,
            filestem=os.path.join(self.log_dir, '%s.%s.log' % (self.title, self.now_str))
        )

    def check_inputs(self):
        """
        Here we can run preliminary checks, e.g. reference file existence
        If these fail, we do not continue with the job.
        Also log details here.
        :return: Boolean
        """
        # no checks - return True by default
        return True

    def set_args(self, *args, **kwargs):
        """
        Set the arguments that will be used later in the command(s) to be run.
        This must be defined in derived classes if more nuanced behaviour is required.
        :param args:
        :param kwargs:
        :return:
        """
        # required args
        self.args = dict([(k, kwargs.pop(k)) for k in self.required_args])
        self.logger.info("Arguments: " + '\n'.join(["%s = %s" % (k, str(self.args[k])) for k in self.args]))
        self.extra_args = args
        if len(args):
            self.logger.info("Additional arguments, passed directly to command: " + '\n'.join(self.extra_args))
        else:
            self.logger.info("No additional arguments")

    def write_conf_file(self):
        """
        Write configuration arguments to a file, so that we have a record of how this job was run.
        This is just a TSV
        """
        # write arguments to a file for reference

        with open(self.conf_fn, 'wb') as f:
            c = csv.writer(f, delimiter='\t')
            c.writerow(['Field', 'Value'])
            c.writerows([(k, str(v)) for k, v in self.args.items()])
        self.logger.info("Wrote config file to %s.", self.conf_fn)

    def prepare_submission(self, *args, **kwargs):
        """
        Here we run any preliminary tasks, e.g. checking for input files, getting number of tasks in array.
        We also need to create a params file if needed to initialise an array of tasks.
        """
        raise NotImplementedError()

    def create_submission_script(self):
        """
        Generate the submission script that we will finally submit to SGE
        :return:
        """

    def submit(self):
        subprocess.call(['qsub', self.script_fn])


class SgeArrayJob(SgeJob):
    def __init__(self, *args, **kwargs):
        self.params = []
        self.params_fn = None
        self.n_tasks = None
        super(SgeArrayJob, self).__init__(*args, **kwargs)
        self.create_params_file()

    def create_params_file(self):
        """
        Write parameters to a file, separated by a comma.
        :return:
        """
        self.params_fn = os.path.join(self.out_dir, "%s.%s.params" % (self.title, self.now_str))
        with open(self.params_fn, 'wb') as f:
            c = csv.writer(f, lineterminator='\n')  # IMPORTANT: the lineterminator command prevents carriage returns
            c.writerows(self.params)


class BamFileIteratorMixin(object):
    """
    Adds a method to build a list of bam files in the input directory.
    We also check for the output directory in each case to avoid overwriting
    """
    def generate_parameters_and_create_subdirs(self, cleanup_regex_arr=None):
        params = []
        seen = []
        rr = re.compile(r'\.bam$', flags=re.IGNORECASE)
        flist = [t for t in os.listdir(self.args['read_dir']) if re.search(rr, t)]
        # check for existing output and identify pairs of files
        for t in flist:
            base = re.sub(r'\.bam', '', t)
            # apply cleanup
            if cleanup_regex_arr is not None:
                for patt, repl in cleanup_regex_arr:
                    base = re.sub(patt, repl, base)
            out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
            self.logger.info("Input file %s. Cleaned filestem %s. Output subdir %s.", t, base, out_subdir)
            # if output file exists, log warning and skip
            if self.require_empty_outdir and os.path.isdir(out_subdir):
                if len(os.listdir(out_subdir)) > 0:
                    self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                    continue
                else:
                    self.logger.info("Using existing empty output subdir %s", out_subdir)

            if not os.path.exists(out_subdir) and self.create_outdir:
                os.makedirs(out_subdir)
                self.logger.info("Created output subdir %s", out_subdir)

            params.append([os.path.abspath(os.path.join(self.args['read_dir'], t)), out_subdir])
            seen.append(base)
        return dict(zip(seen, params))


class PEFastqFileIteratorMixin(object):
    """
    Adds a method to build a list of paired fastq files in the input directory.
    We also check for the output directory in each case to avoid overwriting
    """

    def parse_filename(self, filename, cleanup_regex_arr=None):
        raise NotImplementedError

    def generate_parameters_and_create_subdirs(self, cleanup_regex_arr=None):
        params = []
        seen = []
        rec = {}

        rr = re.compile(r'\.fastq(\.gz)?$', flags=re.IGNORECASE)
        flist = [t for t in os.listdir(self.args['read_dir']) if re.search(rr, t)]
        # check for existing output and identify pairs of files
        for t in flist:
            base, read_num = self.parse_filename(t, cleanup_regex_arr=cleanup_regex_arr)
            rec.setdefault(base, {})

            if 'out_subdir' not in rec[base]:
                out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
                if self.require_empty_outdir and os.path.isdir(out_subdir):
                    if len(os.listdir(out_subdir)) > 0:
                        self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                        continue
                    else:
                        self.logger.info("Using existing empty output subdir %s", out_subdir)
                if self.create_outdir and not os.path.exists(out_subdir):
                    os.makedirs(out_subdir)
                    self.logger.info("Created output subdir %s", out_subdir)
                rec[base]['out_subdir'] = out_subdir

            rec[base][read_num] = os.path.abspath(os.path.join(self.args['read_dir'], t))

        for base, p in rec.items():
            if len(p) == 0:
                # skip
                continue
            params.append([p[1], p[2], p['out_subdir']])
            seen.append(base)

        return dict(zip(seen, params))


class PEFastqIlluminaIteratorMixin(PEFastqFileIteratorMixin):
    def parse_filename(self, filename, cleanup_regex_arr=None):
        base = re.sub(r'\.fastq(\.gz)?$', '', filename)
        # get read number
        read_num = int(base[-1])
        # apply cleanup
        if cleanup_regex_arr is not None:
            for patt, repl in cleanup_regex_arr:
                base = re.sub(patt, repl, base)
        return base, read_num


class SEFastqFileIteratorMixin(object):
    """
    Adds a method to build a list of single fastq files in the input directory.
    We also check for the output directory in each case to avoid overwriting
    """

    def parse_filename(self, filename, cleanup_regex_arr=None):
        base = re.sub(r'\.fastq(\.gz)?$', '', filename)
        # apply cleanup
        if cleanup_regex_arr is not None:
            for patt, repl in cleanup_regex_arr:
                base = re.sub(patt, repl, base)
        return base

    def generate_parameters_and_create_subdirs(self, cleanup_regex_arr=None):
        params = []
        seen = []

        rr = re.compile(r'\.fastq(\.gz)?$', flags=re.IGNORECASE)
        flist = [t for t in os.listdir(self.args['read_dir']) if re.search(rr, t)]
        # check for existing output and identify pairs of files
        for t in flist:
            base = self.parse_filename(t, cleanup_regex_arr=cleanup_regex_arr)
            seen.append(base)
            out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
            if self.require_empty_outdir and os.path.isdir(out_subdir):
                if len(os.listdir(out_subdir)) > 0:
                    self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                    continue
                else:
                    self.logger.info("Using existing empty output subdir %s", out_subdir)
            if self.create_outdir and not os.path.exists(out_subdir):
                os.makedirs(out_subdir)
                self.logger.info("Created output subdir %s", out_subdir)
            params.append([os.path.abspath(os.path.join(self.args['read_dir'], t)), out_subdir])

        return dict(zip(seen, params))


class CufflinksSgeJob(SgeArrayJob, BamFileIteratorMixin):
    title = 'cufflinks'
    required_args = ['read_dir', 'threads', 'library_type', 'GTF']
    require_empty_outdir = False
    create_outdir = True

    def check_inputs(self):
        if not os.path.isfile(self.args['GTF']):
            self.logger.error("Unable to find specified GTF file %s", self.args['GTF'])
            return False
        return True

    def prepare_submission(self):
        cleanup_regex_arr = [
            (r'Aligned\.sortedByCoord\.out', ''),
            (r'Aligned\.out', ''),
        ]
        res = self.generate_parameters_and_create_subdirs(cleanup_regex_arr)
        self.params = res.values()

        # log the filelist
        self.logger.info("Found %d BAM files: %s.", len(self.params), ', '.join(res.keys()))
        self.n_tasks = len(self.params)

    def create_submission_script(self):
        # parameter names as they will appear in the bash script
        param_names = ['BAM', 'SUBDIR']

        # generate the main command
        cmd = "cufflinks --no-update-check -G {GTF} -p {threads} -o $SUBDIR --library-type {library_type} {extra} $BAM"\
            .format(extra=' '.join(self.extra_args), **self.args)

        sh = []

        sh.append(
            sge_submission_header(
                work_dir=self.out_dir,
                threads=self.args['threads'],
                ram_per_core='512M',
                runtime="0:30:0",
                arr_size=self.n_tasks
            )
        )
        sh.append(sge_array_params_boilerplate(self.params_fn, param_names))

        submit, complete = sge_tracking_files_boilerplate(self.submitted_fn, self.completed_fn)
        sh.append(submit)

        sh.append("""
        if [[ -f $BAM && ! -z $SUBDIR ]]; then
            {cmd}
            STATUS=$?
        else
            echo "Unable to execute run ${{SGE_TASK_ID}} as the read file did not exist or the output dir variable is empty."
            echo "Read file: $BAM"
            echo "Output dir: $SUBDIR"
            STATUS=1  # set this so that the task is not marked as completed
        fi
        """.format(cmd=cmd))

        sh.append(complete)

        # create the actual submission script
        with open(self.script_fn, 'wb') as f:
            f.write('\n'.join(sh))

        self.logger.info("Cluster submission script written to %s: \n%s", self.script_fn, sh)


class SalmonIlluminaPESgeJob(SgeArrayJob, PEFastqIlluminaIteratorMixin):
    title = 'salmon'
    required_args = ['read_dir', 'threads', 'library_type', 'index_dir']
    require_empty_outdir = False
    create_outdir = False

    def check_inputs(self):
        if not os.path.isdir(self.args['index_dir']):
            self.logger.error("Unable to find specified index directory %s", self.args['index_dir'])
            return False
        return True

    def prepare_submission(self):
        cleanup_regex_arr = [
            (r'_[12]$', ''),
        ]
        res = self.generate_parameters_and_create_subdirs(cleanup_regex_arr)
        self.params = res.values()

        # log the filelist
        self.logger.info("Found %d fastq pairs: %s.", len(self.params), ', '.join(res.keys()))
        self.n_tasks = len(self.params)

    def create_submission_script(self):
        # parameter names as they will appear in the bash script
        param_names = ['READ1', 'READ2', 'SUBDIR']

        # generate the main command
        cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -1 $READ1 -2 $READ2 -o $SUBDIR".format(**self.args)
        if len(self.extra_args):
            cmd += " {extra}".format(extra=' '.join(self.extra_args))

        sh = []

        # NB: uses one more core than the number we request (if it can)
        # Aim to provide 4Gb overall
        eff_threads = int(self.args['threads']) + 1
        ram_per_core = 1.
        if (ram_per_core * eff_threads) < 4:
            ram_per_core = int(math.ceil(4. / float(eff_threads)))

        sh.append(
            sge_submission_header(
                work_dir=self.out_dir,
                threads=eff_threads,
                ram_per_core='%dG' % ram_per_core,
                runtime="0:20:0",
                arr_size=self.n_tasks
            )
        )
        sh.append(sge_array_params_boilerplate(self.params_fn, param_names))

        submit, complete = sge_tracking_files_boilerplate(self.submitted_fn, self.completed_fn)
        sh.append(submit)

        sh.append("""
        if [[ -f $READ1 && -f $READ2 && ! -z $SUBDIR ]]; then
            {cmd}
            STATUS=$?
        else
            echo "Unable to execute run ${{SGE_TASK_ID}} as the read file did not exist or the output dir variable is empty."
            echo "Read files: $READ1 $READ2"
            echo "Output dir: $SUBDIR"
            STATUS=1  # set this so that the task is not masked as completed
        fi
        """.format(cmd=cmd))

        sh.append(complete)

        # create the actual submission script
        with open(self.script_fn, 'wb') as f:
            f.write('\n'.join(sh))

        self.logger.info("Cluster submission script written to %s: \n%s", self.script_fn, sh)


class SalmonIlluminaSESgeJob(SgeArrayJob, SEFastqFileIteratorMixin):
    title = 'salmon'
    required_args = ['read_dir', 'threads', 'library_type', 'index_dir']
    require_empty_outdir = False
    create_outdir = False

    def check_inputs(self):
        if not os.path.isdir(self.args['index_dir']):
            self.logger.error("Unable to find specified index directory %s", self.args['index_dir'])
            return False
        return True

    def prepare_submission(self):
        res = self.generate_parameters_and_create_subdirs()
        self.params = res.values()

        # log the filelist
        self.logger.info("Found %d fastq files: %s.", len(self.params), ', '.join(res.keys()))
        self.n_tasks = len(self.params)

    def create_submission_script(self):
        # parameter names as they will appear in the bash script
        param_names = ['READ', 'SUBDIR']

        # generate the main command
        cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -r $READ -o $SUBDIR".format(**self.args)
        if len(self.extra_args):
            cmd += " {extra}".format(extra=' '.join(self.extra_args))

        sh = []

        # NB: uses one more core than the number we request (if it can)
        # Aim to provide 4Gb overall
        eff_threads = int(self.args['threads']) + 1
        ram_per_core = 1.
        if (ram_per_core * eff_threads) < 4:
            ram_per_core = int(math.ceil(4. / float(eff_threads)))

        sh.append(
            sge_submission_header(
                work_dir=self.out_dir,
                threads=eff_threads,
                ram_per_core='%dG' % ram_per_core,
                runtime="0:20:0",
                arr_size=self.n_tasks
            )
        )
        sh.append(sge_array_params_boilerplate(self.params_fn, param_names))

        submit, complete = sge_tracking_files_boilerplate(self.submitted_fn, self.completed_fn)
        sh.append(submit)

        sh.append("""
        if [[ -f $READ && ! -z $SUBDIR ]]; then
            {cmd}
            STATUS=$?
        else
            echo "Unable to execute run ${{SGE_TASK_ID}} as the read file did not exist or the output dir variable is empty."
            echo "Read file: $READ"
            echo "Output dir: $SUBDIR"
            STATUS=1  # set this so that the task is not masked as completed
        fi
        """.format(cmd=cmd))

        sh.append(complete)

        # create the actual submission script
        with open(self.script_fn, 'wb') as f:
            f.write('\n'.join(sh))

        self.logger.info("Cluster submission script written to %s: \n%s", self.script_fn, sh)

import os
from utils.log import get_file_logger, get_console_logger
import datetime
import csv
import subprocess
import re


def array_params_boilerplate(params_fn, argnames, sep=','):
    """
    Produces the script lines necessary to read multiple parameters from an input file
    :param params_fn: Path to the file containing params (one per line)
    :param argnames:
    :param sep: Optionally specify a different delimiter to the default (comma)
    :return:
    """
    sh = ['INPUTS=$(sed -n "${{SGE_TASK_ID}}p" {params_fn})'.format(params_fn=params_fn)]
    for i, nm in enumerate(argnames):
        sh.append('%s=$(echo $INPUTS | cut -d %s -f %d)' % (nm, sep, i + 1))
    return '\n'.join(sh)


def tracking_files_boilerplate(submitted_fn, completed_fn):
    """
    Generate two pieces of code:
    - One to output a task number to the 'submitted tasks' file
    - One to check the exit code ($STATUS) and mark as complete if it is OK
    We assume two variables are available: $SRR_TASK_ID (integer giving the array index) and $ID (string giving the
    name of this run).
    :param submitted_fn:
    :param completed_fn:
    :return:
    """
    sh_submit = 'printf "${SGE_TASK_ID}, $ID\\n" >> %s' % submitted_fn
    sh_complete = 'if [ $STATUS == 0 ]; then printf "${SGE_TASK_ID}, $ID\\n" >> %s; fi' % completed_fn
    return sh_submit, sh_complete


class BashJobMixin(object):
    def shebang(self):
        """
        Generate the shebang line(s) for the top of the script
        :return:
        """
        # Default behaviour: run with bash
        return "#!/bin/bash\n\n"

    def submit(self):
        print "Bash submit()"
        subprocess.call(['bash', self.script_fn])


class BashArrayJobMixin(BashJobMixin):
    def generate_script(self):
        """
        Generate the full script object (in a list, to be joined with newlines)
        Set the self.sh variable
        We can make use of self.shebang(), self.script_body()
        """
        sh = [
            self.shebang(),
            "for SGE_TASK_ID in $(seq {n}); do".format(n=len(self.params)),
        ] + self.script_body() + ['done']
        self.sh = sh


class Job(object):
    title = None
    core_cmd = None
    log_dir = os.path.join(os.environ['HOME'], 'log')
    required_args = []
    create_outdir = False

    def __init__(
            self,
            out_dir=os.path.abspath('.'),
            extra_args=tuple(),
            include=None,
            exclude=None,
            **kwargs
    ):
        self.out_dir = os.path.abspath(out_dir)
        self.now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.include = include
        self.exclude = exclude

        # initiate logger
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

        self.sh = []

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
        created_dir = False
        if not os.path.isdir(self.log_dir):
            try:
                os.makedirs(self.log_dir)
                created_dir = True
            except Exception:
                raise AttributeError("Failed to create log output directory %s" % self.log_dir)
        self.logger = get_file_logger(
            name=self.title,
            filestem=os.path.join(self.log_dir, '%s.%s.log' % (self.title, self.now_str))
        )
        if created_dir:
            self.logger.info("Created log output directory %s.", self.log_dir)

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
        self.extra_args = list(args)
        if len(args):
            self.logger.info("Additional arguments, passed directly to command: " + '\n'.join(self.extra_args))
        else:
            self.logger.info("No additional arguments")
        self.set_default_extras()

    def set_default_extras(self):
        """
        Here we can optionally append to or modify the contents of self.extra_args.
        :return:
        """

    def write_conf_file(self):
        """
        Write configuration arguments to a file, so that we have a record of how this job was run.
        This is just a TSV
        """
        # write arguments to a file for reference

        with open(self.conf_fn, 'wb') as f:
            c = csv.writer(f, delimiter='\t')
            c.writerow(['Field', 'Value'])
            c.writerow(['out_dir', self.out_dir])
            c.writerows([(k, str(v)) for k, v in self.args.items()])
        self.logger.info("Wrote config file to %s.", self.conf_fn)

    def prepare_submission(self, *args, **kwargs):
        """
        Here we run any preliminary tasks, e.g. checking for input files, getting number of tasks in array.
        We also need to create a params file if needed to initialise an array of tasks.
        Basic operation requires nothing here
        """
        pass

    def shebang(self):
        """
        Implement this in a mixin.
        :return:
        """
        raise NotImplementedError()

    def core_command(self):
        """
        Create the core lines of code necessary to run the script.
        Boiler plate will be attached around around this.
        Assume we have access to the variables as usual ($VAR_NAME). Leave any args in their `format` token ({varA}).
        Extra arguments are all in {extra}.
        :return: string
        """
        cmd = self.core_cmd
        subs = {}
        if '{extra}' in cmd:
            subs['extra'] = ' '.join(self.extra_args) if len(self.extra_args) else ''
        if '{out_dir}' in cmd:
            subs['out_dir'] = self.out_dir
        subs.update(self.args)
        return cmd.format(**subs)

    def generate_script(self):
        """
        Generate the submission script that we will finally submit to SGE
        Store this in self.sh as a list
        :return:
        """
        cmd = self.core_command()
        self.sh = [self.shebang(), cmd, "STATUS=$?"]

    def create_script(self):
        """
        Write self.sh to self.script_fn
        :return:
        """
        self.generate_script()
        s = '\n'.join(self.sh)
        with open(self.script_fn, 'wb') as f:
            f.write(s)

        self.logger.info("Cluster submission script written to %s: \n%s", self.script_fn, s)


class ArrayJob(Job):
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        # ('$param_name', '! -z')
    ]
    param_delim = ','

    def __init__(self, *args, **kwargs):
        self.params = []
        self.params_fn = None
        self.run_names = []
        super(ArrayJob, self).__init__(*args, **kwargs)
        self.create_params_file()

    def create_params_file(self):
        """
        Write parameters to a file, separated by a comma.
        Params have already been set by this stage
        :return:
        """
        self.params_fn = os.path.join(self.out_dir, "%s.%s.params" % (self.title, self.now_str))
        self.logger.info("Generating parameters file for %d runs: %s", len(self.run_names), ', '.join(self.run_names))
        with open(self.params_fn, 'wb') as f:
            c = csv.writer(f, delimiter=self.param_delim, lineterminator='\n')  # IMPORTANT: the lineterminator command prevents carriage returns
            c.writerows(self.params)

    def script_body(self):
        # generate core
        cmd = self.core_command()

        # submission and completion lines
        submit, complete = tracking_files_boilerplate(self.submitted_fn, self.completed_fn)

        # reading parameters from file
        argnames = [t[0].replace('$', '') for t in self.parameters]
        read_params_str = array_params_boilerplate(self.params_fn, argnames, sep=self.param_delim)

        sh = [
            # self.shebang(),
            # "for SGE_TASK_ID in $(seq {n}); do".format(n=len(self.params)),
            read_params_str,
            submit,
            '\n'
        ]

        # checking parameters (if necessary)
        check = False
        # if any parameters have associated checking criteria, we need to use an if statement
        if len(self.parameters):
            checking = []
            failing = []
            for p, chk in self.parameters:
                if chk is not None:
                    checking.append('%s %s' % (chk, p))
                    failing.append('echo "{varname}": {var}'.format(varname=p.replace('$', ''), var=p))
            if len(checking):
                check = True
                failing = ['echo "Unable to execute run ${SGE_TASK_ID} as a validation check failed"'] + failing

        if check:
            sh.append('if [[ {check_str} ]]; then'.format(check_str=' && '.join(checking)))
        sh.append(cmd)
        sh.append("STATUS=$?")
        if check:
            sh.append("else")
            sh.append('\n'.join(failing))
            sh.append("STATUS=1  # set this so that the task is not marked as completed")
            sh.append("fi")
        sh.append(complete)
        # sh.append('done')

        return sh

    def generate_script(self):
        """
        Implement this in a mixin
        :return:
        """
        raise NotImplementedError()


def filename_to_name(fname, ext, cleanup_regex_arr=None):
    base = re.sub(r'\.{ext}'.format(ext=ext), '', fname)
    # apply cleanup
    if cleanup_regex_arr is not None:
        for patt, repl in cleanup_regex_arr:
            base = re.sub(patt, repl, base)

    return base


def filename_to_name_and_read_num(fname, ext, cleanup_regex_arr=None):
    base = re.sub(r'\.{ext}'.format(ext=ext), '', fname)
    # get read number - this is the last remaining character
    read_num = int(base[-1])
    # apply cleanup
    if cleanup_regex_arr is not None:
        for patt, repl in cleanup_regex_arr:
            base = re.sub(patt, repl, base)
    return base, read_num


class FileIteratorMixin(object):
    cleanup_regex = None
    skip_non_empty = False
    create_subdirs = False
    ext = ''

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
            base = filename_to_name(t, self.ext, cleanup_regex_arr=self.cleanup_regex)

            if self.include is not None:
                if base not in self.include:
                    self.logger.info("Skipping file %s as it is not in the list of included files", base)
                    continue
            if self.exclude is not None:
                if base in self.exclude:
                    self.logger.info("Skipping file %s as it is in the list of excluded files", base)
                    continue

            out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
            # if output file exists, log warning and skip
            if self.skip_non_empty and os.path.isdir(out_subdir):
                if len(os.listdir(out_subdir)) > 0:
                    self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                    continue
                else:
                    self.logger.info("Using existing empty output subdir %s", out_subdir)

            if not os.path.exists(out_subdir) and self.create_subdirs:
                os.makedirs(out_subdir)
                self.logger.info("Created output subdir %s", out_subdir)

            self.params.append([base, os.path.abspath(os.path.join(read_dir, t)), out_subdir])
            self.run_names.append(base)


class PairedFileIteratorMixin(FileIteratorMixin):

    def setup_params(self, read_dir, *args, **kwargs):
        self.params = []
        self.run_names = []
        rec = {}

        rr = re.compile(r'\.{ext}$'.format(ext=self.ext), flags=re.IGNORECASE)
        flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
        # check for existing output and identify pairs of files
        for t in flist:
            base, read_num = filename_to_name_and_read_num(t, self.ext, cleanup_regex_arr=self.cleanup_regex)

            if self.include is not None:
                if base not in self.include:
                    self.logger.info("Skipping file %s as it is not in the list of included files", base)
                    continue
            if self.exclude is not None:
                if base in self.exclude:
                    self.logger.info("Skipping file %s as it is in the list of excluded files", base)
                    continue

            rec.setdefault(base, {})

            if 'out_subdir' not in rec[base]:
                out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
                if self.skip_non_empty and os.path.isdir(out_subdir):
                    if len(os.listdir(out_subdir)) > 0:
                        self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                        continue
                    else:
                        self.logger.info("Using existing empty output subdir %s", out_subdir)

                if self.create_subdirs and not os.path.exists(out_subdir):
                    os.makedirs(out_subdir)
                    self.logger.info("Created output subdir %s", out_subdir)
                rec[base]['out_subdir'] = out_subdir

            rec[base][read_num] = os.path.abspath(os.path.join(read_dir, t))

        for base, p in rec.items():
            if len(p) == 0:
                # skip
                continue
            self.params.append([base, p[1], p[2], p['out_subdir']])
            self.run_names.append(base)


class BamFileIteratorMixin(FileIteratorMixin):
    ext = 'bam'


class PEFastqFileIteratorMixin(PairedFileIteratorMixin):
    ext = r'fastq(\.gz)?'
    cleanup_regex = [
        (r'_[12]$', ''),
    ]

    def setup_params(self, *args, **kwargs):
        super(PEFastqFileIteratorMixin, self).setup_params(*args, **kwargs)
        print self.params
        for p_arr in self.params:
            if re.search(r'\.gz$', p_arr[1], flags=re.IGNORECASE):
                p_arr.append('gz')
            else:
                p_arr.append('fastq')


class SEFastqFileIteratorMixin(FileIteratorMixin):
    ext = r'fastq(\.gz)?'
    cleanup_regex = [
        (r'_[12]$', ''),
    ]

    def setup_params(self, *args, **kwargs):
        super(SEFastqFileIteratorMixin, self).setup_params(*args, **kwargs)
        for p_arr in self.params:
            if re.search(r'\.gz$', p_arr[1], flags=re.IGNORECASE):
                p_arr.append('gz')
            else:
                p_arr.append('fastq')


class PEFastqIlluminaMultiLaneMixin(PairedFileIteratorMixin):
    ext = r'fastq(\.gz)?'
    cleanup_regex = [
        (r'_[12]$', ''),
        ('^WTCHG_[0-9]+_', ''),  # specific to WTCHG, but otherwise ignored
    ]
    file_sep = ','  # the character used to separate files of the same read number in different lanes

    def setup_params(self, read_dir, *args, **kwargs):
        self.params = []
        self.run_names = []

        rr = re.compile(r'\.{ext}$'.format(ext=self.ext), flags=re.IGNORECASE)
        root_dir = os.path.abspath(read_dir)
        dlist = [os.path.join(root_dir, t) for t in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, t))]

        to_join = {}
        rec = {}

        for d in dlist:
            flist = [t for t in os.listdir(d) if re.search(rr, t)]
            for t in flist:
                base, read_num = filename_to_name_and_read_num(t, self.ext, cleanup_regex_arr=self.cleanup_regex)

                if self.include is not None:
                    if base not in self.include:
                        self.logger.info("Skipping file %s as it is not in the list of included files", base)
                        continue
                if self.exclude is not None:
                    if base in self.exclude:
                        self.logger.info("Skipping file %s as it is in the list of excluded files", base)
                        continue

                to_join.setdefault(base, {}).setdefault(d, {})[read_num] = os.path.join(d, t)

        n = None
        self.logger.info(
            "Identified %d lanes that can be combined: %s", len(to_join), ', '.join(to_join.keys())
        )

        for base, d in to_join.items():
            # seen.append(read_id)
            x = d.values()[0]
            # get all valid pairs / singles
            if len(x) != 2:
                self.logger.error("Found %d corresponding reads - expected 2.", len(x))
                raise ValueError("Incorrect number of corresponding reads: %d" % len(x))

            self.logger.info("Read group %s.", base)

            # check the number of directories - this should be consistent across all reads
            if n is None:
                n = len(d)
            else:
                if len(d) != n:
                    raise AttributeError(
                        "Previous read group had %d matching directories; this one has %d." % (n, len(d)))

            # check output subdirectory and create if necessary

            out_subdir = os.path.abspath(os.path.join(self.out_dir, base))
            if self.skip_non_empty and os.path.isdir(out_subdir):
                if len(os.listdir(out_subdir)) > 0:
                    self.logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                    continue
                else:
                    self.logger.info("Using existing empty output subdir %s", out_subdir)

            if self.create_subdirs and not os.path.exists(out_subdir):
                os.makedirs(out_subdir)
                self.logger.info("Created output subdir %s", out_subdir)
            rec.setdefault(base, {})
            rec[base]['out_subdir'] = out_subdir

            for the_dir, the_dict in d.items():
                for read_num, the_file in the_dict.items():
                    rec[base].setdefault(read_num, []).append(the_file)

            this_param = [base]
            for i in [1, 2]:
                # join equivalent read files with a space
                this_param.append(self.file_sep.join(rec[base][i]))
            this_param.append(out_subdir)
            if re.search(r'\.gz$', rec[base][1][0], flags=re.IGNORECASE):
                this_param.append('gz')
            else:
                this_param.append('fastq')
            self.params.append(this_param)
            self.run_names.append(base)


class SraRunIteratorMixin(object):
    def setup_params(self, project_id, *args, **kwargs):
        """
        Generate the parameters array and run names
        :return:
        """
        CMD = "esearch -db sra -query {pid} | efetch --format runinfo | cut -d ',' -f 1 | grep 'SRR'"
        srr_list = subprocess.check_output(CMD.format(pid=project_id), shell=True)
        self.run_names = [t for t in srr_list.split('\n') if len(t)]
        self.params = [[t, t] for t in self.run_names]


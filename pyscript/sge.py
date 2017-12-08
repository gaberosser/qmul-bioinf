import subprocess


def apocrita_submission_header(work_dir=None, threads=1, ram_per_core="1G", runtime_mins=120, arr_size=1):

    # format the runtime variable
    rt_min = int(runtime_mins)
    if rt_min > 59:
        rt_hr = rt_min / 60  # modulo division
        rt_min = rt_min % 60
    else:
        rt_hr = 0

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
    header.append("#$ -l h_rt={rt_hr}:{rt_min}:0  # Runtime".format(rt_hr=rt_hr, rt_min=rt_min))
    header.append("#$ -l h_vmem={ram} # Request RAM / core".format(ram=ram_per_core))

    if arr_size is not None and arr_size > 1:
        header.append("#$ -t 1-{nfile}".format(nfile=arr_size))

    header.append(
        "# source bashrc to get correct path\n"
        ". $HOME/.bashrc\n"
        "# load modules (this script should be in the path)\n"
        ". load_modules.sh"
    )

    return '\n'.join(header)


class ApocritaJobMixin(object):

    def shebang(self):
        """
        Generate the shebang line(s) for the top of the script
        :return:
        """
        # Default behaviour: run with bash
        return apocrita_submission_header(
            work_dir=self.out_dir,
            threads=self.args['threads'],
            ram_per_core=self.ram_per_core,
            runtime_mins=self.runtime_mins,
            arr_size=len(self.params)
        )

    @property
    def ram_per_core(self):
        """
        Here we can perform more sophisticated estimates of the RAM required.
        For example, we might wish to maintain a minimum amount of total RAM available.
        :return:
        """
        return '1G'

    @property
    def runtime_mins(self):
        """
        Here we can perform more sophisticated estimates of the time required.
        :return:
        """
        return 120

    def submit(self):
        print "Apocrita submit()"
        subprocess.call(['qsub', self.script_fn])


class ApocritaArrayJobMixin(ApocritaJobMixin):
    def generate_script(self):
        """
        Generate the full script object (in a list, to be joined with newlines)
        Set the self.sh variable
        We can make use of self.shebang(), self.script_body()
        """
        self.sh = [self.shebang()] + self.script_body()

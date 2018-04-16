import math
import csv
import os
from pyscript import jobs, sge


class MACS2SgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "2G"

    @property
    def runtime_mins(self):
        return 20


class MACS2Base(jobs.ArrayJob):
    title = 'MACS2'
    required_args = ['config_file']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$TARGET', '! -z'),
        ('$CONTROL', None),
    ]
    param_delim = ':'
    core_cmd = """
    macs2 callpeak -t $TARGET -n $ID {extra} --outdir {out_dir} $([ ! -z "$CONTROL" ] && echo "-c $CONTROL")
    """

    def prepare_submission(self, *args, **kwargs):
        self.params = []
        self.run_names = []
        with open(self.args['config_file'], 'rb') as f:
            c = csv.DictReader(f)
            for fld in ['name', 'target']:
                if fld not in c.fieldnames:
                    raise ValueError("Config file must contain at least the fields 'name' and 'target'.")
                if 'control' not in c.fieldnames:
                    self.logger.warn("No control samples supplied in the config; running without.")
            run_info = list(c)
        self.logger.info("Found %d runs in the supplied config file %s", len(run_info), self.args['config_file'])
        for row in run_info:
            if 'control' in row and row['control'] != '':
                ctrl = os.path.abspath(row['control'])
            else:
                ctrl = ''
            self.params.append(
                [row['name'], os.path.abspath(row['target']), ctrl]
            )
            self.run_names.append(row['name'])


class MACS2Apocrita(MACS2SgeRequirements, MACS2Base):
    pass


class MACS2Bash(jobs.BashArrayJobMixin, MACS2Base):
    pass


def multilane_run(run_type):
    import argparse
    import sys

    if run_type == 'bash':
        cls = MACS2Bash
    elif run_type == 'apocrita':
        cls = MACS2Apocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("-p", "--threads", help="Number of threads", default='1')

    parser.add_argument("-c", "--config_file", help="Path to config. file", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()
    print args
    print extra

    if not os.path.exists(args.config_file):
        raise ValueError("No config file found at %s" % args.config_file)

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(os.path.abspath(os.getcwd()), 'macs2')
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    obj = cls(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()
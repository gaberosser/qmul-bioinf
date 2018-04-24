import math
import csv
import os
import re
from pyscript import jobs, sge


class MACS2CallPeaksSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "8G"

    @property
    def runtime_mins(self):
        return 240


class MACS2BdgCmpSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "4G"

    @property
    def runtime_mins(self):
        return 120


class MACS2CallPeaksBase(jobs.ArrayJob):
    title = 'macs2_callpeaks'
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


class MACS2CallPeaksApocrita(MACS2CallPeaksSgeRequirements, MACS2CallPeaksBase):
    core_cmd = """
    cp $TARGET $TMPDIR
    TARGET="$TMPDIR/$TARGET"
    TARGET="$TMPDIR/$(basename $TARGET)"
    if [ ! -z $CONTROL ]; then
        cp $CONTROL $TMPDIR
        CONTROL="$TMPDIR/$(basename $CONTROL)"
    fi
    """ + MACS2CallPeaksBase.core_cmd


class MACS2CallPeaksBash(jobs.BashArrayJobMixin, MACS2CallPeaksBase):
    pass


class MACS2BdgCmpBase(jobs.ArrayJob):
    title = 'macs2_bdgcmp'
    required_args = ['macs_dir']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$TARGET', '! -z'),
        ('$CONTROL', '! -z'),
        ('$GZ', None)
    ]
    param_delim = ':'
    core_cmd = """
    if [ $GZ == "gz" ]; then
        # decompress the input files, keeping the existing zipped version
        # we'll delete these uncompressed files later
        echo "Decompressing target file $TARGET"
        gzip -kd $TARGET
        TARGET=$(echo $TARGET | sed 's/\.gz$//')
        echo "New target file is at $TARGET"

        echo "Decompressing control file $CONTROL"
        gzip -kd $CONTROL
        CONTROL=$(echo $CONTROL | sed 's/\.gz$//')
        echo "New control file is at $CONTROL"
    fi
    macs2 bdgcmp -t $TARGET -c $CONTROL {extra} --outdir {out_dir} --o-prefix $ID

    if [ $GZ == "gz" ]; then
        # delete uncompressed files
        echo "Deleting uncompressed files $TARGET and $CONTROL"
        rm $TARGET
        rm $CONTROL
    fi
    """

    def set_default_extras(self):
        if '-m' not in self.extra_args:
            self.extra_args.extend(
                ['-m', 'qpois']
            )

    def prepare_submission(self, *args, **kwargs):
        self.params = []
        self.run_names = []
        bdg_dir = self.args['macs_dir']

        ext = r'_treat_pileup.bdg(\.gz)?'
        flist = [t for t in os.listdir(bdg_dir) if re.search(ext, t)]

        for trt_fn in flist:
            base = re.sub(r'_treat_pileup.*', '', os.path.split(trt_fn)[-1])
            ctr_fn = trt_fn.replace('_treat_pileup', '_control_lambda')
            gz = bool(re.search(r'.gz$', trt_fn, flags=re.IGNORECASE))
            if os.path.exists(ctr_fn):
                self.params.append(
                    [
                        base,
                        os.path.abspath(os.path.join(bdg_dir, trt_fn)),
                        os.path.abspath(os.path.join(bdg_dir, ctr_fn)),
                        "gz" if gz else "bdg"
                     ]
                )
                self.run_names.append(base)
            else:
                self.logger.info(
                    "Skipping treatment file %s as there is no matching control file %s",
                    trt_fn,
                    ctr_fn
                )


class MACS2BdgCmpApocrita(MACS2BdgCmpSgeRequirements, MACS2BdgCmpBase):
    core_cmd = """
    cp $TARGET $TMPDIR
    TARGET="$TMPDIR/$(basename $TARGET)"
    cp $CONTROL $TMPDIR
    CONTROL="$TMPDIR/$(basename $CONTROL)"
    fi
    """ + MACS2BdgCmpBase.core_cmd


class MACS2BdgCmpBash(jobs.BashArrayJobMixin, MACS2BdgCmpBase):
    pass


def bdgcmp_run(run_type):
    import argparse
    import sys

    if run_type == 'bash':
        cls = MACS2BdgCmpBash
    elif run_type == 'apocrita':
        cls = MACS2BdgCmpApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("--macs_dir", help="Path to bdg files", default='./')

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, use the
        args.out_dir = os.path.abspath(args.macs_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    obj = cls(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()


def callpeaks_multilane_run(run_type):
    import argparse
    import sys

    if run_type == 'bash':
        cls = MACS2CallPeaksBash
    elif run_type == 'apocrita':
        cls = MACS2CallPeaksApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("-p", "--threads", help="Number of threads", default='1')

    parser.add_argument("-c", "--config_file", help="Path to config. file", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

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

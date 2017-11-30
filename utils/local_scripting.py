import sge
import subprocess


class LocalScriptConversionMixin(object):
    """
    Converts the object from a SGE job into a local job
    """
    def create_submission_script(self):
        super(LocalScriptConversionMixin, self).create_submission_script()
        # remove header and add enclosing for loop
        self.sh[0] = "#!/bin/bash\n" + \
        "for SGE_TASK_ID in $(seq {n}); do".format(n=len(self.params))
        self.sh.append('done')

    def submit(self):
        subprocess.call(['bash', self.script_fn])


class SalmonIlluminaSEJob(LocalScriptConversionMixin, sge.SalmonIlluminaSESgeJob):
    pass

class SalmonIlluminaPEJob(LocalScriptConversionMixin, sge.SalmonIlluminaPESgeJob):
    pass

class SalmonIlluminaMultiLanePEJob(LocalScriptConversionMixin, sge.SalmonIlluminaMultiLanePESgeJob):
    pass

class CufflinksJob(LocalScriptConversionMixin, sge.CufflinksSgeJob):
    pass

class TrimgalorePEJob(LocalScriptConversionMixin, sge.TrimgalorePESgeJob):
    pass

class TrimgaloreSEJob(LocalScriptConversionMixin, sge.TrimgaloreSESgeJob):
    pass
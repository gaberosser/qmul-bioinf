try:
    from rpy2 import robjects, rinterface
    from rpy2.robjects import FloatVector, Formula
    import rpy2.robjects.packages as rpackages
    RFUNCTIONS_PRESENT = True
except ImportError:
    RFUNCTIONS_PRESENT = False

try:
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    RPANDAS_PRESENT = True
except Exception:
    RPANDAS_PRESENT = False


class RFunctionDeferred(object):
    """
    Used for R functions that also require library imports.
    Rather than importing all these libraries when the Python module is imported, we defer it until the function is
    actually required
    """
    def __init__(self, func, imports=None, redirect_stdout=False, redirect_stderr=False):
        self.imports = imports or []
        self.func = func
        self.ready = False
        self.stdout = self.stderr = None
        self._redirect_stderr = redirect_stderr
        self._redirect_stdout = redirect_stdout
        self.reset_required = redirect_stderr or redirect_stdout

    def import_packages(self):
        for im in self.imports:
            if not rpackages.isinstalled(im):
                utils = rpackages.importr('utils')
                utils.chooseCRANmirror(ind=1)
                ## FIXME: this returns rpy2.rinterface.NULL if package is not found
                utils.install_packages(im)
            rpackages.importr(im)
        self.ready = True

    def redirect_stdout(self):
        self.stdout = []

        def capture_stdout(x):
            self.stdout.append(x)

        rinterface.set_writeconsole_regular(capture_stdout)

    def redirect_stderr(self):
        self.stderr = []

        def capture_stderr(x):
            self.stderr.append(x)

        rinterface.set_writeconsole_warnerror(capture_stderr)

    def reset_redirects(self):
        rinterface.set_writeconsole_warnerror(rinterface.consolePrint)
        rinterface.set_writeconsole_regular(rinterface.consolePrint)

    def __call__(self, *args, **kwargs):
        if not self.ready:
            self.import_packages()

        if self._redirect_stdout:
            self.redirect_stdout()

        if self._redirect_stderr:
            self.redirect_stderr()

        try:
            return self.func(*args, **kwargs)
        finally:
            if self.reset_required:
                self.reset_redirects()

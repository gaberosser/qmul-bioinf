try:
    from rpy2 import robjects
    from rpy2.robjects import FloatVector, Formula
    import rpy2.robjects.packages as rpackages
except ImportError:
    RFUNCTIONS_PRESENT = False


class RFunctionDeferred(object):
    """
    Used for R functions that also require library imports.
    Rather than importing all these libraries when the Python module is imported, we defer it until the function is
    actually required
    """
    def __init__(self, func, imports=None):
        self.imports = imports or []
        self.func = func
        self.ready = False

    def import_packages(self):
        for im in self.imports:
            if not rpackages.isinstalled(im):
                utils = rpackages.importr('utils')
                utils.chooseCRANmirror(ind=1)
                utils.install_packages(im)
            rpackages.importr(im)
        self.ready = True

    def __call__(self, *args, **kwargs):
        if not self.ready:
            self.import_packages()
        return self.func(*args, **kwargs)
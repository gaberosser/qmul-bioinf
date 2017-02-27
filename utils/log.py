import logging
import os


LOG_FMT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


def get_console_logger(name):
    """
    Get the logger with the supplied name, resetting the handlers to include only a single console logger.
    :return: logging.Logger object
    """
    logger = logging.getLogger(name)

    # reset handlers
    logger.handlers = []
    sh = logging.StreamHandler()
    fmt = logging.Formatter(LOG_FMT)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    logger.setLevel(logging.INFO)

    return logger


def get_file_logger(name, filestem=None):
    if filestem is None:
        filestem = 'output.log'
    i = 0
    fn = "%s.%d" % (filestem, i)
    while os.path.exists(fn):
        i += 1
        fn = "%s.%d" % (filestem, i)
    logger = logging.getLogger(name)

    # reset handlers
    logger.handlers = []
    fh = logging.FileHandler(fn)
    fmt = logging.Formatter(LOG_FMT)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)

    return logger
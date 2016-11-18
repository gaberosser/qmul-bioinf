import logging


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
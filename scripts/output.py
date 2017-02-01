import os
from settings import OUTPUT_DIR
from log import get_console_logger
logger = get_console_logger(__name__)


def unique_output_dir(title, reuse_empty=False):

    if not os.path.exists(OUTPUT_DIR):
        logger.info("Creating global output directory %s", OUTPUT_DIR)
        os.makedirs(OUTPUT_DIR)
    i = 1
    reused = False
    outdir = os.path.join(OUTPUT_DIR, "%s.0" % title)
    while os.path.exists(outdir):
        if reuse_empty and len(os.listdir(outdir)) == 0:
            # directory is empty, so let's reuse it
            logger.info("Reusing empty output directory %s", outdir)
            reused = True
            break
        outdir = os.path.join(OUTPUT_DIR, "%s.%d" % (title, i))
        i += 1

    if not reused:
        logger.info("Creating output directory %s", outdir)
        os.makedirs(outdir)

    return outdir

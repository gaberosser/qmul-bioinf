import os, sys
import datetime
from settings import OUTPUT_DIR
from log import get_console_logger
logger = get_console_logger(__name__)


def unique_output_dir(title=None, root_output_dir=OUTPUT_DIR, reuse_empty=True, add_date=True):
    if title is None:
        # title is set to the name of the calling script
        title = os.path.splitext(os.path.basename(sys.argv[0]))[0]

    if add_date:
        title = "%s.%s" % (title, datetime.date.today().strftime('%Y-%m-%d'))

    if not os.path.exists(root_output_dir):
        logger.info("Creating global output directory %s", root_output_dir)
        os.makedirs(root_output_dir)
    i = 1
    reused = False
    outdir = os.path.join(root_output_dir, "%s.0" % title)
    while os.path.exists(outdir):
        if reuse_empty and len(os.listdir(outdir)) == 0:
            # directory is empty, so let's reuse it
            logger.info("Reusing empty output directory %s", outdir)
            reused = True
            break
        outdir = os.path.join(root_output_dir, "%s.%d" % (title, i))
        i += 1

    if not reused:
        logger.info("Creating output directory %s", outdir)
        os.makedirs(outdir)

    return outdir

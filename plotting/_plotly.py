from plotly import graph_objs as go
import plotly

from settings import PLOTLY_API_CREDENTIALS
from utils import log
logger = log.get_console_logger()


MATPLOTLIB_TO_PLOTLY_MARKERS = {
    'o': 'circle',
    's': 'square'
}


# check plotly credentials file and create if required
try:
    cred = plotly.tools.get_credentials_file()
    if cred['username'] == '':
        logger.info(
            "Plotly credentials are missing, generating these now for user %s",
            PLOTLY_API_CREDENTIALS['username']
        )
        plotly.tools.set_credentials_file(**PLOTLY_API_CREDENTIALS)
    else:
        logger.info(
            "Using existing Plotly credentials for user %s. Run set_plotly_credentials() to regenerate.",
            PLOTLY_API_CREDENTIALS['username']
        )
    ## TODO: is there any way to test this config beforehand?
    _PLOTLY_PRESENT = True
except Exception:
    logger.exception("Something went wrong importing plotly or setting the credentials.")
    _PLOTLY_PRESENT = False


def set_plotly_credentials():
    plotly.tools.set_credentials_file(**PLOTLY_API_CREDENTIALS)

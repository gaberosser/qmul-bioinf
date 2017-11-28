import pandas as pd
import requests
from StringIO import StringIO


SRA_CGI = 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi'


def get_runinfo_table(srp_id):
    """
    Get the runinfo table for the supplied project code (SRPxxxxxx)
    :param srp_id:
    :return:
    """
    params = {
        'save': 'efetch',
        'db': 'sra',
        'rettype': 'runinfo',
        'term': srp_id
    }
    r = requests.get(SRA_CGI, params=params)
    if r.status_code != 200:
        raise ValueError("Failed to get data")
    s = StringIO(r.content)
    return pd.read_csv(s, header=0, index_col=0)
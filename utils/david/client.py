from settings import DAVID_WEB_SERVICES_CONFIG
from suds.client import Client
from utils import log


logger = log.get_console_logger('DAVID_web_services')

WSDL_URL = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
SOAP_ENDPOINT = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'


class WSDLApi(object):

    def __init__(self, url=WSDL_URL, user=DAVID_WEB_SERVICES_CONFIG['email']):
        self.user = user
        self.url = url
        self.client = None
        self.connect(endpoint=SOAP_ENDPOINT)

    def connect(self, endpoint):
        self.client = Client(self.url)
        self.client.wsdl.services[0].setlocation(
            endpoint
        )
        # authenticate user email
        self.client.service.authenticate(self.user)

    def introspection(self):
        """
        Print the service (introspection)
        :return:
        """
        print self.client


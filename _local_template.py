# This is a template for a local config file that needs to be named _local.py
# Make a copy and modify with your own configuration variables
# Many of these are only needed for specific use cases. It is fine to leave them out - errors won't appear until
# they are needed (which may be never!)

import os

# this data dir contains any external files NOT stored in the Git repo itself
DATA_DIR = '/path/to/data/'

# location for intermediate results
INTERMEDIATE_DIR = '//path/to/intermediate_results'

# these data directories are really just for convenience, but you might wish to change them if for some reason
# different modalities are stored in different locations
RNASEQ_DIR = os.path.join(DATA_DIR, 'rnaseq')
CHIPSEQ_DIR = os.path.join(DATA_DIR, 'chipseq')
METHYLATION_ARRAY_DIR = os.path.join(DATA_DIR, 'methylation')

# this data dir contains 'local' data - typically reference genome files
# no reason why it can't be in the same DATA_DIR if you like - I have an SSD so it's quicker to put local files there
# LOCAL_DATA_DIR = DATA_DIR
LOCAL_DATA_DIR = '/path/to/local_data'

# all outputs go here in their own subdirectories
OUTPUT_DIR = '/path/to/python_outputs/'

# This is the path to the 'final' hGIC project folder
# It contains data and plots that others use, so shouldn't be changed without warning
# Generally, it isn't advised to write to this directory programmatically. However, it can be useful to read in
# pre-existing results directly.
HGIC_LOCAL_DIR = '~/Dropbox/research/qmul/hGIC_project/' # NB I don't think ~ works here so replace with the actual path

# credentials for the Heidelberg classifier (https://www.molecularneuropathology.org/mnp)
# I don't mind if others in the group want to use this account
HEIDELBERG_CLASSIFIER_CONFIG = {
    'email': "g.rosser@qmul.ac.uk",
    'password': "IStillWantToClassify"
}

PLOTLY_API_CREDENTIALS = {
    'username': '',
    'api_key': ''
}
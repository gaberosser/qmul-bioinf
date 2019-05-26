import os
from settings import LOCAL_DATA_DIR


NH_ID_TO_PATIENT_ID_MAP = {
    'NH15-1661': '017',
    'NH15-1877': '018',
    'NH15-2101': '019',
    'NH16-270': '026',
    'NH16-616': '030',
    'NH16-677': '031',
    'NH16-1574': '044',
    'NH16-1976': '049',
    'NH16-2063': '050',
    'NH16-2214': '052',
    'NH16-2255': '054',
    'NH16-2806': '061',
}

PATIENT_ID_TO_NH_ID_MAP = dict([t[::1] for t in NH_ID_TO_PATIENT_ID_MAP.items()])

REFERENCE_GENOMES_BASEDIR = os.path.join(LOCAL_DATA_DIR, 'reference_genomes')

REFERENCE_GENOME_DIRECTORIES = {
    9606: {
        'default': os.path.join(REFERENCE_GENOMES_BASEDIR, 'human', 'ensembl', 'GRCh38.p10.release90'),
        'GRCh37': os.path.join(REFERENCE_GENOMES_BASEDIR, 'human', 'ensembl', 'GRCh37'),
        'GRCh38': os.path.join(REFERENCE_GENOMES_BASEDIR, 'human', 'ensembl', 'GRCh38.p10.release90'),
    },
    10090: {
        'default': os.path.join(REFERENCE_GENOMES_BASEDIR, 'mouse', 'ensembl', 'GRCm38.p5.r88'),
        'GRCm38': os.path.join(REFERENCE_GENOMES_BASEDIR, 'mouse', 'ensembl', 'GRCm38.p5.r88'),
    },
}

REFERENCE_GENOME_GTFS = {
    9606: {
        'default': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['default'], 'gtf', 'Homo_sapiens.GRCh38.90'),
        'GRCh37': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['GRCh37'], 'gtf', 'Homo_sapiens.GRCh37.87'),
        'GRCh38': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['GRCh38'], 'gtf', 'Homo_sapiens.GRCh38.90'),
    },
    10090: {
        'default': os.path.join(REFERENCE_GENOME_DIRECTORIES[10090]['default'], 'gtf', 'Mus_musculus.GRCm38.88'),
        'GRCm38': os.path.join(REFERENCE_GENOME_DIRECTORIES[10090]['GRCm38'], 'gtf', 'Mus_musculus.GRCm38.88'),
    },
}

REFERENCE_GENOME_FA = {
    9606: {
        'default': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['default'], 'fa', 'Homo_sapiens.GRCh38.dna.primary_assembly'),
        'GRCh37': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['GRCh37'], 'fa', 'Homo_sapiens.GRCh37.dna.primary_assembly'),
        'GRCh38': os.path.join(REFERENCE_GENOME_DIRECTORIES[9606]['GRCh38'], 'fa', 'Homo_sapiens.GRCh38.dna.primary_assembly'),
    },
    10090: {
        'default': os.path.join(REFERENCE_GENOME_DIRECTORIES[10090]['default'], 'fa', 'Mus_musculus.GRCm38.dna.primary_assembly'),
        'GRCm38': os.path.join(REFERENCE_GENOME_DIRECTORIES[10090]['GRCm38'], 'fa', 'Mus_musculus.GRCm38.dna.primary_assembly'),
    },
}
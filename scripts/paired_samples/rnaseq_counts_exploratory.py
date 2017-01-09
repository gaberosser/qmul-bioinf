from load_data import rnaseq_data
from scripts.comparison_rnaseq_microarray import consts


if __name__ == '__main__':
    res = rnaseq_data.paired_samples_gene_counts(annotate_by='Approved Symbol')
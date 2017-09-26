import sys
import gzip
import itertools
import io
import collections
import os
import re


def fastq_reader(fstream):
    """
    :param fstream: Iterable filestream
    :return: Generator of dict objects
    """
    n = 1
    this_block = {}
    for l in fstream:
        if n % 4 == 0:
            yield this_block
            this_block = {}
        elif n % 4 == 1:
            this_block['id'] = l.split(' ')[0].strip('@')
        elif n % 4 == 2:
            this_block['seq'] = l.strip()
        n += 1


def fastq_paired_reader(fstream1, fstream2):
    """
    :param fstream1, fstream2: Iterable filestreams
    :return: Generator of dict objects
    """
    n = 1
    this_block = {}
    for l1, l2 in itertools.izip(fstream1, fstream2):  # generator version of zip
        if n % 4 == 0:
            yield this_block
            this_block = {}
        elif n % 4 == 1:
            this_id1 = l1.split(' ')[0].strip('@')
            this_id2 = l2.split(' ')[0].strip('@')
            if this_id1 != this_id2:
                raise ValueError("fastq files are not synced: IDs unpaired")
            this_block['id'] = this_id1
        elif n % 4 == 2:
            this_block['seq_1'] = l1.strip()
            this_block['seq_2'] = l2.strip()
        n += 1


def main(fq1, fq2, nseq=100000, seq_length=50, skip_first_N=1, verbose=True):
    is_gz = (fq1[-3:].lower() == '.gz') & (fq2[-3:].lower() == '.gz')

    if is_gz:
        f1 = io.TextIOWrapper(io.BufferedReader(gzip.open(fq1)))
        f2 = io.TextIOWrapper(io.BufferedReader(gzip.open(fq2)))
    else:
        f1 = open(fq1, 'rb')
        f2 = open(fq2, 'rb')

    num_read = 0
    n = 0
    try:
        g = fastq_paired_reader(f1, f2)
        # get candidate sequences
        candidates = {}
        candidate_count = collections.Counter()
        candidate_by_seq = {}
        sequence_list = []
        seq_set = set()
        i = 0
        skipped = 0
        dupes = 0

        while i < nseq:
            the_block = g.next()
            the_id = the_block.pop('id')
            seq_1 = the_block['seq_1'][:seq_length]
            seq_2 = the_block['seq_2'][:seq_length]

            the_seq = seq_1[skip_first_N:] + seq_2[skip_first_N:]
            the_N_count = len(the_seq) - len(the_seq.replace('N', ''))

            if the_N_count > 0:
                skipped += 1
                continue

            candidates[the_id] = {
                'seq_1': seq_1,
                'seq_2': seq_2,
                'hash_seq': the_seq,
            }
            i += 1

            if the_seq in seq_set:
                candidate_count[candidate_by_seq[the_seq]] += 1
                dupes += 1
            else:
                seq_set.add(the_seq)
                candidate_by_seq[the_seq] = the_id
                candidate_count[the_id] += 1

            sequence_list.append((the_seq, the_id))

        if verbose:
            print "Finished recording %d candidate pairs. " \
                  "Skipped %d as they contained undetermined bases. " \
                  "Identified %d duplicates" % (
                nseq,
                skipped,
                dupes
            )
        skipped = 0

        # run through the remainder
        for i, the_block in enumerate(g):
            if verbose and i % 500000 == 0 and i != 0:
                print "%d lines read" % i
            seq_1 = the_block.pop('seq_1')[:seq_length]
            seq_2 = the_block.pop('seq_2')[:seq_length]

            the_seq = seq_1[skip_first_N:] + seq_2[skip_first_N:]

            if 'N' in the_seq:
                skipped += 1
                continue

            if the_seq in seq_set:
                candidate_count[candidate_by_seq[the_seq]] += 1

    finally:
        f1.close()
        f2.close()

    return candidate_count, candidates


def quantity_after_dedupe(candidate_count):
    n_read = float(sum(candidate_count.values()))
    n_cand = float(len(candidate_count))
    return n_cand / n_read * 100.


def duplication_values(candidate_count, bins=None):
    if bins is None:
        bins = range(1, 11) + [100, 1000, 10000, 1e12]  # upper limit is essentially inf

    n_read = float(sum(candidate_count.values()))
    n_cand = float(len(candidate_count))

    binned_data = []

    for i in range(len(bins) - 1):
        x0 = bins[i]
        x1 = bins[i + 1]
        binned_data.append([t for t in candidate_count.values() if x0 <= t < x1])

    # proportion of total reads
    pct_total = [sum(t) / n_read * 100. for t in binned_data]

    # proportion of deduped reads
    pct_dedupe = [len(t) / n_cand * 100. for t in binned_data]

    return pct_total, pct_dedupe


def duplication_plot(candidate_count):
    from matplotlib import pyplot as plt

    lbl = [str(i) for i in range(1, 10)] + ['<100', '<1k', '<10k', '>=10k']
    pct_total, pct_dedupe = duplication_values(candidate_count)

    fig = plt.figure(figsize=(7.2, 3.8))
    ax = fig.add_subplot(111)
    ax.plot(pct_total, label='Percent of total reads')
    ax.plot(pct_dedupe, label='Percent of deduped reads')
    ax.set_xticks(range(len(pct_total)))
    ax.set_xticklabels(lbl)
    ax.legend(loc='upper right')
    fig.tight_layout()

    return ax


if __name__ == "__main__":
    """
    Usage: paired_end_duplicates.py read_1.fastq read2.fastq
    Idea: reproduce FastQC's method for checking for duplication, BUT taking paired ends into account
    The quantity reported represents the % of reads that would remain IF all duplicates (both ends matching) were
    removed.
    """
    # TODO: allow varying seq length and number of sequences to follow

    if len(sys.argv) != 3:
        print "Usage: paired_end_duplicates.py read_1.fastq read2.fastq"
        sys.exit(1)

    # fq1 = '/home/gabriel/data/rnaseq/wtchg_p160704/161219_K00198_0151_BHGYHTBBXX/WTCHG_338493_201101_1.fastq.gz'
    # fq2 = '/home/gabriel/data/rnaseq/wtchg_p160704/161219_K00198_0151_BHGYHTBBXX/WTCHG_338493_201101_2.fastq.gz'

    fq1 = sys.argv[1]
    fq2 = sys.argv[2]

    # get root name
    root = os.path.split(fq1)[1].replace('fastq', '').replace('gz', '').strip('.')
    root = re.sub(r'_1$', '', root)

    candidate_count, candidates = main(fq1, fq2, verbose=False)

    # amount remaining after deduplication
    print "%s, %.2f" % (root, quantity_after_dedupe(candidate_count))

import os
import collections
from settings import DATA_DIR_NON_GIT
import pysam


from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


if __name__ == "__main__":
    indir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163/mouse/bismark/GC-CV-7163-1_S1')
    bam_fn = os.path.join(indir, 'GC-CV-7163-1_S1_pe.sorted.bam')
    s = pysam.AlignmentFile(bam_fn, 'rb')
    chroms = [str(t) for t in range(1, 20)]
    regions = {}
    # for c in chroms:
    c = '1'
    regions[c] = []
    ## TODO: check this actually fetches from the specified chromosome
    it = s.fetch(c)
    for seg in it:
        if seg.flag == 99 or seg.flag == 163:
            # 99: properly mapped in pair, first in pair, forward strand
            # 163: properly mapped in pair, second in pair, forward strand
            if seg.tlen > 0:
                regions[c].append((seg.mpos, seg.tlen))

    # get the sequences from the reference

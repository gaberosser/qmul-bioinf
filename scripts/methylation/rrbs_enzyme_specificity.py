import os
import collections
from settings import DATA_DIR_NON_GIT
import pysam


if __name__ == "__main__":
    indir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163/mouse/bismark/GC-CV-7163-1_S1')
    bam_fn = os.path.join(indir, 'GC-CV-7163-1_S1_pe.sorted.bam')
    s = pysam.AlignmentFile(bam_fn, 'rb')
    chroms = [str(t) for t in range(1, 20)]
    # for c in chroms:
    c = '1'
    it = s.fetch(c)
    for seg in it:
        if seg.flag == 83:
            pass
        elif seg.flag == 99:
            pass
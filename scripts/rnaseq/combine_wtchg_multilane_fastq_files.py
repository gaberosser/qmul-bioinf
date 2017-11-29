import os
import subprocess
import argparse
import sys
import re


if __name__ == "__main__":
    """
    Designed to concatenate fastq files correspondoing to multiple lanes.
    This is useful for Salmon or similar alignment-free approaches, for which it isn't obvious how to join results
    AFTER running them.

    NB I am assuming the naming convention here is the WTCHG standard
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--root_dir', '-r', default='./')

    args = parser.parse_args()
    root_dir = args.root_dir
    print "Starting in root directory %s" % root_dir

    # now we carry out a search for compatible concatenatable files
    dlist = [os.path.join(root_dir, t) for t in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, t))]
    print "Candidate directories: %s" % ', '.join(dlist)

    to_join = {}

    for d in dlist:
        flist = [t for t in os.listdir(d) if '.fastq' in t]
        for fn in flist:
            base = re.sub(r'\.fastq(\.gz)?$', '', fn)
            # get read number
            read_num = int(base[-1])
            # get ID
            read_id = re.sub(r'_[12]$', '', base)
            read_id = re.sub('^WTCHG_[0-9]+_', '', read_id)
            to_join.setdefault(read_id, {}).setdefault(d, {})[read_num] = os.path.join(d, fn)

    n = None
    print "Identified %d fastq groups that can be combined: %s" % (len(to_join), ', '.join(to_join.keys()))
    for read_id, d in to_join.items():
        # get all valid pairs / singles
        if len(d.values()[0]) == 1:
            typ = 'single-end'
        elif len(d.values()[0]) == 2:
            typ = 'paired-end'
        else:
            raise ValueError("Incorrect number of corresponding reads?")

        print "Read group %s. Type: %s" % (read_id, typ)

        if n is None:
            n = len(d)
        else:
            if len(d) != n:
                raise AttributeError("Previous read group had %d matching directories; this one has %d." % (n, len(d)))

        this_join = {}
        for the_dir, the_dict in d.items():
            for read_num, the_file in the_dict.items():
                this_join.setdefault(read_num, []).append(the_file)

        for read_num, arr in this_join.items():
            if 'gz' in arr[0].lower():
                ext = 'fastq.gz'
            else:
                ext = 'fastq'
            out_file = os.path.join(root_dir, '%s_%d.%s' % (read_id, read_num, ext))
            print "Concatenating files (%s) -> %s" % (', '.join(arr), out_file)

            cmd = ["cat"] + arr
            # print ' '.join(cmd)
            with open(out_file, 'wb') as f:
                subprocess.call(cmd, stdout=f)


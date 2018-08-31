#!/usr/bin/env python2
import sys
import os
import shutil
import tempfile
import subprocess
import argparse

__version__ = "0.6.3"

"""
Modified from original code: https://github.com/rvalieris/parallel-fastq-dump
Dependencies: SRA toolkit (fastq-dump, sra-stat)
"""


def pfd(srr_id, extra_args=None, outdir=os.curdir, threads=1):
    if extra_args is None:
        extra_args = []
    tmp_dir = tempfile.mkdtemp(dir=outdir)
    sys.stderr.write("Temporary directory: {}\n".format(tmp_dir))

    try:

        n_spots = get_spot_count(srr_id)
        sys.stderr.write("{} spots: {}\n".format(srr_id,n_spots))
        blocks = split_blocks(n_spots, threads)
        sys.stderr.write("blocks: {}\n".format(blocks))

        ps = []
        for i in range(0, threads):
            d = os.path.join(tmp_dir, str(i))
            os.mkdir(d)
            p = subprocess.Popen(["fastq-dump", "-N", str(blocks[i][0]), "-X", str(blocks[i][1]), "-O", d]+extra_args+[srr_id])
            ps.append(p)

        wfd = {}
        for i in range(0, threads):
            exit_code = ps[i].wait()
            if exit_code != 0:
                sys.stderr.write("fastq-dump error! exit code: {}\n".format(exit_code))
                sys.exit(1)
            else:
                sys.stderr.write("Completed block {blck}/{nblck}.\n".format(blck=i+1, nblck=threads))

            # for each file in the temporary subdir (e.g. 2 if PE sequencing), copy/append to the final output dir.
            # we process this in thread order to ensure ordering is maintained
            tmp_path = os.path.join(tmp_dir, str(i))
            for fo in os.listdir(tmp_path):
                if fo not in wfd:
                    # open a new destination file object
                    wfd[fo] = open(os.path.join(outdir, fo), "wb")
                with open(os.path.join(tmp_path,fo), "rb") as fd:
                    # copy this portion to the final output object
                    shutil.copyfileobj(fd, wfd[fo])

        sys.stderr.write("Completed all blocks. Output files are located at: {outlist}.\n".format(
            outlist='\n'.join([t.name for t in wfd.values()])
        ))

    finally:

        sys.stderr.write("Deleting temporary directory {}\n".format(tmp_dir))
        shutil.rmtree(tmp_dir)


def split_blocks(total, n_pieces):
    avg = int(total / n_pieces)
    out = []
    last = 1
    for i in range(0,n_pieces):
        out.append([last,last + avg-1])
        last += avg
        if i == n_pieces-1: out[i][1] += total % n_pieces
    return out


def get_spot_count(sra_id):
    p = subprocess.Popen(["sra-stat", "--meta", "--quick", sra_id], stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    txt = stdout.decode().rstrip().split("\n")
    total = 0
    for l in txt:
        total += int(l.split("|")[2].split(":")[0])
    return total


def partition(f, l):
    r = ([],[])
    for i in l:
        if f(i):
            r[0].append(i)
        else:
            r[1].append(i)
    return r


def main():
    parser = argparse.ArgumentParser(description="parallel fastq-dump wrapper, extra args will be passed through")
    parser.add_argument("-s","--sra-id", help="SRA id", action="append")
    parser.add_argument("-t","--threads", help="number of threads", default=1, type=int)
    parser.add_argument("-O","--outdir", help="output directory", default=".")
    parser.add_argument("--tmpdir", help="temporary directory", default=None)
    parser.add_argument("-V", "--version", help="shows version", action="store_true")
    args, extra = parser.parse_known_args()

    if args.version:
        print("parallel-fastq-dump : {}".format(__version__))
        subprocess.Popen(["fastq-dump", "-V"]).wait()
        sys.exit(0)

    elif args.sra_id:
        extra_srrs, extra_args = partition(
            lambda s: "SRR" in s.upper() or s.lower().endswith('.sra'),
            extra)
        args.sra_id.extend(extra_srrs)
        sys.stderr.write("SRR ids: {}\n".format(args.sra_id))
        sys.stderr.write("extra args: {}\n".format(extra_args))

        if args.outdir:
            if not os.path.isdir(args.outdir):
                os.mkdir(args.outdir)

        for si in args.sra_id:
            pfd(si, extra_args=extra_args, outdir=args.outdir, threads=args.threads)

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()


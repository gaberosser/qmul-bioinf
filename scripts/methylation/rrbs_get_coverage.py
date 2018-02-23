import pysam


def get_one_coverage(bam_fn, region, region_pad=500):
    chrom, start, end = region
    start = max(start - region_pad, 0)
    end += region_pad  # no problem if we run off the end

    # get the coverage over this region
    region_depth = pysam.depth(bam_fn, "-a", "-r", "%s:%d-%d" % (chrom, start, end))
    return region_depth
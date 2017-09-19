import sys
import zipfile
import re

THE_REGEX = re.compile(r'^#Total Deduplicated Percentage')


def main(fn):
    with zipfile.ZipFile(fn, 'r') as z:
        the_file = [t.filename for t in z.filelist if 'fastqc_data.txt' in t.filename][0]  # should only be 1
        f = z.open(the_file)
        for line in f:
            if THE_REGEX.match(line):
                the_val = re.search(r'.*\t(?P<x>[0-9\.]*)$', line)
                return float(the_val.group('x'))




if __name__ == "__main__":
    fn = sys.argv[1]
    n = main(fn)
    print "%s, %.2f" % (fn, n)
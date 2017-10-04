import gzip
import os
import re
import csv


if __name__ == "__main__":
    LIMIT = int(1e12)
    SOURCES = {'ensembl', 'havana', 'ensembl_havana'}
    distance = 2000 # distance from TSS to include
    fn = '/home/gabriel/Documents/qmul_data/Homo_sapiens.GRCh38.87.gtf.gz'
    fn_out = '/home/gabriel/Documents/qmul_data/tss_pm_%d.bed' % distance

    res = []
    bed = []

    with gzip.open(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        line_no = 1
        for row in c:
            if len(row) == 1:
                continue
            if line_no > LIMIT:
                break
            chr = row[0]
            src = row[1]
            typ = row[2]
            start = int(row[3])
            stop = int(row[4])
            strand = row[6]
            blk = row[8].split('; ')

            if strand == '+':
                start = int(row[3])
                stop = int(row[4])
            elif strand == '-':
                start = int(row[4])
                stop = int(row[3])

            attr = dict([t.split(' ') for t in blk])

            if src in SOURCES \
                    and typ == 'exon' \
                    and attr.get('exon_number') == '"1"'\
                    and attr.get('gene_biotype') == '"protein_coding"':
                attr['chr'] = chr
                attr['source'] = src
                attr['type'] = typ
                attr['start'] = start
                attr['stop'] = stop
                attr['strand'] = strand
                res.append(attr)

                # BED format is 0-based
                # GTF format is 1-based
                # therefore subtract one from all base counts for the BED format
                tss_minus = start - distance - 1
                tss_plus = start + distance - 1

                bed.append([chr, tss_minus, tss_plus, attr['transcript_name'].strip('"')])

            line_no += 1

    with open(fn_out, 'wb') as f:
        c = csv.writer(f, delimiter='\t')
        c.writerows(bed)

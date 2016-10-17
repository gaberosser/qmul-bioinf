# Change the naming convention of chromosomes from '1', '2', ..., 'MT' to 'chr1', 'chr2', ..., 'chrM'
# Also drop any unmapped contigs
# Important to maintain the first few lines beginning with #

sed -nr '/^#/ p; s/^MT/chrM/p; s/(^[XY0-9]+)/chr\1/p' Homo_sapiens.GRCh37.85.gtf > Homo_sapiens.GRCh37.85.chr.gtf
# test it
sed -n '/chr/ !p' Homo_sapiens.GRCh37.85.chr.gtf

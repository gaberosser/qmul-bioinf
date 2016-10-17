source("http://www.bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
# biocLite("GenomicFeatures")

library("Rsamtools")

filenames <- c(
  "/home/gabriel/data/Xinyu_Zhang_NEB_mRNASeq_GC-XZ-2499_270115-19825815/XZ-1-21077637/alignments/XZ-1.alignments.bam",
  "/home/gabriel/data/Xinyu_Zhang_NEB_mRNASeq_GC-XZ-2499_270115-19825815/XZ-2-21099482/alignments/XZ-2.alignments.bam"
)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)

library("GenomicFeatures")
gtffile <- "/home/gabriel/data/Human_Genome/gtf/Homo_sapiens.GRCh38.85.gtf.gz"
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))

# need to rename the chromosomes in the txdb
# also drop the unplaced contigs, since these are not represented in the RNA-Seq data
library("GenomeInfoDb")
keepSeqlevels(txdb, c(paste0("",1:22), c("X", "Y", "MT")))
newnames <- c(paste0("chr", 1:22), c("chrX", "chrY", "chrM"))
names(newnames) <- seqlevels(txdb)
renameSeqlevels(txdb, newnames)

# perform counting
library("GenomicAlignments")
# library("BiocParallel")
# 2 cores in parallel
# register(Serial(workers = 2))

(ebg <- exonsBy(txdb, by="gene"))
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode='Union',
                        singleEnd = FALSE,
                        ignore.strand = FALSE,
                        fragments = TRUE)


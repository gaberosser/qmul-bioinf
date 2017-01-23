source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)

source('io/microarray.R')
source('_settings.R')

in.dirs <- c(
  file.path(
  data.dir.raid, 
  'rnaseq',
  'wtchg_p160704',
  '161219_K00198_0151_BHGYHTBBXX'
  ),
  file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    '161222_K00198_0152_AHGYG3BBXX'
  )
)

samples <- c(
  'GBM018',
  'GBM019',
  'GBM024',
  'GBM026',
  'GBM031',
  'DURA018N2_NSC',
  'DURA019N8C_NSC',
  'DURA024N28_NSC',
  'DURA026N31D_NSC',
  'DURA031N44B_NSC'
)

units <- 'fpkm'

dat <- NULL
nreads <- NULL
first_run <- T
for (d in in.dirs) {
  meta.fn <- file.path(d, 'sources.csv')
  fn = file.path(d, 'featureCounts', 'counts.txt')
  meta <- read.csv(meta.fn, sep=',', header=T, row.names=1)
  meta <- meta[meta$sample %in% samples,]
  
  codes <- paste0(rownames(meta), '.bam')
  x <- read.csv(fn, sep='\t', skip=1, header=T, row.names=1)
  l <- as.data.frame(x$Length)
  nr <- as.data.frame(read_count=meta$read_count)
  rownames(nr) <- samples
  rownames(l) <- rownames(x)
  
  x <- x[, c(codes)]
  colnames(x) <- samples
  if (first_run) {
    dat <- x
    nreads <- nr
    first_run <- F
  } else {
    dat <- dat + x
    nreads <- nreads + nr
  }
  
}

if (units == 'fpkm') {
  dat <- dat / nreads / l * 1e6
}

filenames <- paste0('XZ-', 1:8, '.count')
in.files <- sapply(filenames, function (x) file.path(in.dir, x))

tbl <- lapply(in.files, function(x) read.csv(x, sep='\t', header = F, row.names = 1))
res.xz <- do.call(cbind, tbl)
colnames(res.xz) <- paste0('XZ', 1:8)

# annotate by symbol
annot <- read.csv(file.path(in.dir, 'annotation.csv.gz'), header=1, row.names = NULL)
annot <- annot[!duplicated(annot[,'query']),]
rownames(annot) <- annot[,'query']
annot$query <- NULL

sym <- as.character(annot[rownames(res.xz), 'symbol'])
# discard entries with no symbol
res.xz <- res.xz[!is.na(sym),]
sym <- na.omit(sym)

# aggregate
res.xz.agg <- aggregate(res.xz, by=list(symbol=sym), FUN=sum)
rownames(res.xz.agg) <- res.xz.agg$symbol
res.xz.agg$symbol <- NULL

# allen counts

in.dir <- '../data/allen_human_brain_atlas/rnaseq/'
in.file <- file.path(in.dir, 'cerebellum.counts.csv.gz')

res.he <- read.csv(in.file, row.names = 1)

# combine, keeping only matching rows
match = intersect(rownames(res.he), rownames(res.xz.agg))
res <- cbind(res.he[match,], res.xz.agg[match,])

# filter rows where the read count <2
res <- res[rowSums(res) > 1,]
res <- as.matrix(res)
storage.mode(res) <- "integer"

# metadata
meta <- data.frame(condition = c(
  as.vector(matrix("control", 9)), as.vector(matrix("mb", 8))
), row.names = colnames(res))


# DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(res), colData = meta, design=~condition)
dds <- DESeq(dds)

# don't need to specify contrasts here (same result without) but it helps to ensure the correct ordering
des = results(dds, contrast=c("condition", "mb", "control"))
des = des[order(des$pvalue),]

# find MB-specific genes
genes.wnt <- c("WIF1", "TNC", "GAD1", "DKK2", "EMX2")
genes.shh <- c("PDLIM3", "EYA1", "HHIP", "ATOH1", "SFRP1")
genes.c <- c("IMPG2", "GABRA5", "EYS", "NRL", "MAB21L2", "NPR3")
genes.d <- c("KCNA1", "EOMES", "KHDRBS2", "RBM24", "UNC5D", "OAS1", "OTX2")

genes.all <- c(genes.wnt, genes.shh, genes.c, genes.d)
genes.group <- as.factor(c(
  as.vector(matrix('WNT', length(genes.wnt))),
  as.vector(matrix('SHH', length(genes.shh))),
  as.vector(matrix('C', length(genes.c))),
  as.vector(matrix('D', length(genes.d)))
))

des.ncott <- as.data.frame(des[genes.all,])
des.ncott$subgroup <- genes.group

# repeat with only xz1 and xz2
res.red <- res[,1:11]
meta.red <- data.frame(condition = c(
  as.vector(matrix("control", 9)), as.vector(matrix("mb", 2))
), row.names = colnames(res.red))

dds.red <- DESeqDataSetFromMatrix(countData = as.matrix(res.red), colData = meta.red, design=~condition)
dds.red <- DESeq(dds.red)

des.red = results(dds.red, contrast=c("condition", "mb", "control"))
des.red = des.red[order(des.red$pvalue),]

des.red.ncott <- as.data.frame(des.red[genes.all,])
des.red.ncott$subgroup <- genes.group
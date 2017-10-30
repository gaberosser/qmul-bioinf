source("http://www.bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library(dplyr)
library(DESeq2)

# load data

# XZ counts

in.dir <- '../data/rnaseq_GSE83696/htseq-count'
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

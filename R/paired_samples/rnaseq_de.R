source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')

source('io/microarray.R')
source('_settings.R')

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

# map from Ensembl to gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens.map <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id"), mart=mart)
ens.map <- ens.map[ens.map$hgnc_symbol != '',]
ens.map <- ens.map[isUnique(ens.map$ensembl_gene_id),]
rownames(ens.map) <- ens.map$ensembl_gene_id
ens.map$ensembl_gene_id <- NULL

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
  # 'GBM024',  # remove as it is quite weird
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
  nr <- data.frame(read_count=meta$read_count)
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

dat.fpkm <- dat / rep.col(l[[1]], ncol(dat)) * 1e9
dat.fpkm <- dat.fpkm / rep.row(nreads[[1]], nrow(dat))

dat.tpm <- dat / rep.col(l[[1]], ncol(dat))
dat.tpm <- dat.tpm / rep.row(colSums(dat.tpm), nrow(dat)) * 1e6

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData = meta, design=~type + disease_subgroup)
dds <- DESeq(dds)

# test 1: GBM vs iNSC 
des = results(dds, contrast=c("type", "GBM", "iNSC"))
des = des[order(des$padj),]
# add gene symbol column
des$gene_symbol <- ens.map[rownames(des), 1]

# volcano plot
padj.threshold <- 1e-10
log2fc.threshold <- 5

# Make a basic volcano plot
with(des, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15, 15)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(des, padj < padj.threshold ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(des, abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(des, padj < padj.threshold & abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(des, padj < padj.threshold & abs(log2FoldChange) > log2fc.threshold), textxy(log2FoldChange, -log10(padj), labs=gene_symbol, cex=.8))


# test 2: GBM vs iNSC (TODO)
des = results(dds, contrast=c("disease_subgroup", "GBM", "iNSC"))
des = des[order(des$padj),]
# add gene symbol column
des$gene_symbol <- ens.map[rownames(des), 1]

# volcano plot
padj.threshold <- 1e-10
log2fc.threshold <- 5

# Make a basic volcano plot
with(des, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15, 15)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(des, padj < padj.threshold ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(des, abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(des, padj < padj.threshold & abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(des, padj < padj.threshold & abs(log2FoldChange) > log2fc.threshold), textxy(log2FoldChange, -log10(padj), labs=gene_symbol, cex=.8))

source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')
library(calibrate)
library("pheatmap")

source('io/microarray.R')
source('io/output.R')
source('_settings.R')

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

get_de_genes <- function(res, p.col="padj", p.threshold=1e-10, log2fc.threshold=5) {
  # we'll always need to avoid NAs
  na_idx = !is.na(res[[p.col]])
  subres <- res[res[[p.col]] < p.threshold & abs(res$log2FoldChange) > log2fc.threshold & na_idx,]
  return(subres)
}

volcano <- function(res, p.col="padj", p.threshold=1e-10, log2fc.threshold=5, xlim=NULL, label.field=NULL, title=NULL) {
  
  # Make a basic volcano plot
  plot(res$log2FoldChange, -log10(res[[p.col]]), pch=20, main=title, xlim=xlim)
  
  # we'll always need to avoid NAs
  na_idx = !is.na(res[[p.col]])

  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  subres <- res[res[[p.col]] < p.threshold & na_idx,]
  points(subres$log2FoldChange, -log10(subres[[p.col]]), pch=20, col="red")
  
  subres <- res[abs(res$log2FoldChange) > log2fc.threshold & na_idx,]
  points(subres$log2FoldChange, -log10(subres[[p.col]]), pch=20, col="orange")
  
  subres <- res[res[[p.col]] < p.threshold & abs(res$log2FoldChange) > log2fc.threshold & na_idx,]
  points(subres$log2FoldChange, -log10(subres[[p.col]]), pch=20, col="green", cex=1.1)
  
  # with(subset(res, padj < p.threshold ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  # with(subset(res, abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
  # with(subset(res, padj < p.threshold & abs(log2FoldChange) > log2fc.threshold), points(log2FoldChange, -log10(padj), pch=20, col="green"))
  
  if (!is.null(label.field)) {
    # Label points with the textxy function from the calibrate plot
    c <- res[na_idx & res[[p.col]] < p.threshold & abs(res$log2FoldChange) > log2fc.threshold,]
    textxy(c$log2FoldChange, -log10(c[[p.col]]), labs=c[[label.field]], cex=.8)
  }
  
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
  # 'DURA024N28_NSC',  # remove as it is unpaired
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
  rownames(nr) <- samples  # sample names
  rownames(l) <- rownames(x)  # genes
  
  x <- x[, c(codes)]
  colnames(x) <- samples
  # set sample name as meta rownames instead
  rownames(meta) <- meta$sample
  meta$sample <- NULL
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

# output directory
out.subdir <- getOutputDir("deseq2")

# test 1: GBM (all) vs iNSC (all)

dds.1 <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData = meta, design=~type + disease_subgroup)
dds.1 <- DESeq(dds.1)

des.1 = results(dds.1, contrast=c("type", "GBM", "iNSC"))
des.1 = des.1[order(des.1$padj),]
# add gene symbol column
des.1$gene_symbol <- ens.map[rownames(des.1), 1]
# save
write.csv(des.1, file=file.path(out.subdir, "GBM_vs_iNSC_all.csv"))

# volcano plot
padj.threshold <- 1e-10
log2fc.threshold <- 5
volcano(des.1, label.field = "gene_symbol", xlim=c(-15, 15), title="GBM (all) vs iNSC (all)")

# test 2: GBM vs iNSC (RTK I)

# Add group data to meta
# This is a simple way to achieve "GMB RTK I vs iNSC RTK I"
meta$group <- factor(paste(meta$type, meta$disease_subgroup))

dds.2 <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData = meta, design=~group)
dds.2 <- DESeq(dds.2)
des.2 = results(dds.2, contrast = c("group", "GBM RTK I", "iNSC RTK I"))
des.2 = des.2[order(des.2$padj),]
# add gene symbol column
des.2$gene_symbol <- ens.map[rownames(des.2), 1]
# save
write.csv(des.2, file=file.path(out.subdir, "GBM_vs_iNSC_rtki.csv"))

volcano(des.2, label.field = "gene_symbol", xlim=c(-15, 15), title="RTK I cohort GBM vs iNSC")

# MA plot: scatter log2 FC vs mean of normalized counts, coloured based on padj
plotMA(des.2, alpha=padj.threshold, ylim=c(-5, 5))
if (FALSE) {
  # this is how to find gene names interactively
  idx <- identify(des.2$baseMean, des.2$log2FoldChange)
  des.2$gene_symbol[idx]
}

# heatmap of expression of 30 most DE genes
mostDE <- rownames(des.2)[1:30]
# get var stab transform using the original dispersion estimates to speed things up
vsd <- varianceStabilizingTransformation(dds.2, blind=FALSE)
grouping <- meta[, c("type", "disease_subgroup")]
pheatmap(assay(vsd)[mostDE,], cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = grouping, labels_row = des.2[mostDE, "gene_symbol"])

# test 3: GBM MES vs paired iNSC
des.3 <- results(dds.2, contrast = c("group", "GBM MES", "iNSC MES"))
des.3 <- des.3[order(des.3$padj),]
des.3$gene_symbol <- ens.map[rownames(des.3), 1]
# save
write.csv(des.3, file=file.path(out.subdir, "GBM_vs_iNSC_mes.csv"))

volcano(des.3, label.field = "gene_symbol", xlim=c(-15, 15), title="MES cohort GBM vs iNSC")
mostDE <- rownames(des.3)[1:30]
pheatmap(assay(vsd)[mostDE,], cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = grouping, labels_row = des.3[mostDE, "gene_symbol"])

# test 4: individual GBM vs iNSC for RTK I
# add parent ID to meta
parent.id <- sub('.*(0[1-9]+).*', '\\1', rownames(meta))
meta$parent <- parent.id
meta$groupb <- factor(paste(meta$type, meta$parent))

dds.3 <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData = meta, design=~groupb)
dds.3 <- DESeq(dds.3)

patient_comparison <- function(patient.id="018") {
  des = results(dds.3, contrast = c("groupb", paste0("GBM ", patient.id), paste0("iNSC ", patient.id)))
  des = des[order(abs(des$log2FoldChange), decreasing = T),]
  des <- des[!is.na(des$log2FoldChange),]
  des$gene_symbol <- ens.map[rownames(des), 1]
  des <- des[!is.na(des$gene_symbol),]
}

des.018 <- patient_comparison(patient.id="018")
write.csv(des.018, file=file.path(out.subdir, "GBM_vs_iNSC_018.csv"))
des.019 <- patient_comparison(patient.id="019")
write.csv(des.019, file=file.path(out.subdir, "GBM_vs_iNSC_019.csv"))
des.026 <- patient_comparison(patient.id="026")
write.csv(des.026, file=file.path(out.subdir, "GBM_vs_iNSC_026.csv"))
des.031 <- patient_comparison(patient.id="031")
write.csv(des.031, file=file.path(out.subdir, "GBM_vs_iNSC_031.csv"))



des.018 = results(dds.3, contrast = c("groupb", "GBM 018", "iNSC 018"))
des.018 = des.018[order(abs(des.018$log2FoldChange), decreasing = T),]
des.018 <- des.018[!is.na(des.018$log2FoldChange),]
des.018$gene_symbol <- ens.map[rownames(des.018), 1]
des.018 <- des.018[!is.na(des.018$gene_symbol),]

des.019 = results(dds.3, contrast = c("groupb", "GBM 019", "iNSC 019"))
des.019 = des.4[order(des.019$padj),]
des.019$gene_symbol <- ens.map[rownames(des.019), 1]

des.026 = results(dds.3, contrast = c("groupb", "GBM 026", "iNSC 026"))
des.026 = des.4[order(des.026$padj),]
des.026$gene_symbol <- ens.map[rownames(des.026), 1]

des.031 = results(dds.3, contrast = c("groupb", "GBM 031", "iNSC 031"))
des.031 = des.4[order(des.031$padj),]
des.031$gene_symbol <- ens.map[rownames(des.031), 1]

# test 4: (GBM RTK I vs paired iNSC) vs (GBM MES vs paired iNSC)
# THIS IS POINTLESS but shows how to run this kind of contrast test
# des.4 <- results(
  # dds.2, 
  # contrast = list(c("groupGBM.RTK.I", "groupiNSC.RTK.I"), c("groupGBM.MES", "groupiNSC.MES"))
# )

source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')
library(calibrate)
library("pheatmap")
library("edgeR")
library("gridExtra")

source('io/rnaseq.R')
source('io/output.R')
source('_settings.R')
source("utils.R")


prepare_de_table <- function(lrt, fdr=0.05) {
  de <- as.data.frame(topTags(lrt, p.value = fdr, n = Inf))
  de$ensembl <- rownames(lrt$table)[as.integer(rownames(de))]
  de$direction <- ifelse(de$logFC > 0, 'U', 'D')
  de <- de[, c("genes", "logFC", "ensembl", "direction", "FDR", "logCPM")]
  de
}


#' Export to CSV lists
#' We have two LRT objects (result of glmLRT). We first extract the DE genes for a given FDR in each list.
#' Then we match the DE genes and generate a CSV in blocks: (A and B), A only, B only
export_de_list <- function(lrt1, lrt2, outfile, fdr=0.05) {
  de1 <- prepare_de_table(lrt1, fdr=fdr)
  de2 <- prepare_de_table(lrt2, fdr=fdr)

  # block 1: intersection
  
  x1 <- de1[de1$ensembl %in% de2$ensembl,]
  y1 <- de2[de2$ensembl %in% de1$ensembl,]
  # order by FDR in iNSC (arbitrary)
  ord <- order(-x1$FDR)
  x1 <- x1[ord,]
  y1 <- y1[rownames(x1),]
  
  # block 2: iNSC only
  x2 <- de1[setdiff(rownames(de1), rownames(de2)),]
  y2 <- data.frame(row.names = rownames(x2))
  y2[,colnames(x2)] <- NA
  
  # block 3: eNSC only
  y3 <- y[setdiff(rownames(y), rownames(x)),]
  x3 <- data.frame(row.names = rownames(y3))
  x3[,colnames(y3)] <- NA
  
  xt <- rbind(x1, x2, x3)
  yt <- rbind(y1, y2, y3)
  
  xy <- cbind(xt, yt)
  colnames(xy) <- rep(c("HGNC Symbol", "logFC", "Ensembl ID", "Direction"), 2)
  rownames(xy) <- 1:nrow(xy)
  # replace NA with empty string
  xy[is.na(xy)] <- ''
  
  write.csv(xy, outfile, row.names = F)
  
}


loaded.wtchg <- paired_gbm_nsc_data()
dat.wtchg <- loaded.wtchg$data
meta.wtchg <- loaded.wtchg$meta

loaded.tcga <- tcga_gbm_data()
dat.tcga <- loaded.tcga$data
meta.tcga <- loaded.tcga$meta

loaded.h9 <- duan_nsc_data()
dat.h9 <- loaded.h9$data
meta.h9 <- loaded.h9$meta
n.h9 <- ncol(dat.h9)

loaded.ip <- pollard_nsc_data()
dat.ip <- loaded.ip$data
meta.ip <- loaded.ip$meta
n.ip <- ncol(dat.ip)

#' FILTER
#' The smallest library is ~10mi, the mean lib size is 45mi. 
#' We only keep genes that are expressed at CPM > 1 (i.e. >~5 counts for the avg library) in >=3 samples
#' We use non-TCGA libraries for this purpose, as they are of main interest
#' The TCGA dataset retains around 5000 more genes based on this cutoff criterion. This is probably in part down to the increased sample size.
#' It could also be down to a different protocol?
#' Only 23 genes in our gene list are not in TCGA, so we remove those to harmonise the datasets

# non-TCGA samples - no need to limit genes as they are all aligned alike
dat.non_tcga <- bind_cols(
  dat.wtchg,
  dat.ip,
  dat.h9
)
rownames(dat.non_tcga) <- rownames(dat.wtchg)

y.non_tcga <- DGEList(counts=dat.non_tcga)
keep.non_tcga <- rowSums(cpm(y.non_tcga) > 1) >= 3
print("summary(keep.non_tcga)")
print(summary(keep.non_tcga))

genes <- intersect(rownames(dat.wtchg)[keep.non_tcga], rownames(dat.tcga))
genes.no_tcga <- rownames(dat.wtchg)[keep.non_tcga]

# Restrict data to this gene list in-place
dat.wtchg <- dat.wtchg[genes,]
dat.ip <- dat.ip[genes,]
dat.h9 <- data.frame(H9NSC=dat.h9[genes,], row.names = genes)  # required for any dataset with only 1 column
dat.tcga <- dat.tcga[genes,]

#' Load annotations from biomart. 
#' We will use this to annotate DE results.
#' 
ens.map <- biomart_annotation(index.by='ensembl_gene_id')

# recreate the counts matrix with the new intersecting, filtered gene list

dat.all <- bind_cols(
  dat.wtchg,
  dat.ip,
  dat.h9,
  dat.tcga
)
rownames(dat.all) <- genes


#' Lump iNSC, eNSC, GBM together and use them to estimate dispersion

dat.lumped <- bind_cols(
  dat.wtchg,
  dat.h9
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', n.h9)
)
y.lumped <- DGEList(counts=dat.lumped, group=groups.lumped)
y.lumped <- calcNormFactors(y.lumped)

plotMDS(y.lumped)
title("MDS plot for data without IP")

design <- model.matrix(~as.factor(groups.lumped))

# this estimates tagwise dispersion
y.lumped <- estimateDisp(y.lumped, design)

#' Store the values for later use
#' This is one of the recommended approaches from the authors of edgeR in the situation where no replicates are available
dispersion.trended.lumped <- y.lumped$trended.dispersion
dispersion.common.lumped <- y.lumped$common.dispersion
dispersion.tagwise.lumped <- y.lumped$tagwise.dispersion


#' Now run the analysis with all data included, separated into the correct groups (GBM.018, etc.)


filt = meta.wtchg$disease_subgroup == 'RTK I'

dat <- bind_cols(
  dat.wtchg[, filt],
  dat.h9
)
rownames(dat) <- genes
grp = data.frame(
  cell_type=c(rep(c('GBM', 'iNSC'), each=3), 'eNSC'),
  patient=c(rep(c('018', '019', '031'), 2), ''), 
  row.names = colnames(dat)
)
grp <- factor(paste(grp$cell_type, grp$patient, sep="."))

y <- DGEList(dat, genes = ens.map[rownames(dat), "hgnc_symbol"])
y <- calcNormFactors(y)

design <- model.matrix(~0 + grp)
colnames(design) <- levels(grp)

# in this case, we need to use the dispersion estimated earlier
y$common.dispersion <- dispersion.common.lumped
y$trended.dispersion <- dispersion.trended.lumped
y$tagwise.dispersion <- dispersion.tagwise.lumped

fit.glm <- glmFit(y, design)
my.contrasts <- makeContrasts(
  GBMvsiNSC=(GBM.018+GBM.019+GBM.031)-(iNSC.018+iNSC.019+iNSC.031), 
  GBM018vsiNSC018=GBM.018-iNSC.018,
  GBM019vsiNSC019=GBM.019-iNSC.019,
  GBM031vsiNSC031=GBM.031-iNSC.031,
  GBMvseNSC=(GBM.018+GBM.019+GBM.031)/3 - eNSC.,
  GBM018vseNSC=GBM.018-eNSC.,
  GBM019vseNSC=GBM.019-eNSC.,
  GBM031vseNSC=GBM.031-eNSC.,
  levels=design
)
# we need to take the negative effect to get GBM relative to iNSC
lrt <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvsiNSC"])
lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vsiNSC018"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vseNSC"])
lrt <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vsiNSC019"])
lrt <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vsiNSC031"])

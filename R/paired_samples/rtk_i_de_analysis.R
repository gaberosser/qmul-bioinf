source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')
library(calibrate)
library("pheatmap")
library("edgeR")

source('io/rnaseq.R')
source('io/output.R')
source('_settings.R')
source("utils.R")

loaded.wtchg <- paired_gbm_nsc_data()
dat.wtchg <- loaded.wtchg$data
meta.wtchg <- loaded.wtchg$meta

loaded.tcga <- tcga_gbm_data()
dat.tcga <- loaded.tcga$data
meta.tcga <- loaded.tcga$meta

loaded.h9 <- duan_nsc_data()
dat.h9 <- loaded.h9$data
meta.h9 <- loaded.h9$meta

loaded.ip <- pollard_nsc_data()
dat.ip <- loaded.ip$data
meta.ip <- loaded.ip$meta

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
summary(keep.non_tcga)

genes <- intersect(rownames(dat.wtchg)[keep.non_tcga], rownames(dat.tcga))

# recreate the counts matrix with the new intersecting, filtered gene list

dat.all <- bind_cols(
  dat.wtchg[genes,],
  dat.ip[genes,],
  dat.h9[genes,],
  dat.tcga[genes,]
)
rownames(dat.all) <- genes

#' MDS plot for all data
#' Strangely, this suggests that the Pollard data are very different from the other NSC samples
#' As a result, perhaps we shouldn't use them?
y.all <- DGEList(dat.all)
y.all <- calcNormFactors(y.all)
plotMDS(y.all)
title("MDS plot for all data")

#' Lump iNSC, eNSC, GBM together and use them to estimate dispersion

dat.lumped <- bind_cols(
  dat.wtchg[genes,],
  dat.h9[genes,]
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', 2)
)
y.lumped <- DGEList(counts=dat.lumped, group=groups.lumped)
y.lumped <- calcNormFactors(y.lumped)

plotMDS(y.lumped)
title("MDS plot for used data")

design <- model.matrix(~as.factor(groups.lumped))

# this estimates tagwise dispersion
y.lumped <- estimateDisp(y.lumped, design)
plotBCV(y.lumped)
title("BCV lumped data")

#' Store the values for later use
#' This is one of the recommended approaches from the authors of edgeR in the situation where no replicates are available
dispersion.trended.lumped <- y.lumped$trended.dispersion
dispersion.common.lumped <- y.lumped$common.dispersion
dispersion.tagwise.lumped <- y.lumped$tagwise.dispersion

run_one_comparison <- function(gbm.sample_name, insc.sample_name, include.ip=F, p.value=0.05) {
  if (include.ip) {
    this.dat <- bind_cols(
      dat.wtchg[genes, c(gbm.sample_name, insc.sample_name)],
      dat.h9[genes,],
      dat.ip[genes,]
    )
    this.groups <- c(gbm.sample_name, insc.sample_name, rep('NSC', 4))    
  } else {
    this.dat <- bind_cols(
      dat.wtchg[genes, c(gbm.sample_name, insc.sample_name)],
      dat.h9[genes,]
    )
    this.groups <- c('GBM', 'iNSC', rep('eNSC', 2))
  }
  rownames(this.dat) <- genes
  this.y <- DGEList(this.dat, group = this.groups)
  this.y <- calcNormFactors(this.y)
  this.y$common.dispersion <- dispersion.common.lumped
  this.y$trended.dispersion <- dispersion.trended.lumped
  # THIS LINE HAS A VERY SIGNIFICANT EFFECT:
  # this.y$tagwise.dispersion <- dispersion.tagwise.lumped
  
  design <- model.matrix(~0+group, data=this.y$samples)
  colnames(design) <- levels(this.y$samples$group)
  contrasts = makeContrasts(GBMvsiNSC=GBM-iNSC, GBMvseNSC=GBM-eNSC, levels=design)
  
  fit <- glmFit(this.y, design)
  lrt.paired <- glmLRT(fit, contrast=contrasts[,"GBMvsiNSC"])
  print(
    paste0("LR test finds ", dim(topTags(lrt.paired, n=Inf, p.value=p.value))[1], " significantly DE genes between ", gbm.sample_name, " and ", insc.sample_name)
  )
  toptags.gbm_insc <- topTags(lrt.paired, n=Inf, p.value=p.value)
  
  lrt.ensc <- glmLRT(fit, contrast=contrasts[,"GBMvseNSC"])
  print(
    paste0("LR test finds ", dim(topTags(lrt.ensc, n=Inf, p.value=p.value))[1], " significantly DE genes between ", gbm.sample_name, " and endogenous NSC")
  )
  toptags.gbm_ensc <- topTags(lrt.ensc, n=Inf, p.value=p.value)
  
  return(list(gbm_insc=toptags.gbm_insc, gbm_ensc=toptags.gbm_ensc))
}


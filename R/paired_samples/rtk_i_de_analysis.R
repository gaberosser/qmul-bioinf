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
print(summary(keep.non_tcga))

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
title("MDS plot for data without IP")

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

#' run it again with the IP data included and store the dispersion estimates again
dat.lumped <- bind_cols(
dat.wtchg[genes,],
dat.h9[genes,],
dat.ip[genes,]
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', 4)
)
y.lumped <- DGEList(counts=dat.lumped, group=groups.lumped)
y.lumped <- calcNormFactors(y.lumped)
design <- model.matrix(~as.factor(groups.lumped))

# this estimates tagwise dispersion
y.lumped <- estimateDisp(y.lumped, design)
dispersion.trended.lumped_ip <- y.lumped$trended.dispersion
dispersion.common.lumped_ip <- y.lumped$common.dispersion
dispersion.tagwise.lumped_ip <- y.lumped$tagwise.dispersion


#' Load annotations from biomart. 
#' We will use this to annotate DE results.
ens.map <- biomart_annotation(index.by='ensembl_gene_id')

setup_comparison_data <- function(gbm.sample_name, insc.sample_name, include.ip=F) {
  if (include.ip) {
    this.dat <- bind_cols(
      dat.wtchg[genes, c(gbm.sample_name, insc.sample_name)],
      dat.h9[genes,],
      dat.ip[genes,]
    )
    this.groups <- c('GBM', 'iNSC', rep('eNSC', 4))
  } else {
    this.dat <- bind_cols(
      dat.wtchg[genes, c(gbm.sample_name, insc.sample_name)],
      dat.h9[genes,]
    )
    this.groups <- c('GBM', 'iNSC', rep('eNSC', 2))
  }
  rownames(this.dat) <- genes
  return(list(
    dat=this.dat,
    groups=this.groups
  ))
}


run_one_comparison <- function(
  gbm.sample_name, 
  insc.sample_name, 
  include.ip=F, 
  p.value=0.05
) {
  loaded <- setup_comparison_data(gbm.sample_name, insc.sample_name, include.ip = include.ip)
  this.y <- DGEList(loaded$dat, group = loaded$groups, genes = ens.map[rownames(loaded$dat), "hgnc_symbol"])
  this.y <- calcNormFactors(this.y)
  if (include.ip) {
    this.y$common.dispersion <- dispersion.common.lumped_ip
    this.y$trended.dispersion <- dispersion.trended.lumped_ip   
    # THIS LINE HAS A VERY SIGNIFICANT EFFECT:
    # this.y$tagwise.dispersion <- dispersion.tagwise.lumped_ip
  } else {
    this.y$common.dispersion <- dispersion.common.lumped
    this.y$trended.dispersion <- dispersion.trended.lumped
    # THIS LINE HAS A VERY SIGNIFICANT EFFECT:
    # this.y$tagwise.dispersion <- dispersion.tagwise.lumped 
  }
  
  design <- model.matrix(~0+group, data=this.y$samples)
  colnames(design) <- levels(this.y$samples$group)
  contrasts = makeContrasts(GBMvsiNSC=GBM-iNSC, GBMvseNSC=GBM-eNSC, levels=design)
  
  fit <- glmFit(this.y, design)
  lrt.paired <- glmLRT(fit, contrast=contrasts[,"GBMvsiNSC"])
  print(
    paste0("LR test finds ", dim(topTags(lrt.paired, n=Inf, p.value=p.value))[1], " significantly DE genes between ", gbm.sample_name, " and ", insc.sample_name)
  )
  toptags.gbm_insc <- as.data.frame(topTags(lrt.paired, n=Inf, p.value=p.value))
  # fix rownames
  rownames(toptags.gbm_insc) <- rownames(this.y$counts)[as.integer(rownames(toptags.gbm_insc))]

  lrt.ensc <- glmLRT(fit, contrast=contrasts[,"GBMvseNSC"])
  print(
    paste0("LR test finds ", dim(topTags(lrt.ensc, n=Inf, p.value=p.value))[1], " significantly DE genes between ", gbm.sample_name, " and endogenous NSC")
  )
  toptags.gbm_ensc <- as.data.frame(topTags(lrt.ensc, n=Inf, p.value=p.value))
  # fix rownames
  rownames(toptags.gbm_ensc) <- rownames(this.y$counts)[as.integer(rownames(toptags.gbm_ensc))]

  return(list(y=this.y, gbm_insc=toptags.gbm_insc, gbm_ensc=toptags.gbm_ensc))
}

res.018 <- run_one_comparison('GBM018', 'DURA018N2_NSC', include.ip = F)
res.ip.018 <- run_one_comparison('GBM018', 'DURA018N2_NSC', include.ip = T)

res.019 <- run_one_comparison('GBM019', 'DURA019N8C_NSC', include.ip = F)
res.ip.019 <- run_one_comparison('GBM019', 'DURA019N8C_NSC', include.ip = T)

res.031 <- run_one_comparison('GBM031', 'DURA031N44B_NSC', include.ip = F)
res.ip.031 <- run_one_comparison('GBM031', 'DURA031N44B_NSC', include.ip = T)

#' Export to CSV lists
export_de_list <- function(res, outfile) {
  x <- res$gbm_insc[, c("genes", "logFC")]
  x[, 'ensembl'] <- rownames(x)
  x[, 'direction'] <- ifelse(x$logFC > 0, 'U', 'D')

  y <- res$gbm_ensc[, c("genes", "logFC")]
  y[, 'ensembl'] <- rownames(y)
  y[, 'direction'] <- ifelse(y$logFC > 0, 'U', 'D')
  
  # block 1: intersection
  
  x1 <- x[x$ensembl %in% y$ensembl,]
  y1 <- y[y$ensembl %in% x$ensembl,]
  # order by logFC in iNSC (arbitrary)
  ord <- order(-abs(x1$logFC))
  x1 <- x1[ord,]
  y1 <- y1[rownames(x1),]
  
  # block 2: iNSC only
  x2 <- x[setdiff(rownames(x), rownames(y)),]
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
list.outdir <- getOutputDir(name = "paired_analysis_de_gene_lists")

export_de_list(res.018, file.path(list.outdir, "gbm018_nsc_de_gene.csv"))
export_de_list(res.019, file.path(list.outdir, "gbm019_nsc_de_gene.csv"))
export_de_list(res.031, file.path(list.outdir, "gbm031_nsc_de_gene.csv"))

# plot venn diagrams showing number of genes overlapping
library(VennDiagram)
plot_venn <- function(res, obj.title) {
  area1 = nrow(res$gbm_insc)
  area2 = nrow(res$gbm_ensc)
  cross_area = length(intersect(
    rownames(res$gbm_insc), rownames(res$gbm_ensc)
  ))
  grid.newpage()
  venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c(paste0(obj.title, " vs paired NSC"), paste0(obj.title, " vs reference NSC")),
                     lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                     cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
  grid.draw(venn.plot)
}

plot_venn(res.018, "GBM018")
plot_venn(res.019, "GBM019")
plot_venn(res.031, "GBM031")

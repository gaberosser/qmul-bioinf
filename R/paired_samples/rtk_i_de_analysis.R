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

#' MDS plot for all data
#' Strangely, this suggests that the Pollard data are very different from the other NSC samples
#' As a result, perhaps we shouldn't use them?
y.all <- DGEList(dat.all)
y.all <- calcNormFactors(y.all)
plotMDS(y.all)
title("MDS plot for all data")

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
plotBCV(y.lumped)
title("BCV lumped data")

#' Store the values for later use
#' This is one of the recommended approaches from the authors of edgeR in the situation where no replicates are available
dispersion.trended.lumped <- y.lumped$trended.dispersion
dispersion.common.lumped <- y.lumped$common.dispersion
dispersion.tagwise.lumped <- y.lumped$tagwise.dispersion

#' run it again with the IP data included and store the dispersion estimates again
dat.lumped <- bind_cols(
dat.wtchg,
dat.h9,
dat.ip
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', n.h9 + n.ip)
)
y.lumped <- DGEList(counts=dat.lumped, group=groups.lumped)
y.lumped <- calcNormFactors(y.lumped)
design <- model.matrix(~as.factor(groups.lumped))

# this estimates tagwise dispersion
y.lumped <- estimateDisp(y.lumped, design)
dispersion.trended.lumped_ip <- y.lumped$trended.dispersion
dispersion.common.lumped_ip <- y.lumped$common.dispersion
dispersion.tagwise.lumped_ip <- y.lumped$tagwise.dispersion

setup_comparison_data <- function(gbm.sample_name, insc.sample_name, include.ip=F) {
  if (include.ip) {
    this.dat <- bind_cols(
      dat.wtchg[, c(gbm.sample_name, insc.sample_name)],
      dat.h9,
      dat.ip
    )
    this.groups <- c('GBM', 'iNSC', rep('eNSC', n.h9 + n.ip))
  } else {
    this.dat <- bind_cols(
      dat.wtchg[, c(gbm.sample_name, insc.sample_name)],
      dat.h9
    )
    this.groups <- c('GBM', 'iNSC', rep('eNSC', n.h9))
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

  lrt.ref <- glmLRT(fit, contrast=contrasts[,"GBMvseNSC"])
  print(
    paste0("LR test finds ", dim(topTags(lrt.ref, n=Inf, p.value=p.value))[1], " significantly DE genes between ", gbm.sample_name, " and endogenous NSC")
  )
  toptags.gbm_ensc <- as.data.frame(topTags(lrt.ref, n=Inf, p.value=p.value))
  # fix rownames
  rownames(toptags.gbm_ensc) <- rownames(this.y$counts)[as.integer(rownames(toptags.gbm_ensc))]

  return(list(
    y=this.y, 
    gbm_insc=toptags.gbm_insc, 
    gbm_ensc=toptags.gbm_ensc,
    lrt.paired=lrt.paired,
    lrt.ref=lrt.ref
    ))
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
list.outdir <- getOutputDir(name = "paired_analysis_de_rtkI")

export_de_list(res.018, file.path(list.outdir, "gbm018_nsc_h9.csv"))
export_de_list(res.ip.018, file.path(list.outdir, "gbm018_nsc_all.csv"))

export_de_list(res.019, file.path(list.outdir, "gbm019_nsc_h9.csv"))
export_de_list(res.ip.019, file.path(list.outdir, "gbm019_nsc_all.csv"))

export_de_list(res.031, file.path(list.outdir, "gbm031_nsc_h9.csv"))
export_de_list(res.ip.031, file.path(list.outdir, "gbm031_nsc_all.csv"))

#' Run the paired sample analysis.
#' Here we're looking for common effects present across all 3 pairs.

y.paired <- DGEList(dat.wtchg[, meta.wtchg$disease_subgroup == 'RTK I'], genes = ens.map[rownames(dat.wtchg), "hgnc_symbol"])
y.paired <- calcNormFactors(y.paired)
groups.paired <- data.frame(row.names = rownames(meta.wtchg[meta.wtchg$disease_subgroup == 'RTK I',]))
groups.paired$cell_type <- c(rep('GBM', 3), rep('iNSC', 3))
groups.paired$patient <- c(rep('018', 2), rep('019', 2), rep('031', 2))  
design <- model.matrix(~0+patient+cell_type, data=groups.paired)
y.paired <- estimateDisp(y.paired, design)
fit.glm <- glmFit(y.paired, design)
# we need to take the negative effect to get GBM relative to iNSC
lrt <- glmLRT(fit.glm, contrast=c(0, 0, 0, -1))
de <- as.data.frame(topTags(lrt, n=Inf, p.value=0.05))

de$ensembl <- rownames(y.paired$counts)[as.integer(rownames(de))]
de$direction <- ifelse(de$logFC > 0, 'U', 'D')
de <- de[, c("genes", "logFC", "ensembl", "direction")]
colnames(de) <- c("HGNC Symbol", "logFC", "Ensembl ID", "Direction")
# save
write.csv(de, file.path(list.outdir, "all_gbm_paired.csv"), row.names = F)

#' Run for all GBM (lumped) vs all reference NSC
#' Here we're looking for common effects present across all 3 pairs.

filt <- meta.wtchg$disease_subgroup == 'RTK I' & meta.wtchg$type == 'GBM'
dat.lumped <- bind_cols(
  dat.wtchg[, filt],
  dat.h9,
  dat.ip
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg[filt, 'type']),
  rep('eNSC', n.h9 + n.ip)
)
y.lumped <- DGEList(dat.lumped, genes = ens.map[rownames(dat.lumped), "hgnc_symbol"])
y.lumped <- calcNormFactors(y.lumped)
design <- model.matrix(~0+groups.lumped)
y.lumped <- estimateDisp(y.lumped, design)
fit.glm.lumped <- glmFit(y.lumped, design)
lrt.lumped <- glmLRT(fit.glm.lumped, contrast=c(-1, 1))  # GBM rel to eNSC
de <- as.data.frame(topTags(lrt.lumped, n=Inf, p.value=0.05))

print(paste0("When we compare all our RTK I GBM samples against all ref NSC samples, we find ", nrow(de), " DE genes."))

de$ensembl <- rownames(y.paired$counts)[as.integer(rownames(de))]
de$direction <- ifelse(de$logFC > 0, 'U', 'D')
de <- de[, c("genes", "logFC", "ensembl", "direction")]
colnames(de) <- c("HGNC Symbol", "logFC", "Ensembl ID", "Direction")
# save
write.csv(de, file.path(list.outdir, "all_gbm_vs_all_reference.csv"), row.names = F)

# repeat with H9 only

dat.lumped <- bind_cols(
  dat.wtchg[, filt],
  dat.h9
)
rownames(dat.lumped) <- genes
groups.lumped <- c(
  as.vector(meta.wtchg[filt, 'type']),
  rep('eNSC', n.h9)
)
y.lumped <- DGEList(dat.lumped, genes = ens.map[rownames(dat.lumped), "hgnc_symbol"])
y.lumped <- calcNormFactors(y.lumped)
design <- model.matrix(~0+groups.lumped)
y.lumped <- estimateDisp(y.lumped, design)
fit.glm.lumped <- glmFit(y.lumped, design)
lrt.lumped <- glmLRT(fit.glm.lumped, contrast=c(-1, 1)) # GBM rel to eNSC
de <- as.data.frame(topTags(lrt.lumped, n=Inf, p.value=0.05))

print(paste0("When we compare all our RTK I GBM samples against only the H9 sample, we find ", nrow(de), " DE genes."))

de$ensembl <- rownames(y.lumped$counts)[as.integer(rownames(de))]
de$direction <- ifelse(de$logFC > 0, 'U', 'D')
de <- de[, c("genes", "logFC", "ensembl", "direction")]
colnames(de) <- c("HGNC Symbol", "logFC", "Ensembl ID", "Direction")
# save
write.csv(de, file.path(list.outdir, "all_gbm_vs_h9_reference.csv"), row.names = F)


#' Run GO (and KEGG) analysis
#' This will highlight the pathways that are enriched in the gene lists. Indexing is by Entrez ID
#' Requires GO.db: biocLite("GO.db")
# go.018.paired <- goana(res.018$lrt.paired, geneid = ens.map[rownames(res.018$lrt.paired), 'entrezgene'])
# go.018.ref <- goana(res.018$lrt.ref, geneid = ens.map[rownames(res.018$lrt.ref), 'entrezgene'])
# go.019.paired <- goana(res.019$lrt.paired, geneid = ens.map[rownames(res.019$lrt.paired), 'entrezgene'])
# go.019.ref <- goana(res.019$lrt.ref, geneid = ens.map[rownames(res.019$lrt.ref), 'entrezgene'])
# go.031.paired <- goana(res.031$lrt.paired, geneid = ens.map[rownames(res.031$lrt.paired), 'entrezgene'])
# go.031.ref <- goana(res.031$lrt.ref, geneid = ens.map[rownames(res.031$lrt.ref), 'entrezgene'])

# plot venn diagrams showing number of genes overlapping
library(VennDiagram)
plot_venn <- function(res, obj.title=NULL, plot.direction=NULL, png.file=NULL, pdf.file=NULL) {
  #' plot.direction: NULL, 'up' or 'down'. Controls which DE genes are included (NULL means all)
  if (is.null(plot.direction)) {
    x1 = res$gbm_insc
    x2 = res$gbm_ensc
  } else if (tolower(plot.direction) == 'up') {
    x1 = res$gbm_insc[res$gbm_insc$logFC > 0,]
    x2 = res$gbm_ensc[res$gbm_ensc$logFC > 0,]
  } else if (tolower(plot.direction) == 'down') {
    x1 = res$gbm_insc[res$gbm_insc$logFC < 0,]
    x2 = res$gbm_ensc[res$gbm_ensc$logFC < 0,]
  }
  area1 = nrow(x1)
  area2 = nrow(x2)
  cross_area = length(intersect(
    rownames(x1), rownames(x2)
  ))

  if (!is.null(png.file)) {
    png(png.file)
    plot.new()
    venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
                                    lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                    cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
    title(obj.title)
    dev.off()
  }
  if (!is.null(pdf.file)) {
    pdf(pdf.file)
    plot.new()
    venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
                                    lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                    cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
    title(obj.title)
    dev.off()
  }
  
  plot.new()
  venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
                     lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                     cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
  title(obj.title)
  return(venn.plot)
}

plot_multivenn <- function(res, png.file=NULL) {
  x.1 = res$gbm_insc
  x.2 = res$gbm_ensc

  x.up.1 = res$gbm_insc[res$gbm_insc$logFC > 0,]
  x.up.2 = res$gbm_ensc[res$gbm_ensc$logFC > 0,]

  x.down.1 = res$gbm_insc[res$gbm_insc$logFC < 0,]
  x.down.2 = res$gbm_ensc[res$gbm_ensc$logFC < 0,]
    
  area.1 = nrow(x.1)
  area.2 = nrow(x.2)
  cross_area = length(intersect(
    rownames(x.1), rownames(x.2)
  ))
  
  area.up.1 <- nrow(x.up.1)
  area.up.2 <- nrow(x.up.2)
  cross_area.up <- length(intersect(
    rownames(x.up.1), rownames(x.up.2)
  ))
  
  area.down.1 <- nrow(x.down.1)
  area.down.2 <- nrow(x.down.2)
  cross_area.down <- length(intersect(
    rownames(x.down.1), rownames(x.down.2)
  ))
  
  if (!is.null(png.file)) {
    png(png.file, width=600, height=1200, units='px', res=120)
  }
  
  venn <- draw.pairwise.venn(area.1, area.2, cross_area, category = c("GBM - paired NSC", "GBM - reference NSC"),
                             lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                             cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
  venn.up <- draw.pairwise.venn(area.up.1, area.up.2, cross_area.up, category = c("GBM - paired NSC (UP)", "GBM - reference NSC (UP)"),
                             lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                             cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
  venn.down <- draw.pairwise.venn(area.down.1, area.down.2, cross_area.down, category = c("GBM - paired NSC (DOWN)", "GBM - reference NSC (DOWN)"),
                             lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                             cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
  grid.arrange(gTree(children=venn), gTree(children=venn.up), gTree(children=venn.down), nrow=3)

  if (!is.null(png.file)) {
    dev.off()
  }
}


plot_multivenn(res.018, png.file=file.path(list.outdir, "gbm018.png"))
plot_multivenn(res.ip.018, png.file=file.path(list.outdir, "gbm018_ip.png"))

plot_multivenn(res.019, png.file=file.path(list.outdir, "gbm019.png"))
plot_multivenn(res.ip.019, png.file=file.path(list.outdir, "gbm019_ip.png"))

plot_multivenn(res.031, png.file=file.path(list.outdir, "gbm031.png"))
plot_multivenn(res.ip.031, png.file=file.path(list.outdir, "gbm031_ip.png"))


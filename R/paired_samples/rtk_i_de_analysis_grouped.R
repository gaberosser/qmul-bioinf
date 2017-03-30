source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')
library(calibrate)
library("pheatmap")
library("edgeR")
library("gridExtra")
library(VennDiagram)
library(ggplot2)

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

#' Get intersecting and unique results
decompose_de_lists <- function(lrt1, lrt2, fdr=0.05) {
  de1 <- prepare_de_table(lrt1, fdr=fdr)
  de2 <- prepare_de_table(lrt2, fdr=fdr)
  
  # block 1: intersection
  
  x1 <- de1[de1$ensembl %in% de2$ensembl,]
  y1 <- de2[de2$ensembl %in% de1$ensembl,]
  # order by FDR in iNSC (arbitrary)
  ord <- order(x1$FDR)
  x1 <- x1[ord,]
  y1 <- y1[rownames(x1),]
  
  # block 2: iNSC only
  x2 <- de1[setdiff(rownames(de1), rownames(de2)),]
  y2 <- data.frame(row.names = rownames(x2))
  y2[,colnames(x2)] <- NA
  
  # block 3: eNSC only
  y3 <- de2[setdiff(rownames(de2), rownames(de1)),]
  x3 <- data.frame(row.names = rownames(y3))
  x3[,colnames(y3)] <- NA
  
  return(list(x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3))
}


#' Get intersecting and unique results for 3 components
decompose_de_lists.3 <- function(lrt1, lrt2, lrt3, fdr=0.05) {
  de1 <- prepare_de_table(lrt1, fdr=fdr)
  de2 <- prepare_de_table(lrt2, fdr=fdr)
  de3 <- prepare_de_table(lrt3, fdr=fdr)
  de <- list(
    de1, de2, de3
  )
  ens <- list(
    de1$ensembl,
    de2$ensembl,
    de3$ensembl
  )
  
  blocks <- list()
  
  for (i in seq(1, 7)) {
    comb <- as.integer(intToBits(i))[1:3]
    idx.in <- which(comb == 1)
    idx.out <- which(comb == 0)
    
    ens.in <- Reduce(intersect, lapply(idx.in, function(x){ ens[[x]] }))
    ens.out <- Reduce(union, lapply(idx.out, function(x){ ens[[x]] }))
    ens.this <- setdiff(ens.in, ens.out)
    
    get_de <- function(x) {
      tmp <- de[[x]][de[[x]]$ensembl %in% ens.this,]
      tmp <- tmp[order(tmp$ensembl),]
    }
    
    de.in <- lapply(idx.in, get_de)
    
    blocks[[paste0(comb, collapse = '')]] <- do.call(cbind, de.in)
  }
  
  blocks
}


#' Export to CSV lists
#' We have two LRT objects (result of glmLRT). We first extract the DE genes for a given FDR in each list.
#' Then we match the DE genes and generate a CSV in blocks: (A and B), A only, B only
export_de_list <- function(lrt1, lrt2, outfile, fdr=0.05) {

  res <- decompose_de_lists(lrt1, lrt2, fdr)
  x1 <- res$x1
  x2 <- res$x2
  x3 <- res$x3
  y1 <- res$y1
  y2 <- res$y2
  y3 <- res$y3
  
  xt <- rbind(x1, x2, x3)
  yt <- rbind(y1, y2, y3)
  
  xy <- cbind(xt, yt)
  # colnames(xy) <- rep(c("HGNC Symbol", "logFC", "Ensembl ID", "Direction"), 2)
  rownames(xy) <- 1:nrow(xy)
  # replace NA with empty string
  xy[is.na(xy)] <- ''
  
  write.csv(xy, outfile, row.names = F)
  return(list(de1=prepare_de_table(lrt1, fdr=fdr), de2=prepare_de_table(lrt2, fdr=fdr)))
}

# plot venn diagrams showing number of genes overlapping in 2 DE results
plot_multivenn <- function(de1, de2, png.file=NULL) {

  de1.up = de1[de1$logFC > 0,]
  de2.up = de2[de2$logFC > 0,]
  
  de1.down = de1[de1$logFC < 0,]
  de2.down = de2[de2$logFC < 0,]
  
  area.1 = nrow(de1)
  area.2 = nrow(de2)
  cross_area = length(intersect(
    rownames(de1), rownames(de2)
  ))
  
  area.up.1 <- nrow(de1.up)
  area.up.2 <- nrow(de2.up)
  cross_area.up <- length(intersect(
    rownames(de1.up), rownames(de2.up)
  ))
  
  area.down.1 <- nrow(de1.down)
  area.down.2 <- nrow(de2.down)
  cross_area.down <- length(intersect(
    rownames(de1.down), rownames(de2.down)
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


# plot venn diagrams showing number of genes overlapping in 3 DE results
plot_multivenn.3 <- function(blocks, png.file=NULL) {
  
  n123 = nrow(blocks[["111"]])
  
  n1 = nrow(blocks[["100"]])
  n2 = nrow(blocks[["010"]])
  n3 = nrow(blocks[["001"]])
  
  n12 = nrow(blocks[["110"]])
  n23 = nrow(blocks[["011"]])
  n13 = nrow(blocks[["101"]])
  
  area1 = n1 + n12 + n13 + n123
  area2 = n2 + n12 + n23 + n123
  area3 = n3 + n13 + n23 + n123
  
  n12 = n12 + n123
  n23 = n23 + n123
  n13 = n13 + n123
  
  category = c("GBM - paired NSC", "GBM - H9 NSC", "GBM - IP NSC")

  if (!is.null(png.file)) {
    png(png.file, width=800, height=800, units='px', res=120)
  }
  
  draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = category, scaled = T)

  if (!is.null(png.file)) {
    dev.off()
  }
}


loaded.wtchg <- paired_gbm_nsc_data()
dat.wtchg <- loaded.wtchg$data
meta.wtchg <- loaded.wtchg$meta

loaded.tcga <- tcga_gbm_data()
dat.tcga <- loaded.tcga$data
meta.tcga <- loaded.tcga$meta

loaded.h9.sep <- duan_nsc_data(collapse.replicates = F)
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
dat.h9 <- dat.h9[genes, , drop=F]
dat.tcga <- dat.tcga[genes,]

#' Load annotations from biomart. 
#' We will use this to annotate DE results.
#' 
ens.map <- biomart_annotation(index.by='ensembl_gene_id')

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

#' QQ plot of deviance values, used to show the utility of tagwise estimates
#' This will not work for a 'saturated' design matrix because the dof = 0
#' Basically, the results show that you have a lot of poor fits to the model if you don't use genewise dispersion estimates
qqplot_chisq_fit <- function(y, design, dispersion, outfile=NULL, title=NULL) {

  fit <- glmFit(y, design, dispersion=dispersion)
  p <- pchisq(fit$deviance, fit$df.residual, lower.tail = F)
  p.holm <- p.adjust(p, method='holm')
  z <- zscoreGamma(fit$deviance, shape=fit$df.residual/2, scale=2)
  z[abs(z) == Inf] <- NA
  col <- ifelse(p.holm < 0.05, 'blue', 'black')

  plot.new()
  if (!is.null(outfile)){
    png(outfile)
  }
  
  qqnorm(z, pch=16, col=col, main=title)
  qqline(z)

  if (!is.null(outfile)) {
    dev.off()
  }
  
}

list.outdir <- getOutputDir(name = "paired_analysis_de_rtkI")

qqplot_chisq_fit(y.lumped, design, dispersion.common.lumped, outfile=file.path(list.outdir, "qqplot_fit_common.png"), title="Common dispersion")
qqplot_chisq_fit(y.lumped, design, dispersion.trended.lumped, outfile=file.path(list.outdir, "qqplot_fit_trended.png"), title="Trended dispersion")
qqplot_chisq_fit(y.lumped, design, dispersion.tagwise.lumped, outfile=file.path(list.outdir, "qqplot_fit_tagwise.png"), title="Genewise dispersion")

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
  GBMvsiNSC=(GBM.018+GBM.019+GBM.031)/3-(iNSC.018+iNSC.019+iNSC.031)/3, 
  GBM018vsiNSC018=GBM.018-iNSC.018,
  GBM019vsiNSC019=GBM.019-iNSC.019,
  GBM031vsiNSC031=GBM.031-iNSC.031,
  GBMvseNSC=(GBM.018+GBM.019+GBM.031)/3 - eNSC.,
  GBM018vseNSC=GBM.018-eNSC.,
  GBM019vseNSC=GBM.019-eNSC.,
  GBM031vseNSC=GBM.031-eNSC.,
  levels=design
)

lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvsiNSC"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvseNSC"])
res.tot <- export_de_list(lrt1, lrt2, file.path(list.outdir, "gbm-insc-ensc-all.csv"))

lrt1.018 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vsiNSC018"])
lrt2.018 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vseNSC"])
res.018 <- export_de_list(lrt1.018, lrt2.018, file.path(list.outdir, "gbm-insc-ensc-018.csv"))

lrt1.019 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vsiNSC019"])
lrt2.019 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vseNSC"])
res.019 <- export_de_list(lrt1.019, lrt2.019, file.path(list.outdir, "gbm-insc-ensc-019.csv"))

lrt1.031 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vsiNSC031"])
lrt2.031 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vseNSC"])
res.031 <- export_de_list(lrt1.031, lrt2.031, file.path(list.outdir, "gbm-insc-ensc-031.csv"))

plot_multivenn(res.tot$de1, res.tot$de2, png.file=file.path(list.outdir, "all.png"))
plot_multivenn(res.018$de1, res.018$de2, png.file=file.path(list.outdir, "018.png"))
plot_multivenn(res.019$de1, res.019$de2, png.file=file.path(list.outdir, "019.png"))
plot_multivenn(res.031$de1, res.031$de2, png.file=file.path(list.outdir, "031.png"))

#' Repeat but include TWO reference datasets

filt = meta.wtchg$disease_subgroup == 'RTK I'

dat <- bind_cols(
  dat.wtchg[, filt],
  dat.h9,
  dat.ip
)
rownames(dat) <- genes
grp = data.frame(
  cell_type=c(rep(c('GBM', 'iNSC'), each=3), 'eNSC.H9', rep('eNSC.IP', n.ip)),
  patient=c(rep(c('018', '019', '031'), 2), rep('', n.h9 + n.ip)), 
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
  GBMvsiNSC=(GBM.018+GBM.019+GBM.031)/3-(iNSC.018+iNSC.019+iNSC.031)/3, 
  GBM018vsiNSC018=GBM.018-iNSC.018,
  GBM019vsiNSC019=GBM.019-iNSC.019,
  GBM031vsiNSC031=GBM.031-iNSC.031,
  GBMvseNSCH9=(GBM.018+GBM.019+GBM.031)/3 - eNSC.H9.,
  GBMvseNSCIP=(GBM.018+GBM.019+GBM.031)/3 - eNSC.IP.,
  GBM018vseNSCH9=GBM.018-eNSC.H9.,
  GBM019vseNSCH9=GBM.019-eNSC.H9.,
  GBM031vseNSCH9=GBM.031-eNSC.H9.,
  GBM018vseNSCIP=GBM.018-eNSC.IP.,
  GBM019vseNSCIP=GBM.019-eNSC.IP.,
  GBM031vseNSCIP=GBM.031-eNSC.IP.,
  levels=design
)


lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvsiNSC"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvseNSCH9"])
lrt3 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBMvseNSCIP"])
blocks <- decompose_de_lists.3(lrt1, lrt2, lrt3)
plot_multivenn.3(blocks, png.file = file.path(list.outdir, "gbm-insc-ensc_all.png"))

lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vsiNSC018"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vseNSCH9"])
lrt3 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vseNSCIP"])
blocks <- decompose_de_lists.3(lrt1, lrt2, lrt3)
plot_multivenn.3(blocks, png.file = file.path(list.outdir, "gbm-insc-ensc_018.png"))

lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vsiNSC019"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vseNSCH9"])
lrt3 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vseNSCIP"])
blocks <- decompose_de_lists.3(lrt1, lrt2, lrt3)
plot_multivenn.3(blocks, png.file = file.path(list.outdir, "gbm-insc-ensc_019.png"))

lrt1 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vsiNSC031"])
lrt2 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vseNSCH9"])
lrt3 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vseNSCIP"])
blocks <- decompose_de_lists.3(lrt1, lrt2, lrt3)
plot_multivenn.3(blocks, png.file = file.path(list.outdir, "gbm-insc-ensc_031.png"))


res.tot <- export_de_list(lrt1, lrt2, file.path(list.outdir, "gbm-insc-ensc-all.csv"))

lrt1.018 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vsiNSC018"])
lrt2.018 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM018vseNSC"])
res.018 <- export_de_list(lrt1.018, lrt2.018, file.path(list.outdir, "gbm-insc-ensc-018.csv"))

lrt1.019 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vsiNSC019"])
lrt2.019 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM019vseNSC"])
res.019 <- export_de_list(lrt1.019, lrt2.019, file.path(list.outdir, "gbm-insc-ensc-019.csv"))

lrt1.031 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vsiNSC031"])
lrt2.031 <- glmLRT(fit.glm, contrast=my.contrasts[, "GBM031vseNSC"])
res.031 <- export_de_list(lrt1.031, lrt2.031, file.path(list.outdir, "gbm-insc-ensc-031.csv"))

run_one_go <- function(ens.all, ens.up, ens.down, p.value=0.05) {
  go <- goana(lapply(list(ens.all, ens.up, ens.down), FUN = function(x) {ens.map[x, 'entrezgene']}))
  colnames(go)[4:9] <- c('N.all', 'N.up', 'N.down', 'P.all', 'P.up', 'P.down')
  go <- go[go$P.all < p.value,]
  go <- go[order(go$P.all),]
  go
}


#' Run GO (and KEGG) analysis
#' This will highlight the pathways that are enriched in the gene lists. Indexing is by Entrez ID
#' Requires GO.db: biocLite("GO.db")
run_go <- function(lrt1, lrt2, fdr=0.05) {
  res <- decompose_de_lists(lrt1, lrt2, fdr = fdr)
  
  x.all <- rbind(res$x1, res$x2)
  ens1.all <- x.all$ensembl
  ens1.up <- x.all[x.all$direction == 'U', 'ensembl']
  ens1.down <- x.all[x.all$direction == 'D', 'ensembl']
  go1 <- run_one_go(ens1.all, ens1.up, ens1.down)

  ens1.only.all <- res$x2$ensembl
  ens1.only.up <- res$x2[res$x2$direction == 'U', 'ensembl']
  ens1.only.down <- res$x2[res$x2$direction == 'D', 'ensembl']
  
  go1.only <- run_one_go(ens1.only.all, ens1.only.up, ens1.only.down)

  y.all <- rbind(res$y1, res$y3)
  ens2.all <- y.all$ensembl
  ens2.up <- y.all[y.all$direction == 'U', 'ensembl']
  ens2.down <- y.all[y.all$direction == 'D', 'ensembl']
  
  go2 <- run_one_go(ens2.all, ens2.up, ens2.down)

  ens2.only.all <- res$y3$ensembl
  ens2.only.up <- res$y3[res$y3$direction == 'U', 'ensembl']
  ens2.only.down <- res$y3[res$y3$direction == 'D', 'ensembl']
  
  go2.only <- run_one_go(ens2.only.all, ens2.only.up, ens2.only.down)

  ens.both.all <- res$x1[res$x1$direction == res$y1$direction, 'ensembl']
  ens.both.up <- res$x1[(res$x1$direction == 'U') & (res$y1$direction == 'U'), 'ensembl']
  ens.both.down <- res$x1[(res$x1$direction == 'D') & (res$y1$direction == 'D'), 'ensembl']
  
  go.both <- run_one_go(ens.both.all, ens.both.up, ens.both.down)
  
  return(list(
    go1=go1,
    go1.only=go1.only,
    go2=go2,
    go2.only=go2.only,
    go.both=go.both
  ))

}

go.018 <- run_go(lrt1.018, lrt2.018, fdr = 0.05)
go.019 <- run_go(lrt1.019, lrt2.019, fdr = 0.05)
go.031 <- run_go(lrt1.031, lrt2.031, fdr = 0.05)

source("http://www.bioconductor.org/biocLite.R")

library(dplyr)
library(DESeq2)
library('biomaRt')
library(calibrate)
library("edgeR")
library("gridExtra")

source('io/rnaseq.R')
source('differential_expression/edger_de.R')
source('io/output.R')
source('_settings.R')
source("utils.R")

fdr <- 0.05
log2FC.min <- 1.

output.file <- getOutputDir('RTKI.de')

REFERENCE <- 'gibco'
LUMPED_LBL <- 'groups.lumped'
GROUP_LBL <- 'groups'

loaded.wtchg <- paired_rtki_data()
dat.wtchg <- loaded.wtchg$data
meta.wtchg <- loaded.wtchg$meta

# add group info
meta.wtchg[, LUMPED_LBL] <- meta.wtchg$type
meta.wtchg[, GROUP_LBL] <- rownames(meta.wtchg)

if (REFERENCE == 'h9') {
  loaded.ref <- duan_nsc_data()
  dat.ref <- loaded.ref$data
  meta.ref <- loaded.ref$meta
  
} else if (REFERENCE == 'gibco') {
  dat.ref <- dat.wtchg[, grep('GIBCO', colnames(dat.wtchg)), drop=F]
  meta.ref <- meta.wtchg[grep('GIBCO', rownames(meta.wtchg)), , drop=F]
}

# mark the reference group (whatever it is)
meta.ref[, LUMPED_LBL] <- 'NSC_control'
meta.ref[, GROUP_LBL] <- 'NSC_control'

# remove Gibco from WTCHG
dat.wtchg <- dat.wtchg[, -grep('GIBCO', colnames(dat.wtchg))]
meta.wtchg <- meta.wtchg[-grep('GIBCO', rownames(meta.wtchg)), ]

n.ref <- ncol(dat.ref)

genes <- intersect(rownames(dat.wtchg), rownames(dat.ref))
dat.wtchg <- dat.wtchg[genes,]
dat.ref <- dat.ref[genes, , drop=F]

dat <- cbind.outer(
  dat.wtchg,
  dat.ref
)

meta <- rbind.outer(meta.wtchg, meta.ref)

#' FILTER
#' The smallest library is ~10mi, the mean lib size is 45mi. 
#' We only keep genes that are expressed at CPM > 1 (i.e. >~5 counts for the avg library) in >=3 samples
#' The exception: any gene that has a single CPM value >= 10 is retained

dat <- filter_genes(dat, cpm.min = 1., nsamples.min = 3, unless.cpm.gte = 10.)

#' Define contrasts of interest
#' Important: these must be aligned in their order as this is used for downstream comparison
contrasts = list(
  GBM.vs.iNSC="(GBM018_P10+GBM018_P12+GBM019_P4+GBM031_P4)/4-(DURA018_NSC_N4_P4+DURA018_NSC_N2_P6+DURA019_NSC_N8C_P2+DURA031_NSC_N44B_P2)/4",
  GBM018.vs.iNSC018="(GBM018_P10+GBM018_P12)/2-(DURA018_NSC_N4_P4+DURA018_NSC_N2_P6)/2",
  GBM019.vs.iNSC019="GBM019_P4-DURA019_NSC_N8C_P2",
  GBM030.vs.iNSC030="GBM030_P5-DURA030_NSC_N16B6_P1",
  GBM031.vs.iNSC031="GBM031_P4-DURA031_NSC_N44B_P2"
)

contrasts.ref <- list(
  GBM.vs.refNSC="(GBM018_P10+GBM018_P12+GBM019_P4+GBM031_P4)/4-NSC_control",
  GBM018.vs.refNSC="(GBM018_P10+GBM018_P12)/2-NSC_control",
  GBM019.vs.refNSC="GBM019_P4-NSC_control",
  GBM030.vs.refNSC="GBM030_P5-NSC_control",
  GBM031.vs.refNSC="GBM031_P4-NSC_control"
)

res <- grouped_analysis(dat, meta[,GROUP_LBL], meta[,LUMPED_LBL], c(contrasts, contrasts.ref), gene.symbols=NULL, output.dir=NULL)

#' Export the individual lists to CSV

for (i in seq(1, length(res))) {
  this.name <- names(res)[i]
  this.de <- prepare_de_table(res[[i]], fdr = fdr, log2FC.min = log2FC.min)
  write.csv(this.de, file.path(output.file, sprintf("%s.csv", this.name)), row.names = F)
}

#' For each contrast, compare with the matching reference contrast
#' For example: GBM018.vs.iNSC018 vs GBM018.vs.refNSC
#' This leads to a Venn set.
#' We also export the results to a CSV

output <- list()

for (i in seq(1, length(contrasts))) {
  c1 <- names(contrasts)[i]
  c2 <- names(contrasts.ref)[i]
  this.name <- paste(c1, c2, sep='-')
  output[[this.name]] <- venn_edger_de_lists(res[[c1]], res[[c2]], fdr = fdr, log2FC.min = log2FC.min)
  export_de_list(output[[this.name]], file.path(output.file, sprintf("%s.csv", this.name)))
}

ens.notinref <- lapply(output, function(x) {x$`10`$ensembl})
venn_ens <- do.call(venn_sets, ens.notinref)
venn_ens$contrasts <- names(contrasts)
venn_diagram.from_blocks(venn_ens, count_func = length)

#' Let's look at the core genes (shared by all, but not in the reference comparisons)
core_ens <- venn_ens$`11111`

#' And the core genes shared by 4/5
almost_core_ens <- as.vector(unlist(
  sapply(
    names(venn_ens)[sapply(names(venn_ens)[1:31], function(x){sum(as.integer(strsplit(x, "")[[1]]))}) == 4],
    function(x) {venn_ens[[x]]}
  )
))

ens.map <- biomart_annotation(index.by='ensembl_gene_id')
writeLines(as.vector(ens.map[core_ens, 'hgnc_symbol']))
writeLines(as.vector(ens.map[c(core_ens, almost_core_ens), 'hgnc_symbol']))

#' #' Run GO (and KEGG) analysis
#' #' This will highlight the pathways that are enriched in the gene lists. Indexing is by Entrez ID
#' #' Requires GO.db: biocLite("GO.db")
#' # go.018.paired <- goana(res.018$lrt.paired, geneid = ens.map[rownames(res.018$lrt.paired), 'entrezgene'])
#' # go.018.ref <- goana(res.018$lrt.ref, geneid = ens.map[rownames(res.018$lrt.ref), 'entrezgene'])
#' # go.019.paired <- goana(res.019$lrt.paired, geneid = ens.map[rownames(res.019$lrt.paired), 'entrezgene'])
#' # go.019.ref <- goana(res.019$lrt.ref, geneid = ens.map[rownames(res.019$lrt.ref), 'entrezgene'])
#' # go.031.paired <- goana(res.031$lrt.paired, geneid = ens.map[rownames(res.031$lrt.paired), 'entrezgene'])
#' # go.031.ref <- goana(res.031$lrt.ref, geneid = ens.map[rownames(res.031$lrt.ref), 'entrezgene'])
#' 
#' # plot venn diagrams showing number of genes overlapping
#' library(VennDiagram)
#' plot_venn <- function(res, obj.title=NULL, plot.direction=NULL, png.file=NULL, pdf.file=NULL) {
#'   #' plot.direction: NULL, 'up' or 'down'. Controls which DE genes are included (NULL means all)
#'   if (is.null(plot.direction)) {
#'     x1 = res$gbm_insc
#'     x2 = res$gbm_ensc
#'   } else if (tolower(plot.direction) == 'up') {
#'     x1 = res$gbm_insc[res$gbm_insc$logFC > 0,]
#'     x2 = res$gbm_ensc[res$gbm_ensc$logFC > 0,]
#'   } else if (tolower(plot.direction) == 'down') {
#'     x1 = res$gbm_insc[res$gbm_insc$logFC < 0,]
#'     x2 = res$gbm_ensc[res$gbm_ensc$logFC < 0,]
#'   }
#'   area1 = nrow(x1)
#'   area2 = nrow(x2)
#'   cross_area = length(intersect(
#'     rownames(x1), rownames(x2)
#'   ))
#' 
#'   if (!is.null(png.file)) {
#'     png(png.file)
#'     plot.new()
#'     venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
#'                                     lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                                     cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'     title(obj.title)
#'     dev.off()
#'   }
#'   if (!is.null(pdf.file)) {
#'     pdf(pdf.file)
#'     plot.new()
#'     venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
#'                                     lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                                     cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'     title(obj.title)
#'     dev.off()
#'   }
#'   
#'   plot.new()
#'   venn.plot <- draw.pairwise.venn(area1, area2, cross_area, category = c("GBM vs paired NSC", "GBM vs reference NSC"),
#'                      lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                      cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'   title(obj.title)
#'   return(venn.plot)
#' }
#' 
#' plot_multivenn <- function(res, png.file=NULL) {
#'   x.1 = res$gbm_insc
#'   x.2 = res$gbm_ensc
#' 
#'   x.up.1 = res$gbm_insc[res$gbm_insc$logFC > 0,]
#'   x.up.2 = res$gbm_ensc[res$gbm_ensc$logFC > 0,]
#' 
#'   x.down.1 = res$gbm_insc[res$gbm_insc$logFC < 0,]
#'   x.down.2 = res$gbm_ensc[res$gbm_ensc$logFC < 0,]
#'     
#'   area.1 = nrow(x.1)
#'   area.2 = nrow(x.2)
#'   cross_area = length(intersect(
#'     rownames(x.1), rownames(x.2)
#'   ))
#'   
#'   area.up.1 <- nrow(x.up.1)
#'   area.up.2 <- nrow(x.up.2)
#'   cross_area.up <- length(intersect(
#'     rownames(x.up.1), rownames(x.up.2)
#'   ))
#'   
#'   area.down.1 <- nrow(x.down.1)
#'   area.down.2 <- nrow(x.down.2)
#'   cross_area.down <- length(intersect(
#'     rownames(x.down.1), rownames(x.down.2)
#'   ))
#'   
#'   if (!is.null(png.file)) {
#'     png(png.file, width=600, height=1200, units='px', res=120)
#'   }
#'   
#'   venn <- draw.pairwise.venn(area.1, area.2, cross_area, category = c("GBM - paired NSC", "GBM - reference NSC"),
#'                              lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                              cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'   venn.up <- draw.pairwise.venn(area.up.1, area.up.2, cross_area.up, category = c("GBM - paired NSC (UP)", "GBM - reference NSC (UP)"),
#'                              lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                              cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'   venn.down <- draw.pairwise.venn(area.down.1, area.down.2, cross_area.down, category = c("GBM - paired NSC (DOWN)", "GBM - reference NSC (DOWN)"),
#'                              lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
#'                              cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = T)
#'   grid.arrange(gTree(children=venn), gTree(children=venn.up), gTree(children=venn.down), nrow=3)
#' 
#'   if (!is.null(png.file)) {
#'     dev.off()
#'   }
#' }
#' 
#' 
#' plot_multivenn(res.018, png.file=file.path(list.outdir, "gbm018.png"))
#' plot_multivenn(res.ip.018, png.file=file.path(list.outdir, "gbm018_ip.png"))
#' 
#' plot_multivenn(res.019, png.file=file.path(list.outdir, "gbm019.png"))
#' plot_multivenn(res.ip.019, png.file=file.path(list.outdir, "gbm019_ip.png"))
#' 
#' plot_multivenn(res.031, png.file=file.path(list.outdir, "gbm031.png"))
#' plot_multivenn(res.ip.031, png.file=file.path(list.outdir, "gbm031_ip.png"))


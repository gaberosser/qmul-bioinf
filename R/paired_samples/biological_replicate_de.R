
source('io/rnaseq.R')
source('plotting/venn.R')
source('differential_expression/edger_de.R')
source('io/output.R')
source('_settings.R')
source("utils.R")

fdr = 0.05
output.dir <- getOutputDir("replicates.de")

loaded.wtchg <- paired_gbm_nsc_data()

samples.rtki <- c(
  'GBM018_P10',
  'GBM018_P12',
  'GBM019_P4',
  'GBM030_P5',
  'GBM031_P4',
  'DURA018_NSC_N4_P4',
  'DURA018_NSC_N2_P6',
  'DURA019_NSC_N8C_P2',
  'DURA030_NSC_N16B6_P1',
  'DURA031_NSC_N44B_P2',
  'GIBCO_NSC_P4'
)
# data <- loaded.wtchg$data
# meta <- loaded.wtchg$meta

dat.wtchg.rtki <- loaded.wtchg$data[,samples.rtki]
meta.wtchg.rtki <- loaded.wtchg$meta[samples.rtki,]

loaded.h9 <- duan_nsc_data(collapse.replicates = T)
dat.h9 <- loaded.h9$data
meta.h9 <- loaded.h9$meta

data <- cbind.outer(
 dat.wtchg.rtki,
 dat.h9
)
meta <- rbind.outer(
 meta.wtchg.rtki,
 meta.h9
)

#' FILTER
#' The smallest library is ~10mi, the mean lib size is 45mi. 
#' We only keep genes that are expressed at CPM > 1 (i.e. >~5 counts for the avg library) in >=3 samples

data <- filter_genes(data)

#' Defining groups
meta$groups.lumped <- as.vector(meta$type)
meta['H9_NSC', 'groups.lumped'] <- 'control'
meta['GIBCO_NSC_P4', 'groups.lumped'] <- 'control'

meta$groups <- rownames(meta)

#' CONTRASTS
#' We want to include every possible combination: 2 x 2 individual comparison, pooled n=2 comparison
#' We'll carry out all comparisons against the general background of the GBM vs iNSC samples
contrasts.rtki = list(
  RTKI.vs.H9="(GBM018_P10+GBM018_P12+GBM019_P4+GBM031_P4)/4-H9_NSC",
  RTKI.vs.GIBCO="(GBM018_P10+GBM018_P12+GBM019_P4+GBM031_P4)/4-GIBCO_NSC_P4",
  RTKI.vs.iNSC="(GBM018_P10+GBM018_P12+GBM019_P4+GBM031_P4)/4-(DURA018_NSC_N4_P4+DURA018_NSC_N2_P6+DURA019_NSC_N8C_P2+DURA031_NSC_N44B_P2)/4",
  GBM018A.vs.iNSC018A="GBM018_P10-DURA018_NSC_N4_P4",
  GBM018B.vs.iNSC018B="GBM018_P12-DURA018_NSC_N2_P6",
  GBM018A.vs.iNSC018B="GBM018_P10-DURA018_NSC_N2_P6",
  GBM018B.vs.iNSC018A="GBM018_P12-DURA018_NSC_N4_P4",
  GBM018.vs.iNSC018="(GBM018_P10+GBM018_P12)/2-(DURA018_NSC_N4_P4+DURA018_NSC_N2_P6)/2",
  GBM018.vs.H9="(GBM018_P10+GBM018_P12)/2-H9_NSC",
  GBM018.vs.GIBCO="(GBM018_P10+GBM018_P12)/2-GIBCO_NSC_P4"
)

res.rtki <- grouped_analysis(data, meta$groups, meta$groups.lumped, contrasts.rtki, output.dir=output.dir)

#' Compute Venn blocks based on overlaps of the 4 single sample pair comparisons
blocks.rtki.all <- do.call(venn_edger_de_lists, c(res.rtki[4:7], list(fdr=fdr)))
blocks.rtki <- do.call(venn_edger_de_lists, c(res.rtki[4:7], list(fdr=fdr, background=res.rtki[["GBM018.vs.GIBCO"]])))
counts.rtki.all <- lapply(blocks.rtki.all, nrow)
counts.rtki <- lapply(blocks.rtki, nrow)

png(file.path(output.dir, "replicate_de_venn_018.all.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.rtki.all, print.mode = c('raw', 'percent'))
dev.off()

png(file.path(output.dir, "replicate_de_venn_018.vsGIBCO.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.rtki, print.mode = c('raw', 'percent'))
dev.off()

blocks.018a_018all <- do.call(venn_edger_de_lists, c(res.rtki[c(4, 8)], list(fdr=fdr)))
png(file.path(output.dir, "replicate_de_venn_018AA.vs.018all.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.018a_018all, print.mode = c('raw', 'percent'))
dev.off()

blocks.018a_018b_018all <- do.call(venn_edger_de_lists, c(res.rtki[c(4, 5, 8)], list(fdr=fdr)))
png(file.path(output.dir, "replicate_de_venn_018AA.vs.018BB.vs.018all.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.018a_018b_018all, print.mode = c('raw', 'percent'))
dev.off()

blocks.018fourway_018all <- do.call(venn_edger_de_lists, c(res.rtki[4:8], list(fdr=fdr)))
png(file.path(output.dir, "replicate_de_venn_018fourway.vs.018all.png"), width=1200, height=600)
venn_diagram.from_blocks(blocks.018fourway_018all, print.mode = c('raw', 'percent'), margin = 0.1)
dev.off()

#TODO: repeat for 044 (mesenchymal, but seems to drift)
samples.mes <- c(
  'GBM026_P3n4',
  'GBM026_P8',
  'GBM044_P4',
  'GBM044_P8',
  'DURA026_NSC_N31D_P5',
  'DURA044_NSC_N8_P2',
  'DURA044_NSC_N17_P3',
  'GIBCO_NSC_P4'
)

dat.wtchg.mes <- loaded.wtchg$data[,samples.mes]
meta.wtchg.mes <- loaded.wtchg$meta[samples.mes,]

data <- cbind.outer(
  dat.wtchg.mes,
  dat.h9
)
meta <- rbind.outer(
  meta.wtchg.mes,
  meta.h9
)

data <- filter_genes(data)

#' Defining groups
meta$groups.lumped <- as.vector(meta$type)
meta['H9_NSC', 'groups.lumped'] <- 'control'
meta['GIBCO_NSC_P4', 'groups.lumped'] <- 'control'

meta$groups <- rownames(meta)

#' CONTRASTS
#' We want to include every possible combination: 2 x 2 individual comparison, pooled n=2 comparison
#' We'll carry out all comparisons against the general background of the GBM vs iNSC samples
contrasts.mes = list(
  MES.vs.H9="(GBM026_P3n4+GBM026_P8+GBM044_P4+GBM044_P8)/4-H9_NSC",
  MES.vs.GIBCO="(GBM026_P3n4+GBM026_P8+GBM044_P4+GBM044_P8)/4-GIBCO_NSC_P4",
  MES.vs.iNSC="(GBM026_P3n4+GBM026_P8+GBM044_P4+GBM044_P8)/4-(DURA026_NSC_N31D_P5+DURA044_NSC_N8_P2+DURA044_NSC_N17_P3)/3",
  GBM044A.vs.iNSC044A="GBM044_P4-DURA044_NSC_N8_P2",
  GBM044B.vs.iNSC044B="GBM044_P8-DURA044_NSC_N17_P3",
  GBM044A.vs.iNSC044B="GBM044_P4-DURA044_NSC_N17_P3",
  GBM044B.vs.iNSC044A="GBM044_P8-DURA044_NSC_N8_P2",
  GBM044.vs.iNSC044="(GBM044_P4+GBM044_P8)/2-(DURA044_NSC_N8_P2+DURA044_NSC_N17_P3)/2",
  GBM044.vs.H9="(GBM044_P4+GBM044_P8)/2-H9_NSC",
  GBM044.vs.GIBCO="(GBM044_P4+GBM044_P8)/2-GIBCO_NSC_P4"
)

res.mes <- grouped_analysis(data, meta$groups, meta$groups.lumped, contrasts.mes, output.dir=output.dir)

#' Compute Venn blocks based on overlaps of the 4 single sample pair comparisons
blocks.mes.all <- do.call(venn_edger_de_lists, c(res.mes[4:7], list(fdr=fdr)))
blocks.mes <- do.call(venn_edger_de_lists, c(res.mes[4:7], list(fdr=fdr, background=res.mes[["GBM044.vs.GIBCO"]])))
counts.mes.all <- lapply(blocks.mes.all, nrow)
counts.mes <- lapply(blocks.mes, nrow)

png(file.path(output.dir, "replicate_de_venn_044.all.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.mes.all, print.mode = c('raw', 'percent'))
dev.off()

png(file.path(output.dir, "replicate_de_venn_044.vsGIBCO.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.mes, print.mode = c('raw', 'percent'))
dev.off()

blocks.044a_044all <- do.call(venn_edger_de_lists, c(res.mes[c(4, 8)], list(fdr=fdr)))
png(file.path(output.dir, "replicate_de_venn_044AA.vs.044all.png"), width=800, height=500)
venn_diagram.from_blocks(blocks.044a_044all, print.mode = c('raw', 'percent'))
dev.off()

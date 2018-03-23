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

# load paired data

in.dir <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    'human',
    'star_alignment'
)

meta.file <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    'sources.csv'
)

samples <- c(
  'GBM018_P10',
  'GBM019_P4',
  'GBM026_P8',
  'GBM031_P4',
  'DURA018_NSC_N2_P6',
  'DURA019_NSC_N8C_P2',
  'DURA026_NSC_N31D_P5',
  'DURA031_NSC_N44B_P2'
)

loaded.wtchg <- star.load_all(in.dir, metafile = meta.file, stranded='r')
dat.wtchg <- loaded.wtchg$data[grep("ENSG", rownames(loaded.wtchg$data)),]
colnames(dat.wtchg) <- rownames(loaded.wtchg$meta)
meta.wtchg <- loaded.wtchg$meta[loaded.wtchg$meta$sample %in% samples,]
dat.wtchg <- dat.wtchg[,rownames(meta.wtchg)]

meta.wtchg$type <- as.factor(replace(as.vector(meta.wtchg$type), meta.wtchg$type == 'iNSC', 'NSC'))
meta.wtchg$pid <- as.factor(c('018', '019', '026', '031', '018', '019', '026', '031'))

y <- DGEList(counts=dat.wtchg, group=meta.wtchg$type)
y <- calcNormFactors(y)
des <- model.matrix(~0 + type, data=meta.wtchg)
y <- estimateDisp(y, design = des)

fit <- glmFit(y)
qlfit <- glmQLFit(y)

treat <- glmTreat(fit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des), lfc=1)
qltreat <- glmTreat(qlfit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des), lfc=1)

lrt <- glmLRT(fit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des))
qllrt <- glmQLFTest(qlfit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des))

# these revert to the expected behaviour (glmLRT and glmQLFTest, as appropriate)
# lrt_treat0 <- glmTreat(fit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des), lfc=0)
# qllrt_treat0 <- glmTreat(qlfit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des), lfc=0)

lrt_ql <- glmLRT(qlfit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des))
# this raises an error:
# ql_lrt <- glmQLFTest(fit, contrast=makeContrasts(c("typeGBM", "typeNSC"), levels=des))

# load TCGA-GBM data

in.dir.tcga = file.path(
  data.dir.raid,
  'rnaseq',
  'tcga_gbm',
  'htseq_count',
  'counts'
)

meta.file.tcga = file.path(
  data.dir.raid,
  'rnaseq',
  'tcga_gbm',
  'sources.csv'
)

loaded.tcga <- htseq.load_all(in.dir.tcga, metafile = meta.file.tcga, file.pattern = '.gz')
dat.tcga <- loaded.tcga$data
meta.tcga <- loaded.tcga$meta

# change subgroup names to match
meta.tcga$subgroup <- replace(as.vector(meta.tcga$subgroup), meta.tcga$subgroup == 'GBM_RTK_I', 'RTK I')
meta.tcga$subgroup <- replace(as.vector(meta.tcga$subgroup), meta.tcga$subgroup == 'GBM_RTK_II', 'RTK II')
meta.tcga$subgroup <- replace(as.vector(meta.tcga$subgroup), meta.tcga$subgroup == 'GBM_MES', 'MES')

# Pollard NSC

in.dir.ip = file.path(
  data.dir.raid,
  'rnaseq',
  'E-MTAB-3867',
  'star_alignment'
)

meta.file.ip = file.path(
  data.dir.raid,
  'rnaseq',
  'E-MTAB-3867',
  'sources.csv'
)

loaded.ip <- star.load_all(in.dir.ip, metafile = meta.file.ip, stranded='u')
dat.ip <- loaded.ip$data[grep("ENSG", rownames(loaded.ip$data)), ]
meta.ip <- loaded.ip$meta

# H9 NSC

in.dir.h9 = file.path(
  data.dir.raid,
  'rnaseq',
  'GSE61794',
  'star_alignment'
)

meta.file.h9 = file.path(
  data.dir.raid,
  'rnaseq',
  'GSE61794',
  'sources.csv'
)

loaded.h9 <- star.load_all(in.dir.h9, metafile = meta.file.h9, stranded='u')
dat.h9 <- loaded.h9$data[grep("ENSG", rownames(loaded.h9$data)), ]
meta.h9 <- loaded.h9$meta

# COMBINE

# find the intersecting genes (required because the TCGA dataset has been aligned by a different method)
# might be worth checking whether any of the discarded genes in WTCHG are relevant (TODO)?
genes <- intersect(rownames(dat.tcga), rownames(dat.wtchg))

# combine all data using only intersecting genes

dat.all <- bind_cols(
  dat.wtchg[genes,],
  dat.ip[genes,],
  dat.h9[genes,],
  dat.tcga[genes,]
)
rownames(dat.all) <- genes

# non-TCGA samples - no need to limit genes as they are all aligned alike
dat.non_tcga <- bind_cols(
  dat.wtchg,
  dat.ip,
  dat.h9
)
rownames(dat.non_tcga) <- rownames(dat.wtchg)

#' FILTER
#' The smallest library is ~10mi, the mean lib size is 45mi. 
#' We only keep genes that are expressed at CPM > 1 (i.e. >~5 counts for the avg library) in >=3 samples

y.all <- DGEList(counts=dat.all)
keep.all <- rowSums(cpm(y.all) > 1) >= 3
print("summary(keep.all)")
summary(keep.all)

# repeat but with the non-TCGA samples

y.non_tcga <- DGEList(counts=dat.non_tcga)
keep.non_tcga <- rowSums(cpm(y.non_tcga) > 1) >= 3
print("summary(keep.non_tcga)")
summary(keep.non_tcga)

#' Outcome:
#' The TCGA dataset retains around 5000 more genes based on this cutoff criterion. This is probably in part down to the increased sample size.
#' It could also be down to a different protocol?
#' Only 23 genes in our gene list are not in TCGA, so we remove those to harmonise the datasets

genes <- intersect(rownames(dat.wtchg)[keep.non_tcga], rownames(dat.tcga))
dat.all <- bind_cols(
  dat.wtchg[genes,],
  dat.ip[genes,],
  dat.h9[genes,],
  dat.tcga[genes,]
)
rownames(dat.all) <- genes

dat.non_tcga <- bind_cols(
  dat.wtchg[genes,],
  dat.ip[genes,],
  dat.h9[genes,]
)
rownames(dat.non_tcga) <- genes

#' edgeR
#' Lump iNSC, eNSC, GBM together and use them to estimate dispersion

groups.non_tcga.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', 4)
)
y.non_tcga <- DGEList(counts=dat.non_tcga, group=groups.non_tcga.lumped)
y.non_tcga <- calcNormFactors(y.non_tcga)

#' MDS looks a bit like PCA - it separates the data in 2D
#' Strangely, this suggests that the Pollard data are very different from the other NSC samples
plotMDS(y.non_tcga)

#' Press on. Define the design matrix and estimate dispersion
design <- model.matrix(~as.factor(groups.non_tcga.lumped))
# this estimates tagwise dispersion
y.non_tcga <- estimateDisp(y.non_tcga, design)
plotBCV(y.non_tcga)
title("BCV all non-TCGA data")

dispersion.trended.non_tcga <- y.non_tcga$trended.dispersion
dispersion.common.non_tcga <- y.non_tcga$common.dispersion

#' Motivated by this, let's try this again but without Pollard
dat.non_tcga <- bind_cols(
  dat.wtchg[genes,],
  dat.h9[genes,]
)
rownames(dat.non_tcga) <- genes

groups.non_tcga.lumped <- c(
  as.vector(meta.wtchg$type),
  rep('eNSC', 2)
)
y.non_tcga <- DGEList(counts=dat.non_tcga, group=groups.non_tcga.lumped)
y.non_tcga <- calcNormFactors(y.non_tcga)

#' Replot MDS. This now looks more sensible?
plotMDS(y.non_tcga)

design <- model.matrix(~as.factor(groups.non_tcga.lumped))
# these methods all estimate just the common dispersion
y.non_tcga.robust <- estimateGLMCommonDisp(y.non_tcga, design=design, method='deviance', robust=T, subset=NULL)
print(y.non_tcga.robust$common.dispersion)
y.non_tcga.notrobust <- estimateGLMCommonDisp(y.non_tcga, design=design, method='deviance', robust=F, subset=NULL)
print(y.non_tcga.notrobust$common.dispersion)
y.non_tcga.default <- estimateGLMCommonDisp(y.non_tcga, design=design)
print(y.non_tcga.default$common.dispersion)

# this estimates tagwise dispersion
y.non_tcga <- estimateDisp(y.non_tcga, design)

dispersion.trended.non_tcga_ip <- y.non_tcga$trended.dispersion
dispersion.common.non_tcga_ip <- y.non_tcga$common.dispersion

plotBCV(y.non_tcga)
title("BCV paired samples and H9 NSC")

# RTK I paired analysis
# start with 018

# dat.018 <- bind_cols(
#   dat.wtchg[genes, c('GBM018', 'DURA018N2_NSC')],
#   dat.h9[genes,]
# )
# groups.018 <- c('GBM018', 'iNSC018', rep('NSC', 2))
# dispersion.common <- dispersion.common.non_tcga_ip
# dispersion.trended <- dispersion.trended.non_tcga_ip
# dispersion.common <- NULL
# dispersion.trended <- NULL
# contrast <- c(1, -1, 0)

dat.018 <- bind_cols(
  dat.wtchg[genes, c('GBM018', 'DURA018N2_NSC')],
  dat.ip[genes,],
  dat.h9[genes,]
)
groups.018 <- c('GBM018', 'iNSC018', rep('NSC', 4))
dispersion.common <- dispersion.common.non_tcga_ip
dispersion.trended <- dispersion.trended.non_tcga_ip
# dispersion.common <- NULL
# dispersion.trended <- NULL
contrast <- c(1, -1, 0)

# dat.018 <- bind_cols(
#   dat.wtchg[genes, c('GBM018', 'DURA018N2_NSC')],
#   dat.ip[genes,][1]
# )
# groups.018 <- c('GBM018', 'iNSC018', rep('NSC', 1))
# dispersion.common <- dispersion.common.non_tcga_ip
# dispersion.trended <- dispersion.trended.non_tcga_ip
# contrast <- c(1, -1, 0)

# dat.018 <- dat.wtchg[genes, c('GBM018', 'DURA018N2_NSC')]
# groups.018 <- c('GBM018', 'iNSC018')
# dispersion.common <- dispersion.common.non_tcga_ip
# dispersion.trended <- dispersion.trended.non_tcga_ip
# dispersion.trended <- NULL
# contrast <- c(1, -1)

rownames(dat.018) <- genes
y.018 <- DGEList(dat.018, group = groups.018)
y.018 <- calcNormFactors(y.018)
design <- model.matrix(~0+groups.018)
if (!is.null(dispersion.common)) {
  # manually add dispersion estimates back in
  y.018$common.dispersion <- dispersion.common
  #' Copying across the trended distribution makes for a messed up plotBCV, since the average log CPM (x) values shift, resulting in a jagged looking line.
  #' This is OK if we are sticking to our guns and stating that the previous estimates of dispersion are 'more correct'.
  y.018$trended.dispersion <- dispersion.trended
} else {
  y.018 <- estimateDisp(y.018, design)
}

#' fit model and look for DE genes
#' The QL F test fails if there are no replicates
#' But works if there is no dispersion?!
qlfit <- glmQLFit(y.018, design, robust=T)
qlf.paired <- glmQLFTest(qlfit, contrast = contrast)
print(
  paste0("QL F test finds ", dim(topTags(qlf.paired, n=Inf, p.value=0.05))[1], " significantly DE genes.")
)

#' The likelihood ratio test always works
fit <- glmFit(y.018, design)
# lrt.paired <- glmLRT(fit, contrast = contrast)
lrt.paired <- glmLRT(fit, contrast=contrast)
print(
  paste0("LR test finds ", dim(topTags(lrt.paired, n=Inf, p.value=0.05))[1], " significantly DE genes.")
)
toptags.lr <- topTags(lrt.paired, n=20, p.value=0.05)

# WTCHG and Pollard NSC

dat.wtchg_ip <- bind_cols(
  dat.wtchg,
  dat.ip
)
rownames(dat.wtchg_ip) <- rownames(dat.wtchg)

meta.wtchg_ip <- data.frame(row.names = c(
  as.vector(meta.wtchg[,'sample']), 
  as.vector(meta.ip[,'sample'])
))

# meta.wtchg_ip[grep(pattern = '018', x = rownames(meta.wtchg_ip)), 'pair'] = 1
# meta.wtchg_ip[grep(pattern = '019', x = rownames(meta.wtchg_ip)), 'pair'] = 2
# meta.wtchg_ip[grep(pattern = '026', x = rownames(meta.wtchg_ip)), 'pair'] = 3
# meta.wtchg_ip[grep(pattern = '031', x = rownames(meta.wtchg_ip)), 'pair'] = 4
# meta.wtchg_ip[grep(pattern = 'Pollard', x = rownames(meta.wtchg_ip)), 'pair'] = -1
# meta.wtchg_ip$pair <- as.factor(meta.wtchg_ip$pair)
# 
# meta.wtchg_ip[grep(pattern = 'GBM', x = rownames(meta.wtchg_ip)), 'type'] = 'GBM'
# meta.wtchg_ip[grep(pattern = 'NSC', x = rownames(meta.wtchg_ip)), 'type'] = 'NSC'
# meta.wtchg_ip$type <- as.factor(meta.wtchg_ip$type)

meta.wtchg_ip[, 'group'] = 'other'
meta.wtchg_ip['GBM018', 'group'] = 'GBM018'
meta.wtchg_ip['DURA018N2_NSC', 'group'] = 'iNSC018'
meta.wtchg_ip[grep(pattern = 'Pollard', x = rownames(meta.wtchg_ip)), 'group'] = 'eNSC'
meta.wtchg_ip[, 'group'] = as.factor(meta.wtchg_ip[, 'group'])



# GBM vs NSC for 3 x RTK I
# y <- DGEList(counts = dat.wtchg_ip, group = meta.wtchg_ip$group)
# y <- calcNormFactors(y)
# design <- model.matrix(~meta.wtchg_ip$group)
# y <- estimateDisp(y, design)
# 
# fit <- glmQLFit(y, design)
# qlf <- glmQLFTest(fit, contrast = c(0, 1, -1, 0))


# 
# meta.all <- data.frame(row.names = c(
#   as.vector(meta.wtchg[,'sample']), 
#   as.vector(meta.ip[,'sample']), 
#   as.vector(meta.h9[,'sample']), 
#   as.vector(meta.tcga[,'sample'])
#   ))
# 
# meta.all$cell_type <- as.factor(
#   c(
#     as.vector(meta.wtchg$type),
#     rep('NSC', 2),
#     rep('NSC', 2),
#     rep('GBM', nrow(meta.tcga))
#   )
# )
# 
# meta.all$study <- as.factor(
#   c(
#     rep(1, nrow(meta.wtchg)),
#     rep(2, nrow(meta.ip)),
#     rep(3, nrow(meta.h9)),
#     rep(4, nrow(meta.tcga))
#   )
# )
# 
# meta.all$subtype <- as.factor(
#   c(
#     as.vector(meta.wtchg$disease_subgroup),
#     rep("None", 4),
#     as.vector(meta.tcga$subgroup)
#   )
# )

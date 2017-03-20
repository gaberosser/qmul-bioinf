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

in.dirs <- c(
  file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    '161219_K00198_0151_BHGYHTBBXX',
    'star_alignment'
  ),
  file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    '161222_K00198_0152_AHGYG3BBXX',
    'star_alignment'
  )
)

meta.files <- c(
  file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    '161219_K00198_0151_BHGYHTBBXX',
    'sources.csv'
  ),
  file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704',
    '161222_K00198_0152_AHGYG3BBXX',
    'sources.csv'
  )
)

samples <- c(
  'GBM018',
  'GBM019',
  'GBM026',
  'GBM031',
  'DURA018N2_NSC',
  'DURA019N8C_NSC',
  'DURA026N31D_NSC',
  'DURA031N44B_NSC'
)

loaded.wtchg <- star.combine_lanes(in.dirs, metafiles = meta.files, stranded='r')
dat.wtchg <- loaded.wtchg$data[grep("ENSG", rownames(loaded.wtchg$data)), samples]
meta.wtchg <- loaded.wtchg$meta[loaded.wtchg$meta$sample %in% samples,]

meta.wtchg$type <- as.factor(replace(as.vector(meta.wtchg$type), meta.wtchg$type == 'iNSC', 'NSC'))

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

# combine

# find the intersecting genes (required because the TCGA dataset has been aligned by a different method)
# might be worth checking whether any of the discarded genes in WTCHG are relevant (TODO)?
genes <- intersect(rownames(dat.tcga), rownames(dat.wtchg))

# WTCHG and Pollard NSC

dat.wtchg_ip <- bind_cols(
  dat.wtchg,
  dat.ip
)

meta.wtchg_ip <- data.frame(row.names = c(
  as.vector(meta.wtchg[,'sample']), 
  as.vector(meta.ip[,'sample'])
))

meta.wtchg_ip[grep(pattern = '018', x = rownames(meta.wtchg_ip)), 'pair'] = 1
meta.wtchg_ip[grep(pattern = '019', x = rownames(meta.wtchg_ip)), 'pair'] = 2
meta.wtchg_ip[grep(pattern = '026', x = rownames(meta.wtchg_ip)), 'pair'] = 3
meta.wtchg_ip[grep(pattern = '031', x = rownames(meta.wtchg_ip)), 'pair'] = 4
meta.wtchg_ip[grep(pattern = 'Pollard', x = rownames(meta.wtchg_ip)), 'pair'] = -1
meta.wtchg_ip$pair <- as.factor(meta.wtchg_ip$pair)

meta.wtchg_ip[grep(pattern = 'GBM', x = rownames(meta.wtchg_ip)), 'type'] = 'GBM'
meta.wtchg_ip[grep(pattern = 'NSC', x = rownames(meta.wtchg_ip)), 'type'] = 'NSC'
meta.wtchg_ip$type <- as.factor(meta.wtchg_ip$type)

# edgeR

# GBM vs NSC for 3 x RTK I
y <- DGEList(counts = dat.wtchg_ip)
y <- calcNormFactors(y)
design <- model.matrix(~meta.wtchg_ip$pair + meta.wtchg_ip$type)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast = c(0, 1, 0, 0, 0, -1))


# dat.all <- bind_cols(
#   dat.wtchg[genes,],
#   dat.ip[genes,],
#   dat.h9[genes,],
#   dat.tcga[genes,]
# )
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

#'
#' Investigating the strange results we see with GBM044 paired sampled comparison
#' (2 x GBM) vs (2 x iNSC)
#' 

fdr <- 0.01

library(edgeR)
source("io/rnaseq.R")
source("differential_expression/edger_de.R")

root.dir <- file.path(
  data.dir.raid, 
  'rnaseq',
  'wtchg_p170218'
)

in.dirs <- sapply(
  c('170509_K00150_0192_BHJKCLBBXX', '170515_K00150_0196_BHJKC5BBXX_lane_2', '170515_K00150_0196_BHJKC5BBXX_lane_3'),
  function(x) file.path(root.dir, x)
)

meta.files <- sapply(
  in.dirs,
  function(x) file.path(x, 'sources.csv')
)

alignment.dirs <- sapply(
  in.dirs,
  function(x) file.path(x, 'human', 'star_alignment')
)

samples <- c(
  'GBM044_P4',
  'GBM044_P8',
  'DURA044_NSC_N8_P2',
  'DURA044_NSC_N17_P3'
)

loaded.wtchg <- star.combine_lanes(alignment.dirs, metafiles = meta.files, stranded='r')
dat.wtchg <- loaded.wtchg$data[grep("ENSG", rownames(loaded.wtchg$data)), samples]
the.data <- filter_genes(dat.wtchg, nsamples.min = 1)
meta.wtchg <- loaded.wtchg$meta[loaded.wtchg$meta$sample %in% samples,]
the.groups <- factor(meta.wtchg[, 'type'])
y <- DGEList(counts=the.data, group=the.groups)
y <- calcNormFactors(y)
design <- model.matrix(~0+as.factor(the.groups))
colnames(design) <- levels(the.groups)

y <- estimateDisp(y, design)
print(sprintf("Common dispersion %f", y$common.dispersion))
contrast <- makeContrasts(c("GBM-iNSC"), levels = design)
fit <- glmFit(y, design)
fit <- glmQLFit(y, design)
et <- exactTest(y)
lrt <- glmTreat(fit, contrast = contrast, lfc=1)
nrow(topTags(lrt, n=Inf, p.value=fdr))
nrow(topTags(et, n=Inf, p.value=fdr))

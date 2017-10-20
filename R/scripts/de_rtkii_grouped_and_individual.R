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
source("differential_expression/edger_de.R")

fdr = 0.01

# load paired data

root.dir <- file.path(
  data.dir.raid, 
  'rnaseq',
  'wtchg_p170503'
)

in.dirs <- sapply(
  c('170929_K00150_0250_BHLGNHBBXX', '171003_K00198_0242_AHLGYVBBXX_1', '171003_K00198_0242_AHLGYVBBXX_2'),
  function(x) file.path(root.dir, x)
)

meta.files <- sapply(
  in.dirs,
  function(x) file.path(x, 'sources.csv')
)

alignment.dirs <- sapply(
  in.dirs,
  function(x) file.path(x, 'star_alignment')
)

samples <- c(
  'GBM050_P7',
  'GBM050_P9',
  'GBM054_P4',
  'GBM054_P6',
  'GBM061_P3',
  'GBM061_P5',
  'DURA050_NSC_N12_P3',
  'DURA050_NSC_N16_P4',
  'DURA054_NSC_N3C_P2',
  'DURA054_NSC_N2E_P1',
  'DURA061_NSC_N4_P2',
  'DURA061_NSC_N6_P4'
)
pids = c('050', '054', '061')

loaded.wtchg <- star.combine_lanes(alignment.dirs, metafiles = meta.files, stranded='r')
dat.wtchg <- loaded.wtchg$data[grep("ENSG", rownames(loaded.wtchg$data)), samples]
meta.wtchg <- loaded.wtchg$meta[loaded.wtchg$meta$sample %in% samples,]

# filter
# dat.wtchg <- filter_genes(dat.wtchg, nsamples.min = 2)

individual.comparison <- list()

# individual comparison
for (pid in pids) {
  idx <- grep(pid, meta.wtchg$sample)
  the.data <- dat.wtchg[, idx]
  the.data <- filter_genes(the.data, nsamples.min = 1)
  the.groups <- meta.wtchg[idx, 'type']
  y <- DGEList(counts=the.data, group=the.groups)
  y <- calcNormFactors(y)
  design <- model.matrix(~0+as.factor(the.groups))
  colnames(design) <- levels(the.groups)
  # this estimates tagwise dispersion
  y <- estimateDisp(y, design)
  print(sprintf("GBM%s has common dispersion %f", pid, y$common.dispersion))
  contrast <- makeContrasts(c("GBM-iNSC"), levels = design)
  # fit <- glmFit(y, design)
  # lrt <- glmLRT(fit, contrast = contrast)
  fit <- glmQLFit(y, design)
  lrt <- glmTreat(fit, contrast = contrast, lfc=1)
  individual.comparison[[pid]] <- topTags(lrt, n=Inf, p.value=fdr)
}

#' there's something strange about 061 P4: the library normalisation constants are much lower than the other samples
#' it also has the highest BCV
#' so, let's look at the ECDFs
ecdf <- list()
cols <- rainbow(nrow(meta.wtchg))
dev.new()
plot(0, 0, ylim=c(0, 1), type='n', xlab='Normalised', ylab='CDF')
for (i in seq(nrow(meta.wtchg))) {
  sname <- as.character(meta.wtchg$sample[i])
  the.data <- dat.wtchg[,sname, drop = F]
  the.data <- filter_genes(the.data, nsamples.min = 1)
  ord <- order(the.data)
  xnorm <- the.data[ord] / sum(the.data)
  # pct <- seq(length(the.data)) / length(the.data)
  cdf <- cumsum(xnorm)
  # ecdf[[sname]] = list(pct=pct, ecdf=cdf)
  # lines(cdf, pct, col=cols[i])
  lines(xnorm, cdf, col=cols[])
}


# group comparison

dat.wtchg.filt <- filter_genes(dat.wtchg, nsamples.min = 2)

the.groups <- as.factor(meta.wtchg[, 'type'])
y <- DGEList(counts=dat.wtchg.filt, group=the.groups)
y <- calcNormFactors(y)
design <- model.matrix(~0+as.factor(the.groups))
colnames(design) <- levels(the.groups)
# this estimates tagwise dispersion
y <- estimateDisp(y, design)
contrast <- makeContrasts(c("GBM-iNSC"), levels = design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = contrast)
group.comparison <- topTags(lrt, n=Inf, p.value=fdr)


# individual comparison with grouped dispersion estimate
#' the group structure here leads to a much lower estimate of common dispersion than in the full paired structure
#' I presume that's because the dispersion 'within group' here is within technical replicates, while the full design matrix
#' has multiple patients within each group?

individual_from_group.comparison <- list()

the.pid <- as.factor(sub('[^0-9]*([0-9]{3})_.*', '\\1', meta.wtchg$sample))
the.combined_groups <- as.factor(paste(meta.wtchg$type, the.pid, sep = '_'))
y <- DGEList(counts=dat.wtchg.filt)
y <- calcNormFactors(y)
design <- model.matrix(~0 + the.combined_groups)
colnames(design) <- levels(the.combined_groups)
# this estimates tagwise dispersion
y <- estimateDisp(y, design)
for (pid in pids) {
  contrast <- makeContrasts(c(sprintf("GBM_%s - iNSC_%s", pid, pid)), levels = design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, contrast = contrast)
  individual_from_group.comparison[[pid]] <- topTags(lrt, n=Inf, p.value=fdr)
}

venns <- list()
library(VennDiagram)

for (pid in pids) {
  a <- individual.comparison[[pid]]
  b <- individual_from_group.comparison[[pid]]
  c <- intersect(rownames(a), rownames(b))
  venns[[pid]] <- list(
    ind.not.grp = setdiff(rownames(a), rownames(b)),
    grp.not.ind = setdiff(rownames(b), rownames(a)),
    ind.and.grp = intersect(rownames(a), rownames(b))
  )
  png(filename=paste0("GBM", pid, "_de.png"), width=800, height=800)
  plot.new()
  draw.pairwise.venn(nrow(a), nrow(b), length(c), category = c("Individual only", "Group dispersion"), margin = .1, fill=c("red", "blue"))
  title(main=sprintf("GBM%s", pid))
  dev.off()
}


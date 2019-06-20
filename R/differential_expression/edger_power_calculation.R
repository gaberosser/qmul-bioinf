library(PROPER)
library(readxl)
source('_settings.R')

star.fn <- file.path(
  hgic.dir,
  'current/input_data/our_data/rnaseq/cell_culture/star_pipeline/',
  'cell_culture_star_counts.xlsx'
)
star_counts <- read_excel(star.fn)
star_counts$`Gene Symbol` <- NULL
star_counts <- data.frame(star_counts)
ix <- star_counts[,1]
star_counts <- star_counts[,-1]
rownames(star_counts) <- ix
star_counts <- star_counts[,grep('GBM|NSC', colnames(star_counts))]

hgic.params <- estParam(as.matrix(star_counts), type=2)



# in theory, we can run estParams() here, but in practice estimateTrendedDisp() breaks (in edgeR) - why??
# This is copied from the underlying PROPER code. It turns out that the trended dispersion is not used anyway?

X = DGEList(counts = star_counts, lib.size = colSums(star_counts))
y = estimateCommonDisp(X)
y = estimateTagwiseDisp(y)

lmean = y$AveLogCPM * log(2) + log(y$pseudo.lib.size/1e+06)
phi.g = y$tagwise.dispersion
hgic.params <- list(seqDepth = y$sample$lib.size, lmean = lmean, lOD = log(phi.g))
hgic.baseline <- log(rowMeans(star_counts))

sim.real <- RNAseq.SimOptions.2grp(ngenes=nrow(X), p.DE=0.05, lOD = hgic.params$lOD, lBaselineExpr = hgic.baseline)

simres = runSims(Nreps = c(2, 3), sim.opts=sim.real, DEmethod="edgeR", nsims=20)
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.05, stratify.by="expr", delta=1., target.by = 'lfc')

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05, lOD="cheung", lBaselineExpr="cheung")
simres.Cheung = runSims(Nreps = c(2, 3), sim.opts=sim.opts.Cheung, DEmethod="edgeR", nsims=20)
powers.Cheung = comparePower(simres.Cheung, alpha.type="fdr", alpha.nominal=0.05, stratify.by="expr", delta=2, target.by = 'lfc')

plotPower(
  comparePower(
    simres.Cheung, 
    alpha.type="fdr", 
    alpha.nominal=0.05, 
    stratify.by="expr", 
    delta=1, 
    target.by = 'lfc', 
    strata.filtered = 1, 
    strata = c(0, 10, 2^(1:7) * 10, Inf), 
    filter.by = 'expr'
  )
)

library(AGDEX)
library("RColorBrewer")

res.dir <- '/home/gabriel/Dropbox/research/qmul/results/mb_agdex/by_gene_ncott100'

grp.a <- c("hu.mb", "C", "D", "SHH")
names(grp.a) <- c("All MB", "Group C", "Group D", "SHH")
grp.b <- c("mo.mb", "sb.chd7", "sb.nochd7")
names(grp.b) <- c("All MB", "SB CHD7 insertion", "SB no CHD7 insertion")

ctrl.a <- "hu.control"
ctrl.b <- "mo.control"

grid <- expand.grid(grp.a, grp.b)
grid.labels <- expand.grid(names(grp.a), names(grp.b))

generateFilename <- function(a, b) {
  s <- "agdex_{a}-{act}_{b}-{bct}.res.csv"
  s <- sub('{a}', a, s, fixed = T)
  s <- sub('{act}', ctrl.a, s, fixed = T)
  s <- sub('{b}', b, s, fixed = T)
  s <- sub('{bct}', ctrl.b, s, fixed = T)
  return(s)
}

filenames <- apply(grid, 1, function (x) generateFilename(x[1], x[2]))

loadResults <- function(fn) {
  ff <- file.path(res.dir, fn)
  res <- read.agdex.result(ff)
  return(res$gwide.agdex.res)
}

res <- lapply(filenames, loadResults)

extractMetric <- function(func) {
  m <- as.numeric(lapply(res, func))
  m <- matrix(m, nrow=length(grp.b), ncol=length(grp.a), byrow=T)
  rownames(m) <- names(grp.b)
  colnames(m) <- names(grp.a)
  return(m)
}

cos.vals <- extractMetric(function(x) x$stat.value[x$stat.name == 'cos'])
cos.pvals.a <- extractMetric(function(x) x$A.pval[x$stat.name == 'cos'])
cos.pvals.b <- extractMetric(function(x) x$B.pval[x$stat.name == 'cos'])
cos.pvals.worst <- extractMetric(function(x)  max(x[x$stat.name == 'cos', c("A.pval", "B.pval")]))

dop.vals <- extractMetric(function(x) x$stat.value[x$stat.name == 'dop'])
dop.pvals.a <- extractMetric(function(x) x$A.pval[x$stat.name == 'dop'])
dop.pvals.b <- extractMetric(function(x) x$B.pval[x$stat.name == 'dop'])
dop.pvals.worst <- extractMetric(function(x)  max(x[x$stat.name == 'dop', c("A.pval", "B.pval")]))

plotLabelledHeatmap <- function(
  data, 
  xlab="Mouse (relative to healthy)", 
  ylab="Human (relative to healthy)", 
  clim=NULL,
  cdirection = 1,
  title=NULL,
  subtitle=NULL,
  save.filestem=NULL) {
  dat <- melt(data)
  colnames(dat) <- c(xlab, ylab, "value")
  
  g = ggplot(dat, aes_q(x = as.name(xlab), y = as.name(ylab))) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(dat$value, 3))) +
    scale_fill_distiller(limits=clim, palette = "Reds", direction = cdirection)
  if (!is.null(title)) {
    g = g + ggtitle(title, subtitle = subtitle) + theme(plot.title = element_text(lineheight=.8, hjust=0.5))
    
  }
  if (!is.null(save.filestem)) {
    ggsave(paste0(out.filestem, '.pdf'), plot=g, width=8, height=6, units="in", dpi=200)
    ggsave(paste0(out.filestem, '.png'), plot=g, width=8, height=6, units="in", dpi=200)
  } else {
    print(g)
  }
  
}

# COS
out.filestem = file.path(res.dir, "cos_values")
plotLabelledHeatmap(cos.vals, title = "AGDEX cos value", save.filestem = out.filestem)

out.filestem = file.path(res.dir, "cos_pvalues_human")
plotLabelledHeatmap(cos.pvals.a, title = "AGDEX cos pvalues", subtitle = "Permuting human group labels", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

out.filestem = file.path(res.dir, "cos_pvalues_mouse")
plotLabelledHeatmap(cos.pvals.b, title = "AGDEX cos pvalues", subtitle = "Permuting mouse group labels", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

out.filestem = file.path(res.dir, "cos_pvalues_worst")
plotLabelledHeatmap(cos.pvals.worst, title = "AGDEX cos pvalues", subtitle = "Worst result from permuting group labels on both species", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

# DOP

out.filestem = file.path(res.dir, "dop_values")
plotLabelledHeatmap(dop.vals, title = "AGDEX dop value", save.filestem = out.filestem)

out.filestem = file.path(res.dir, "dop_pvalues_human")
plotLabelledHeatmap(dop.pvals.a, title = "AGDEX dop pvalues", subtitle = "Permuting human group labels", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

out.filestem = file.path(res.dir, "dop_pvalues_mouse")
plotLabelledHeatmap(dop.pvals.b, title = "AGDEX dop pvalues", subtitle = "Permuting mouse group labels", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

out.filestem = file.path(res.dir, "dop_pvalues_worst")
plotLabelledHeatmap(dop.pvals.worst, title = "AGDEX dop pvalues", subtitle = "Worst result from permuting group labels on both species", cdirection = -1, clim=c(0, 1), save.filestem = out.filestem)

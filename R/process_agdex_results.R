library(AGDEX)
library("RColorBrewer")

res.dir <- '/home/gabriel/Dropbox/research/qmul/results/mb_agdex/by_gene_ncott100_yugene'

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
  m <- lapply(res, func)
  m <- matrix(m, nrow=length(grp.b), ncol=length(grp.a))
  rownames(m) <- grp.b
  colnames(m) <- grp.a
  return(m)
}


cos.vals <- lapply(res, function(x) x$stat.value[x$stat.name == 'cos'])
cos.vals <- matrix(cos.vals, nrow=length(grp.b), ncol=length(grp.a), byrow = T)
rownames(cos.vals) <- grp.b
colnames(cos.vals) <- grp.a

heatmap(cos.vals, Rowv = NA, Colv = NA, col = rev(brewer.pal(11, "RdBu")))

cos.pvals.a <- lapply(res, function(x) max(x[x$stat.name == 'cos', "A.pval"]))


cos.pvals.worst <- lapply(res, function(x) max(x[x$stat.name == 'cos', c("A.pval", "B.pval")]))
cos.pvals.worst <- matrix(cos.pvals.worst, nrow=length(grp.b), ncol=length(grp.a))
rownames(cos.pvals.worst) <- grp.b
colnames(cos.pvals.worst) <- grp.a

heatmap(cos.pvals.worst, Rowv = NA, Colv = NA, col = rev(brewer.pal(11, "RdBu")))


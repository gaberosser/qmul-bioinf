source("io/microarray.R")
library(pamr)


# lst <- gse37382(aggr.by = 'SYMBOL', aggr.method='max')
lst <- gse37418(aggr.by = 'SYMBOL', aggr.method='max')

expr <- lst$expr
meta <- lst$meta

# remove the outlier
meta <- meta[meta$subgroup != 'SHH OUTLIER',]
meta <- meta[meta$subgroup != 'U',]
expr <- expr[, rownames(meta)]

dat = list()
dat$x <- data.matrix(expr)
dat$y <- meta$subgroup
dat$geneid <- rownames(expr)

expr.train <- pamr.train(dat)
expr.results <- pamr.cv(expr.train, dat)

print(expr.results)
pamr.plotcv(expr.results)

thresh <- expr.results$threshold[max(which(expr.results$error == min(expr.results$error)))]

pamr.confusion(expr.results, threshold=thresh)
pamr.plotcvprob(expr.results, dat, threshold=thresh)

# plot of LDA components at threshold
pamr.plotcen(expr.train, dat, threshold=15)

fdr.obj <- pamr.fdr(expr.train, dat)
pamr.plotfdr(fdr.obj)

pamr.listgenes(expr.train, dat, threshold=5.0)

threshold.scale <- pamr.adaptthresh(expr.train)
expr.train.adaptive <- pamr.train(dat, threshold.scale = threshold.scale)
expr.results.adaptive <- pamr.cv(expr.train.adaptive, dat)
pamr.plotcv(expr.results.adaptive)

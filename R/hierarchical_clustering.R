source("io/microarray.R")
library(pvclust)
library(parallel)

n.core <- as.integer(16)
n.boot <- 2000

# load Northcott data: max by Entrez ID
lst <- gse37382(aggr.by = 'SYMBOL', aggr.method='max')
expr <- lst$expr
meta <- lst$meta

gene.std <- sort(apply(expr, 1, sd), decreasing = T)
genes <- names(gene.std[1:5000])

res.37382 <- pvclust(expr[genes,], method.dist="cor", method.hclust="average", nboot=n.boot, parallel = n.core)
# res.37382 <- pvclust(expr, method.dist="cor", method.hclust="average", nboot=n.boot, parallel = n.core)

lst <- gse10327(aggr.by = 'SYMBOL', aggr.method='max')
expr <- lst$expr
meta <- lst$meta

gene.std <- sort(apply(expr, 1, sd), decreasing = T)
genes <- names(gene.std[1:5000])

res.10327 <- pvclust(expr[genes,], method.dist="cor", method.hclust="average", nboot=n.boot, parallel = n.core)
# res.10327 <- pvclust(expr, method.dist="cor", method.hclust="average", nboot=n.boot, parallel = n.core)
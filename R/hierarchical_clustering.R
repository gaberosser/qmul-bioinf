source("io/microarray.R")
library(pvclust)
library(parallel)

# load Northcott data: max by Entrez ID
# lst <- gse37382(aggr.by = 'ENTREZID', aggr.method='max')
lst <- gse10327(aggr.by = 'SYMBOL', aggr.method='max')
expr <- lst$expr
meta <- lst$meta

result <- pvclust(expr, method.dist="cor", method.hclust="average", nboot=100, parallel = as.integer(12))

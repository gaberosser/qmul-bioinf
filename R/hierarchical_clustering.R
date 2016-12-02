source("io/microarray.R")
library(pvclust)
library(parallel)

# load Northcott data: max by Entrez ID
lst <- gse37382(aggr.by = 'ENTREZID', aggr.method='max')
expr <- lst$expr
meta <- lst$meta

cl <- makeCluster(12, type = "FORK")
result <- pvclust(cl = cl, expr, method.dist="cor", method.hclust="average", nboot=100, parallel = T)

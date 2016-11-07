#source("http://www.bioconductor.org/biocLite.R")

library(AGDEX)
library(Biobase)
library(vsn)

hu.expr.raw <- read.csv('hu_yg.txt', row.names = 1, header = 1, check.names = TRUE)
# hu data is only BG corrected, so apply a var. stabilising transform
# v <- vsn2(data.matrix(hu.expr.raw))
# hu.expr.vst <- data.frame(v@hx)

mo.expr.raw <- read.csv('mo_yg.txt', row.names = 1, header = 1, check.names = FALSE)

hu.pdata.raw <- read.csv('hu_pdata.txt', row.names=1, header = 1)

# manually alter the names
rownames(hu.pdata.raw) <- colnames(hu.expr.raw)


mo.pdata.raw <- read.csv('mo_pdata.txt', row.names=1, header = 1)

hu.pdata.meta <- data.frame(labelDescription=c("Case/control", "MB Northcott classification"), row.names=colnames(hu.pdata.raw))
mo.pdata.meta <- data.frame(labelDescription=c("Case/control", "CHD7 insertion status"), row.names=colnames(mo.pdata.raw))

hu.phenoData <- new("AnnotatedDataFrame", data=hu.pdata.raw, varMetadata = hu.pdata.meta)
mo.phenoData <- new("AnnotatedDataFrame", data=mo.pdata.raw, varMetadata = mo.pdata.meta)

# hu.expr.eset <- ExpressionSet(assayData = as.matrix(hu.expr.vst), phenoData = hu.phenoData)
hu.expr.eset <- ExpressionSet(assayData = as.matrix(hu.expr.raw), phenoData = hu.phenoData)
mo.expr.eset <- ExpressionSet(assayData = as.matrix(mo.expr.raw), phenoData = mo.phenoData)

#quick consistency check on sample names: these should be TRUE
all(rownames(pData(hu.expr.eset))==colnames(exprs(hu.expr.eset)))
all(rownames(pData(mo.expr.eset))==colnames(exprs(mo.expr.eset)))

# compute dex
dex.set.hu <- make.dex.set.object(Eset.data = hu.expr.eset, comp.var=2, comp.def = "WNT-hu.control")
# dex.set.hu <- make.dex.set.object(Eset.data = hu.expr.eset, comp.var=1, comp.def = "hu.mb-hu.control")
dex.set.mo <- make.dex.set.object(Eset.data = mo.expr.eset, comp.var=1, comp.def = "mo.mb-mo.control")

probe.map <- read.csv('homolog_mapping.txt', header = 1, row.names = 1)
map.data <- list(
  probe.map = probe.map,
  map.Aprobe.col = 3,
  map.Bprobe.col = 4
)

agdex.res <- agdex(dex.setA = dex.set.hu, dex.setB = dex.set.mo, map.data = map.data)

# write main results table
# write.csv(agdex.res$meta.dex.res, "agdex_hu_mo.csv")

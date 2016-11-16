# source("http://www.bioconductor.org/biocLite.R")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("hugene11sttranscriptcluster.db")

library(AGDEX)
library(Biobase)
library(parallel)
library(mogene10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)
source("io/microarray.R")

data.dir <- '../data/'
out.dir <- '/home/gabriel/Dropbox/research/qmul/results/mb_agdex/'

hu_he <- allen_cerebellum(by.gene = T)
hu.he.expr <- hu_he$expr
hu.he.meta <- hu_he$meta

mo.mb.expr <- dubuc_sb_screen(by.gene = T)
mo.he.expr <- gse54650()
hu_mb <- gse37382()
hu.mb.expr <- hu_mb$expr
hu.mb.meta <- hu_mb$meta

# combine samples
hu.expr <- merge(hu.he.expr, hu.mb.expr, by=0)
rownames(hu.expr) <- hu.expr$Row.names
hu.expr$Row.names <- NULL

mo.expr <- merge(mo.he.expr, mo.mb.expr, by=0)
rownames(mo.expr) <- mo.expr$Row.names
mo.expr$Row.names <- NULL

# define phenodata
hu.he.pdata <- data.frame(
  row.names=colnames(hu.he.expr)
)
hu.he.pdata[,'case'] = 'control'
hu.he.pdata[,'subgroup'] = NA

hu.mb.pdata <- data.frame(
  row.names=colnames(hu.mb.expr)
)
hu.mb.pdata[rownames(hu.mb.meta), 'case'] = 'mb'
hu.mb.pdata[rownames(hu.mb.meta), 'subgroup'] = hu.mb.meta$subgroup

mo.he.pdata <- data.frame(
  row.names=colnames(mo.he.expr)
)
mo.he.pdata[,'case'] = 'control'
mo.he.pdata[,'chd7'] = F

mo.mb.pdata <- data.frame(
  row.names=colnames(mo.mb.expr)
)
mo.mb.pdata[,'case'] = 'mb'
mo.mb.pdata[1:3,'chd7'] = T
mo.mb.pdata[4:8,'chd7'] = F

##### FROM HERE ######

# Not needed, but here is how to apply a var. stabilising transform:
# library(vsn)
# v <- vsn2(data.matrix(hu.expr.raw))
# hu.expr.vst <- data.frame(v@hx)

mo.expr.raw <- read.csv(file.path(data.dir, 'mo_yg.txt.gz'), row.names = 1, header = 1, check.names = FALSE)

hu.pdata.raw <- read.csv(file.path(data.dir, 'hu_pdata.txt'), row.names=1, header = 1)

# manually alter the names
rownames(hu.pdata.raw) <- colnames(hu.expr.raw)


mo.pdata.raw <- read.csv(file.path(data.dir, 'mo_pdata.txt'), row.names=1, header = 1)

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

# SHH vs WT
# dex.set.hu <- make.dex.set.object(Eset.data = hu.expr.eset, comp.var=2, comp.def = "WNT-hu.control")

probe.map <- read.csv(file.path(data.dir, 'homolog_mapping.txt'), header = 1, row.names = 1)
map.data <- list(
  probe.map = probe.map,
  map.Aprobe.col = 3,
  map.Bprobe.col = 4
)

runAgdex <- function(comp.def.hu, comp.var.hu, comp.def.mo='mo.mb-mo.control', comp.var.mo=1) {
  filestem <- file.path(out.dir, paste('agdex', comp.def.hu, comp.def.mo, sep='_'))
  # compute dex
  dex.set.mo <- make.dex.set.object(Eset.data = mo.expr.eset, comp.var=comp.var.mo, comp.def = comp.def.mo)
  dex.set.hu <- make.dex.set.object(Eset.data = hu.expr.eset, comp.var=comp.var.hu, comp.def = comp.def.hu)
  
  # run main AGDEX routine
  agdex.res <- agdex(dex.setA = dex.set.hu, dex.setB = dex.set.mo, map.data = map.data)
  
  # write all results (load again with read.agdex.result)
  write.agdex.result(agdex.res, paste0(filestem, '.res.csv'))

  # scatterplot
  png(filename = paste0(filestem, '.png'), width = 1024, height = 768, res=200)
  agdex.scatterplot(agdex.res)
  dev.off()
  
  pdf(file = paste0(filestem, '.pdf'), width = 8, height = 6)
  agdex.scatterplot(agdex.res)
  dev.off()
}

# runAgdex(comp.def.hu = 'hu.mb-hu.control', comp.var.hu = 1)
cl <- makeCluster(8, type = "FORK")
# subruns = c("SHH-hu.control", 
#             "C-hu.control",
#             "D-hu.control")
# parLapply(
#   cl,
#   subruns,
#   function(x) runAgdex(comp.def.hu = x, comp.var.hu = 2)
# )
subruns = list(
  list(comp.def.hu = "SHH-hu.control", comp.def.mo = "mo.chd7-mo.control"),
  list(comp.def.hu = "C-hu.control", comp.def.mo = "mo.chd7-mo.control"),
  list(comp.def.hu = "D-hu.control", comp.def.mo = "mo.chd7-mo.control"),
  list(comp.def.hu = "SHH-hu.control", comp.def.mo = "mo.no_chd7-mo.control"),
  list(comp.def.hu = "C-hu.control", comp.def.mo = "mo.no_chd7-mo.control"),
  list(comp.def.hu = "D-hu.control", comp.def.mo = "mo.no_chd7-mo.control")
)

parLapply(
  cl,
  subruns,
  function(x) runAgdex(comp.def.hu = x$comp.def.hu, comp.var.hu = 2, comp.def.mo = x$comp.def.mo, comp.var.mo = 2)
)

stopCluster(cl)

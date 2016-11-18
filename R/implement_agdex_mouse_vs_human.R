# source("http://www.bioconductor.org/biocLite.R")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("hugene11sttranscriptcluster.db")
# biocLite("YuGene")

library(AGDEX)
library(YuGene)
library(Biobase)
library(parallel)
library(mogene10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)
source("io/microarray.R")

NCOTT_ONLY = T
APPLY_YUGENE = F

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

# reduce to Northcott genes
if (NCOTT_ONLY) {
  ncott <- readRDS(file.path(data.dir, 'agdex_mb_comparison_rdata', 'ncott.rds'))
  hu.expr <- na.omit(hu.expr[names(ncott),])  
}

# YuGene transform
if (APPLY_YUGENE) {
  hu.expr <- YuGene(hu.expr, progressBar = F)
  mo.expr <- YuGene(mo.expr, progressBar = F)
  
  # this is bizarre but don't know how else to do it!?
  class(hu.expr) <- 'matrix'
  class(mo.expr) <- 'matrix'
}

# define phenodata
hu.he.pdata <- data.frame(
  row.names=colnames(hu.he.expr)
)
hu.he.pdata[,'case'] = 'hu.control'
hu.he.pdata[,'subgroup'] = 'hu.control'

hu.mb.pdata <- data.frame(
  row.names=colnames(hu.mb.expr)
)
hu.mb.pdata[rownames(hu.mb.meta), 'case'] = 'hu.mb'
hu.mb.pdata[rownames(hu.mb.meta), 'subgroup'] = hu.mb.meta$subgroup
# rename for brevity: Group 3 -> C, etc.
hu.mb.pdata$subgroup <- as.character(hu.mb.pdata$subgroup)
hu.mb.pdata$subgroup <- replace(hu.mb.pdata$subgroup, hu.mb.pdata$subgroup == 'Group 3', 'C')
hu.mb.pdata$subgroup <- replace(hu.mb.pdata$subgroup, hu.mb.pdata$subgroup == 'Group 4', 'D')
hu.mb.pdata$subgroup <- as.factor(hu.mb.pdata$subgroup)

mo.he.pdata <- data.frame(
  row.names=colnames(mo.he.expr)
)
mo.he.pdata[,'case'] = 'mo.control'
mo.he.pdata[,'sb.chd7'] = 'mo.control'
mo.he.pdata[,'sb.nochd7'] = 'mo.control'

mo.mb.pdata <- data.frame(
  row.names=colnames(mo.mb.expr)
)
mo.mb.pdata[,'case'] = 'mo.mb'
mo.mb.pdata[1:3,'sb.chd7'] = 'sb.chd7'
mo.mb.pdata[4:8,'sb.nochd7'] = 'sb.nochd7'

# For a comparison of CHD7 insertion vs everything else (SB and healthy):
# mo.mb.pdata[4:8,'sb.chd7'] = 'mo.control'
# mo.mb.pdata[1:3,'sb.nochd7'] = 'mo.control'

# For a comparison of CHD7 insertion vs healthy only:
mo.mb.pdata[4:8,'sb.chd7'] = 'sb.nochd7'
mo.mb.pdata[1:3,'sb.nochd7'] = 'sb.chd7'

# merge phenodata
hu.pdata <- rbind(hu.he.pdata, hu.mb.pdata)
mo.pdata <- rbind(mo.he.pdata, mo.mb.pdata)

# create ExpressionSet objects
hu.phenoData <- new("AnnotatedDataFrame", data=hu.pdata)
mo.phenoData <- new("AnnotatedDataFrame", data=mo.pdata)
hu.expr.eset <- ExpressionSet(assayData = as.matrix(hu.expr), phenoData = hu.phenoData)
mo.expr.eset <- ExpressionSet(assayData = as.matrix(mo.expr), phenoData = mo.phenoData)

# set the mapping between genes
probe.map <- read.csv(file.path(data.dir, 'agdex_mb_comparison_rdata', 'homolog_mapping.txt'), header = 1, row.names = 1)
# use entrez ID to map
map.data <- list(
  probe.map = probe.map,
  map.Aprobe.col = 3,
  map.Bprobe.col = 4
)

runAgdex <- function(
  comp.def.hu, 
  comp.var.hu, 
  comp.def.mo='mo.mb-mo.control', 
  comp.var.mo=1,
  min.nperms = 100,
  max.nperms = 2000,
  save.plots = T,
  save.res = T) {
  
  filestem <- file.path(out.dir, paste('agdex', comp.def.hu, comp.def.mo, sep='_'))
  # compute dex
  dex.set.mo <- make.dex.set.object(Eset.data = mo.expr.eset, comp.var=comp.var.mo, comp.def = comp.def.mo)
  dex.set.hu <- make.dex.set.object(Eset.data = hu.expr.eset, comp.var=comp.var.hu, comp.def = comp.def.hu)
  
  # run main AGDEX routine
  agdex.res <- agdex(
    dex.setA = dex.set.hu, 
    dex.setB = dex.set.mo, 
    map.data = map.data,
    min.nperms = min.nperms,
    max.nperms = max.nperms
    )
  
  # write all results (load again with read.agdex.result)
  if (save.res) {
    write.agdex.result(agdex.res, paste0(filestem, '.res.csv'))
  }
  
  # scatterplot
  if (save.plots) {
    png(filename = paste0(filestem, '.png'), width = 1024, height = 768, res=200)
    agdex.scatterplot(agdex.res)
    dev.off()
    pdf(file = paste0(filestem, '.pdf'), width = 8, height = 6)
    agdex.scatterplot(agdex.res)
    dev.off()
  } else {
    agdex.scatterplot(agdex.res)
  }
  
  return(agdex.res)
  
}

# Not needed, but here is how to apply a var. stabilising transform:
# library(vsn)
# v <- vsn2(data.matrix(hu.expr.raw))
# hu.expr.vst <- data.frame(v@hx)

cl <- makeCluster(4, type = "FORK")

subruns = list(
  list(comp.def.hu = "hu.mb-hu.control", comp.def.mo = "mo.mb-mo.control", comp.var.hu=1, comp.var.mo = 1),
  list(comp.def.hu = "hu.mb-hu.control", comp.def.mo = "sb.chd7-mo.control", comp.var.hu=1, comp.var.mo = 2),
  list(comp.def.hu = "hu.mb-hu.control", comp.def.mo = "sb.nochd7-mo.control", comp.var.hu=1, comp.var.mo = 3),
  
  list(comp.def.hu = "SHH-hu.control", comp.def.mo = "mo.mb-mo.control", comp.var.hu=2, comp.var.mo = 1),
  list(comp.def.hu = "C-hu.control", comp.def.mo = "mo.mb-mo.control", comp.var.hu=2, comp.var.mo = 1),
  list(comp.def.hu = "D-hu.control", comp.def.mo = "mo.mb-mo.control", comp.var.hu=2, comp.var.mo = 1),
  
  list(comp.def.hu = "SHH-hu.control", comp.def.mo = "sb.chd7-mo.control", comp.var.hu=2, comp.var.mo = 2),
  list(comp.def.hu = "C-hu.control", comp.def.mo = "sb.chd7-mo.control", comp.var.hu=2, comp.var.mo = 2),
  list(comp.def.hu = "D-hu.control", comp.def.mo = "sb.chd7-mo.control", comp.var.hu=2, comp.var.mo = 2),
  list(comp.def.hu = "SHH-hu.control", comp.def.mo = "sb.nochd7-mo.control", comp.var.hu=2, comp.var.mo = 3),
  list(comp.def.hu = "C-hu.control", comp.def.mo = "sb.nochd7-mo.control", comp.var.hu=2, comp.var.mo = 3),
  list(comp.def.hu = "D-hu.control", comp.def.mo = "sb.nochd7-mo.control", comp.var.hu=2, comp.var.mo = 3)
)

parLapply(
  cl,
  subruns,
  function(x) {
    do.call(runAgdex, x)
  }
)

stopCluster(cl)


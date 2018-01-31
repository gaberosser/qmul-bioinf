source('io/output.R')
source('_settings.R')
source("utils.R")
library("ChAMP")
library("minfi")
library("wateRmelon")

MfromBeta <- function(beta) {
  log2(beta / (1 - beta))
}


get_idat_basenames <- function(idat.dir) {
  #' Get the basename of all idat files found recursively under the provided directory
  #' The basename is the full path, minus the trailing _Red.idat
  flist <- list.files(path = idat.dir, recursive = T)
  flist <- flist[grep('_Red.idat', flist)]
  
  # convert to basenames for loading
  basenames <- file.path(idat.dir, sub(pattern = "_Red.idat", "", flist))
}


process_and_save <- function(
  idat.dir, 
  meta.file, 
  samples=NULL, 
  output.dir=file.path(idat.dir, '..', 'beta'),
  gzip=T,
  arraytype='EPIC'
  ) {
  # load meta
  meta <- read.csv(meta.file)
  # set the rownames as filenames
  rownames(meta) <- paste(meta$Sentrix_ID, meta$Sentrix_Position, sep = '_')

  basenames <- get_idat_basenames(idat.dir)
  rgSet <- read.metharray(basenames, extended = T)
  
  print(rgSet@annotation)
  
  # rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19")

  mset <- preprocessRaw(rgSet)
  detP <- detectionP(rgSet)

  beta.raw <- getBeta(mset, "Illumina")
  
  # ensure that meta has the same order
  meta <- meta[colnames(beta.raw),]
  
  colnames(beta.raw) <- meta[colnames(beta.raw), 'Sample_Name']
  colnames(detP) <- meta[colnames(detP), 'Sample_Name']
  
  if (!is.null(samples)) {
    keep <- meta$Sample_Name %in% samples
    if (sum(keep) != length(samples)) {
      warning(sprintf("The number of samples supplied (%i) != the number of samples kept (%i)", length(samples), sum(keep)))
    }
    beta.raw <- beta.raw[,keep]
    meta <- meta[keep,]
  }
  
  champLoad <- champ.filter(beta.raw, detP = detP, pd = meta, arraytype = arraytype)
  
  beta.raw <- champLoad$beta

  beta.bmiq <- champ.norm(beta = beta.raw, method = 'BMIQ', arraytype = arraytype, cores=4)

  beta.pbc <- champ.norm(beta = beta.raw, method = 'PBC', arraytype = arraytype)

  mset.swan <- preprocessSWAN(rgSet, mSet = mset)
  beta.swan <- getBeta(mset.swan)
  beta.swan <- beta.swan[rownames(beta.raw),]
  colnames(beta.swan) <- meta[colnames(beta.swan), 'Sample_Name']

  grSet.funnorm <- preprocessFunnorm(rgSet)
  beta.funnorm <- getBeta(grSet.funnorm)[rownames(beta.raw),]

  dir.create(output.dir, showWarnings = FALSE)
  if (gzip) {
    write.csv(beta.raw, file=gzfile(file.path(output.dir, "beta_raw.csv.gz")))
    write.csv(beta.bmiq, file=gzfile(file.path(output.dir, "beta_bmiq.csv.gz")))
    write.csv(beta.swan, file=gzfile(file.path(output.dir, "beta_swan.csv.gz")))
    write.csv(beta.pbc, file=gzfile(file.path(output.dir, "beta_pbc.csv.gz")))
    write.csv(beta.funnorm, file=gzfile(file.path(output.dir, "beta_funnorm.csv.gz")))
  } else {
    write.csv(beta.raw, file = file.path(output.dir, "beta_raw.csv"))
    write.csv(beta.bmiq, file = file.path(output.dir, "beta_bmiq.csv"))
    write.csv(beta.swan, file = file.path(output.dir, "beta_swan.csv"))
    write.csv(beta.pbc, file = file.path(output.dir, "beta_pbc.csv"))
    write.csv(beta.funnorm, file = file.path(output.dir, "beta_funnorm.csv"))
  }
  
}


# base.dir <- file.path(data.dir.raid, 'methylation', '2016-12-19_ucl_genomics')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-05-12')
# base.dir <- file.path(data.dir.raid, 'methylation', 'tcga_gbm')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-08-23')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-09-19')
base.dir <- file.path(data.dir.raid, 'methylation', '2018-01-12')


idat.dir <- file.path(base.dir, 'idat')
meta.file <- file.path(base.dir, 'sources.csv')
# process_and_save(idat.dir, meta.file, arraytype = "450K")
process_and_save(idat.dir, meta.file, arraytype = "EPIC")

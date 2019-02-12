require("minfi")

MfromBeta <- function(beta) {
  log2(beta / (1 - beta))
}


split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))


get_idat_basenames <- function(idat.dir) {
  #' Get the basename of all idat files found recursively under the provided directory
  #' The basename is the full path, minus the trailing _Red.idat
  flist <- list.files(path = idat.dir, recursive = T)
  flist <- flist[grep('_Red.idat', flist)]
  
  # convert to basenames for loading
  basenames <- file.path(idat.dir, sub(pattern = "_Red.idat", "", flist))
}


get_filenames_from_batches <- function(
  in.dirs,
  samples,
  meta_fn='sources.csv',
  idat_subdir='idat'
) {
  in.files <- NULL
  snames <- NULL
  batches <- NULL
  
  for (b in in.dirs) {
    meta <- read.csv(file.path(b, meta_fn))
    # set the rownames as filenames
    rownames(meta) <- paste(meta$Sentrix_ID, meta$Sentrix_Position, sep = '_')
    this_files <- get_idat_basenames(file.path(b, idat_subdir))
    
    # reorder meta
    meta <- meta[basename(this_files),]
    
    # filter meta
    idx <- meta[, 'sample'] %in% samples
    
    # define file and sample names and add to list
    this_files <- this_files[idx]
    this_snames <- as.vector(meta[idx, 'sample'])
    this_batches <- as.vector(sapply(this_files, function(x){split_path(x)[4]}))
    
    in.files <- c(in.files, this_files)
    snames <- c(snames, this_snames)
    batches <- c(batches, this_batches)
  }
  
  list(in.files=in.files, snames=snames, batches=batches)
}


process_idats <- function(
  in.files,
  snames,
  norm.fun=c('swan', 'bmiq', 'funnorm', 'noob', 'quantile', 'pbc', 'raw'),
  arraytype='EPIC',
  force=F
) {
  norm.fun = match.arg(norm.fun)
  rgSet <- read.metharray(in.files, extended = T, force = force)
  colnames(rgSet) <- snames
  
  mset <- preprocessRaw(rgSet)
  detP <- detectionP(rgSet)
  
  # Load beta values (raw), then apply default ChAMP filtering
  beta.raw <- getBeta(mset, "Illumina")
  champLoad <- champ.filter(beta.raw, detP = detP, pd = NULL, arraytype = arraytype)
  beta.raw <- champLoad$beta
  
  if (norm.fun == 'raw') beta <- beta.raw
  
  if (norm.fun == 'swan') {
    mset.swan <- preprocessSWAN(rgSet, mSet = mset)
    beta.swan <- getBeta(mset.swan)
    beta <- beta.swan[rownames(beta.raw),]
  }
  
  if (norm.fun == 'bmiq') {
    beta <- champ.norm(beta = beta.raw, method = 'BMIQ', arraytype = arraytype, cores=4)
  }
  
  if (norm.fun == 'pbc') {
    beta <- champ.norm(beta = beta.raw, method = 'PBC', arraytype = arraytype)
  }
  
  if (norm.fun == 'funnorm') {
    grSet.funnorm <- preprocessFunnorm(rgSet)
    beta <- getBeta(grSet.funnorm)[rownames(beta.raw),]
  }
  
  if (norm.fun == 'quantile') {
    grSet.quantile <- preprocessQuantile(rgSet)
    beta <- getBeta(grSet.quantile)[rownames(beta.raw),]
  }
  
  if (norm.fun == 'noob') {
    mset.noob <- preprocessNoob(rgSet)
    beta <- getBeta(mset.noob)[rownames(beta.raw),]
  }
  
  return(list(beta.raw=beta.raw, beta=beta))
  
}

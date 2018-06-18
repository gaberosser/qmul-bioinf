source('io/output.R')
source('_settings.R')
source("utils.R")
library("ChAMP")
library("minfi")
library("wateRmelon")
library("data.table")

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


process_and_save_idat_ucl <- function(
  idat.dir, 
  meta.file, 
  samples=NULL, 
  output.dir=file.path(idat.dir, '..', 'beta'),
  gzip=T,
  arraytype='EPIC',
  name.col='sample'
  ) {
  # load meta
  meta <- read.csv(meta.file)
  # set the rownames as filenames
  rownames(meta) <- paste(meta$Sentrix_ID, meta$Sentrix_Position, sep = '_')

  basenames <- get_idat_basenames(idat.dir)
  rgSet <- read.metharray(basenames, extended = T)
  
  print(rgSet@annotation)
  dir.create(output.dir, showWarnings = FALSE)
  
  if (gzip) {
    file_ext <- ".csv.gz"
    open_func <- gzfile
  } else {
    file_ext <- ".csv"
    open_func <- file
  }
  
  # rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19")

  mset <- preprocessRaw(rgSet)
  detP <- detectionP(rgSet)

  beta.raw <- getBeta(mset, "Illumina")
  
  # ensure that meta and detP have the same order
  meta <- meta[colnames(beta.raw),]
  detP <- detP[, colnames(beta.raw)]
  
  colnames(beta.raw) <- meta[, name.col]
  colnames(detP) <- meta[, name.col]
  colnames(mset) <- meta[, name.col]
  
  if (!is.null(samples)) {
    keep <- meta[, name.col] %in% samples
    if (sum(keep) != length(samples)) {
      warning(sprintf("The number of samples supplied (%i) != the number of samples kept (%i)", length(samples), sum(keep)))
    }
    beta.raw <- beta.raw[,keep]
    meta <- meta[keep,]
  }
  
  champLoad <- champ.filter(beta.raw, detP = detP, pd = meta, arraytype = arraytype)
  
  beta.raw <- champLoad$beta
  write.csv(beta.raw, file = open_func(file.path(output.dir, paste0("beta_raw", file_ext))))

  beta.bmiq <- champ.norm(beta = beta.raw, method = 'BMIQ', arraytype = arraytype, cores=4)
  write.csv(beta.bmiq, file = open_func(file.path(output.dir, paste0("beta_bmiq", file_ext))))
  
  beta.pbc <- champ.norm(beta = beta.raw, method = 'PBC', arraytype = arraytype)
  write.csv(beta.pbc, file = open_func(file.path(output.dir, paste0("beta_pbc", file_ext))))

  mset.swan <- preprocessSWAN(rgSet, mSet = mset)
  beta.swan <- getBeta(mset.swan)
  beta.swan <- beta.swan[rownames(beta.raw),]
  write.csv(beta.swan, file = open_func(file.path(output.dir, paste0("beta_swan", file_ext))))

  if (arraytype == 'EPIC') {
    grSet.funnorm <- preprocessFunnorm(rgSet)
    beta.funnorm <- getBeta(grSet.funnorm)[rownames(beta.raw),]
    colnames(beta.funnorm) <- meta[colnames(beta.funnorm), name.col]
    write.csv(beta.funnorm, file = open_func(file.path(output.dir, paste0("beta_funnorm", file_ext))))
  }

  # if (gzip) {
  #   write.csv(beta.raw, file=gzfile(file.path(output.dir, "beta_raw.csv.gz")))
  #   write.csv(beta.bmiq, file=gzfile(file.path(output.dir, "beta_bmiq.csv.gz")))
  #   write.csv(beta.swan, file=gzfile(file.path(output.dir, "beta_swan.csv.gz")))
  #   write.csv(beta.pbc, file=gzfile(file.path(output.dir, "beta_pbc.csv.gz")))
  #   if (arraytype == 'EPIC') {
  #     write.csv(beta.funnorm, file =gzfile(file.path(output.dir, "beta_funnorm.csv.gz")))
  #   }
  # } else {
  #   write.csv(beta.raw, file = file.path(output.dir, "beta_raw.csv"))
  #   write.csv(beta.bmiq, file = file.path(output.dir, "beta_bmiq.csv"))
  #   write.csv(beta.swan, file = file.path(output.dir, "beta_swan.csv"))
  #   write.csv(beta.pbc, file = file.path(output.dir, "beta_pbc.csv"))
  #   if (arraytype == 'EPIC') {
  #     write.csv(beta.funnorm, file = file.path(output.dir, "beta_funnorm.csv"))
  #   }
  # }
  
}


process_and_save_methylset <- function(
  mset,
  meta=NULL,
  detP=NULL,
  samples=NULL, 
  output.dir='./',
  gzip=T,
  arraytype='EPIC'
) {
  
  beta.raw <- getBeta(mset, "Illumina")
  
  # ensure that meta has the same order
  meta <- meta[colnames(beta.raw),]
  
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
  
  # grSet.funnorm <- preprocessFunnorm(rgSet)
  # beta.funnorm <- getBeta(grSet.funnorm)[rownames(beta.raw),]
  
  dir.create(output.dir, showWarnings = FALSE)
  if (gzip) {
    write.csv(beta.raw, file=gzfile(file.path(output.dir, "beta_raw.csv.gz")))
    write.csv(beta.bmiq, file=gzfile(file.path(output.dir, "beta_bmiq.csv.gz")))
    # write.csv(beta.swan, file=gzfile(file.path(output.dir, "beta_swan.csv.gz")))
    write.csv(beta.pbc, file=gzfile(file.path(output.dir, "beta_pbc.csv.gz")))
    # write.csv(beta.funnorm, file=gzfile(file.path(output.dir, "beta_funnorm.csv.gz")))
  } else {
    write.csv(beta.raw, file = file.path(output.dir, "beta_raw.csv"))
    write.csv(beta.bmiq, file = file.path(output.dir, "beta_bmiq.csv"))
    # write.csv(beta.swan, file = file.path(output.dir, "beta_swan.csv"))
    write.csv(beta.pbc, file = file.path(output.dir, "beta_pbc.csv"))
    # write.csv(beta.funnorm, file = file.path(output.dir, "beta_funnorm.csv"))
  }
  
}


GenomicMethylSetfromGEORaw <- function(
  filename, 
  pData, 
  array_type='450K',
  sep='\t',
  Uname="Unmethylated Signal",
  Mname="Methylated Signal",
  Pname="Detection Pval",
  row.names=1
  ) {
  colnames <- strsplit(readLines(filename, n = 1), sep)[[1]]
  
  if(all(!grepl(Uname, colnames)))
    stop("No columns contain Uname. Use readLines or look at file header to see column names.")
  
  if(all(!grepl(Mname, colnames)))
    stop("No columns contain Mname. Use readLines or look at file header to see column names.")
  
  if (!is.null(Pname) & all(!grepl(Pname, colnames))) {
    stop("No columns contain Pname. Use readLines or look at file header to see column names.")
  } 
  
  if (is.null(Pname)) {
    select <- sort(c(row.names, grep(Uname,colnames), grep(Mname,colnames)))
  } else {
    select <- sort(c(row.names, grep(Uname,colnames), grep(Mname,colnames), grep(Pname, colnames)))
  }
  
  mat <- fread(filename, header = TRUE, sep = sep, select=select)
  
  rowNames <- as.matrix(mat[,1,with=FALSE])
  mat <- as.matrix(mat[,-1,with=FALSE])
  rownames(mat) <- rowNames
  rm(rowNames)
  
  uindex <- grep(Uname,colnames(mat))
  mindex <- grep(Mname,colnames(mat))
  pindex <- grep(Pname,colnames(mat))
  
  trim <- function (x){
    x<-gsub("^\\s+|\\s+$", "", x)
    x<-gsub("^\\.+|\\.+$", "", x)
    x<-gsub("^\\_+|\\_$", "", x)
    return(x)
  }
  
  UsampleNames <- trim(sub(Uname, "", colnames(mat)[uindex]))
  MsampleNames <- trim(sub(Mname, "", colnames(mat)[mindex]))
  PsampleNames <- trim(sub(Pname, "", colnames(mat)[pindex]))
  
  index <- match(UsampleNames,MsampleNames)
  MsampleNames <- MsampleNames[index]
  mindex <- mindex[index]
  
  if(!identical(UsampleNames,MsampleNames))
    stop("Sample names do not match for Meth and Unmeth channels.")
  
  if (!is.null(Pname)) {
    PsampleNames <- PsampleNames[index]
    pindex <- pindex[index]
    if(!identical(UsampleNames,PsampleNames))
      stop("Sample names do not match for Meth/Unmeth channels and Pdet.")
    Pdat <- mat[,pindex]
  }
  
  if(is.data.frame(pData)) pData <- as(pData,"DataFrame")
  
  if(is.null(pData))  pData <- DataFrame(
    X1=seq_along(UsampleNames),
    row.names=UsampleNames)
  
  if (array_type == "450K") {
    array = "IlluminaHumanMethylation450k"
    annotation = "ilmn12.hg19"
  } else {
    stop(sprintf("No annotations currently defined for array type %s", array_type))
  }
  
  ann <- sprintf("%sanno.%s", array, annotation)
  if(!require(ann, character.only = TRUE))
    stop(sprintf("cannot load annotation package %s", ann))
  object <- get(ann)
  
  gr <- getLocations(object, mergeManifest = F, orderByLocation = TRUE)
  
  locusNames <- names(gr)
  
  ##this might return NAs but it's ok
  ###fix this. return only what is sent
  common <- intersect(locusNames,rownames(mat))
  if(length(common)==0)
    stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")
  ##Note we give no warning if some of the rownmaes have no match.
  
  ind1 <- match(common,rownames(mat))
  ind2 <- match(common,locusNames)
  
  preprocessing <- c(rg.norm=paste0("Data read from file ",filename,"."))
  colnames(mat) <- NULL
  
  gms  <- GenomicMethylSet(gr =  gr[ind2,],
                                    Meth = mat[ind1,mindex],
                                    Unmeth = mat[ind1,uindex],
                                    colData = pData,
                                    preprocessMethod = preprocessing,
                                    annotation = c(array=array,annotation=annotation))
  
  if (!is.null(Pname)) {
    colnames(Pdat) <- rownames(pData)
    Pdat <- Pdat[ind1,]
  }
  
  
  if (!is.null(Pname)) {
    return(
      list(
        GenomicMethylSet=gms,
        Pdat=Pdat
        )
    )
  } else {
    return(gms)
  }
}


# base.dir <- file.path(data.dir.raid, 'methylation', '2016-06-10_brandner')
# base.dir <- file.path(data.dir.raid, 'methylation', '2016-12-19_ucl_genomics')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-05-12')
# base.dir <- file.path(data.dir.raid, 'methylation', 'tcga_gbm')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-08-23')
# base.dir <- file.path(data.dir.raid, 'methylation', '2017-09-19')
# base.dir <- file.path(data.dir.raid, 'methylation', '2018-01-12')
# base.dir <- file.path(data.dir.raid, 'methylation', 'ENCODE_EPIC')
# base.dir <- file.path(data.dir.raid, 'methylation', 'ENCODE_450k')
# base.dir <- file.path(data.dir.raid, 'methylation', '2018-04-09')
# base.dir <- file.path(data.dir.raid, 'methylation', '2018-03-19')
base.dir <- file.path(data.dir.raid, 'methylation', 'E-MTAB-6194')

idat.dir <- file.path(base.dir, 'idat')
# raw.file <- file.path(base.dir, 'geo_raw.txt')
meta.file <- file.path(base.dir, 'sources.csv')
# process_and_save(idat.dir, meta.file, arraytype = "450K")
process_and_save_idat_ucl(idat.dir, meta.file, arraytype = "EPIC", name.col = "sample")
# process_and_save_idat_ucl(idat.dir, meta.file, arraytype = "450k", name.col="sample")

# meta <- read.csv(meta.file)
# rownames(meta) <- meta$Sample_Name
# res <- GenomicMethylSetfromGEORaw(raw.file, meta, Uname = "Signal_A", Mname = "Signal_B", Pname = "Detection Pval", sep='\t')

# process_and_save_methylset(res$GenomicMethylSet, meta=meta, detP=res$Pdat, arraytype = "450K", output.dir = file.path(base.dir, "beta"))

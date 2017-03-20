source('_settings.R')
library('biomaRt')
library('Biobase')
source('io/output.R')

bm_attributes.defaults = c(
  "hgnc_symbol",
  "ensembl_gene_id",
  "entrezgene"
)


get_dated_files <- function(indir) {
  flist <- list.files(indir, pattern='\\.csv')
  dates <- list()
  for (f in flist) {
    theDate <- as.Date(f, "%Y-%m-%d.csv")
    if (!is.na(theDate)) {
      dates[f] <- theDate
    }
  }
  if (length(dates) > 0) {
    theFile <- names(sort(unlist(dates), decreasing = T))[1]
    return(file.path(indir, theFile))
  } else {
    return(NULL)
  }
}


download_annotations_from_biomart <- function(
  mart="ensembl", 
  dataset="hsapiens_gene_ensembl") {
  
  mart <- useDataset(dataset, useMart(mart))
  ens.map <- getBM(attributes=bm_attributes.defaults, mart=mart)
  
  return(ens.map)
  
}


biomart_annotation <- function(  
  mart="ensembl", 
  dataset="hsapiens_gene_ensembl", 
  attributes=bm_attributes.defaults, 
  index.by='ensembl_gene_id',
  force_download=F
) {
  
  outdir <- file.path(out.dir, dataset)
  if (!file.exists(outdir)) {
    print(paste0("Creating output directory ", outdir))
    dir.create(outdir)
  }
  
  # attempt to find a suitable input file
  outfile = get_dated_files(outdir)
  
  if (!is.null(outfile) & file.exists(outfile) & !force_download) {
    print(paste0("Found existing file: ", outfile))
    ens.map <- read.csv(outfile, header=TRUE, row.names=1)
  } else {
    print("Downloading files.")
    ens.map <- download_annotations_from_biomart(
      mart = mart,
      dataset = dataset
    )
    print(paste0("Saving result to ", outfile))
    write.csv(ens.map, file=outfile)
  }
  
  # continue with processing: reindex
  if (!is.null(index.by)) {
    ens.map <- ens.map[isUnique(ens.map[,index.by]),]
  }
  
  rownames(ens.map) <- ens.map[,index.by]
  ens.map[,index.by] <- NULL
  
  return(ens.map)
  
}


annotate_ensembl_data <- function(dat, annotation.col='all', keep.summary=T) {
  # keep.summary: if T (default) then maintain the N_unmapped, etc. rows
  if (is.null(annotation.col)) {
    return(dat)
  }
  
  if (keep.summary) {
    summary <- dat[grep('N_', rownames(dat)),]
  }
  
  ens.map <- biomart_annotation()
  if (annotation.col == 'all') {
    # add additional columns: gene symbol, entrez ID
    for (c in colnames(ens.map)) {
      dat[, c] <- ens.map[rownames(dat), c]
    }
  } else {
    both.idx <- rownames(dat[rownames(dat) %in% rownames(ens.map),])
    new_rownames <- ens.map[both.idx, annotation.col]
    uniq.idx <- isUnique(new_rownames)
    # filter out so that it's only the unique ones
    dat <- dat[both.idx[uniq.idx],]
    new_rownames <- new_rownames[uniq.idx]
    rownames(dat) <- new_rownames
  }
  
  
  if (keep.summary) {
    dat[rownames(summary),] <- summary
  }
  
  return(dat)
}


star.load_one <- function(filename, label=NULL, annotation.col=NULL, stranded='u') {
  if (!(stranded %in% c('u', 'f', 'r'))) {
    stop("Permitted values for stranded are 'u', 'r', 'f'")
  }
  dat = read.csv(filename, sep='\t', row.names = 1, header = F)
  if (stranded == 'u') {
    dat = dat[1]
  } else if (stranded == 'f') {
    dat = dat[2]
  } else {
    dat = dat[3]
  }
  
  if (!is.null(label)) {
    colnames(dat) <- label
  }
  
  dat <- annotate_ensembl_data(dat, annotation.col = annotation.col)
  
}

star.load_all <- function(indir, metafile=NULL, annotation.col=NULL, stranded='u') {
  fname.pattern <- 'ReadsPerGene\\.out\\.tab'
  flist <- list.files(indir, pattern=fname.pattern)
  
  meta <- NULL
  if (!is.null(metafile)) {
    meta <- read.csv(metafile, header = 1, row.names = 1)
    # compare file count in dir with meta
    if (nrow(meta) != length(flist)) {
      warning("Number of meta entries(", nrow(meta), ") is not equal to number of matching files (", length(flist), ")")
    }
  }
  
  # load individual files, deferring annotation to later
  dat <- data.frame(lapply(flist, FUN = function (x) {
    star.load_one(file.path(indir, x), label = gsub(fname.pattern, "", x), annotation.col=NULL, stranded=stranded)
  }))
  
  # apply annotation if requested
  dat <- annotate_ensembl_data(dat, annotation.col = annotation.col)
  
  # no need to test for null
  if ('sample' %in% colnames(meta)) {
    colnames(dat) <- meta[colnames(dat), 'sample']
  }
  
  return(list(data=dat, meta=meta))
  
}

star.combine_lanes <- function(indirs, metafiles=NULL, annotation.col=NULL, stranded='u') {
  # combine multiple lanes
  # TODO: test

  if (!is.null(metafiles) & length(indirs) != length(metafiles)) {
    throw("Length of indirs doesn't match length of metafiles")
  }
  
  a <- star.load_all(indirs[1], metafiles[1], annotation.col = annotation.col, stranded = stranded)
  data <- a$data
  meta <- a$meta
  if (length(indirs) > 1) {
    for (i in 2:length(indirs)) {
      a <- star.load_all(indirs[i], metafiles[i], annotation.col = annotation.col, stranded = stranded)
      data <- data + a$data
      if ('read_count' %in% a$meta) {
        meta[,'read_count'] <- meta[,'read_count'] + a$meta[,'read_count']
      }
    }
  }
  
  return(list(data=data, meta=meta))
  
}


htseq.load_one <- function(infile, label=NULL, annotation.col=NULL) {
  dat <- read.csv(infile, header = F, row.names = 1, sep='\t', check.names = F)
  
  if (!is.null(label)) {
    colnames(dat) <- label
  }
  
  # remove accession version numbers if required
  rownames(dat) <- gsub("\\.[[:digit:]]+$", "", rownames(dat))
  
  # annotate if required
  dat <- annotate_ensembl_data(dat, annotation.col = annotation.col)
  
  return(dat)
}


htseq.load_all <- function(indir, metafile=NULL, annotation.col=NULL, file.pattern=NULL) {
  # assuming that all files in the directory are relevant, unless a pattern has been supplied
  flist <- list.files(indir, pattern=file.pattern)
  
  meta <- NULL
  if (!is.null(metafile)) {
    meta <- read.csv(metafile, header = 1, row.names = 1, check.names = F)
    # compare file count in dir with meta
    if (nrow(meta) != length(flist)) {
      warning("Number of meta entries(", nrow(meta), ") is not equal to number of matching files (", length(flist), ")")
    }
  }
  
  # load individual files, deferring annotation to later
  get_label <- function(x) {
    if (!is.null(file.pattern)) {
      y <- gsub(file.pattern, "", x)
    } else {
      y <- x
    }
    return(y)
  }
  dat <- data.frame(lapply(flist, FUN = function (x) {
    htseq.load_one(file.path(indir, x), label = get_label(x), annotation.col=NULL)
  }), check.names = F)
  
  # apply annotation if requested
  dat <- annotate_ensembl_data(dat, annotation.col = annotation.col)
  
  # no need to test for null
  if ('sample' %in% colnames(meta)) {
    # check whether the data column names match and rename if they do
    if (!any(colnames(dat) %in% rownames(meta))) {
      warning(
        "The data column names ",  
        paste0(paste(colnames(dat)[1:3], collapse = ', '), ', ...'),
        " do not match the meta row names ", 
        paste0(paste(rownames(meta)[1:3], collapse = ', '), ', ...'),
        " so we won't rename.")
    } else {
      colnames(dat) <- meta[colnames(dat), 'sample']
    }
  }
  
  return(list(data=dat, meta=meta))  
  
}

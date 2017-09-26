source('_settings.R')
library('biomaRt')
library('Biobase')
source('io/output.R')
source('utils.R')

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
  
  if (!is.null(outfile) & !force_download) {
    if (file.exists(outfile)) {
      print(paste0("Found existing file: ", outfile))
      ens.map <- read.csv(outfile, header=TRUE, row.names=1)
    }
  } else {
    today <- Sys.Date()
    outfile <- file.path(outdir, format(today, "%Y-%m-%d.csv"))
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
  # load data from multiple lanes, combining where sample names match

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
      if ('read_count' %in% colnames(a$meta) & 'read_count' %in% colnames(meta)) {
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


paired_gbm_nsc_data <- function() {
  samples <- c(
    'GBM018_P10',
    'GBM018_P12',
    'GBM019_P4',
    'GBM026_P8',
    'GBM026_P3n4',
    'GBM030_P5',
    'GBM031_P4',
    'GBM044_P4',
    'GBM044_P8',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'DURA019_NSC_N8C_P2',
    'DURA026_NSC_N31D_P5',
    'DURA030_NSC_N16B6_P1',
    'DURA031_NSC_N44B_P2',
    'DURA044_NSC_N17_P3',
    'DURA044_NSC_N8_P2',
    'GIBCO_NSC_P4'
  )
  dir.1 <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p160704'
  )
  dir.2 <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p170218'
  )
  
  in.dirs <- list()
  meta.files <- list()
  
  in.dirs[[1]] <- c(
    file.path(
      dir.1,
      '161219_K00198_0151_BHGYHTBBXX',
      'star_alignment'
    ),
    file.path(
      dir.1,
      '161222_K00198_0152_AHGYG3BBXX',
      'star_alignment'
    )
  )
  
  in.dirs[[2]] <- c(
    file.path(
      dir.2,
      '170509_K00150_0192_BHJKCLBBXX',
      'human',
      'star_alignment'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_2',
      'human',
      'star_alignment'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_3',
      'human',
      'star_alignment'
    )
  )
  
  meta.files[[1]] <- c(
    file.path(
      dir.1,
      '161219_K00198_0151_BHGYHTBBXX',
      'sources.csv'
    ),
    file.path(
      dir.1,
      '161222_K00198_0152_AHGYG3BBXX',
      'sources.csv'
    )
  )
  
  meta.files[[2]] <- c(
    file.path(
      dir.2,
      '170509_K00150_0192_BHJKCLBBXX',
      'sources.csv'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_2',
      'sources.csv'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_3',
      'sources.csv'
    )
  )
  
  dat <- list()
  meta <- list()
  
  for (n in seq(1, length(in.dirs))) {
    loaded <- star.combine_lanes(in.dirs[[n]], metafiles = meta.files[[n]], stranded='r')
    dat[[n]] = loaded$data[grep("ENSG", rownames(loaded$data)), colnames(loaded$data) %in% samples]
    meta[[n]] <- loaded$meta[loaded$meta$sample %in% samples,]
  }
  
  dat <- do.call(cbind.outer, dat)
  meta <- do.call(rbind.outer, meta)
  
  meta$filename <- rownames(meta)
  rownames(meta) <- meta$sample
  meta$sample <- NULL

  return(list(data=dat, meta=meta))
  
}


paired_rtki_data <- function(include_control = T) {
  samples <- c(
    'GBM018_P10',
    'GBM018_P12',
    'GBM019_P4',
    'GBM030_P5',
    'GBM031_P4',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'DURA019_NSC_N8C_P2',
    'DURA030_NSC_N16B6_P1',
    'DURA031_NSC_N44B_P2'
  )
  if (include_control) {
    samples <- c(samples, 'GIBCO_NSC_P4')
  }
  
  loader <- paired_gbm_nsc_data()
  
  keep <- rownames(loader$meta[rownames(loader$meta) %in% samples,])
  
  data <- loader$data[,keep]
  meta <- loader$meta[keep,]
  
  return(list(data=data, meta=meta))
}

gibco_nsc_data <- function() {
  samples = c('GIBCO_NSC_P4')
  dir.2 <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p170218'
  )
  in.dirs <- c(
    file.path(
      dir.2,
      '170509_K00150_0192_BHJKCLBBXX',
      'human',
      'star_alignment'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_2',
      'human',
      'star_alignment'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_3',
      'human',
      'star_alignment'
    )
  )
  meta.files <- c(
    file.path(
      dir.2,
      '170509_K00150_0192_BHJKCLBBXX',
      'sources.csv'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_2',
      'sources.csv'
    ),
    file.path(
      dir.2,
      '170515_K00150_0196_BHJKC5BBXX_lane_3',
      'sources.csv'
    )
  )
  loaded <- star.combine_lanes(in.dirs, metafiles = meta.files, stranded='r')
  dat = loaded$data[grep("ENSG", rownames(loaded$data)), colnames(loaded$data) %in% samples, drop=F]
  meta <- loaded$meta[loaded$meta$sample %in% samples, , drop=F]
  
  meta$filename <- rownames(meta)
  rownames(meta) <- meta$sample
  meta$sample <- NULL
  
  return(list(data=dat, meta=meta))
}


tcga_gbm_data <- function(modify.subgroup.names=T) {
  in.dir.tcga = file.path(
    data.dir.raid,
    'rnaseq',
    'tcga_gbm',
    'htseq_count',
    'counts'
  )
  
  meta.file.tcga = file.path(
    data.dir.raid,
    'rnaseq',
    'tcga_gbm',
    'sources.csv'
  )
  
  loaded <- htseq.load_all(in.dir.tcga, metafile = meta.file.tcga, file.pattern = '.gz')
  dat <- loaded$data
  meta <- loaded$meta
  
  if (modify.subgroup.names) {
    # change subgroup names to match
    meta$subgroup <- replace(as.vector(meta$subgroup), meta$subgroup == 'GBM_RTK_I', 'RTK I')
    meta$subgroup <- replace(as.vector(meta$subgroup), meta$subgroup == 'GBM_RTK_II', 'RTK II')
    meta$subgroup <- replace(as.vector(meta$subgroup), meta$subgroup == 'GBM_MES', 'MES')    
  }

  return(list(data=dat, meta=meta))
}

pollard_nsc_data <- function() {
  in.dir.ip = file.path(
    data.dir.raid,
    'rnaseq',
    'E-MTAB-3867',
    'star_alignment'
  )
  
  meta.file.ip = file.path(
    data.dir.raid,
    'rnaseq',
    'E-MTAB-3867',
    'sources.csv'
  )
  
  loaded.ip <- star.load_all(in.dir.ip, metafile = meta.file.ip, stranded='u')
  dat.ip <- loaded.ip$data[grep("ENSG", rownames(loaded.ip$data)), ]
  meta.ip <- loaded.ip$meta  
  return(list(data=dat.ip, meta=meta.ip))
}


duan_nsc_data <- function(collapse.replicates = T) {
  in.dir.h9 = file.path(
    data.dir.raid,
    'rnaseq',
    'GSE61794',
    'star_alignment'
  )
  
  meta.file.h9 = file.path(
    data.dir.raid,
    'rnaseq',
    'GSE61794',
    'sources.csv'
  )
  
  loaded.h9 <- star.load_all(in.dir.h9, metafile = meta.file.h9, stranded='u')
  dat.h9 <- loaded.h9$data[grep("ENSG", rownames(loaded.h9$data)), ]
  meta.h9 <- loaded.h9$meta
  
  if (collapse.replicates) {
    vals <- rowSums(dat.h9)
    df.data <- list()
    df.data[['H9_NSC']] = vals
    dat.h9 <- data.frame(df.data, row.names = names(vals), check.names = F)
    meta.h9 <- data.frame(type='NSC', sample='H9_NSC', read_count=sum(meta.h9$read_count), row.names = 'H9_NSC')
  }

  return(list(data=dat.h9, meta=meta.h9))
}


mouse_validation_data <- function() {
  samples.1 = c('eNSC3med', 'eNSC5med', 'eNSC6med')
  dir.1 <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p170506'
  )
  in.dir <- file.path(
      dir.1,
      '170829_K00150_0236_AHL5YHBBXX',
      'mouse',
      'star_alignment'
    )
  meta.file <- file.path(
      dir.1,
      '170829_K00150_0236_AHL5YHBBXX',
      'sources.csv'
    )
  loaded.1 <- star.load_all(in.dir, metafile = meta.file, stranded='r')
  dat.1 = loaded.1$data[grep("ENS", rownames(loaded.1$data)), colnames(loaded.1$data) %in% samples.1, drop=F]
  meta.1 <- loaded.1$meta[loaded.1$meta$sample %in% samples.1, , drop=F]
  
  meta.1$filename <- rownames(meta.1)
  rownames(meta.1) <- meta.1$sample
  meta.1$sample <- NULL
  
  dir.2 <- file.path(
    data.dir.raid, 
    'rnaseq',
    'wtchg_p170390'
  )
  in.dirs <- c(
    file.path(dir.2, '170727_K00198_0222_AHKWW5BBXX', 'mouse', 'star_alignment'),
    file.path(dir.2, '170731_K00150_0226_AHL2CJBBXX_1', 'mouse', 'star_alignment'),
    file.path(dir.2, '170731_K00150_0226_AHL2CJBBXX_2', 'mouse', 'star_alignment')
  )
  meta.files <- c(
    file.path(dir.2, '170727_K00198_0222_AHKWW5BBXX', 'sources.csv'),
    file.path(dir.2, '170731_K00150_0226_AHL2CJBBXX_1', 'sources.csv'),
    file.path(dir.2, '170731_K00150_0226_AHL2CJBBXX_2', 'sources.csv')
  )
  samples.2 = c('eNSC3mouse', 'eNSC5mouse', 'eNSC6mouse', 'mDura3N1mouse', 'mDura5N24Amouse', 'mDura6N6mouse', 'mDura3N1human', 'mDura5N24Ahuman', 'mDura6N6human')
  loaded.2 <- star.combine_lanes(in.dirs, metafiles = meta.files, stranded = 'r')
  dat.2 <- loaded.2$data[grep("ENS", rownames(loaded.2$data)), colnames(loaded.2$data) %in% samples.2, drop=F]
  meta.2 <- loaded.2$meta[loaded.2$meta$sample %in% samples.2, , drop=F]
  
  meta.2$filename <- rownames(meta.2)
  rownames(meta.2) <- meta.2$sample
  meta.2$sample <- NULL
  
  return(list(data=cbind.outer(dat.1, dat.2), meta=rbind.outer(meta.1, meta.2)))
}

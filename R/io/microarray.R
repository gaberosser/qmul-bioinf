# source("http://www.bioconductor.org/biocLite.R")
# biocLite("oligo")
# biocLite("AnnotationDbi")

library(oligo)
library(AnnotationDbi)
library(data.table)
source("utils.R")  # contains median_by and other aggregation routines
source("_settings.R")


eset_from_celdir <- function(cel.dir, gzipped = FALSE) {
  celFiles <- oligoClasses::list.celfiles(cel.dir, listGzipped = gzipped, full.names = TRUE)
  affyRaw <- read.celfiles(celFiles)
}

add_annotation_column <- function(expr_df, annotlib, col.name, multiVals = 'asNA') {
  pids <- row.names(expr_df)
  col.data <- mapIds(annotlib, keys = pids, column = col.name, keytype = 'PROBEID', multiVals = multiVals)
  expr_df[names(col.data), col.name] <- col.data
  return(expr_df)
}


annotate_by_entrez <- function(expr_df, annotlib, aggr.method='median') {

  # annotate
  df <- data.frame(expr_df)
  annot_df <- add_annotation_column(df, annotlib, col.name = 'ENTREZID')
  
  # only keep items with an Entrez ID
  df = df[!is.na(annot_df[['ENTREZID']]),]
  annot_df = annot_df[!is.na(annot_df[['ENTREZID']]),]
  
  if (is.null(aggr.method)) {
    return(annot_df)
  } else {
    if (aggr.method == 'median') {
      aggr_df <- median_by(df, annot_df[['ENTREZID']])
    } else if (aggr.method == 'max') {
      aggr_df <- max_by(df, annot_df[['ENTREZID']])
    }
    # to matrix
    aggr_df <- data.matrix(aggr_df)
    return(aggr_df)
  }
  
}


annotated_expr_from_celdir <- function(cel.dir, annotlib, gzipped = T, strip.title = '.CEL.gz', aggr.by = NULL, aggr.method = 'median') {
  # aggr.method is ignored if aggr.by is not supplied
  
  # TODO: permit different pre-processing steps?
  eset <- rma(eset_from_celdir(cel.dir, gzipped = gzipped))
  
  # extract expression values
  expr_df <- data.frame(exprs(eset))
  
  # optionally rename samples
  if (is.character(strip.title)) {
    # remove portion from sample title
    colnames(expr_df) <- sapply(colnames(expr_df), function(x) {gsub(strip.title, '', x)})    
  }
  
  if (!is.null(aggr.by)) {
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = aggr.by)
    x <- expr_df[, colnames(expr_df) != aggr.by]
    labels <- expr_df[, aggr.by]
    if (aggr.method == 'median') {
      expr_df <- median_by(x, labels)
    } else if (aggr.method == 'mean') {
      expr_df <- mean_by(x, labels)
    } else if (aggr.method == 'min') {
      expr_df <- min_by(x, labels)
    } else if (aggr.method == 'max') {
      expr_df <- max_by(x, labels)
    } else {
      stop("aggr.by is specified but aggr.method is not a recognised option")
    }
    
  } else {
    # add all annotation columns
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'ENTREZID')
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'SYMBOL')
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'ENSEMBL')
  }
  
  return(expr_df)
}


preprocessed_filename <- function(aggr.by = NULL, aggr.method = 'median') {
  if (!is.null(aggr.by)) {
    if (is.null(aggr.method)) {
      stop("aggr.by is specified but aggr.method is NULL.")
    }
    fn = 'expr.rma.{aggr.method}_{aggr.by}.csv.gz'
    fn = sub('{aggr.by}', aggr.by, fn, fixed = T)
    fn = sub('{aggr.method}', aggr.method, fn, fixed = T)
  } else {
    fn = "expr.rma.csv.gz"
  }
  return(fn)
}


allen_cerebellum <- function(by.gene = T) {
  if (by.gene) {
    in.file.expr <- file.path(data.dir, 'allen_human_brain_atlas/microarray', 'cerebellum_expression.by_entrez_id.agg_median.csv.gz')
    missing_msg <- "Unable to find input file {FIN}. Run Python code load_data.allen_human_brain_atlas.save_cerebellum_microarray_data_by_entrez_id()."
  }
  else {
    in.file.expr <- file.path(data.dir, 'allen_human_brain_atlas/microarray', 'cerebellum_expression.csv.gz')
    missing_msg <- "Unable to find input file {FIN}. Run Python code load_data.allen_human_brain_atlas.cerebellum_microarray_reference_data()."
  }
  
  in.file.meta <- file.path(data.dir, 'allen_human_brain_atlas/microarray', 'cerebellum_meta.csv')
  if (!file.exists(in.file.expr)) {
    stop(sub('{FIN}', in.file.expr, missing_msg, fixed = T))
  }
  if (!file.exists(in.file.meta)) {
    msg <- "Unable to find input file {FIN}. Run Python code load_data.allen_human_brain_atlas.cerebellum_microarray_reference_data()."
    stop(sub('{FIN}', in.file.meta, msg, fixed = T))
  }
  
  expr <- read.csv(in.file.expr, header=1, row.names = 1, check.names = F)
  meta <- read.csv(in.file.meta, header=1, row.names = 1)
  
  return(list(
    expr=expr,
    meta=meta
  ))
  
}

dubuc_sb_screen <- function(by.gene = T) {
  library(mogene10sttranscriptcluster.db)
  
  in.file <- file.path(data.dir, 'sleeping_beauty_mouse_screen', 'Dubuc_BMi1_Tg Expression profile.csv')
  arr_data <- read.csv(in.file, sep = ',', skip = 1, header = TRUE, row.names=1)
  arr_data <- arr_data[,1:8]
  sampleNames <- c("Wu050", "Wu053", "Wu054", "Wu051", "Wu052", "Wu055", "Wu056", "Wu057")
  colnames(arr_data) <- sampleNames
  
  if (by.gene) {
    # annotate by Entrez ID and take median over probe sets
    arr_data <- annotate_by_entrez(arr_data, mogene10sttranscriptcluster.db)
  }
  
  return(arr_data)
}


gse54650 <- function() {
  pre_saved.file = file.path(data.dir, 'microarray_GSE54650', 'expr.rma.median_entrez_id.rds')
  if (file.exists(pre_saved.file)) {
    expr <- readRDS(pre_saved.file)
  } else {
    library(mogene10sttranscriptcluster.db)
    in.dir = file.path(data.dir, 'microarray_GSE54650', 'raw')
    expr <- annotated_expr_from_celdir(
      in.dir, 
      annotlib = mogene10sttranscriptcluster.db, 
      gzipped = T, 
      strip.title = '_MoGene1.0ST.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  
  return(expr)
}


load_microarray_data_from_raw <- function(in.dir, annotlib, aggr.by = NULL, aggr.method = 'median', strip.title = '.CEL.gz', gzipped = T) {
  fn <- preprocessed_filename(aggr.by = aggr.by, aggr.method = aggr.method)
  pre_saved.file = file.path(in.dir, fn)
  
  if (file.exists(pre_saved.file)) {
    print(paste0("Loading from existing file ", pre_saved.file))
    expr <- read.csv(pre_saved.file, sep='\t', header=1, row.names = 1)
  } else {
    raw.dir = file.path(in.dir, 'raw')
    print(paste0("Loading from raw files in ", raw.dir))
    expr <- annotated_expr_from_celdir(
      raw.dir,
      annotlib = annotlib, 
      gzipped = gzipped, 
      strip.title = strip.title,
      aggr.by = aggr.by,
      aggr.method = aggr.method)
    print(paste0("Saving to file ", pre_saved.file))
    gz1 <- gzfile(pre_saved.file, "w")
    write.table(expr, file=gz1, sep="\t", col.names = NA, row.names = TRUE)
    close(gz1)
  }
  return(expr)
}


gse37382 <- function(aggr.by = NULL, aggr.method = 'median') {
  library(hugene11sttranscriptcluster.db)
  in.dir <- file.path(data.dir.raid, 'GSE37382')
  
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = hugene11sttranscriptcluster.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)
  
  # load meta
  meta.file = file.path(data.dir.raid, 'GSE37382', 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  # arrange so they are sorted in the same order
  if (!is.null(aggr.by)) {
    expr <- expr[, rownames(meta)]
  } else {
    expr <- expr[, c(rownames(meta), 'SYMBOL', 'ENTREZID', 'ENSEMBL')]
  }

  return(list(expr=expr, meta=meta))
}


gse10327 <- function(aggr.by = NULL, aggr.method = 'median') {
  in.dir <- file.path(data.dir.raid, 'GSE10327')
  library(hgu133plus2.db)
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = hgu133plus2.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)

  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  return(list(expr=expr, meta=meta))
}


gse12992 <- function() {
  in.dir <- file.path(data.dir.raid, 'GSE12992')
  in.dir.raw <- file.path(in.dir, 'raw')
  pre_saved.file = file.path(in.dir, 'expr.rma.median_entrez_id.rds')
  if (file.exists(pre_saved.file)) {
    expr <- readRDS(pre_saved.file)
  } else {
    library(hgu133plus2.db)
    expr <- annotated_expr_from_celdir(
      in.dir.raw, 
      annotlib = hgu133plus2.db,
      gzipped = T, 
      strip.title = '.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  return(list(expr=expr, meta=meta))
}


gse37418 <- function(aggr.by = NULL, aggr.method = 'median') {
  in.dir <- file.path(data.dir.raid, 'GSE37418')
  library(hgu133plus2.db)
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = hgu133plus2.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)
  
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 2)

  # arrange so they are sorted in the same order
  if (!is.null(aggr.by)) {
    expr <- expr[, rownames(meta)]
  } else {
    expr <- expr[, c(rownames(meta), 'SYMBOL', 'ENTREZID', 'ENSEMBL')]
  }  
  
  return(list(expr=expr, meta=meta))
}


thompson2006 <- function(aggr.by = NULL, aggr.method = 'median') {
  in.dir <- file.path(data.dir.raid, 'microarray', 'thompson2006')
  library(hgu133plus2.db)
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = hgu133plus2.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)
  
  return(list(expr=expr, meta=NULL))
}


gse50161 <- function(aggr.by = NULL, aggr.method = 'median') {
  in.dir <- file.path(data.dir.raid, 'microarray', 'GSE50161')
  library(hgu133plus2.db)
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = hgu133plus2.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)
  
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 6)
  
  # arrange so they are sorted in the same order
  if (!is.null(aggr.by)) {
    expr <- expr[, rownames(meta)]
  } else {
    expr <- expr[, c(rownames(meta), 'SYMBOL', 'ENTREZID', 'ENSEMBL')]
  }  
  
  return(list(expr=expr, meta=meta))
}



gse33199 <- function(aggr.by = NULL, aggr.method = 'median') {
  in.dir <- file.path(data.dir.raid, 'microarray', 'GSE33199')
  library(mouse4302.db)
  expr <- load_microarray_data_from_raw(
    in.dir = in.dir,
    annotlib = mouse4302.db,
    aggr.by = aggr.by,
    aggr.method = aggr.method,
    gzipped = T)
  
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  
  # arrange so they are sorted in the same order
  if (!is.null(aggr.by)) {
    expr <- expr[, rownames(meta)]
  } else {
    expr <- expr[, c(rownames(meta), 'SYMBOL', 'ENTREZID', 'ENSEMBL')]
  }  
  
  return(list(expr=expr, meta=meta))
}


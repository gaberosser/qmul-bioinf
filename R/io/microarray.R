# source("http://www.bioconductor.org/biocLite.R")
# biocLite("oligo")
# biocLite("AnnotationDbi")

library(oligo)
library(AnnotationDbi)
library(data.table)
source("utils.R")  # contains median_by and other aggregation routines
data.dir <- '../data/'
data.dir.raid <- '/media/gabriel/raid1_4tb/data/microarray/'


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
  annot_df <- data.frame(expr_df)
  add_annotation_column(annot_df, annotlib, col.name = 'ENTREZID')
  
  # pids <- row.names(expr_df)
  # entrez <- mapIds(annotlib, keys = pids, column = 'ENTREZID', keytype = 'PROBEID', multiVals = 'asNA')
  # entrez <- entrez[!is.na(entrez)]  # only keep probes that match a single gene
  # symbol <- mapIds(annotlib, keys = pids, column = 'SYMBOL', keytype = 'PROBEID', multiVals = 'asNA')
  # symbol <- symbol[!is.na(symbol)]  # only keep probes that match a single gene
  # 
  # annot_df <- data.frame(expr_df)
  # annot_df[names(entrez), 'ENTREZID'] = entrez
  # annot_df[names(symbol), 'SYMBOL'] = symbol
  
  # only keep items with an Entrez ID
  annot_df = annot_df[!is.na(annot_df[['ENTREZID']]),]
  
  if (is.null(aggr.method)) {
    return(annot_df)
  }
  if (aggr.method == 'median') {
    # median aggregation
    aggr_df <- median_by(annot_df[,1:(length(annot_df) - 2)], annot_df[['ENTREZID']])
    
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
    colnames(expr_df) <- sapply(colnames(expr_df), function(x) {gsub(strip_title, '', x)})    
  }
  
  if (!is.null(aggr.by)) {
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = aggr.by)
    
    if (aggr.method == 'median') {
      x <- expr_df[, colnames(expr_df) != aggr.by]
      labels <- expr_df[, aggr.by]
      expr_df <- median_by(x, labels)
    } else if (aggr.method == 'mean') {
      x <- expr_df[, colnames(expr_df) != aggr.by]
      labels <- expr_df[, aggr.by]
      expr_df <- mean_by(x, labels)
    } else if (aggr.method == 'min') {
      x <- expr_df[, colnames(expr_df) != aggr.by]
      labels <- expr_df[, aggr.by]
      expr_df <- min_by(x, labels)
    } else if (aggr.method == 'max') {
      x <- expr_df[, colnames(expr_df) != aggr.by]
      labels <- expr_df[, aggr.by]
      expr_df <- max_by(x, labels)
    } else {
      stop("aggr.by is specified but aggr.method is not a recognised option")
    }
    
  } else {
    # add all annotation columns
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'ENTREZID')
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'SYMBOL')
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'GENENAME')
    expr_df <- add_annotation_column(expr_df, annotlib, col.name = 'ENSEMBL')
  }
  
  return(expr_df)
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
      strip_title = '_MoGene1.0ST.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  
  return(expr)
}


gse37382 <- function() {
  pre_saved.file = file.path(data.dir, 'microarray_GSE37382', 'expr.rma.median_entrez_id.rds')
  if (file.exists(pre_saved.file)) {
    expr <- readRDS(pre_saved.file)
  } else {
    library(hugene11sttranscriptcluster.db)
    in.dir = file.path(data.dir, 'microarray_GSE37382', 'raw')
    expr <- annotated_expr_from_celdir(
      in.dir, 
      annotlib = hugene11sttranscriptcluster.db, 
      gzipped = T, 
      strip_title = '.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  # load meta
  meta.file = file.path(data.dir, 'microarray_GSE37382', 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  return(list(expr=expr, meta=meta))
}


gse10327 <- function() {
  in.dir <- file.path(data.dir.raid, 'GSE10327')
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
      strip_title = '.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
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
      strip_title = '.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  return(list(expr=expr, meta=meta))
}


thompson2006 <- function() {
  in.dir <- file.path(data.dir.raid, 'thompson2006')
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
      strip_title = '.CEL.gz')
    saveRDS(expr, pre_saved.file)
  }
  # load meta
  meta.file = file.path(in.dir, 'sources.csv')
  meta <- read.csv(meta.file, header=1, row.names = 1)
  rownames(meta) <- sapply(rownames(meta), function(x) sub('_', '.', x))
  
  return(list(expr=expr, meta=meta))
}


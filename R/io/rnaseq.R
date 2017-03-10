source('_settings.R')
library('biomaRt')
library('Biobase')
source('io/output.R')

bm_attributes.defaults = c(
  "hgnc_symbol",
  "ensembl_gene_id",
  "entrezgene"
)


download_annotations_from_biomart <- function(
  mart="ensembl", 
  dataset="hsapiens_gene_ensembl", 
  attributes=bm_attributes.defaults, 
  index.by='ensemble_gene_id') {
  
  outdir = getOutputDir(dataset)
  outfile = file.path(outdir, paste0(format(Sys.Date(), "%Y-%m-%d"), '.txt'))
  mart <- useDataset(dataset, useMart(mart))
  ens.map <- getBM(attributes=attributes, mart=mart)
  if (!is.null(index.by)) {
    ens.map <- ens.map[isUnique(ens.map[index.by]),]
  }
  
  rownames(ens.map) <- ens.map$ensembl_gene_id
  ens.map$ensembl_gene_id <- NULL
}


load_from_star <- function(filename, metafile=NULL, annotation.source=NULL, annotation.dest='all', stranded='u') {
  if (!(stranded %in% c('u', 'f', 'r'))) {
    stop("Permitted values for stranded are 'u', 'r', 'f'")
  }
  dat = read.csv(filename, sep='\t', row.names = 1, col.names = 0)
  if (stranded == 'u') {
    dat = dat[1]
  } else if (stranded == 'f') {
    dat = dat[2]
  } else {
    dat = dat[3]
  }
  
  meta <- NULL
  if (!is.null(metafile)) {
    meta <- read.csv(metafile, header = 1, row.names = 1)
  }
  
  if (!is.null(annotation.dest)) {
    
  }
  
  return list(data=dat, meta=meta)
  
}
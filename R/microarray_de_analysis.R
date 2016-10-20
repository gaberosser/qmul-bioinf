# source("http://www.bioconductor.org/biocLite.R")
# biocLite("GEOquery")
library(Biobase)
library(GEOquery)
library("limma")

dataDir <- '../data/'
microarrayDir <- file.path(dataDir, 'microarray_GSE28192')
microarrayFile <- file.path(microarrayDir, 'microarray.RData')

loadCsv <- function (x) {
  read.csv(x, sep='\t', header = 1, row.names = 1)
}


load_from_raw <- function() {
  # load the sample CSV table
  samples <- read.csv(file.path(microarrayDir, "sources.csv"), stringsAsFactors = FALSE)
  
  # load the expression data
  fileList <- apply(samples[2], 1, function(x) file.path(microarrayDir, paste(x, '.gz', sep='')))
  # load first file
  d <- loadCsv(fileList[1])
  
  expr <- data.frame(row.names = rownames(d))
  pVal <- data.frame(row.names = rownames(d))
  expr[samples[1, 2]] <- d[, 1]
  pVal[samples[1, 2]] <- d[["Detection.Pval"]]
  
  for (i in 2:length(fileList)) {
    d <- loadCsv(fileList[i])
    expr[samples[i, 2]] <- d[, 1]
    pVal[samples[i, 2]] <- d[["Detection.Pval"]]
  }
  
  save(list=c("expr", "pVal", "samples"), file=microarrayFile, compress=TRUE)
  return(
    list("expr" = expr, "pVal" = pVal, "samples" = samples)
  )
}


load(microarrayFile)


healthySampleNames = c(
  'NT1197',
  'NCb1',
  'NCb2',
  'A911105',
  'A508112',
  'A508285'
)

# setup design matrix
es1 <- factor(samples$northcott.classification)
es2 <- factor(sapply(samples[,2], function(x) substr(x, nchar(x) - 1, nchar(x)) == '-R'))

X <- model.matrix(~0 + es1 + es2)
# Does this do what we want??
lm.fit(expr, X)

# gse <- getGEO('GSE28192', destdir='./tmp-download')
gse <- getGEO(filename='tmp-download/GSE28192_series_matrix.txt.gz')
# get the expression matrix
mat <- exprs(gse)
# lookup by name match
mat[1:10, samples[colnames(mat), "name"] == "Pt1299"]

# load the probe set definition library
anno <- getGEO(filename = 'tmp-download/GPL6102.soft')


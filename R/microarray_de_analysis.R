# source("http://www.bioconductor.org/biocLite.R")
# biocLite("GEOquery")
library(Biobase)
library(GEOquery)

# load the sample CSV table
samples <- read.csv("../data/microarray_GSE28192/sources.csv")

# gse <- getGEO('GSE28192', destdir='./tmp-download')
gse <- getGEO(filename='tmp-download/GSE28192_series_matrix.txt.gz')
# get the expression matrix
mat <- exprs(gse)
# lookup by name match
mat[1:10, samples[colnames(mat), "name"] == "Pt1299"]

# load the probe set definition library
anno <- getGEO(filename = 'tmp-download/GPL6102.soft')


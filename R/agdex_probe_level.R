# source("http://bioconductor.org/biocLite.R")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("hugene11sttranscriptcluster.db")
# biocLite("annotationTools")
# biocLite("oligo")

library(mogene10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)
library(oligo)
library(Biobase)

dataDir <- '../data/'

mo_data_mb_file = file.path(dataDir, "sleeping_beauty_mouse_screen/Dubuc_BMi1_Tg Expression profile.csv")
mo_data_he_dir = file.path(dataDir, "microarray_GSE54650", "raw")
# hu_data_file = file.path(dataDir, )

celFiles <- list.celfiles(mo_data_he_dir, listGzipped = TRUE, full.names = TRUE)
affyRaw <- read.celfiles(celFiles)

# should auto-load stmogene definitions
eset <- rma(affyRaw)

mo.he.expr <- data.frame(exprs(eset))
mo.he.pData <- data.frame(row.names=colnames(mo.he.expr))
mo.he.phenoData <- AnnotatedDataFrame(data=mo.he.pData)
mo.he.eset <- ExpressionSet(assayData = as.matrix(mo.he.expr), phenoData = mo.he.phenoData)

# remove .CEL.gz from title
colnames(mo.he.expr) <- sapply(colnames(mo.he.expr), function(x) {gsub('_MoGene1.0ST.CEL.gz', '', x)})

mo.mb.raw <- read.csv(mo_data_mb_file, header=TRUE, skip=1, row.names=1)[, 1:8]

# convert MB version to ExpressionSet

mo.mb.pData <- data.frame(row.names=colnames(mo.mb.raw))
mo.mb.phenoData <- AnnotatedDataFrame(data=mo.mb.pData)
mo.mb.eset <- ExpressionSet(assayData = as.matrix(mo.mb.raw), phenoData = mo.mb.phenoData)


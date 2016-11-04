source("http://bioconductor.org/biocLite.R")
biocLite("oligo")
biocLite("limma")
biocLite("mogene10sttranscriptcluster.db")

library(oligo)
library(limma)
library(mogene10sttranscriptcluster.db)

dataDir <- '../data/'
celDir <- file.path(dataDir, 'microarray_GSE54650', 'raw')
celFiles <- list.celfiles(celDir, listGzipped = TRUE, full.names = TRUE)
affyRaw <- read.celfiles(celFiles)

# should auto-load stmogene definitions
eset <- rma(affyRaw)

# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="data.txt")

# annotate
my_frame <- data.frame(exprs(eset))

# Put annotation information in a data frame.  To get specific fields, use packageNameSYMBOL, where the caps part names the type of data you're after
# To get a list of available annotation information, run the packagename with () at the end, i.e. mogene20sttranscriptcluster()
Annot <- data.frame(
  ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), 
  SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), 
  DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "),
  ENSEMBL=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=", "),
  ENTREZ=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", ")
  )

# Merge data frames together (like a database table join)
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)

# remove probesets that have no listed target gene
all = all[(all$ENTREZ != 'NA'),]

# Write out to a file:
write.table(all,file="data.ann.txt",sep="\t")

# Optionally generate a histogram
hist(as.matrix(my_frame), 100)

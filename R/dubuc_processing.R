
library(mogene10sttranscriptcluster.db)
dataDir <- '../data/'
inFile <- file.path(dataDir, 'sleeping_beauty_mouse_screen', 'Dubuc_BMi1_Tg Expression profile.csv')

arr_data <- read.csv(inFile, sep = ',', skip = 1, header = TRUE)
sampleNames <- c("Wu050", "Wu053", "Wu054", "Wu051", "Wu052", "Wu055", "Wu056", "Wu057")
colnames(arr_data) <- c("probeset_id", "Wu050", "Wu053", "Wu054", "Wu051", "Wu052", "Wu055", "Wu056", "Wu057", "gene_symbol", "mrna_accession")
rownames(arr_data) <- arr_data$probeset_id
arr_data$probeset_id <- NULL

# Put annotation information in a data frame.  To get specific fields, use packageNameSYMBOL, where the caps part names the type of data you're after
# To get a list of available annotation information, run the packagename with () at the end, i.e. mogene20sttranscriptcluster()
Annot <- data.frame(
  ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), 
  SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), 
  DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "),
  ENSEMBL=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=", "),
  ENTREZ=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", ")
)

all <- merge(Annot, arr_data, by.x=0, by.y=0, all=T)

# useful to check the gene_symbol and mrna_accession columns - verify they are the same (they are)
# then drop them
all$gene_symbol <- NULL
all$mrna_accession <- NULL

# remove probesets that have no listed target gene
all = all[(all$ENTREZ != 'NA'),]

# Write out to a file:
write.table(all,file="data.ann.txt",sep="\t")

# Optionally generate a histogram
hist(as.matrix(arr_data[, sampleNames]), 100)

# write annotation table on its own, too
annot <- data.frame(
  accession_id=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), 
  gene_symbol=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), 
  gene_description=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "),
  ensembl_id=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=", "),
  entrez_id=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", ")
)

# remove rows that are all NA
annot <- annot[!apply(annot, 1, function(x){all(x == "NA")}),]

write.table(data.frame("probeset_id"=rownames(annot), annot), file="mogene10sttranscriptcluster.tsv",sep="\t", row.names = FALSE)

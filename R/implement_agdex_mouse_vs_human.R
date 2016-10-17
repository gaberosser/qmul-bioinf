#source("http://www.bioconductor.org/biocLite.R")
#biocLite("pd.mogene.1.1.st.v1")
#biocLite("pd.hugene.1.1.st.v1")
#library("pd.mogene.1.1.st.v1")
#library("pd.hugene.1.1.st.v1")

library(annotationTools)
library(biomaRt)

# *Should* be able to automate the process of finding homologous probesets in the two microarrays
# using the code found here: https://mikedewar.wordpress.com/2010/05/14/generating-homologues-using-biomart/
# Unfortunately, the probeset is not defined ni biomaRt.

mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)
martAttrs <- listAttributes(mart)
martAttrs[grep('affy', martAttrs$name),]  # see - doesn't have our microarray

# gen_hs2mm <- function(affyids){
#   ensembl_hs <- useMart(
#     "ensembl",
#     dataset = "hsapiens_gene_ensembl"
#   )
#   hs2mm_filters <- c(
#     "affy_hugene_1_1_st_v1",
#     "with_mmusculus_homolog"
#   )
#   hs2mm_gene_atts <- c(
#     "affy_hugene_1_1_st_v1",  ## NUH-UH
#     "ensembl_gene_id"
#   )
#   hs2mm_homo_atts <- c(
#     "ensembl_gene_id",
#     "mouse_ensembl_gene"
#   )
#   # the names in these lists are arbitrary
#   hs2mm_value = list(
#     affyid=affyids,
#     with_homolog=TRUE
#   )
#   # get the human genes and mouse orthologues
#   hs2mm_gene <- getBM(
#     attributes = hs2mm_gene_atts,
#     filters = hs2mm_filters,
#     value = hs2mm_value,
#     mart = ensembl_hs
#   )
#   
#   hs2mm_homo <- getBM(
#     attributes = hs2mm_homo_atts,
#     filters = hs2mm_filters,
#     value = hs2mm_value,
#     mart = ensembl_hs
#   )
#   # merge the two lists!
#   hs2mm <- merge(hs2mm_gene,hs2mm_homo)
# }
# 
# gen_mm2hs <- function(affyids){
#   ensembl_mm <- useMart("ensembl",
#                         dataset = "mmusculus_gene_ensembl")
#   mm2hs_filters <- c(
#     "affy_mogene_1_1_st_v1",  ## NUH-UH
#     "with_hsapiens_homolog"
#   )
#   mm2hs_gene_atts <- c(
#     "affy_mogene_1_0_st_v1",
#     "ensembl_gene_id"
#   )
#   mm2hs_homo_atts <- c(
#     "ensembl_gene_id",
#     "human_ensembl_gene"
#   )
#   # the names in these lists are arbitrary
#   mm2hs_value = list(
#     affyids=affyids,
#     with_homolog=TRUE
#   )
#   # get the mouse genes and human orthologues
#   mm2hs_gene <- getBM(
#     attributes = mm2hs_gene_atts ,
#     filters = mm2hs_filters,
#     value = mm2hs_value,
#     mart = ensembl_mm
#   )
#   mm2hs_homo <- getBM(
#     attributes = mm2hs_homo_atts,
#     filters = mm2hs_filters,
#     value = mm2hs_value,
#     mart = ensembl_mm
#   )
#   mm2hs <- merge(mm2hs_gene,mm2hs_homo)
#   
# }

muID = 10090
huID = 9606

dataDir <- "../data/sleeping_beauty_mouse_screen/"
#annoDir <- "/home/gabriel/Dropbox/research/qmul/data/microarray_annotation/"
annoDir <- "../## NUH-UHdata/microarray_annotation/"

homologeneFile <- file.path("../data/homologene", "homologene.data")
moAnnoFile <- file.path(annoDir, "MoGene-1_1-st-v1.na36.mm10.transcript.csv.gz")
huAnnoFile <- file.path(annoDir, "HuGene-1_1-st-v1.na36.hg19.transcript.csv.gz")
exprFile <- file.path(dataDir, "Dubuc_BMi1_Tg Expression profile.csv")

# huAnno <- read.csv(huAnnoFile, colClasses='character', comment.char='#')
# moAnno <- read.csv(moAnnoFile, colClasses='character', comment.char='#')
# homologene <- read.delim(homologeneFile, header=FALSE)

# load raw expression data
expr <- read.csv(exprFile, header = TRUE, skip = 1, row.names = 1)

# assay data: the raw expression counts, as a matrix F (# probe sets) x S (# samples)
assayData <- as.matrix(expr[setdiff(colnames(expr), c('Gene.Symbol', 'mRNA.Accession'))])

# pheno data: the conditions, as a table S x V (# covariates)
chd7Status <- t(read.csv(exprFile, nrows = 1, header = FALSE))[2:9]
chd7Status <- data.frame(chd7Status, row.names = colnames(assayData))
pData <- table(rownames(chd7Status), chd7Status$chd7Status)

# orthoTable<-ps2ps(huAnno, moAnno, homologene, 10090)  # 10090: mus musculus


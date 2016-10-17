#source("http://www.bioconductor.org/biocLite.R")
#biocLite("pd.mogene.1.1.st.v1")
#biocLite("pd.hugene.1.1.st.v1")
#library("pd.mogene.1.1.st.v1")
#library("pd.hugene.1.1.st.v1")

library(annotationTools)

musID = 10090

dataDir <- "/home/gabriel/Dropbox/research/qmul/data/sleeping_beauty_mouse_screen/"
annoDir <- "/home/gabriel/Dropbox/research/qmul/data/microarray_annotation/"

homologeneFile <- file.path(dataDir, "homologene.data")
moAnnoFile <- file.path(annoDir, "MoGene-1_1-st-v1.na36.mm10.probeset.csv.gz")
huAnnoFile <- file.path(annoDir, "HuGene-1_1-st-v1.na36.hg19.probeset.csv.gz")

huAnno <- read.csv(huAnnoFile, colClasses='character', comment.char='#')
moAnno <- read.csv(moAnnoFile, colClasses='character', comment.char='#')
homologene <- read.delim(homologeneFile, header=FALSE)

orthoTable<-ps2ps(huAnno, moAnno, homologene, 10090)  # 10090: mus musculus

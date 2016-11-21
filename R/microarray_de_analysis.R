

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("GEOquery")
# library(Biobase)
# library(GEOquery)
library("limma")
library("vsn")
library('illuminaHumanv2.db')
source('io/microarray.R')

dataDir <- '../data/'
data.dir.raid <- '/media/gabriel/raid1_4tb/data/microarray/'



microarrayDir <- file.path(dataDir, 'microarray_GSE28192')
microarrayFile <- file.path(microarrayDir, 'microarray.RData')

loadCsv <- function (x) {
  read.csv(x, sep='\t', header = 1, row.names = 1)
}

load_from_raw <- function(avg.tech_rpts = T) {
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
  
  if (avg.tech_rpts) {
    # average over technical repeats
    expt.runs <- colnames(expr)[!grepl('-R', colnames(expr))]
    expr <- sapply(expt.runs, function(e) rowMeans(expr[,c(e, paste0(e, '-R'))]))
    samples <- samples[!grepl('-R', samples$name),]    
  }

  save(list=c("expr", "pVal", "samples"), file=microarrayFile, compress=TRUE)
  return(
    list("expr" = expr, "pVal" = pVal, "samples" = samples)
  )
}

# load_from_raw(F)
load(microarrayFile)
expr.entrez <- annotate_by_entrez(expr, illuminaHumanv2.db)
colnames(expr.entrez) <- colnames(expr)

# transform using the Variance Stabilised Transform (VST)
v <- vsn2(data.matrix(expr))
v.entrez <- vsn2(data.matrix(expr.entrez))

# These plots show the effect of the VST
if (F) {
  meanSdPlot(v)
  meanSdPlot(expr)
}

# return to a data.frame and reset sample names
expr.vst <- data.frame(v@hx)
colnames(expr.vst) <- colnames(expr)
expr.entrez.vst <- data.frame(v.entrez@hx)
colnames(expr.entrez.vst) <- colnames(expr)

healthy.samples = c(
  'NT1197',
  'NCb1',
  'NCb2',
  'A911105',
  'A508112',
  'A508285'
)
healthy.samples <- c(healthy.samples, paste0(healthy.samples, '-R'))

mb.samples <- c(
  "Pt1299",
  "ICb1299-I",
  "ICb1299-III",
  "ICb1299-IV"
)
mb.samples = c(mb.samples, paste0(mb.samples, '-R'))

# subsamples <- samples[grepl(paste(c(healthy.samples, mb.samples), collapse = '|'), samples$name),]
subsamples <- subset(
  samples,
  grepl(paste(c(healthy.samples, mb.samples), collapse = '|'), name)
)


subexpr <- expr[,subsamples$name]
subexpr.entrez <- expr.entrez[,subsamples$name]
subexpr.v <- expr.vst[,subsamples$name]
subexpr.entrez.v <- expr.entrez.vst[,subsamples$name]

# setup design matrix
Group <- factor(subsamples$northcott.classification)
# also check that no technical repeat effects exist
TechRpt <- factor(sapply(subsamples[,2], function(x) if (substr(x, nchar(x) - 1, nchar(x)) == '-R') {'B'} else {'A'}))

# the design matrix
X <- model.matrix(~Group + TechRpt)

# fit data to design matrix
fit <- lmFit(subexpr, X)
fit.v <- lmFit(subexpr.v, X)
fit.entrez <- lmFit(subexpr.entrez, X)
fit.entrez.v <- lmFit(subexpr.entrez.v, X)

# if we need contrasts, make them like this:
# contrasts.matrix <- makeContrasts(
#   tech_rpt = es2B, 
#   D_healthy = es1D - es1, 
#   levels=X
# )
# 
# ebfit <- eBayes(contrasts.fit(fit, contrasts.matrix))
# ebfit.v <- eBayes(contrasts.fit(fit.v, contrasts.matrix))
# ebfit.entrez <- eBayes(contrasts.fit(fit.entrez, contrasts.matrix))
# ebfit.entrez.v <- eBayes(contrasts.fit(fit.entrez.v, contrasts.matrix))

ebfit <- eBayes(fit)
ebfit.v <- eBayes(fit.v)
ebfit.entrez <- eBayes(fit.entrez)
ebfit.entrez.v <- eBayes(fit.entrez.v)

# extract top N DE genes / probes above a certain P value
number <- Inf
pval_max = 1.
dehits <- topTable(ebfit, coef='GroupD', number=number, p.value=pval_max)
dehits.v <- topTable(ebfit.v, coef='GroupD', number=number, p.value=pval_max)
dehits.entrez <- topTable(ebfit.entrez, coef='GroupD', number=number, p.value=pval_max)
dehits.entrez.v <- topTable(ebfit.entrez.v, coef='GroupD', number=number, p.value=pval_max)

# change Entrez ID to gene symbol for clarity
eids <- keys(x = illuminaHumanv2.db, keytype = "ENTREZID")
map <- mapIds(x = illuminaHumanv2.db, keys=eids, column='SYMBOL', keytype = 'ENTREZID')

dehits.entrez$symbol <- map[rownames(dehits.entrez)]
dehits.entrez.v$symbol <- map[rownames(dehits.entrez.v)]

# switch probe set ID to gene symbol
pids <- keys(x = illuminaHumanv2.db, keytype='PROBEID')
map.pid <- mapIds(x = illuminaHumanv2.db, keys=pids, column='SYMBOL', keytype = 'PROBEID')

dehits$symbol <- map.pid[rownames(dehits)]
dehits.v$symbol <- map.pid[rownames(dehits.v)]

# remove non-genes
dehits <- na.omit(dehits)
dehits.v <- na.omit(dehits.v)

# find MB-specific genes
genes.wnt <- c("WIF1", "TNC", "GAD1", "DKK2", "EMX2")
genes.shh <- c("PDLIM3", "EYA1", "HHIP", "ATOH1", "SFRP1")
genes.c <- c("IMPG2", "GABRA5", "EYS", "NRL", "MAB21L2", "NPR3")
genes.d <- c("KCNA1", "EOMES", "KHDRBS2", "RBM24", "UNC5D", "OAS1", "OTX2")

# check they are all in the mapper
all(sapply(c(genes.wnt, genes.shh, genes.c, genes.d), function (x) as.logical(length(grep(x, map)))))
all(sapply(c(genes.wnt, genes.shh, genes.c, genes.d), function (x) as.logical(length(grep(x, map.pid)))))

# look for them in the top hits symbol column
dehits.v[melt(lapply(genes.wnt, function(x) which(dehits.v$symbol == x)))$value,]
dehits.v[melt(lapply(genes.shh, function(x) which(dehits.v$symbol == x)))$value,]
dehits.v[melt(lapply(genes.c, function(x) which(dehits.v$symbol == x)))$value,]
dehits.v[melt(lapply(genes.d, function(x) which(dehits.v$symbol == x)))$value,]

dehits.entrez.v[melt(lapply(genes.wnt, function(x) which(dehits.entrez.v$symbol == x)))$value,]
dehits.entrez.v[melt(lapply(genes.shh, function(x) which(dehits.entrez.v$symbol == x)))$value,]
dehits.entrez.v[melt(lapply(genes.c, function(x) which(dehits.entrez.v$symbol == x)))$value,]
dehits.entrez.v[melt(lapply(genes.d, function(x) which(dehits.entrez.v$symbol == x)))$value,]

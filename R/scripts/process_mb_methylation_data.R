source('io/output.R')
source('_settings.R')
source("utils.R")
library("ChAMP")
library("minfi")
library("wateRmelon")
library("data.table")
library(RColorBrewer)
library('openxlsx')

MfromBeta <- function(beta) {
  log2(beta / (1 - beta))
}


split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))


get_idat_basenames <- function(idat.dir) {
  #' Get the basename of all idat files found recursively under the provided directory
  #' The basename is the full path, minus the trailing _Red.idat
  flist <- list.files(path = idat.dir, recursive = T)
  flist <- flist[grep('_Red.idat', flist)]
  
  # convert to basenames for loading
  basenames <- file.path(idat.dir, sub(pattern = "_Red.idat", "", flist))
}


process_idats <- function(
  in.files,
  snames,
  norm.fun=c('swan', 'bmiq', 'funnorm', 'quantile', 'pbc', 'raw'),
  arraytype='EPIC'
) {
  norm.fun = match.arg(norm.fun)
  rgSet <- read.metharray(in.files, extended = T)
  colnames(rgSet) <- snames
  
  mset <- preprocessRaw(rgSet)
  detP <- detectionP(rgSet)
  
  # Load beta values (raw), then apply default ChAMP filtering
  beta.raw <- getBeta(mset, "Illumina")
  champLoad <- champ.filter(beta.raw, detP = detP, pd = NULL, arraytype = arraytype)
  beta.raw <- champLoad$beta
  
  if (norm.fun == 'raw') beta <- beta.raw
  
  if (norm.fun == 'swan') {
    mset.swan <- preprocessSWAN(rgSet, mSet = mset)
    beta.swan <- getBeta(mset.swan)
    beta <- beta.swan[rownames(beta.raw),]
  }
  
  if (norm.fun == 'bmiq') {
    beta <- champ.norm(beta = beta.raw, method = 'BMIQ', arraytype = arraytype, cores=4)
  }
  
  if (norm.fun == 'pbc') {
    beta <- champ.norm(beta = beta.raw, method = 'PBC', arraytype = arraytype)
  }
  
  if (norm.fun == 'funnorm') {
    grSet.funnorm <- preprocessFunnorm(rgSet)
    beta <- getBeta(grSet.funnorm)[rownames(beta.raw),]
  }
  
  if (norm.fun == 'quantile') {
    grSet.quantile <- preprocessQuantile(rgSet)
    beta <- getBeta(grSet.quantile)[rownames(beta.raw),]
  }
  
  return(list(beta.raw=beta.raw, beta=beta))
  
}

# main script starts here
norm.fun <- 'quantile'

samples <- c(
  'ICb1299_Scr', 
  'ICb1299_shBMI1', 
  'ICb1299_shCHD7', 
  'ICb1299_shBMI1CHD7', 
  'p62_3_shBmi1', 
  'p62_3_shChd7', 
  'p62_3_shB+C',
  'p62_3_Scr', 
  '3021_1_Scr', 
  '3021_1_shB', 
  '3021_1_shC', 
  '3021_1_shB+C',
  'S', 
  'B', 
  'C', 
  'B+C'
)

base.dirs <- c(
  file.path(data.dir.raid, 'methylation', '2017-09-19'),
  file.path(data.dir.raid, 'methylation', '2018-01-12'),
  file.path(data.dir.raid, 'methylation', '2018-03-19'),
  file.path(data.dir.raid, 'methylation', '2018-03-26'),
  file.path(data.dir.raid, 'methylation', '2018-04-09')
)

in.files <- NULL
snames <- NULL
batches <- NULL

for (b in base.dirs) {
  meta <- read.csv(file.path(b, 'sources.csv'))
  # set the rownames as filenames
  rownames(meta) <- paste(meta$Sentrix_ID, meta$Sentrix_Position, sep = '_')
  this_files <- get_idat_basenames(file.path(b, 'idat'))
  
  # reorder meta
  meta <- meta[basename(this_files),]
  
  # filter meta
  idx <- meta[, 'sample'] %in% samples
  
  # define file and sample names and add to list
  this_files <- this_files[idx]
  this_snames <- as.vector(meta[idx, 'sample'])
  this_batches <- as.vector(sapply(this_files, function(x){split_path(x)[4]}))
  
  in.files <- c(in.files, this_files)
  snames <- c(snames, this_snames)
  batches <- c(batches, this_batches)
}

# manually remove a failed sample
idx <- in.files != "/home/gabriel/data/methylation/2018-03-26/idat/202081130238/202081130238_R01C01"
in.files <- in.files[idx]
snames <- snames[idx]
batches <- batches[idx]

res <- process_idats(in.files, snames, norm.fun=norm.fun)
beta <- res$beta
beta.raw <- res$beta.raw

# simple metadata
meta <- data.frame(row.names = snames)

meta$cell_line <- '3021'
meta[grep('1299', rownames(meta)), 'cell_line'] <- '1299'
meta[grep('p62', rownames(meta)), 'cell_line'] <- '1299'

meta$batch <- batches

meta$condition <- 'Scr'
# meta[grep('Scr', rownames(meta)), 'condition'] <- 'Scr'
meta[grep('bmi1', rownames(meta), ignore.case = T), 'condition'] <- 'shBMI1'
meta[grep('shb', rownames(meta), ignore.case = T), 'condition'] <- 'shBMI1'

meta[grep('chd', rownames(meta), ignore.case = T), 'condition'] <- 'shCHD7'
meta[grep('shc', rownames(meta), ignore.case = T), 'condition'] <- 'shCHD7'

meta[grep('\\+', rownames(meta)), 'condition'] <- 'shBMI1shCHD7'
meta[grep('bmi1chd7', rownames(meta), ignore.case = T), 'condition'] <- 'shBMI1shCHD7'
meta[rownames(meta) == 'B', 'condition'] <- 'shBMI1'
meta[rownames(meta) == 'C', 'condition'] <- 'shCHD7'

densityPlot(beta, sampGroups = meta$batch, ylim = c(0,7))
densityPlot(beta, sampGroups = meta$cell_line, ylim = c(0,7))

m <- MfromBeta(beta)

pal <- brewer.pal(8,"Dark2")
plotMDS(m, top=NULL, col=pal[factor(meta$batch)])
legend("topright", legend=levels(factor(meta$batch)), text.col=pal, bg="white")

# Run DMP analysis, using batch as a blocking variable
#' We'll need to filter some probes, based on the range of values exhibited. I don't know if this is valid, but without it our statistics
#' get swamped in the MHT correction step.

min_range <- 1.

# 1. Lumping the two cell lines together
rg <- apply(m, 1, range)
rg <- rg[2,] - rg[1,]
ix <- rg > min_range
this.m <- m[ix,]

condition <- factor(meta$condition)
cell_line <- factor(meta$cell_line)
batch <- make.names(factor(meta$batch))

design <- model.matrix(~0 + condition + batch)
fit <- lmFit(this.m, design)

contrasts <- makeContrasts(
  conditionshBMI1-conditionScr, 
  conditionshCHD7-conditionScr, 
  conditionshBMI1shCHD7-conditionScr, 
  levels=design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

#' Outcome: no DMPs

# 2. Running separately for each cell line 
dmps <- list()

for (cl in c('3021', '1299')) {
  dmps[[cl]] <- list()
  print(paste0("Cell line ", cl))
  this.meta <- meta[meta$cell_line == cl,]
  this.m <- m[, rownames(this.meta)]
  
  # filter by range: require a minimum change over the samples
  rg <- apply(this.m, 1, range)
  rg <- rg[2,] - rg[1,]
  hist(rg, 40, xlab='M value range', main=paste0(cl, ' M value range'))
  abline(v=min_range, col='red', lty=2)
  
  ix <- rg > min_range
  print(paste0("Of the ", length(ix), " probes, ", sum(ix), " are retained for DMP testing."))
  
  this.m <- this.m[ix,]
  
  condition <- factor(this.meta$condition)
  batch <- make.names(factor(this.meta$batch))
  
  design <- model.matrix(~0 + condition + batch)
  fit <- lmFit(this.m, design)
  
  contrasts <- makeContrasts(
    conditionshBMI1-conditionScr, 
    conditionshCHD7-conditionScr, 
    conditionshBMI1shCHD7-conditionScr, 
    levels=design
  )
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  print(summary(decideTests(fit2)))
  
  for (i in seq(ncol(contrasts))) {
    ttl <- colnames(contrasts)[i]
    ttl <- gsub(pattern = "condition", replacement = '', x = ttl)
    dmps[[cl]][[ttl]] <- topTable(fit2, coef = i, number = Inf, p.value = 0.05)
    ## TODO: export results to xlsx or similar
    # write.csv(dmps[[cl]][[colnames(contrasts)[i]]])
  }
  write.xlsx(dmps[[cl]], paste0("dmps_", cl, ".xlsx"))

}

for (cl in c('3021', '1299')) {
  dmps[[cl]] <- list()
  print(paste0("Cell line ", cl))
  for (i in seq(ncol(contrasts))) {
    ttl <- colnames(contrasts)[i]
    ttl <- gsub(pattern = "condition", replacement = '', x = ttl)
    dmps[[cl]][[ttl]] <- dmps2[[cl]][[colnames(contrasts)[i]]]
  }
}

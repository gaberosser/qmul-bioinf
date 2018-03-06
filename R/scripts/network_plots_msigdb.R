source("_settings.R")

indir <- file.path(data.dir, "msigdb")
fn <- file.path(indir, "h.all.v6.1.symbols.gmt")

msig <- read.csv(fn, sep='\t', header = F, row.names = 1)

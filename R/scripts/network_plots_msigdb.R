source("_settings.R")

indir <- file.path(data.dir, "msigdb")
fn <- file.path(indir, "h.all.v6.1.symbols.gmt")

msig <- read.csv(fn, sep='\t', header = F, row.names = 1)
all_genes <- lapply(msig[, 2:ncol(msig)], function (x) {as.vector(x[x != ""])})
all_genes <- as.vector(unlist(all_genes))
pathways <- rownames(msig)


sprintf("Loaded %d genes (%d unique) involved in %d pathways.", length(all_genes), length(unique(all_genes)), nrow(msig))

vertices <- data.frame(node=unique(all_genes))

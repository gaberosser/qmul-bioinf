library(edgeR)
source('io/rnaseq.R')


#' Prepare a table of DE genes with the given FDR cutoff.
#' lrt: must be compatible with topTags, for example the output of glmLRT
prepare_de_table <- function(lrt, fdr=0.05) {
  de <- as.data.frame(topTags(lrt, p.value = fdr, n = Inf))
  de$ensembl <- rownames(lrt$table)[as.integer(rownames(de))]
  de$direction <- ifelse(de$logFC > 0, 'U', 'D')
  de <- de[, c("genes", "logFC", "ensembl", "direction", "FDR", "logCPM")]
  de
}

#' Get lists of DE genes for all possible Venn segments for an arbitrary number of comparisons.
#' The input arguments are DGELRT objects or anything else that can be passed into prepare_de_table.
venn_edger_de_lists <- function(..., fdr=0.05, id.key='ensembl') {

  de = lapply(list(...), function (x) {prepare_de_table(x, fdr=fdr)})
  ids <- lapply(de, function(x){x[[id.key]]})

  blocks <- list()
  
  nbit <- length(ids)
  n <- 2 ** nbit - 1
  
  
  for (i in seq(1, n)) {
    comb <- as.integer(intToBits(i))[1:nbit]
    idx.in <- which(comb == 1)
    idx.out <- which(comb == 0)
    
    ids.in <- Reduce(intersect, lapply(idx.in, function(x){ ids[[x]] }))
    ids.out <- Reduce(union, lapply(idx.out, function(x){ ids[[x]] }))
    ids.this <- setdiff(ids.in, ids.out)
    
    get_de <- function(x) {
      tmp <- de[[x]][ids[[x]] %in% ids.this,]
      tmp <- tmp[order(tmp[[id.key]]),]
    }
    
    de.in <- lapply(idx.in, get_de)
    
    blocks[[paste0(comb, collapse = '')]] <- do.call(cbind, de.in)
  }
  
  blocks
}

grouped_analysis <- function(data, groups, groups.lumped, contrasts, gene.symbols=NULL, output.dir=NULL) {
  #' contrasts: list of characters containing valid formulae based on design matrix ~0 + groups

  groups <- as.factor(groups)
  groups.lumped <- as.factor(groups.lumped)
  
  if (is.null(gene.symbols)) {
    gene.symbols <- ens.map[rownames(data), "hgnc_symbol"]
  }
  
  y.lumped <- DGEList(counts=data, group=groups.lumped)
  y.lumped <- calcNormFactors(y.lumped)

  design <- model.matrix(~groups.lumped)
  
  # estimate dispersion of lumped groups
  y.lumped <- estimateDisp(y.lumped, design)
  
  #' Store the values for later use
  #' This is one of the recommended approaches from the authors of edgeR in the situation where no replicates are available
  dispersion.trended.lumped <- y.lumped$trended.dispersion
  dispersion.common.lumped <- y.lumped$common.dispersion
  dispersion.tagwise.lumped <- y.lumped$tagwise.dispersion
  
  # now re-initialise with the correct groupings
  y <- DGEList(data, genes = gene.symbols)
  y <- calcNormFactors(y)
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  
  # in this case, we need to use the dispersion estimated earlier
  y$common.dispersion <- dispersion.common.lumped
  y$trended.dispersion <- dispersion.trended.lumped
  y$tagwise.dispersion <- dispersion.tagwise.lumped
  
  fit.glm <- glmFit(y, design)
  
  contrasts[['levels']] <- design
  contrasts.made <- do.call(makeContrasts, contrasts)
  
  lrt <- list()
  for (t in colnames(contrasts.made)) {
    lrt[[t]] <- glmLRT(fit.glm, contrast=contrasts.made[, t])
    if (!is.null(output.dir)) {
      #' Save the results
      output.fn = file.path(output.dir, paste0(t, '.csv'))
      write.csv(prepare_de_table(lrt[[t]], fdr=fdr), file=output.fn)
    }
  }
  
  lrt
  
}


filter_genes <- function(data, cpm.min = 1, nsamples.min = 3) {
  #' Filter the genes (rows in data) based on prevalence, in order to remove genes with consistently low expression. 
  #' cpm.min: The minimum CPM required to 'pass'
  #' nsamples.min: The minimum number of samples that must pass in order to keep this gene
  #' Returns a reduced data frame
  
  y <- DGEList(counts=data)
  keep <- rowSums(cpm(y) > cpm.min) >= nsamples.min
  data <- data[keep,]
}

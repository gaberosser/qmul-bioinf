library(edgeR)
source('io/rnaseq.R')


#' Prepare a table of DE genes with the given FDR cutoff.
#' lrt: must be compatible with topTags, for example the output of glmLRT
prepare_de_table <- function(lrt, fdr = 0.05, log2FC.min = NULL) {
  de <- as.data.frame(topTags(lrt, p.value = fdr, n = Inf))
  if (!is.null(log2FC.min)) {
    de <- de[abs(de$logFC) >= log2FC.min,]
  }
  de$ensembl <- rownames(lrt$table)[as.integer(rownames(de))]
  de$direction <- ifelse(de$logFC > 0, 'U', 'D')
  de <- de[, c("genes", "logFC", "ensembl", "direction", "FDR", "logCPM")]
  de
}


#' Get lists of DE genes for all possible Venn segments for an arbitrary number of comparisons.
#' The input arguments are DGELRT objects or anything else that can be passed into prepare_de_table.
#' If `background` is supplied, only features that are NOT included in the background are counted
venn_edger_de_lists <- function(..., background=NULL, fdr=0.05, log2FC.min = NULL, id.key='ensembl') {

  de <- lapply(list(...), function (x) {prepare_de_table(x, fdr=fdr, log2FC.min = log2FC.min)})
  
  if (!is.null(background)) {
    print("BACKGROUND")
    de.ref <- prepare_de_table(background, fdr=fdr)
    de <- lapply(de, function(x) {
      these.ids <- setdiff(x[[id.key]], de.ref[[id.key]])
      x[x[[id.key]] %in% these.ids,]
    }) 
  } 
  
  ids <- lapply(de, function(x){x[[id.key]]})
  
  blocks.ids <- do.call(venn_sets, ids)
  
  #' Function to get the original topTags data corresponding to each Venn block
  #' The input is a string representing the binary index, e.g. '0011' and a list of the ids
  get_de <- function(bin) {
    these.ids <- blocks.ids[[bin]]
    idx.in <- which(strsplit(bin, "")[[1]] == 1)
    tmp <- lapply(
      idx.in, 
      function (x) {
        aa <- de[[x]][de[[x]][[id.key]] %in% these.ids,]
        aa[order(aa[, 'FDR']),]
      }
    )
    do.call(cbind, tmp)
  }
  
  blocks <- lapply(names(blocks.ids), get_de)
  names(blocks) <- names(blocks.ids)

  blocks$contrasts <- names(de)
  blocks
}

grouped_analysis <- function(data, groups, groups.lumped, contrasts, gene.symbols=NULL, output.dir=NULL, fdr = 0.05) {
  #' Carry out analysis based on a paired study, where there are no replicates available, e.g. cancer vs healthy tissue in individual patients.
  #' In this case, we need to 'lump' individuals together based on a simpler grouping (cancer vs healthy) to estimate dispersion, then transfer 
  #' this estimate across for the purpose of computing differential expression whilst honouring the paired structure.
  #' contrasts: list of characters containing valid formulae based on design matrix ~0 + groups. makeContrasts will be called on each element of the list.
  #' data: Numeric matrix containing gene counts.
  #' groups: vector or factor giving the group each sample belongs to.
  #' groups.lumped: vector or factor giving the lumped groups for each sample.

  groups <- factor(groups)
  groups.lumped <- factor(groups.lumped)
  
  if (is.null(gene.symbols)) {
    ens.map <- biomart_annotation(index.by='ensembl_gene_id')
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

filter_genes <- function(data, cpm.min = 1, nsamples.min = 3, unless.cpm.gte = NULL) {
  #' Filter the genes (rows in data) based on prevalence, in order to remove genes with consistently low expression. 
  #' cpm.min: The minimum CPM required to 'pass'
  #' nsamples.min: The minimum number of samples that must pass in order to keep this gene
  #' unless.cpm.gte: If supplied, this acts as an override; if any single CPM value is >= this value, the gene is retained, even if
  #' it would otherwise be earmarked for removal.
  #' Returns a reduced data frame
  
  y <- DGEList(counts=data)
  keep <- rowSums(cpm(y) > cpm.min) >= nsamples.min
  if (!is.null(unless.cpm.gte)) {
    keep <- keep | (rowSums(cpm(y) >= unless.cpm.gte) > 0)
  }
  data <- data[keep,]
}

export_de_list <- function(blocks, outfile) {
  #' Generate a Venn-like DE CSV and write to disk
  #' Blocks is a list with names as a binary representation of the Venn region, e.g. 1101 means "in groups 1, 2 and 4 but not 3".
  #'
  
  # only keep names that are in binary format
  idx <- names(blocks)[grep('^[01]+$', names(blocks))]
  idx <- idx[order(idx, decreasing = T)]
  
  # check the format: all should have the same length
  ns <- sapply(idx, nchar)
  if (!all(ns == ns[1])) {
    stop("Unequal block names. Expecting them to have the same format, e.g. `011`.")
  }
  
  n <- ns[[1]]
  message(sprintf("Exporting %i way DE comparison to %s.", n, outfile))
  
  
  csv.data <- data.frame(blocks[[idx[1]]])
  if (ncol(csv.data) %% n != 0) {
    stop(sprintf("Unequal number of rows detected (%i / %i)", ncol(csv.data), n))
  } else {
    csv.ncol <- as.integer(ncol(csv.data) / n)
    message(sprintf("Detected %i columns per block.", csv.ncol))
  }
  block.colnames <- colnames(csv.data)[1:csv.ncol]

  for (i in seq(2, 2^n - 1)) {
    # build this block
    this.data <- blocks[[idx[i]]]
    this.nrow <- nrow(this.data)
    this.block <- list()
    k <- 1
    l <- 1
    for (j in strsplit(idx[i], "")[[1]]) {
      if (j == '1') {
        this.block[[k]] <- this.data[, (csv.ncol * (l - 1) + 1):(csv.ncol * l)]
        l <- l + 1
      } else {
        this.block[[k]] <- data.frame(rep.col(rep.row("", this.nrow), csv.ncol))
        colnames(this.block[[k]]) <- block.colnames
      }
      k <- k + 1
    }
    this.csvdata <- do.call(cbind, c(this.block, list(deparse.level=0)))
    csv.data <- rbind(csv.data, this.csvdata)
  }
  
  write.csv(csv.data, file = outfile, row.names = F)
  
  
  
}

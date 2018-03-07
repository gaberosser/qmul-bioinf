source("_settings.R")

indir <- file.path(data.dir, "msigdb")
fn <- file.path(indir, "h.all.v6.1.symbols.gmt")

msig <- read.csv(fn, sep='\t', header = F, row.names = 1)
all_genes <- lapply(msig[, 2:ncol(msig)], function (x) {as.vector(x[x != ""])})
all_genes <- as.vector(unlist(all_genes))
pathways <- rownames(msig)

sprintf("Loaded %d genes (%d unique) involved in %d pathways.", length(all_genes), length(unique(all_genes)), nrow(msig))

#' create a data.frame containing all the edges
#' every node is linked to the origin (which won't be plotted)
#' pathways are represented at the first layer, genes at the leaf nodes
d1 = data.frame(from="origin", to=pathways)

a <- apply(msig[, 2:ncol(msig)], 1, function (x) {as.vector(x[x != ""])})
b <- lapply(seq_along(a), function(i) data.frame(from=names(a)[[i]], to=a[[i]]))
d2 <- do.call("rbind", b)

edges=rbind(d1, d2)

# compute the number of pathways shared by each gene
conn_map <- sapply(all_genes, function (x) {
  sapply(a, function(t) x %in% t)
})
conn_count <- sapply(all_genes, function (x) {
  sum(sapply(a, function(t) x %in% t))
})
conn_count <- conn_count[unique(names(conn_count))]

# create a dataframe with connection between leaves (individuals)

all_leaves=paste("subgroup", seq(1,100), sep="_")
connect=rbind( data.frame( from=sample(all_leaves, 100, replace=T) , to=sample(all_leaves, 100, replace=T)), data.frame( from=sample(head(all_leaves), 30, replace=T) , to=sample( tail(all_leaves), 30, replace=T)), data.frame( from=sample(all_leaves[25:30], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)), data.frame( from=sample(all_leaves[75:80], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)) )
connect$value=runif(nrow(connect))
  
# create a vertices data.frame. One line per object of our hierarchy
all_nodes <- unique(c(as.character(edges$from), as.character(edges$to)))
vertices = data.frame(
  name = all_nodes,
  value = c(rep(0, length(pathways) + 1), conn_count)
)
vertices$group = edges$from[ match( vertices$name, edges$to ) ]

# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# connection object - must refer to the ids of the leaves
from = match( connect$from, vertices$name)
to = match( connect$to, vertices$name)

## TODO: we are going to have to make the leaf labels UNIQUE (i.e. not just gene symbols)
# Otherwise we won't be able to make the connections correctly.
# Suggest PATHWAY_genesymbol
# Then we need to specify the data for the labels manually (or as a separate column in the existing data structure??)
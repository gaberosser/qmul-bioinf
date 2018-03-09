source("_settings.R")

indir <- file.path(data.dir, "msigdb")
fn <- file.path(indir, "h.all.v6.1.symbols.gmt")

msig <- read.csv(fn, sep='\t', header = F, row.names = 1)
all_genes <- lapply(msig[, 2:ncol(msig)], function (x) {as.vector(x[x != ""])})
all_genes <- as.vector(unlist(all_genes))
pathways <- rownames(msig)

sprintf("Loaded %d genes (%d unique) involved in %d pathways.", length(all_genes), length(unique(all_genes)), nrow(msig))
all_genes <- unique(all_genes)

#' create a data.frame containing all the edges
#' every node is linked to the origin (which won't be plotted)
#' pathways are represented at the first layer, genes at the leaf nodes
d1 = data.frame(from="origin", to=pathways)

a <- apply(msig[, 2:ncol(msig)], 1, function (x) {as.vector(x[x != ""])})
b <- lapply(seq_along(a), function(i) data.frame(from=names(a)[[i]], to=paste(names(a)[[i]], a[[i]], sep = '|')))
d2 <- do.call("rbind", b)

edges=rbind(d1, d2)

# compute the number of pathways shared by each gene
conn_map <- sapply(all_genes, function (x) {
  sapply(a, function(t) x %in% t)
})
conn_count <- colSums(conn_map)

# create a dataframe with connection between leaves (individuals)
conns = list()
for (n in names(conn_count)[conn_count > 1]) {
  the_pathways <- names(which(conn_map[,n]))
  the_genes <- paste(the_pathways, n, sep='|')
  this <- data.frame(t(combn(the_genes, 2)))
  colnames(this) <- c('from', 'to')
  this$values <- conn_count[[n]]
  conns[[n]] <- this
}

connect = do.call("rbind", conns)
rownames(connect) <- seq(1, nrow(connect))

# create a vertices data.frame. One line per object of our hierarchy
all_nodes <- unique(c(as.character(edges$from), as.character(edges$to)))
gene_nodes <- sub('[^|]*\\|', '', all_nodes[grep('\\|', all_nodes)])
node_conn_num <- c(rep(0, length(pathways) + 1), conn_count[gene_nodes])
node_labels <- sub('[^|]*\\|', '', all_nodes)

vertices = data.frame(
  name = all_nodes,
  value = node_conn_num ** 2,
  labels = node_labels
)
vertices$group = edges$from[ match( vertices$name, edges$to ) ]

# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# connection object - must refer to the ids of the leaves
from = match( connect$from, vertices$name)
to = match( connect$to, vertices$name)

# node labels 


ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="black", width=0.9) +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=labels), size=1.5, alpha=1) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))


## TODO: we are going to have to make the leaf labels UNIQUE (i.e. not just gene symbols)
# Otherwise we won't be able to make the connections correctly.
# Suggest PATHWAY_genesymbol
# Then we need to specify the data for the labels manually (or as a separate column in the existing data structure??)
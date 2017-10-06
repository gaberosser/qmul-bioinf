source('differential_expression/edger_de.R')
source('io/rnaseq.R')
source('io/output.R')
source('_settings.R')
source("utils.R")


loaded <- mouse_validation_data()
data <- filter_genes(loaded$data)

groups = as.factor(c(
  rep('eNSCmed', 3),
  rep('eNSCmouse', 3),
  rep('iNSCmouse', 3),
  rep('iNSChuman', 3)
))
design <- model.matrix(~0 + groups)
colnames(design) <- sub('groups', '', colnames(design))

y <- DGEList(counts=data, group=groups)
y <- calcNormFactors(y)

# estimate dispersion of lumped groups
y <- estimateDisp(y, design)
fit.glm <- glmFit(y, design)

# CPM
cpm <- cpm(y, log=T, prior.count = 1)

res.imouse_emed <- topTags(glmLRT(fit.glm, contrast=makeContrasts("iNSCmouse-eNSCmed", levels=groups)), n=Inf, p.value = 0.01)
res.imouse_emed$title <- "iNSC (mouse protocol) vs eNSC"
res.emouse_emed <- topTags(glmLRT(fit.glm, contrast=makeContrasts("eNSCmouse-eNSCmed", levels=groups)), n=Inf, p.value = 0.01)
res.emouse_emed$title <- "eNSC (mouse induction media) vs eNSC"
res.ihuman_imouse <- topTags(glmLRT(fit.glm, contrast=makeContrasts("iNSChuman-iNSCmouse", levels=groups)), n=Inf, p.value = 0.01)
res.ihuman_imouse$title <- "iNSC (human protocol) vs iNSC (mouse protocol)"

# get gene symbols using biomaRt
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset = "mmusculus_gene_ensembl")

for (t in list(res.imouse_emed, res.emouse_emed, res.ihuman_imouse)) {
  mgi <- getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), filters='ensembl_gene_id', values=rownames(t), mart=ensembl)
  t$table[mgi$ensembl_gene_id, "gene_symbol"] <- mgi$mgi_symbol

  # select the CPM from DE genes only
  sel_cpm <- cpm[
    rownames(t)[t$table$FDR < 1e-6 & abs(t$table$logFC) > 1],
    c(grep(pattern = 'eNSC[0-9]med', colnames(data)), grep(pattern = 'mDura.*mouse', colnames(data)))
    ]
  
  # clustering the DE genes
  # c <- cor(t(sel_cpm))
  # fill in NA values
  # c[is.na(c)] <- 0.
  # a <- as.dist(1 - c)
  
  a <- dist(sel_cpm, method = 'euclidean')
  b <- hclust(a, method = "average")
  plot(b, labels=F, xlab = "", sub = "", ylab = "Euclidean distance", main = t$title)
  
}



library(gplots)
png("heatmap", width=5, height=10)
heatmap.2(sel_cpm, trace="none", ylab=NULL)
dev.off()

# volcano plots
volcanoData <- cbind(res.imouse_emed$table$logFC, -log10(res.imouse_emed$table$FDR))
colnames(volcanoData) <- c("logFC", "-log10(FDR)")
plot(volcanoData, pch=19)

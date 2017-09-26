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
glmLRT(fit.glm, contrast="iNSCmouse-eNSCmed")

res.imouse_emed <- topTags(glmLRT(fit.glm, contrast=makeContrasts("iNSCmouse-eNSCmed", levels=groups)), n=Inf, p.value = 0.01)
res.emouse_emed <- topTags(glmLRT(fit.glm, contrast=makeContrasts("eNSCmouse-eNSCmed", levels=groups)), n=Inf, p.value = 0.01)
res.ihuman_imouse <- topTags(glmLRT(fit.glm, contrast=makeContrasts("iNSChuman-iNSCmouse", levels=groups)), n=Inf, p.value = 0.01)

# volcano plots
volcanoData <- cbind(res.imouse_emed$table$logFC, -log10(res.imouse_emed$table$FDR))
colnames(volcanoData) <- c("logFC", "-log10(FDR)")
plot(volcanoData, pch=19)

cpm <- cpm(y, log=T, prior.count = 1)
sel_cpm <- cpm[
  rownames(res.imouse_emed)[res.imouse_emed$table$FDR < 1e-6 & abs(res.imouse_emed$table$logFC) > 1],
  c(grep(pattern = 'eNSC[0-9]med', colnames(data)), grep(pattern = 'mDura.*mouse', colnames(data)))
]

library(gplots)
png("heatmap", width=5, height=10)
heatmap.2(sel_cpm, trace="none", ylab=NULL)
dev.off()

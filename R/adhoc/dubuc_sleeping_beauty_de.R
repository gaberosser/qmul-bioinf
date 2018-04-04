library(limma)
library(Biobase)

source("_settings.R")

in.dir <- file.path(data.dir, 'sleeping_beauty_mouse_screen')

# phenoData
tmp <- read.csv(file.path(in.dir, "pheno_data.csv"), row.names = 1)
groups <- factor(as.matrix(tmp))
pdata <- AnnotatedDataFrame(tmp)

# featureData
tmp <- read.csv(file.path(in.dir, "feature_data.csv"), row.names = 1)
fdata <- AnnotatedDataFrame(tmp)

# expressionData
tmp <- read.csv(file.path(in.dir, "expression_data.csv"), row.names = 1)
m <- as.matrix(tmp)

eset <- new("ExpressionSet", exprs = m, phenoData = pdata, featureData = fdata)

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(eset, design)
fit <- eBayes(fit)

contrasts <- makeContrasts(ins="Ins-WT", levels=design)

fit  <- contrasts.fit(fit, contrasts)
fit  <- eBayes(fit)

data("water", package = "HSAUR")
source("io/microarray.R")
library("limma")
data.dir <- '../data/'

probes <- c(
  letters,
  'aa', 'bb', 'cc', 'dd'
)

dat <- dubuc_sb_screen(by.gene = F)
dat <- add_annotation_column(dat, mogene10sttranscriptcluster.db, col.name = 'SYMBOL')
dat <- na.omit(dat)
# dat <- dat[dat$SYMBOL == 'Chd7', colnames(dat) != 'SYMBOL']
dat <- dat[, colnames(dat) != 'SYMBOL']
# we want each probe to occupy a column
# dat <- t(dat)
# colnames(dat) <- probes
# probes <- cbind(colnames(dat))
# arr <- as.data.frame(dat)

chd7.samples <- c(
  "Wu050",
  "Wu053",
  "Wu054"
)
# arr$chd7 <- rownames(arr) %in% chd7.samples

# form = as.formula(paste('cbind(', paste(probes, collapse = ','), ')', ' ~ chd7'))
# summary(manova(form, data=arr), test='Hotelling-Lawley')
## ERROR: doesn't support p > n
## TODO: switch to DE analysis with limma

chd7 <- factor(colnames(dat) %in% chd7.samples)
X <- model.matrix(~chd7)
fit <- lmFit(dat, X)
ebfit <- eBayes(fit)

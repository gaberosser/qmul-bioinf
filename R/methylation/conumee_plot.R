library("minfi")
library("conumee")

source("_settings.R")


ffpe.fn <- file.path(data.dir.raid, "methylation", "2016-06-10_brandner", "idat", "200512330080", "200512330080_R04C01")
gic.fn <- file.path(data.dir.raid, "methylation", "2016-12-19_ucl_genomics", "idat", "200788220017", "200788220017_R08C01")

rgSet <- read.metharray(c(ffpe.fn, gic.fn), extended = T)
colnames(rgSet) <- c('FFPE 019', 'GIC 019 P4')

mSet <- preprocessSWAN(rgSet)

data("detail_regions")  # supplied by conumee
data("exclude_regions")  # supplied by conumee
anno <- CNV.create_anno(array_type = 'EPIC', exclude_regions=exclude_regions, detail_regions = detail_regions)
common_probes <- rownames(mSet)[is.element(rownames(mSet), names(anno@probes))]
anno@probes <- anno@probes[common_probes]

computeCNV <- function(query, control, anno, name=NULL) {
  x <- CNV.fit(query, control, anno, name = name)
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
}

minfi.data <- CNV.load(mSet)

# What shall we use as a control here??

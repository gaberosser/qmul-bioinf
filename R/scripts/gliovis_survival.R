#' For each of these packages, you need to install it using install.packages("PACKAGE")

library(survival)
library(gridExtra)
library(GGally)
library(DT)
library(shinyBS)
# library(limma)
library(survminer)
library(weights)
library(plotly)


survivalPlot <- function (df, gene, cutoff, censor=F, risk.table=F, conf.int=T, font.legend=12) {
  x <- df[, c(gene, "status", "survival")]
  mRNA <- df[, gene]
  surv.status <- x[ ,"status"]
  surv.time <- x[ ,"survival"]
  x$cutoff_group <- ifelse(mRNA >= cutoff, c("high"), c("low"))
  
  my.Surv <- Surv(time = surv.time, event = surv.status== 1)
  smax <- max(surv.time, na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  expr.surv <- survfit(my.Surv ~ cutoff_group, data=x)
  log.rank <- survdiff(my.Surv ~ cutoff_group, rho = 0, data=x)
  mantle.cox <- survdiff(my.Surv ~ cutoff_group, rho = 1, data=x)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ cutoff_group, data=x))
  
  # get data for the plot legend
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  
  legend.labs = c(sprintf("%s High (n=%s, events=%s, median=%s)", gene, surv$records[1], surv$events[1], surv$median[1]),
                  sprintf("%s Low (n=%s, events=%s, median=%s)", gene, surv$records[2], surv$events[2], surv$median[2]))
  xlegend <- 0.60

  # main = "TITLE"
  p <- survminer::ggsurvplot(fit = expr.surv, censor = censor, conf.int = conf.int, legend = c(xlegend,0.9), surv.scale = "percent", ylab = "Surviving", xlab = "Survival time (Months)",
                             xlim = c(0,smax), legend.labs = legend.labs, legend.title = "", font.legend = font.legend, risk.table = risk.table,
                             risk.table.y.text = F, risk.table.y.text.col = T, risk.table.height = 0.4)
  plot <- p$plot
  plot <- plot + annotate("text", x = xmax, y = c(0.725,0.65,0.575), size = font.legend/3,
                          label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                    sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                    sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p)))
  p$plot <- plot
  p$logrank.p <- log.rank.p
  p$wilcoxon.p <- mantle.cox.p
  p$hr <- HR
  return(p)
}


##### CHANGE THESE PARAMETERS

# the list of genes of interest
gois <- c(
  "ANO1",
  "CTPS1",
  "CCKAR",
  "SLC1A4"
)


indir <- "./"
outdir <- "survival_plots/017"

in.expr.rna <- "2018-02-19_TCGA_GBM_expression_RNA.txt"
in.expr.agilent <- "2018-02-19_TCGA_GBM_expression_Agilent.txt"
in.expr.u133a <- "2018-02-19_TCGA_GBM_expression_HGU133A.txt"
in.pheno.rna <- "2018-02-19_TCGA_GBM_pheno_RNA.txt"
in.pheno.agilent <- "2018-02-19_TCGA_GBM_pheno_Agilent.txt"
in.pheno.u133a <- "2018-02-19_TCGA_GBM_pheno_HGU133A.txt"

censor=F
risk.table=F
conf.int=T
font.legend=12

##### DON'T CHANGE ANYTHING ELSE!!

dir.create(outdir, recursive = T)

pdata <- read.csv(file.path(indir, in.pheno.rna), sep = "\t", row.names = 1)

expr.rna <- read.csv(in.expr <- file.path(indir, in.expr.rna), sep = "\t", row.names = 1)
expr.agilent <- read.csv(file.path(indir, in.expr.agilent), sep = "\t", row.names = 1)
expr.u133a <- read.csv(file.path(indir, in.expr.u133a), sep = "\t", row.names = 1)

# filter: primary tumour, non-GCIMP, IDH1 WT
idx <- (pdata$Recurrence == "Primary") & (pdata$CIMP_status == "NON G-CIMP")
# drop anything that doesn't have a recurrence label
idx[is.na(idx)] <- FALSE

pdata <- pdata[idx,]
expr <- list(
  u133a=expr.u133a[rownames(pdata),],
  agilent=expr.agilent[rownames(pdata),],
  rna=expr.rna[rownames(pdata),]
)

####### TO SPEED THINGS UP, RE-RUN FROM HERE

# construct the input for getting the cutoff
cutoff <- data.frame(row.names = gois)
surv.logrankp <- data.frame(row.names = gois)
surv.wilcoxonp <- data.frame(row.names = gois)
surv.hr <- data.frame(row.names = gois)

for (t in names(expr)) {
  this_expr <- expr[[t]]
  
  row.has.na <- apply(this_expr, 1, function(x){any(is.na(x))})
  this_expr <- this_expr[!row.has.na,]
  this_pdata <- pdata[rownames(this_expr),]
  
  print(paste0("Data source: ", t, ". Samples remaining after filtering: ", nrow(this_expr), ". Genes: ", ncol(this_expr), "."))
  
  X <- data.frame(
    this_pdata[,c("survival", "status")],
    this_expr
  )
  
  the_genes <- gois[gois %in% colnames(this_expr)]
  if (length(the_genes) != length(gois)) {
    print(paste0("WARNING: some GOIs were not found in the data (", t, ")"))
  }
  
  this_cutoff <- surv_cutpoint(X, time="survival", event="status", variables=the_genes)
  this_cutpoint <- this_cutoff$cutpoint[the_genes, 'cutpoint']
  cutoff[the_genes, t] <- this_cutoff$cutpoint[the_genes, 'cutpoint']
  
  for (g in the_genes) {
    
    ## it should be possible to run this in the survival_plots() function, but something in the implementation isn't playing nice...
    
    x <- X[, c(g, "status", "survival")]
    mRNA <- X[, g]
    surv.status <- x[ ,"status"]
    surv.time <- x[ ,"survival"]
    x$cutoff_group <- ifelse(mRNA >= cutoff[g, t], c("high"), c("low"))
    
    my.Surv <- Surv(time = surv.time, event = surv.status== 1)
    smax <- max(surv.time, na.rm = TRUE)
    tmax <- smax-(25*smax)/100
    xmax <- (90*tmax)/100
    expr.surv <- survfit(my.Surv ~ cutoff_group, data=x)
    log.rank <- survdiff(my.Surv ~ cutoff_group, rho = 0, data=x)
    mantle.cox <- survdiff(my.Surv ~ cutoff_group, rho = 1, data=x)
    surv <- data.frame(summary(expr.surv)$table)
    model <- summary(coxph(my.Surv ~ cutoff_group, data=x))
    
    # get data for the plot legend
    HR <- round(model$conf.int[1],2)
    HR.lower <- round(model$conf.int[3],2)
    HR.upper <- round(model$conf.int[4],2)
    log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
    mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
    star.log <- weights::starmaker(log.rank.p)
    star.mcox <- weights::starmaker(mantle.cox.p)
    
    legend.labs = c(sprintf("%s High (n=%s, events=%s, median=%s)", g, surv$records[1], surv$events[1], surv$median[1]),
                    sprintf("%s Low (n=%s, events=%s, median=%s)", g, surv$records[2], surv$events[2], surv$median[2]))
    xlegend <- 0.60
    
    # main = "TITLE"
    this_surv <- survminer::ggsurvplot(fit = expr.surv, censor = censor, conf.int = conf.int, legend = c(xlegend,0.9), surv.scale = "percent", ylab = "Surviving", xlab = "Survival time (Months)",
                               xlim = c(0,smax), legend.labs = legend.labs, legend.title = "", font.legend = font.legend, risk.table = risk.table,
                               risk.table.y.text = F, risk.table.y.text.col = T, risk.table.height = 0.4)
    plot <- this_surv$plot
    plot <- plot + annotate("text", x = xmax, y = c(0.725,0.65,0.575), size = font.legend/3,
                            label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                      sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                      sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p)))
    this_surv$plot <- plot
    this_surv$logrank.p <- log.rank.p
    this_surv$wilcoxon.p <- mantle.cox.p
    this_surv$hr <- HR
    
    # this_surv <- survivalPlot(X, g, cutoff[g, t])
    filestem <- file.path(outdir, paste0(t, "_", g, "."))
    ggsave(paste0(filestem, "png"), plot=this_surv$plot, dpi = 300)
    ggsave(paste0(filestem, "tiff"), plot=this_surv$plot, dpi = 300)
    surv.logrankp[g, t] <- this_surv$logrank.p
    surv.wilcoxonp[g, t] <- this_surv$wilcoxon.p
    surv.hr[g, t] <- this_surv$hr
  }
  
}

print("Log rank P values:")
print(surv.logrankp)
print("Wilcoxon P values:")
print(surv.wilcoxonp)
print("Hazard ratios:")
print(surv.hr)

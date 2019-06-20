source('_settings.R')
source('io/output.R')

nanMedian <- function(x) {
  x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
  x[is.na(x)] = median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}


imputeMissing <- function(data) {
  data.frame(apply(data, 2, nanMedian))
}


predictivePerformance <- function(fit, truth) {
  pred <- predict(fit, type='response') > 0.5
  right <- pred == truth
  tp <- right & truth
  tn <- right & (!truth)
  fp <- (!right) & (!truth)
  fn <- (!right) & (truth)
  list(tp=sum(tp), fp=sum(fp), tn=sum(tn), fn=sum(fn))
}


plotDecisionBoundary <- function(fit, coeffx, coeffy, values, outcomes, xlab=NULL, ylab=NULL) {
  # let's plot the decision boundary for this model
  colSd <- apply(values[,c(coeffx, coeffy)], 2, function (x) sd(na.omit(x)))
  colMn <- apply(values[,c(coeffx, coeffy)], 2, function (x) mean(na.omit(x)))
  
  b0 <- coef(fit)[['(Intercept)']]
  b_peak <- coef(fit)[[coeffy]]
  b_trough <- coef(fit)[[coeffx]]
  s_peak <- colSd[[coeffy]]
  s_trough <- colSd[[coeffx]]
  m_peak <- colMn[[coeffy]]
  m_trough <- colMn[[coeffx]]
  
  incpt <- s_peak / b_peak * (-b0 + b_trough / s_trough * m_trough + b_peak / s_peak * m_peak)
  slope <- -(s_peak * b_trough) / (s_trough * b_peak)
  
  fit_x <- c(0, 300)
  fit_y = incpt + slope * fit_x
  plot(fit_x, fit_y, type = "l", pch=22, lty=2, col="red", lwd=2, ylim=c(0, 900), xlim=c(0, 300),
       xlab=xlab, ylab=ylab)
  points(values[!outcomes, coeffx], values[!outcomes, coeffy], col='dodgerblue2', pch=19)
  points(values[outcomes, coeffx], values[outcomes, coeffy], col='green', pch=19)
}

prepareData <- function(dat, value_cols) {
  batch <- as.factor(substr(as.vector(dat$Study.No), 1, 2))
  
  y <- dat$Outcome
  y[y == 1] <- 'Unfav'
  y[y == 2] <- 'Fav'
  y <- as.factor(y)
  
  plt.given <- data.frame(Platelets=rep(F, nrow(dat)), row.names = rownames(dat))
  plt.given[grep('Platelets', dat$Blood.product), 'Platelets'] <- T
  
  mec <- dat$Meconium.Aspiration == 'Y'
  
  culture <- dat$Culture == 'Y'
  
  values_missing <- dat[, value_cols]
  colSd <- apply(values_missing, 2, function (x) sd(na.omit(x)))
  colMn <- apply(values_missing, 2, function (x) mean(na.omit(x)))
  values <- imputeMissing(values_missing)
  values_standard <- data.frame(t(apply(values, 1, function (x) (x - colMn) / colSd)))
  
  X <- data.frame(values, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec, culture=culture)
  Xs <- data.frame(values_standard, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec, culture=culture)
  X_missing <- data.frame(values_missing, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec, culture=culture)
  
  list(
    X=X,
    Xs=Xs,
    X_missing=X_missing
  )
}

output.dir <- getOutputDir('divyen_hie_summary_table')

fn <- file.path(data.dir, 'divyen_shah', 'cleaned_data_feb_2018.csv')
dat <- read.csv(fn)

fn.full <- file.path(data.dir, 'divyen_shah', 'cleaned_data_full_cohort_feb_2018.csv')
dat.full <- read.csv(fn.full)

biomarkers <- c(    
  'Hb',
  'Plt',
  'Neutrophil',
  'Lymphocyte',
  # 'PT',  # removed as it is highly correlated with INR
  'INR',
  # 'APTT',  # removed as it is missing in one batch
  'Urea',
  'Creatinine',
  'ALT',
  'ALP',
  'CRP'
)

peak_cols <- paste(biomarkers, 'peak', sep='.')
trough_cols <- c("Plt.trough")
all_cols <- c(peak_cols, trough_cols)

prep_dat <- prepareData(dat, all_cols)
X <- prep_dat[['X']]

prep_dat.full <- prepareData(dat.full, all_cols)
X.full <- prep_dat.full[['X']]

get_binary_summary <- function(dat, colname, value) {
  summ <- table(dat[[colname]], dat$Outcome)
  colnames(summ) <- c('unfav', 'fav')
  summ_rowsum <- rowSums(summ)
  summ_colsum <- colSums(summ)
  a <- summ[value, 'unfav']
  b <- summ[value, 'fav']
  c <- summ_colsum['unfav'] - a
  d <- summ_colsum['fav'] - b
  str_fav <- sprintf("%d/%d (%.1f)", b, summ_colsum['fav'], b / summ_colsum['fav'] * 100)
  str_unfav <- sprintf("%d/%d (%.1f)", a, summ_colsum['unfav'], a / summ_colsum['unfav'] * 100)
  ctg <- matrix(
    c(a, b, c, d), 
    2,
    2
  )
  ft <- fisher.test(ctg)
  list(
    fav=str_fav,
    unfav=str_unfav,
    fisherp=ft$p.value
  )
}


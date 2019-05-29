source('_settings.R')
source('io/output.R')
source('plotting/heatmap.R')
require(ggplot2)
require(glmnet)
require(ggfortify)
require(NMF)


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
  
  values_missing <- dat[, value_cols]
  colSd <- apply(values_missing, 2, function (x) sd(na.omit(x)))
  colMn <- apply(values_missing, 2, function (x) mean(na.omit(x)))
  values <- imputeMissing(values_missing)
  values_standard <- data.frame(t(apply(values, 1, function (x) (x - colMn) / colSd)))
  
  X <- data.frame(values, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec)
  Xs <- data.frame(values_standard, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec)
  X_missing <- data.frame(values_missing, plt.given, Unit=batch, outcome=y == 'Unfav', meconium=mec)
  
  list(
    X=X,
    Xs=Xs,
    X_missing=X_missing
  )
}

output.dir <- getOutputDir('divyen_hie')

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

# values_missing <- dat[, c(peak_cols, trough_cols)]
# colSd <- apply(values_missing, 2, function (x) sd(na.omit(x)))
# colMn <- apply(values_missing, 2, function (x) mean(na.omit(x)))
# values <- imputeMissing(values_missing)
# values_standard <- data.frame(t(apply(values, 1, function (x) (x - colMn) / colSd)))

# X <- data.frame(values, plt.given, Unit=batch, outcome=y == 'Unfav')
# Xs <- data.frame(values_standard, plt.given, Unit=batch, outcome=y == 'Unfav')

#' Do this on data without imputation of missing values
# X_missing <- data.frame(values_missing, Unit=batch)

prep_dat <- prepareData(dat, all_cols)
X <- prep_dat[['X']]
Xs <- prep_dat[['Xs']]
X_missing <- prep_dat[['X_missing']]

prep_dat.full <- prepareData(dat.full, all_cols)
X.full <- prep_dat.full[['X']]
Xs.full <- prep_dat.full[['Xs']]
X_missing.full <- prep_dat.full[['X_missing']]

#' Preamble: check for batch effects between the four units
p_anova <- list()
p_anova_post_imputation <- list()
for (col in c(peak_cols, trough_cols)) {
  mod <- aov(X_missing[,col] ~ X_missing[,'Unit'])
  p_anova[[col]] <- summary(mod)[[1]][["Pr(>F)"]][1]
  mod <- aov(X[,col] ~ X[,'Unit'])
  p_anova_post_imputation[[col]] <- summary(mod)[[1]][["Pr(>F)"]][1]
}

p_anova <- p.adjust(p_anova, method='BH')
p_anova_post_imputation <- p.adjust(p_anova_post_imputation, method='BH')

#' Preamble: PCA plot to test for batch effect
res_pca <- prcomp(X[,c(peak_cols, trough_cols)], scale=T)
theme_update(text = element_text(size=12))
autoplot(res_pca, data=X, colour='Unit', label.colour='Unit', size=2)
ggsave(file.path(output.dir, "batch_pca.png"), width=5, height = 3)

#' 1: Decide on exclusions
#' Candidates are: blood products (platelets), meconium aspiration
#' Let's use a heatmap first

# aheatmap from the NMF package
outcome_row <- ifelse(X.full$outcome, 'Unfavourable', 'Favourable')
platelet_row <- ifelse(X.full$Platelets, 'Y', 'N')
mec_row <- ifelse(X.full$meconium, 'Y', 'N')
ann_col <- data.frame(cbind(outcome_row, platelet_row, mec_row))
colnames(ann_col) <-c("Outcome group", "Platelets given?", "Meconium aspiration")
ann_colours <- list()
ann_colours[["Outcome group"]] <- c('grey80', 'black')
ann_colours[["Platelets given?"]] <- c('cadetblue2', 'darkblue')
ann_colours[["Meconium aspiration"]] <- c('springgreen', 'darkgreen')

png(filename = file.path(output.dir, "hclust_full.png"), width=7, height=4, units = "in", res = 300)
aheatmap(
  t(Xs.full[, all_cols]),
  distfun = 'pearson',
  hclustfun = 'average',
  Rowv = F,
  scale = 'none',
  annCol = ann_col,
  annColors = ann_colours,
  labCol = NA
)
dev.off()

# custom heatmap function by Obi Griffith (heatmap.3)

# outcome_row <- ifelse(X.full$outcome, 'green4', 'palegreen')
# platelet_row <- ifelse(X.full$Platelets, 'red3', 'salmon')
# mec_row <- ifelse(X.full$meconium, 'navyblue', 'lightskyblue')
# rlab <- t(cbind(outcome_row, platelet_row, mec_row))
# rownames(rlab) <- c("Outcome group", "Platelets given?", "Meconium aspiration")
# mydist = function(c) {dist(c,method="euclidean")}
# myclust = function(c) {hclust(c, method="average")}
# 
# heatmap.3(
#   t(Xs.full[, all_cols]), 
#   ColSideColors = t(rlab), 
#   ColSideColorsSize = 5,
#   keysize = 1.2,
#   KeyValueName = 'Standardized value',
#   symbreaks = F,
#   dendrogram = "col",
#   hclustfun = myclust,
#   distfun = mydist
# )

#' ANOVA for each of these
anova.mec <- list()
anova.blood <- list()

for (col in c(peak_cols, trough_cols)) {
  mod <- aov(X.full[,col] ~ X.full[,'meconium'])
  anova.mec[[col]] <- summary(mod)[[1]][["Pr(>F)"]][1]
  mod <- aov(X.full[,col] ~ X.full[,'Platelets'])
  anova.blood[[col]] <- summary(mod)[[1]][["Pr(>F)"]][1]
}

anova.mec <- p.adjust(anova.mec, method='BH')
anova.blood <- p.adjust(anova.blood, method='BH')

png(file.path(output.dir, "meconium_crp.png"), width = 3, height= 4, units = 'in', res = 300)
ggplot(X.full, aes(x=meconium, y=CRP.peak)) + geom_boxplot() + xlab('Meconium aspiration?') + ylab('CRP peak')
dev.off()

# 1:  variable selection
# We have too many variables to model them all without overfitting. Let's run lasso regression to shrink the number involved
xx <- model.matrix(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                     Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                     ALT.peak + CRP.peak + Plt.trough, data=Xs)
yy <- as.integer(Xs$outcome)

# since this is stochastic, run it several times
n.lasso <- 100
lambda.mins <- list()
lasso.coefs <- list()
for (i in seq(1, n.lasso)) {
  lasso.cv <- cv.glmnet(xx, yy, alpha=1, family="binomial", nfolds = 10)
  lambda.mins[[i]] <- lasso.cv$lambda.min
  lasso.coefs[[i]] <- coef(lasso.cv, s=lasso.cv$lambda.min)
}
# count number of times each variable appears
a <- lasso.coefs[[1]]
coef.count <- data.frame(count=rep(0, 10), row.names = unique(rownames(a))[2:11])
coef.values <- list()

for (i in seq(1, n.lasso)) {
  a <- lasso.coefs[[i]]
  ix <- a@i[a@i > 1] - 1
  coef.count[ix,] <- coef.count[ix,] + 1
  for (i in a@i[a@i > 1]) {
    rn <- rownames(a)[i + 1]
    coef.values[[rn]] <- c(coef.values[[rn]], a[i + 1])
  }
}

# only keep those biomarkers present in 1/2 of the iterations or more
coef.values <- coef.values[rownames(coef.count)[coef.count > n.lasso / 2]]

longform <- do.call(rbind, lapply(names(coef.values), function (x) {data.frame(coef=coef.values[[x]], name=x)}))
png(file.path(output.dir, 'coefficient_value_histogram.png'), width=4, height = 3.5, units='in', res = 300)
ggplot(longform, aes(coef, fill=name)) + geom_histogram(alpha=0.4, position = 'identity') + xlab('Coefficient value') + ylab('Frequency') + labs(fill='Biomarker')
dev.off()

lasso.cv <- cv.glmnet(xx, yy, alpha=1, family="binomial", nfolds = 10)
plot(lasso.cv)
coef(lasso.cv, s=lasso.cv$lambda.min)

# Compare two simple models (plus the null)
fit.null <- glm(outcome ~ 1, data=Xs, family=binomial())
fit.simple1 <- glm(outcome ~ Plt.peak + Plt.trough, data=Xs, family=binomial())
fit.simple2 <- glm(outcome ~ Plt.peak + Plt.trough + CRP.peak, data=Xs, family=binomial())

summary(fit.simple1)
summary(fit.simple2)
anova(fit.null, fit.simple1, test='Chisq')
anova(fit.null, fit.simple2, test='Chisq')
anova(fit.simple1, fit.simple2, test='Chisq')

# let's plot the decision boundary for the simplest (2 factor) model
b0 <- coef(fit.simple1)[['(Intercept)']]
b_peak <- coef(fit.simple1)[['Plt.peak']]
b_trough <- coef(fit.simple1)[['Plt.trough']]
s_peak <- colSd[['Plt.peak']]
s_trough <- colSd[['Plt.trough']]
m_peak <- colMn[['Plt.peak']]
m_trough <- colMn[['Plt.trough']]

incpt <- s_peak / b_peak * (-b0 + b_trough / s_trough * m_trough + b_peak / s_peak * m_peak)
slope <- -(s_peak * b_trough) / (s_trough * b_peak)

fit_x <- c(0, 300)
fit_y = incpt + slope * fit_x
plot(fit_x, fit_y, type = "l", pch=22, lty=2, col="red", lwd=2, ylim=c(0, 650), xlim=c(0, 300),
     xlab="Plt trough", ylab="Plt peak")
points(values$Plt.trough[y == 'Fav'], values$Plt.peak[y == 'Fav'], col='dodgerblue2', pch=19)
points(values$Plt.trough[y == 'Unfav'], values$Plt.peak[y == 'Unfav'], col='green', pch=19)
# highlight cases where platelets were given
points(values$Plt.trough[X$Platelets], values$Plt.peak[X$Platelets], col='black', pch=1)
legend("topleft", legend = c("Decision boundary", "Favourable", "Unfavourable", "Given platelets"), lty=c(2, 0, 0, 0), pch=c(-1, 19, 19, 1), col=c("red", "dodgerblue2", "green", "black"), bty='n')

# now test whether the cohort that were given platelets should be excluded
# scatterplot with linear regression fit
fit.lm <- lm(Plt.peak ~ Plt.trough, data=X[!X$Platelets,])
fit_x <- c(0, 300)
fit_y <- coef(fit.lm)[[1]] + coef(fit.lm)[[2]] * fit_x
plot.new()
plot(fit_x, fit_y, type="l", pch=22, lty=2, col="red", lwd=2, xlim=c(0, 300), ylim=c(0, 650), xlab="Plt trough", ylab="Plt peak")
points(values$Plt.trough[y == 'Fav'], values$Plt.peak[y == 'Fav'], col='dodgerblue2', pch=19)
points(values$Plt.trough[y == 'Unfav'], values$Plt.peak[y == 'Unfav'], col='green', pch=19)
# highlight cases where platelets were given
points(values$Plt.trough[X$Platelets], values$Plt.peak[X$Platelets], col='black', pch=1)
legend("topleft", legend = c("Favourable", "Unfavourable", "Given platelets"), lty=c(0, 0, 0), pch=c(19, 19, 1), col=c("dodgerblue2", "green", "black"), bty='n')

# box and scatter
plot.new()
boxplot(Plt.trough~Platelets, data=X, lwd=2, xlab="Platelets given?", ylab="Plt trough value", ylim=c(0, 300))
stripchart(Plt.trough ~ Platelets, data=X, vertical=T, method="jitter", add=T, pch=20, col='red')

plot.new()
boxplot(Plt.peak~Platelets, data=X, lwd=2, xlab="Platelets given?", ylab="Plt peak value")
stripchart(Plt.peak ~ Platelets, data=X, vertical=T, method="jitter", add=T, pch=20, col='darkgrey')

fit.no_platelets <- glm(outcome ~ Plt.peak + Plt.trough, data=Xs[!Xs$Platelets,], family=binomial())
plotDecisionBoundary(
  fit.no_platelets, 
  'Plt.trough', 
  'Plt.peak', 
  X[!X$Platelets,], 
  Xs[!Xs$Platelets, 'outcome'],
  xlab="Plt trough", ylab="Plt peak"
)
legend("topleft", legend = c("Favourable", "Unfavourable"), lty=c(0, 0), pch=c(19, 19), col=c("dodgerblue2", "green"), bty='n')

fit.plt1 <- glm(outcome ~ Plt.peak + Plt.trough + Platelets, data=Xs, family=binomial())
fit.plt2 <- glm(outcome ~ Plt.peak + Plt.trough + CRP.peak + Platelets, data=Xs, family=binomial())


# Initial model: all biomarkers, dept = outcome
# Not sure this is remotely valid!!

fit.all <- glm(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                 Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                 ALT.peak + CRP.peak + Plt.trough,
               data=Xs, family = binomial(link = "logit"))
summary(fit.all)

# we can't trust the P values from the previous approach as they are based on the Wald test, so we should run a likelihood based test now.
# https://stats.stackexchange.com/questions/315684/i-think-my-logistic-model-is-overfitted-even-with-lasso-r-gives-me-a-perfect-se?rq=1
anova(fit.all, test='Chi')
anova(fit.null, fit.all, test="Chisq")

# coefficients
odds <- exp(coef(fit.all))
odds.ci <- exp(confint.default(fit.all))
odds_table <- data.frame(odds=odds, odds.ci)

# 3: Now let's load the full cohort
fn.full <- file.path(data.dir, 'divyen_shah', 'cleaned_data_full_cohort_feb_2018.csv')
dat.full <- read.csv(fn.full)
batch.full <- as.factor(substr(as.vector(dat.full$Study.No), 1, 2))

y.full <- dat.full$Outcome
y.full[y.full == 1] <- 'Unfav'
y.full[y.full == 2] <- 'Fav'
y.full <- as.factor(y.full)

values.full <- dat.full[, c(peak_cols, trough_cols)]

colSd <- apply(values.full, 2, function (x) sd(na.omit(x)))
colMn <- apply(values.full, 2, function (x) mean(na.omit(x)))
values.full <- imputeMissing(values.full)
values_standard.full <- data.frame(t(apply(values.full, 1, function (x) (x - colMn) / colSd)))

mec.or.cult <- dat.full$Meconium.Aspiration == 'Y'
mec.or.cult[dat.full$Culture == 'Y'] <- T

X.full <- data.frame(values.full, batch=batch.full, mec.or.cult=mec.or.cult, outcome=y.full == 'Unfav')
Xs.full <- data.frame(values_standard.full, batch=batch.full, mec.or.cult=mec.or.cult, outcome=y.full == 'Unfav')

# what is the outcome when we include all data in the basic model?
fit.null.full <- glm(outcome ~ 1, data=Xs.full, family=binomial())
fit.full1 <- glm(outcome ~ Plt.peak + Plt.trough, data=Xs.full, family=binomial())

predictivePerformance(fit.full1, y.full == 'Unfav')

# let's plot the decision boundary for this model
b0 <- coef(fit.full1)[['(Intercept)']]
b_peak <- coef(fit.full1)[['Plt.peak']]
b_trough <- coef(fit.full1)[['Plt.trough']]
s_peak <- colSd[['Plt.peak']]
s_trough <- colSd[['Plt.trough']]
m_peak <- colMn[['Plt.peak']]
m_trough <- colMn[['Plt.trough']]

incpt <- s_peak / b_peak * (-b0 + b_trough / s_trough * m_trough + b_peak / s_peak * m_peak)
slope <- -(s_peak * b_trough) / (s_trough * b_peak)

fit_x <- c(0, 300)
fit_y = incpt + slope * fit_x
plot(fit_x, fit_y, type = "l", pch=22, lty=2, col="red", lwd=2, ylim=c(0, 900), xlim=c(0, 300),
     xlab="Plt trough", ylab="Plt peak")
points(values.full$Plt.trough[y.full == 'Fav'], values.full$Plt.peak[y.full == 'Fav'], col='dodgerblue2', pch=19)
points(values.full$Plt.trough[y.full == 'Unfav'], values.full$Plt.peak[y.full == 'Unfav'], col='green', pch=19)
# highlight cases where platelets were given
points(values.full$Plt.trough[X.full$mec.or.cult], values.full$Plt.peak[X.full$mec.or.cult], col='black', pch=1)
legend("topright", legend = c("Decision boundary", "Favourable", "Unfavourable", "Meconium/culture"), lty=c(2, 0, 0, 0), pch=c(-1, 19, 19, 1), col=c("red", "dodgerblue2", "green", "black"), bty='n')

  

fit.mec <- glm(mec.or.cult ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                 Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                 ALT.peak + CRP.peak + Plt.trough, data = X.full)
fit.mec.null <- glm(mec.or.cult ~ 1, data = X.full)

anova(fit.mec.null, fit.mec, test='Chisq')

odds.mec <- exp(coef(fit.mec))
odds.ci.mec <- exp(confint.default(fit.mec))

# conclusion: yes, the meconium or culture status is predictable from the biomarkers.
# but can we still use these?

fit.outcome.full <- glm(outcome ~ mec.or.cult + Hb.peak + Plt.peak + Neutrophil.peak + 
                          Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                          ALT.peak + CRP.peak + Plt.trough, data = X.full)

fit.outcome.full.hier <- glm(outcome ~ mec.or.cult * (Hb.peak + Plt.peak + Neutrophil.peak + 
                          Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                          ALT.peak + CRP.peak + Plt.trough), data = X.full)

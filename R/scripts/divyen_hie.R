source('_settings.R')
require(ggplot2)
require(glmnet)

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


fn <- file.path(data.dir, 'divyen_shah', 'cleaned_data_feb_2018.csv')
dat <- read.csv(fn)
batch <- as.factor(substr(as.vector(dat$Study.No), 1, 2))

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

y <- dat$Outcome
y[y == 1] <- 'Unfav'
y[y == 2] <- 'Fav'
y <- as.factor(y)

values <- dat[, c(peak_cols, trough_cols)]
colSd <- apply(values, 2, function (x) sd(na.omit(x)))
colMn <- apply(values, 2, function (x) mean(na.omit(x)))
values <- imputeMissing(values)
values_standard <- data.frame(t(apply(values, 1, function (x) (x - colMn) / colSd)))

plt.given <- data.frame(Platelets=rep(F, nrow(dat)), row.names = rownames(dat))
plt.given[grep('Platelets', dat$Blood.product), 'Platelets'] <- T

X <- data.frame(values, plt.given, batch=batch, outcome=y == 'Unfav')
Xs <- data.frame(values_standard, plt.given, batch=batch, outcome=y == 'Unfav')

# 1:  variable selection
# We have too many variables to model them all without overfitting. Let's run lasso regression to shrink the  number involved
xx <- model.matrix(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                     Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                     ALT.peak + CRP.peak + Plt.trough, data=Xs)
yy <- as.integer(y == 'Unfav')

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
m_peak <- colM[['Plt.peak']]
m_trough <- colM[['Plt.trough']]

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

# TODO: code up the performance evaluation

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
m_peak <- colM[['Plt.peak']]
m_trough <- colM[['Plt.trough']]

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
                 ALT.peak + CRP.peak + APTT.peak + Plt.trough, data = X.full)
fit.mec.null <- glm(mec.or.cult ~ 1, data = X.full)

anova(fit.mec.null, fit.mec, test='Chisq')

odds.mec <- exp(coef(fit.mec))
odds.ci.mec <- exp(confint.default(fit.mec))

# conclusion: yes, the meconium or culture status is predictable from the biomarkers.
# but can we still use these?

fit.outcome.full <- glm(outcome ~ mec.or.cult + Hb.peak + Plt.peak + Neutrophil.peak + 
                          Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                          ALT.peak + CRP.peak + APTT.peak + Plt.trough, data = X.full)

fit.outcome.full.hier <- glm(outcome ~ mec.or.cult * (Hb.peak + Plt.peak + Neutrophil.peak + 
                          Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                          ALT.peak + CRP.peak + APTT.peak + Plt.trough), data = X.full)

source('_settings.R')
require(ggplot2)
require(GGally)
require(nnet)
require(pscl)

nanMedian <- function(x) {
  x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
  x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}


imputeMissing <- function(data) {
  data.frame(apply(data, 2, nanMedian))
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
values <- imputeMissing(values)

plt.given <- data.frame(Platelets=rep(F, nrow(dat)), row.names = rownames(dat))
plt.given[grep('Platelets', dat$Blood.product), 'Platelets'] <- T

X <- data.frame(values, plt.given, batch=batch, outcome=y == 'Unfav')

# Initial model: all biomarkers, dept = outcome
fit.glm <- glm(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                 Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                 ALT.peak + CRP.peak + Plt.trough,
               data=X, family = binomial(link = "logit"))

summary(fit.glm)
# we can't trust the P values from the previous approach as they are based on the Wald test, so we should run a likelihood based test now.
# https://stats.stackexchange.com/questions/315684/i-think-my-logistic-model-is-overfitted-even-with-lasso-r-gives-me-a-perfect-se?rq=1
anova(fit.glm, test='Chi')

fit.null <- glm(outcome ~ 1, data=X, family=binomial())

anova(fit.null, fit.glm, test="Chisq")

# fit a reduced model to improve overfitting
fit.reduced <- glm(outcome ~ Plt.peak + Plt.trough, data=X, family=binomial())


# coefficients
odds <- exp(coef(fit.glm))
odds.ci <- exp(confint.default(fit.glm))
odds_table <- data.frame(odds=odds, odds.ci)

# 2: Is the platelet product predictive of outcome?
fit.glm.plt <- glm(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                 Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                 ALT.peak + CRP.peak + Plt.trough + Platelets,
               data=X, family=binomial())
summary(fit.glm.plt)

anova(fit.glm, fit.glm.plt, test="Chisq")

# 3: Now let's load the full cohort
fn.full <- file.path(data.dir, 'divyen_shah', 'cleaned_data_full_cohort_feb_2018.csv')
dat.full <- read.csv(fn.full)

y.full <- dat.full$Outcome
y.full[y.full == 1] <- 'Unfav'
y.full[y.full == 2] <- 'Fav'
y.full <- as.factor(y.full)

values.full <- dat.full[, c(peak_cols, trough_cols)]
values.full <- imputeMissing(values.full)
mec.or.cult <- dat.full$Meconium.Aspiration == 'Y'
mec.or.cult[dat.full$Culture == 'Y'] <- T

X.full <- data.frame(values.full, mec.or.cult=mec.or.cult, outcome=y.full == 'Unfav')

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

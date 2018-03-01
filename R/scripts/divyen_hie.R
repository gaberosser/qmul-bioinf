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
# 'PT',
'INR',
'APTT',
'Urea',
'Creatinine',
'ALT',
'ALP',
'CRP'
)

biomarkers.shortened <- c(
  'Hb', 'Plt', 'Nphil', 'Lcyte', 'PT', 'INR', 'APTT', 'Urea', 'Cnine', 'ALT', 'ALP', 'CRP'
)

peak_cols <- paste(biomarkers, 'peak', sep='.')
trough_cols <- c("Plt.trough")

y <- dat$Outcome
y[y == 1] <- 'Unfav'
y[y == 2] <- 'Fav'
y <- as.factor(y)

values <- dat[, c(peak_cols, trough_cols)]
values <- imputeMissing(values)

X <- data.frame(values, outcome=y == 'Unfav')

fit.glm <- glm(outcome ~ Hb.peak + Plt.peak + Neutrophil.peak + 
                 Lymphocyte.peak + INR.peak + Urea.peak + Creatinine.peak + 
                 ALT.peak + CRP.peak + APTT.peak + Plt.trough,
               data=X)

summary(fit.glm)
fit.null <- glm(outcome ~ 1, data=X)

anova(fit.null, fit.glm, test="Chisq")

# coefficients
odds <- exp(coef(fit.glm))
odds.ci <- exp(confint.default(fit.glm))

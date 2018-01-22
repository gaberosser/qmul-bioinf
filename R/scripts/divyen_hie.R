source('_settings.R')
require(ggplot2)
require(GGally)
require(nnet)
require(pscl)

fn <- file.path(data.dir, 'divyen_shah', 'cleaned_data_sep_2017.csv')
dat <- read.csv(fn)
batch <- as.factor(substr(as.vector(dat$Study.No), 1, 2))

biomarkers <- c(    
'Hb',
'Plt',
'Neutrophil',
'Lymphocyte',
'PT',
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
colnames(values) <- c(paste(biomarkers.shortened, 'hi', sep='.'), 'Plt.lo')

peak_age_cols <- paste(biomarkers, 'peak', 'age', sep='.')
trough_age_cols <- c("Plt.trough.age")

obs_ages = dat[, c(peak_age_cols, trough_age_cols)]
colnames(obs_ages) <- c(paste(biomarkers.shortened, 'hi', sep='.'), 'Plt.lo')

## Some variables are systematically absent from some batches
# Shall we exclude these?
print(summary(batch))
print(summary(batch[rowSums(is.na(values)) > 0]))
print("The following variables are missing from the SU studies")
print(colnames(values)[colSums(is.na(values[batch == 'SU',])) > 0])


# check for significant batch effects using a multinomial regression model
# TODO: how do we test the significance of this model fit?
X <- data.frame(values)
X$batch <- batch
# we don't test for APTT.hi, PT.hi as these are missing in all SU
# also remove INR.hi and Plt.lo since these are missing in some cases
frm <- batch ~ Hb.hi + Plt.hi + Nphil.hi + Lcyte.hi + Urea.hi + Cnine.hi + ALT.hi + ALP.hi + CRP.hi
test <- multinom(frm, X)
print(pR2(test))
paste("Number of batches correctly predicted is", sum(predict(test) == batch), "/", length(batch), collapse = " ")

X <- data.frame(values)
X$outcome <- y


# ggpairs plot
# TODO: reduce fontsize
png("value_scatterplot.png", width=1280, height=1280, res=150, pointsize=8)
ggpairs(X, columns=1:ncol(values), ggplot2::aes(colour=outcome))
dev.off()

# matrix of scatterplots showing observation time - value for each biomarker
png("value_vs_observation_time.png", width=1024, height=1024, res=150)
par(mfrow=c(4, 4), mai=c(0.6, 0.6, 0.1, 0.1))
for (i in seq(ncol(values))) {
  plot(obs_ages[,i], values[,i], xlab=colnames(obs_ages)[i], ylab=colnames(values)[i], col=ifelse(y=='Fav', "blue", "red"), pch=16)
}
plot(0, 1, col='white', frame.plot = F, xlab="", ylab="", xaxt="n", yaxt="n")
legend("center", c("Favourable", "Unfavourable"), col=c("blue", "red"), pch=c(16, 16))
dev.off()

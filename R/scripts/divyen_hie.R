source('_settings.R')


fn <- file.path(data.dir, 'divyen_shah', 'cleaned_data_sep_2017.csv')
dat <- read.csv(fn)

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
peak_cols <- paste(biomarkers, 'peak', sep='.')
trough_cols <- c("Plt.trough")

y <- dat$Outcome
y[y == 1] <- 'Unfav'
y[y == 2] <- 'Fav'
y <- as.factor(y)

X <- dat[, c(peak_cols, trough_cols)]

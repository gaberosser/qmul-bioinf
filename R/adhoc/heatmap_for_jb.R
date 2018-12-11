library(reshape2)
library(ggplot2)

# load the data
ipa_fn <- "/home/gabriel/hGIC_project/archived/2018-11-16/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways/full_de_syngeneic_z.txt"
ipa_data <- read.csv(ipa_fn, header=T, sep='\t', skip = 2, na.strings = 'N/A')

# let's only use a few for easy plotting
# in practice, prioiritise based on some meaningful metric
ipa_data <- ipa_data[1:50,]

# rename columns (they're messed up in my case)
ix <- c(1, grep("X([0-9]{3}).*", colnames(ipa_data)))
ipa_data <- ipa_data[,ix]
pids <- gsub("X([0-9]{3}).*", "\\1", colnames(ipa_data)[-1])
colnames(ipa_data)[-1] <- pids
colnames(ipa_data)[1] <- 'pathway'

# melt
df <- melt(ipa_data, id.vars = 'pathway')
colnames(df)[2:3] <- c("patient", "z")

# plot
ggplot(df, aes(patient, pathway)) + 
  geom_tile(aes(fill=z), color='white') + 
  scale_fill_gradient(low = "red", high = "steelblue") + 
  xlab("Patient") + 
  ylab("Pathway") + 
  theme(legend.title=element_text(size=12)) + 
  labs(fill="Z")

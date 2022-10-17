#### Plot to explore gene count data from the GBE project using PCA! Yeehaw
### Amanda Neves, May 10, 2022

############################################################################

library(ggplot2)
library(DESeq2)
library(ggpubr)
library(viridis)

# Load in data

sample_info <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/plot_sample_info.csv")

# remove the outlier to match the counts data
sample_info <- sample_info[-5,]
# create allele type variable for plotting
allele_type <- substr(sample_info$allele, start = 1, stop = 2)

salmon_counts <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/salmon_counts_raw_all.csv", row.names = 1)

salmon_matrix <- t(salmon_counts) # transpose so that rows = samples and columns = variables

# apply the variance stabilizing transformation
vst_salmon <- varianceStabilizingTransformation(round(t(salmon_matrix)), blind = FALSE)

salmon_pca <- prcomp(t(vst_salmon))
summary(salmon_pca)
salmon_vals <- salmon_pca$x

# create plotting dataframe
s_sub_plot <- data.frame(background = sample_info$background, allele = sample_info$allele, allele_type = allele_type, perturbation = sample_info$perturbation, PC1 = salmon_vals[,1], PC2 = salmon_vals[,2], PC3 = salmon_vals[,3], PC4 = salmon_vals[,4])

pca_plot <- ggplot(data = s_sub_plot, aes(x = PC1, y = PC2, colour = perturbation, shape = background)) + geom_point(size = 3) + xlab("PC1 (25.4%)") + ylab("PC2 (21.8%)") + theme_classic() + scale_shape(solid = TRUE) + theme(legend.position="bottom") + scale_colour_viridis()

# save to output

ggsave("gbe_pca_plot.png", width = 8, height = 6)

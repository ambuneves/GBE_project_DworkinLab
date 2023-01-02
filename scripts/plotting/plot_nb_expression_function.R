# Function to plot expression by perturbation using the salmon nb estimates (for the regression line)
# Amanda Neves Oct 17, 2022

library(lme4)
library(glmmTMB)
library(emmeans)
library(parallel)
library(stringr)
library(effects)
library(ggplot2)

# Read in data needed for the function to work

perturbation <- c(7.434782609, 7.434782609,
                  3.541666667, 3.541666667,
                  4.564102564, 4.564102564,
                  1.733333333, 1.733333333,
                  3.692307692, 3.692307692,
                  5.846153846, 5.846153846,
                  10, 10,
                  9.708333333, 9.708333333,
                  6.739130435, 6.739130435,
                  8.275862069, 8.275862069,
                  1.217391304, 1.217391304,
                  3.695652174, 3.695652174,
                  9.714285714, 9.714285714,
                  10, 10)

perturbation <- abs(perturbation - max(perturbation))

setwd("/Users/amandan/Desktop/Dworkin/background_effects/data")
# read in the coutns without the outlier (ORE sdETX4 and without the hybrid samples)
salmon_counts <- read.csv("salmon_counts_raw_no_outlier.csv", row.names = 1)
# read in the normalization factors calculated using DESeq2 (normalizationFactors)
normFactors <- read.csv("salmon_nb_normfactors.csv", row.names = 1)
rna.design <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/rna_design.csv")

gene_names <- read.delim("FlyGeneDictionary.txt",
                         header = TRUE,
                         sep = "\t",
                         dec = ".") # for adding gene name to plot

plot_PE_gene_nb <- function(gene){
  g <- gene
  the_counts <- as.numeric(salmon_counts[g,])
  the_counts <- 1 + the_counts
  the_counts <- round(the_counts, 0)
  the_normFactors <- as.numeric(normFactors[g,])
  
  # Run the nb model using the norm factors as the offset
  nbinom.model <- glmmTMB(the_counts ~ background + 
                            perturbation + 
                            background:perturbation + 
                            (1 | lane),
                          data = rna.design,
                          offset = log(the_normFactors),
                          family = nbinom2(),
                          #control = glmmTMBControl(parallel = nt)
  )
  
  plotting_df <- as.data.frame(predictorEffect("perturbation", 
                                               nbinom.model, 
                                               xlevels = 13))
  
  plotting_df[,3:6] <- lapply(plotting_df[,3:6], log )
  base_allele <- c("sd[1]", "sd[1]", "sd[E3]", "sd[E3]", "sd[ETX4]", "sd[ETX4]", "vg[1]", "vg[1]", "vg[21-3]", "vg[21-3]", "vg[2a33]", "vg[2a33]", "wild type", "wild type", "sd[1]", "sd[1]", "sd[E3]", "sd[E3]", "sd[ETX4]", "sd[ETX4]", "vg[1]", "vg[1]", "vg[21-3]", "vg[21-3]", "vg[2a33]", "vg[2a33]", "wild type", "wild type")
  
  
  # and i need to add the perturbation of each... 
  
  base_perturbation <- perturbation
  
  # for plotting 
  base_background <- c(rep("ORE", times = 14), 
                       rep("SAM", times = 14))
  
  # add 0.001 to expression values to avoid complete separation issue
  point_df <- data.frame(background = base_background,
                         allele = base_allele,
                         expression = log(the_counts/the_normFactors),
                         expression_old = log(as.numeric(the_counts) + 0.001),
                         perturbation = base_perturbation)
  
  point_df$allele <- factor(point_df$allele, levels = rev(c("wild type", "sd[1]", "vg[2a33]", "sd[ETX4]", "sd[E3]", "vg[21-3]", "vg[1]")))
  
  plot_title <- ifelse(gene %in% gene_names$X.submitted_item,
                       gene_names$current_symbol[match(gene, gene_names$X.submitted_item)],
                       gene)
  
  #ylab_plot <- ylab(expression(log[2](CPM)))
  
  ggplot() + geom_line(data = plotting_df, aes(x = perturbation, y = fit, colour = background)) + geom_ribbon(data=plotting_df,aes(ymin=lower,ymax=upper,fill = background, x = perturbation),alpha=0.3) + geom_point(data = point_df, size = 4, alpha = 0.65, aes(x = perturbation, y = expression, shape = allele, colour = background)) + scale_shape_manual(values=1:7) +  ylab(expression(log[2](CPM))) + xlab("Perturbation") + theme_classic() + ggtitle(substitute(italic(x), list(x=plot_title)))
  
}

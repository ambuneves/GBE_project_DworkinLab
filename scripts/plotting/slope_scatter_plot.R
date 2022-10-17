#### Amanda Neves October 11, 2022 ######

# libraries
library(ggplot2)
library(viridisLite)
library(ggrepel)

# read in data for all of the slopes (salmon counts, nb model)
slopes_df <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/all_slopes_nb.csv")
colnames(slopes_df) <- c("gene_id", "ore_slope", "sam_slope")

# read in estimate data for classifying genes
# salmon negative binomial
salmon_nb_df <- get(load(file = "/Users/amandan/Desktop/Dworkin/background_effects/data/anova_pvalues_salmon_nb.Rdata")) # load in data
salmon_nb_df <- as.data.frame(salmon_nb_df)
dim(salmon_nb_df)

# function to categorize by p values
cat_pvals <- function(dataset){
  pval_df <- dataset
  
  # background
  bg_vals <- p.adjust(pval_df$background, method = "BH")
  #hist(bg_vals)
  #which(bg_vals < 0.1)
  pval_df$b.adj <- bg_vals
  
  # perturbation
  pb_vals <- p.adjust(pval_df$perturbation, method = "fdr")
  #hist(pb_vals)
  #which(pb_vals < 0.1)
  pval_df$p.adj <- pb_vals
  
  # background:perturbation
  bg_pb_vals <- p.adjust(pval_df$`background:perturbation`,  method = "fdr")
  #hist(bg_pb_vals)
  #which(bg_pb_vals < 0.1)
  pval_df$bp.adj <- bg_pb_vals
  
  
  #### Categorize genes ####
  
  # Perturbation Effects only 
  p1 <- rownames(pval_df)[which(pval_df$bp.adj > 0.1 & pval_df$b.adj < 0.1 & pval_df$p.adj < 0.1)]
  
  p2 <- rownames(pval_df)[which(pval_df$bp.adj > 0.1 & pval_df$b.adj > 0.1 & pval_df$p.adj < 0.1)]
  
  p <- c(p1,p2)
  
  # Interaction effects only
  
  pb1 <- rownames(pval_df)[which(pval_df$bp.adj < 0.1 & pval_df$b.adj < 0.1 & pval_df$p.adj > 0.1)]
  
  pb2 <- rownames(pval_df)[which(pval_df$bp.adj < 0.1 & pval_df$b.adj > 0.1 & pval_df$p.adj > 0.1)]
  
  pb <- c(pb1,pb2)
  
  # Interaction and perturbation
  
  p_pb1 <- rownames(pval_df)[which(pval_df$bp.adj < 0.1 & pval_df$b.adj > 0.1 & pval_df$p.adj < 0.1)]
  
  p_pb2 <- rownames(pval_df)[which(pval_df$bp.adj < 0.1 & pval_df$b.adj < 0.1 & pval_df$p.adj < 0.1)]
  
  p_pb <- c(p_pb1, p_pb2)
  
  outlist <- list(p, pb, p_pb)
  names(outlist) <- c("perturbation_only", "interaction_only", "perturbation_and_interaction")
  
  outlist
}

# get the categorizations from all of the genes
all_classes <- cat_pvals(salmon_nb_df)

## create dataframes for each category I am interested and add variables for plotting (category and alpha)

# not significant genes
sig_genes <- unlist(all_classes, use.names = FALSE)
length(sig_genes)
length(unique(sig_genes)) # just checking

unsig_df <- slopes_df[which(slopes_df$gene_id %in% sig_genes == FALSE),]
unsig_df$category <- "Not significant"
unsig_df$alpha = 0.2

# other genes
# make one for each and then combine together?
df_1 <- slopes_df[slopes_df$gene_id %in% all_classes[["interaction_only"]],]
df_1$category <- "Interaction only"
dim(df_1)
df_1$alpha <- 0.3

df_2 <- slopes_df[slopes_df$gene_id %in% all_classes[["perturbation_and_interaction"]],]
df_2$category <- "Perturbation and interaction"
dim(df_2)
df_2$alpha <- 0.3

df_pert <- slopes_df[slopes_df$gene_id %in% all_classes[["perturbation_only"]],]
df_pert$category <- "Perturbation only"
dim(df_pert)
df_pert$alpha <- 0.3


# combine into one dataframe for plotting
plot_df <- rbind.data.frame(unsig_df, df_1, df_2, df_pert)

# same alpha for all
# relevel category
plot_df$category <- factor(plot_df$category, levels = c("Not significant", "Perturbation only", "Interaction only", "Perturbation and interaction"))

scatter_plot <- ggplot() + 
geom_point(data = plot_df, aes(x = ore_slope, y = sam_slope, shape = category, colour = category), alpha = 0.4, size = 3) + 
geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) + theme_classic() + 
theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=10)) + 
xlab("ORE effect") + ylab("SAM effect") + 
xlim(-1,1) + ylim(-1,1) + 
scale_shape_manual(values = c(3, 16, 15,17)) + 
scale_color_manual(values = c("grey", "#440154FF", "#FDE725FF", "#31688EFF"))

# save to output
ggsave("slope_categories_scatterplot.png", width = 7, height = 7)

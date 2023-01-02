#### Amanda Neves October 11, 2022 ######

# libraries
library(ggplot2)
library(viridisLite)
library(ggrepel)
library(ggarrange)
library(gridExtra)

# read in data for all of the slopes (salmon counts, nb model)
slopes_df <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/all_slopes_nb.csv")
colnames(slopes_df) <- c("gene_id", "ore_slope", "sam_slope")
dim(slopes_df)

# read in estimate data for classifying genes
# salmon negative binomial
salmon_nb_df <- get(load(file = "/Users/amandan/Desktop/Dworkin/background_effects/data/anova_pvalues_salmon_nb_add1.Rdata")) # load in data
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

# save this one to use as a control for vector correlations
#write.csv(x = all_classes[["perturbation_only"]], file = "/Users/amandan/Desktop/Dworkin/background_effects/data/perturbation_only_genes.csv")

#### Making plotting dataframes
# Since my idea is to make two separate plotting DFs so that I am able to set two different opacities, I will need to make one for the genes not significant for any term, and one for the genes in the all_classes object

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

#sig_df <- rbind.data.frame(df_1, df_2)
plot_df <- rbind.data.frame(unsig_df, df_1, df_2, df_pert)

unsig_pert <- rbind.data.frame(unsig_df, df_pert)

# same alpha for all

plot_df$category <- factor(plot_df$category, levels = c("Not significant", "Perturbation only", "Interaction only", "Perturbation and interaction"))

gene_name <- list()
# add a column that contains the gene name instead of gene_id
gene_names <- read.delim("/Users/amandan/Desktop/Dworkin/background_effects/data/FlyGeneDictionary.txt",
                         header = TRUE,
                         sep = "\t",
                         dec = ".") # for adding gene name to plot


plot_df$gene_name <- gene_names$current_symbol[match(plot_df$gene_id, gene_names$validated_id)]

# maually adjust some missing values
plot_df$gene_name[which(plot_df$gene_id == "FBgn0032538")] <- "Vajk2"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0040813")] <- "Nplp2"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0034480")] <- "CG16898"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0030157")] <- "CG1468"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0034290")] <- "CG5773"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0035186")] <- "CG13912"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0085320")] <- "CG34291"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0030294")] <- "Pa1"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0033821")] <- "CG10799"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0050196")] <- "CG30196"

plot_df$gene_name[which(plot_df$gene_id == "FBgn0053333")] <- "CG33333"


plot_df$category <- factor(plot_df$category, levels = c("Not significant", "Interaction only", "Perturbation only", "Perturbation and interaction"))

# points in one dataframe

# Create scatterplot with all slopes, no labels 

big_plot <- ggplot(data = plot_df, aes(x = ore_slope, y = sam_slope, shape = category, colour = category, label = gene_name)) + 
  geom_point(size = 5, aes(alpha = category)) + 
  scale_alpha_manual(values = c(0.5,0.9,0.1,0.9)) +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) + 
  theme_classic() + 
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14)) + 
  xlab("ORE effect") + 
  ylab("SAM effect") + 
  xlim(-1,1) + 
  ylim(-1,1) + 
  scale_shape_manual(values = c(3, 16, 15,17)) + 
  scale_color_manual(values = c("grey", "#FDE725FF", "#440154FF", "#31688EFF")) + 
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), alpha = guide_legend(nrow = 2, byrow = TRUE))

big_plot

ggsave("slope_categories_scatterplot.png", width = 7, height = 7, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

### I have looked over the plots and edited the all_interesting_genes .csv file to demarcate which genes are interesting. I would like to now plot the scatter plot but with these genes labeled

gene_table <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/GBE_candidate_genes.csv")

# get the FBgns and fix the NAs

gene_table$gene_id <- gene_names$validated_id[match(gene_table$gene_name, gene_names$current_symbol)]

# interesting but noisy
goi_label_0 <- gene_table$gene_name[gene_table$id == 0]

big_plot + geom_label_repel(data = subset(plot_df, gene_name%in%goi_label_0), size = 5, box.padding = 0.5, show.legend=FALSE, nudge_y = 0.2, force = 10, fontface = "italic")

ggsave("slope_categories_scatterplot_w_labels_0.png", width = 7, height = 7, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")


# interesting 
goi_label_1 <- gene_table$gene_name[gene_table$id == 1]

big_plot + geom_label_repel(data = subset(plot_df, gene_name%in%goi_label_1), size = 5, box.padding = 0.5, show.legend=FALSE, force = 50, fontface = "italic")

ggsave("slope_categories_scatterplot_w_labels_1.png", width = 7, height = 7, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

# Final genes to reccomend follow-up
# interesting but noisy
goi_label_all <- gene_table$gene_name[gene_table$final_id == 1]
goi_label_all <- goi_label_all[!is.na(goi_label_all)]

goi_label_all <- goi_label_all[-c(1,3,5,6,7,8,9)]

# remove a couple I decided I actually don't want

big_plot + geom_label_repel(data = subset(plot_df, gene_name%in%goi_label_all), size = 5, box.padding = 0.5, show.legend=FALSE, nudge_y = 2, nudge_x = 0.5, force = 30, fontface = "italic") 

ggsave("slope_categories_scatterplot_w_labels_final.png", width = 7, height = 7, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

#######

# Next, make a cool plot with all these significant genes.

rec_gene_ids <- gene_names$validated_id[match(goi_label_all, gene_names$current_symbol)]

match(rec_gene_ids, gene_names$validated_id)

add_these <- data.frame(X.submitted_item = c("FBgn0032538", "FBgn0034480", "FBgn0035186"), validated_id = c("FBgn0032538", "FBgn0034480", "FBgn0035186"), current_symbol = c("Vajk2", "CG16898", "CG13912"))

gene_names <- rbind.data.frame(gene_names, add_these)

gene_1 <- plot_PE_gene_nb(rec_gene_ids[1]) 
gene_2 <- plot_PE_gene_nb(rec_gene_ids[2])
gene_3 <- plot_PE_gene_nb(rec_gene_ids[3])
gene_4 <- plot_PE_gene_nb(rec_gene_ids[4])
gene_5 <- plot_PE_gene_nb(rec_gene_ids[5])

ggarrange(gene_1,gene_2,gene_3,gene_4,gene_5, ncol=2, nrow=3, common.legend = TRUE, legend="bottom", labels = "AUTO") + bgcolor("white")

ggsave("all_rec_genes.png", width = 10, height = 10, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

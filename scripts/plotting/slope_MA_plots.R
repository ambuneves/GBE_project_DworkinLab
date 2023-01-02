### GBE MA Plot type thing, one for SAM and one for ORE
### the x-axis will be normalized wild-type expression, and the
### y-axis will be the slope of a gene

#########

library(ggplot2)

# Let's see if this works!

setwd("/Users/amandan/Desktop/Dworkin/background_effects/data")
# Load in gene expression data... normalized, right?
salmon_counts <- read.csv("salmon_counts_raw_no_outlier.csv", row.names = 1)
# read in the normalization factors calculated using DESeq2 (normalizationFactors)
normFactors <- read.csv("salmon_nb_normfactors.csv", row.names = 1)

# the_counts <- 1 + the_counts
# log(the_counts/the_normFactors)

up_counts <- round(salmon_counts + 1)
norm_counts <- log(up_counts/normFactors)

#only want wildtype expression
norm_counts <- norm_counts[,grep("w", colnames(norm_counts))]

# Load in slope data
slopes_df <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/all_slopes_nb.csv")
colnames(slopes_df) <- c("gene_id", "ore_slope", "sam_slope")
dim(slopes_df)


#### Plot???
norm_counts_sub <- norm_counts[match(slopes_df$gene_id, rownames(norm_counts)),]

# average the ORE and SAM counts
norm_counts_sub$ore_avg <- rowMeans(norm_counts_sub[c("ORE_w_r1_ATCACG_L001_quant","ORE_w_r2_GTTTCG_L002_quant")])

norm_counts_sub$sam_avg <- rowMeans(norm_counts_sub[c("SAM_w_r1_ACTTGA_L003_quant","SAM_w_r2_GAGTGG_L004_quant")])


# ORE df
ore_df <- data.frame(expression = norm_counts_sub$ore_avg, slope = slopes_df$ore_slope)

# SAM df
sam_df <- data.frame(expression = norm_counts_sub$sam_avg, slope = slopes_df$sam_slope)

# plots
ggplot(data = ore_df, aes(x = expression, y = slope))+ 
  geom_point(colour = "black", alpha = 0.3) + 
  theme_classic() + 
  xlab("Mean wildtype normalized expression") + 
  ylab("ORE Slope") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

ggsave("ore_MA.png", 
       width = 8, 
       height = 5, 
       path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

ggplot(data = sam_df, aes(x = expression, y = slope)) + 
  geom_point(colour = "black", alpha = 0.3) + 
  theme_classic() + 
  xlab("Mean wildtype normalized expression") + 
  ylab("SAM Slope") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

ggsave("sam_MA.png", 
       width = 8, 
       height = 5, 
       path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")


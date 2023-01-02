#### Gaussian Estimates: Code to create a plotting DF for each gene, and also to categorize each gene based on ANOVA p-values. Need to do another similar script but with the negative bionomial estimates. Then compare and contrast the two
### Amanda Neves August 10
#############################################################################

library(tximport)
library(readr)
library("RColorBrewer")
library(stringr)
library(edgeR)
library(limma)
library(lme4)
library(glmmTMB)
library(emmeans)
library(parallel)
library(car)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(effects)

#Import
setwd("/Users/amandan/Desktop/Dworkin/background_effects/data/")
#setwd("/home/abneves/projects/def-idworkin/arteen/GBE/data")

#I renamed the folder containing all the counts (salmon_counts) to quants
quant_files <- file.path("quants", list.files("quants"), "quant.sf")
quant_files <- quant_files[ -c(14:20, 35:41) ]
file.exists(quant_files)

samples <- c("ORE_sd1_r1_GTGGCC_L002_quant",
             "ORE_sd1_r2_TTAGGC_L003_quant",
             "ORE_sdE3_r1_GTGGCC_L004_quant",
             "ORE_sdE3_r2_TGACCA_L005_quant",
             "ORE_sdETX4_r1_ATCACG_L003_quant",
             "ORE_sdETX4_r2_GTTTCG_L004_quant",
             "ORE_vg1_r1_TGACCA_L007_quant",
             "ORE_vg1_r2_TTAGGC_L001_quant",
             "ORE_vg213_r1_AGTCAA_L006_quant",
             "ORE_vg213_r2_CGATGT_L007_quant",
             "ORE_vg2a33_r1_CGATGT_L005_quant",
             "ORE_vg2a33_r2_AGTTCC_L006_quant",
             "ORE_w_r1_ATCACG_L001_quant",
             "ORE_w_r2_GTTTCG_L002_quant",
             "SAM_sd1_r1_CGTACG_L004_quant",
             "SAM_sd1_r2_GCCAAT_L005_quant",
             "SAM_sdE3_r1_ATGTCA_L006_quant",
             "SAM_sdE3_r2_GCCAAT_L007_quant",
             "SAM_sdETX4_r1_ACAGTG_L005_quant",
             "SAM_sdETX4_r2_CCGTCC_L006_quant",
             "SAM_vg1_r1_CGTACG_L002_quant",
             "SAM_vg1_r2_GATCAG_L003_quant",
             "SAM_vg213_r1_ACTTGA_L001_quant",
             "SAM_vg213_r2_GAGTGG_L002_quant",
             "SAM_vg2a33_r1_ACAGTG_L007_quant",
             "SAM_vg2a33_r2_GATCAG_L001_quant",
             "SAM_w_r1_ACTTGA_L003_quant",
             "SAM_w_r2_GAGTGG_L004_quant")


names(quant_files) <- samples

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(x = txdb, keys = k, "GENEID", "TXNAME")

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)



background <- c(rep("ORE", 14), rep("SAM", 14))

background <- as.factor(background)
background <- relevel(background, "ORE")

mutation_type <- c( rep("sd", 6), rep("vg", 6),
                    rep("w", 2), rep ("sd", 6),
                    rep("vg", 6), rep("w", 2))

mutation_type <- as.factor(mutation_type)
mutation_type <- relevel(mutation_type, "w")



lane <- str_sub(samples, - 10, - 1) 
lane <- str_sub(lane, 1, 4)
lane <- as.factor(lane)

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

rna.design <- data.frame(sample = samples,
                         file = quant_files,
                         background = background,
                         mutation_type = mutation_type,
                         perturbation = perturbation,
                         lane = lane)




design <- model.matrix(~ lane + background + perturbation + perturbation:background,
                       data = rna.design)

dgeList <- DGEList(txi$counts)

keep2 <- filterByExpr(dgeList, design)
dgeList <- dgeList[keep2, ]

# normalize and run voom transformation
dgeList <- calcNormFactors(dgeList)

v <- voom(dgeList, design,plot=FALSE)
logCPM_norm <- v$E

nt <- min(parallel::detectCores(),5)

#Fitting
fit <- lmFit(v, design)
fit <- eBayes(fit)



target.limma <- topTable(fit, coef = 1, sort.by = "p", n = Inf)




target_genes <- row.names(target.limma)

target_genes <- target_genes[1:5]


tmb.neb.genes.check.outlier <- list()

# setting up a list to put the PredictorEffect estimates for plotting
PE_plot_df <- list()

# setting up list to save the pvalues

p_values <- list()

# Run the model for all of the genes!

count <- 1
Sys.time()
for (g in target_genes) {
  count <- count + 1
  if (count%%100 == 0){
    print(count)
    print(Sys.time())
  }
  geneID <- g
  single_cpm <- logCPM_norm[g,]
  tmb.data <- data.frame(single_cpm, background, perturbation, lane)
  tmb.mod.outlier <- glmmTMB(single_cpm ~ background + 
                               perturbation + 
                               background:perturbation + 
                               (1 | lane),
                             data = tmb.data,
                             control = glmmTMBControl(parallel = nt))
  if (anyNA(as.vector(summary(tmb.mod.outlier)$coef$cond[,2:4])) == T) {
    tmb.mod.outlier.reduced <- glmmTMB(single_cpm ~ background + 
                                         perturbation +
                                         background:perturbation,
                                       data = tmb.data,
                                       control = glmmTMBControl(parallel = nt))
    tmb.neb.genes.check.outlier[[g]] <- list(summary(tmb.mod.outlier.reduced)$coef$cond, "Used Reduced Model")
    
    
  } else {
    
    tmb.neb.genes.check.outlier[[g]] <- list(summary(tmb.mod.outlier)$coef$cond, "Used Full Model")
    
    # Compute estimates with predictor effect
    
    plot_df <- as.data.frame(predictorEffect("perturbation", tmb.mod.outlier, xlevels = 13))
    
    # Add to list
    # grab only every other data point for ploting the line
    PE_plot_df[[g]] <- plot_df[seq(from = 1, to = 100, by = 2),]
    
    # Extract p values and save for multiple testing and categorizing genes later 
    
    p.vals <- Anova(tmb.mod.outlier)
    # background, perturbation, background:perturbation
    p_values[[g]] <- p.vals$`Pr(>Chisq)`

  }
}

#save(tmb.neb.genes.check.outlier, file = "/home/arteen/projects/def-idworkin/arteen/GBE/r/full_pert/tmb.neb.genes.check.outlier.Rdata")


# save the p.values into one dataframe

p_values <- do.call(rbind, p_values)
colnames(p_values) <- c("background", "perturbation", "background:perturbation")

#save(p_values, file = "/home/abneves/projects/def-idworkin/abneves/GBE/outputs/anova_pvalues.Rdata")
save(p_values, file = "/Users/amandan/Desktop/anova_pvalues.Rdata")


# save plotting df list into an R data file

#save(PE_plot_df, file = "/Users/amandan/Desktop/predictoreffect_plotting_df.Rdata")
#save(PE_plot_df, file = "/home/abneves/projects/def-idworkin/abneves/GBE/outputs/predictoreffect_plotting_df.Rdata")

# Collapse the list into a dataframe for potentially better saving
 
PE_plot_df <- do.call(rbind, PE_plot_df)
save(PE_plot_df, file = "/home/abneves/projects/def-idworkin/abneves/GBE/outputs/predictoreffect_plotting_df.Rdata")



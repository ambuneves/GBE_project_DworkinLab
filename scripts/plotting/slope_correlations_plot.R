### Script for vector correlations for the GBE stuff!
### First, using the models that are global (not specific to alleles)

######### Amanda Neves November 1, 2022 ########################

# Read in the pathways I chose and set up a dataframe with the gene names and pathway name, ensure genes only appear once in each pathway 
# For each background, grab the vector of genes (in the right order) from these pathways, this vector should contain the slopes from the models (I'll have to re-run these as I'm going so it might take a couple of minutes for each pathway)
# Correlate that vector between the two backgrounds and save this
# Plot!

################################

# libraries
library(ggplot2)
library(viridis)

#### Set up pathways ####
# Here I am just reading in all of the text files that are in a directory and saving them to a list.
# Each text file contains the genes in the GO term group I am looking at

pathway_info <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/pathway_correlations/vec_cor_pathways.csv")

#If the path is different than your working directory
# you'll need to set full.names = TRUE to get the full
# paths.
my_files <- list.files("/Users/amandan/Desktop/Dworkin/background_effects/data/pathway_correlations", pattern = ".txt", full.names = TRUE)

# read in each file and save into a list 
all_csv <- lapply(my_files,read.table)

for(i in 1:length(all_csv)){
  all_csv[[i]] <- unique(all_csv[[i]]$V1) # take the first column, and make sure there only unique values are taken
}

lengths(all_csv)

# rename the lists based on what the actual path names are
file_names <- list.files("/Users/amandan/Desktop/Dworkin/background_effects/data/pathway_correlations", pattern = ".txt")
file_names <- gsub(".txt", "", file_names) # remove the .txt suffix 

pathway_info$ID[match(file_names, pathway_info$ID)] # make sure it has been indexed right

# set the proper name
names(all_csv) <- pathway_info$pathway[match(file_names, pathway_info$ID)]

# some of the pathways can be merged (ex. positive/negative regulators), so do that here

hippo_all <- c(all_csv[["Hippo signalling pathway core components"]],
       all_csv[["Negative regulators of hippo signalling pathway"]],
       all_csv[["Positive regulators of hippo signalling pathway"]])

hippo_all <- unique(hippo_all)

cell_pop_pro <- c(all_csv[["Negative regulation of cell population proliferation"]],
                  all_csv[["Positive regulation of cell population proliferation"]])

cell_pop_pro <- unique(cell_pop_pro)

cell_death <- c(all_csv[["Positive regulation of cell death"]],
                all_csv[["Negative regulation of cell death"]])

cell_death <- unique(cell_death)

cell_diff <- c(all_csv[["Negative regulation of cell differentiation "]],
               all_csv[["Positive regulation of cell differentiation "]])

cell_diff <- unique(cell_diff)

path_list <- all_csv[c("Positive regulation of cell cycle ",
                             "Regulation of intracellular mRNA organization ",
                             "Regulation of cellular response to growth factor stimulus ",
                             "Apoptotic signalling pathway")]

length(path_list)
path_list[["Hippo pathway regulators"]] <- hippo_all
path_list[["Regulation of cell population proliferation"]] <- cell_pop_pro
path_list[["Regulation of cell death"]] <- cell_death
path_list[["Regulation of cell differentiation"]] <- cell_diff

# add perturbation genes as a positive control, we expect these to have a relatively high correlation ####
perturbation_genes <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/perturbation_only_genes.csv")
perturbation_genes <- perturbation_genes$x
path_list[["Perturbation significant genes"]] <- perturbation_genes

length(path_list)



#### grab slope vector for each background ####
# read in necessary data
slopes_df <- read.csv("/Users/amandan/Desktop/Dworkin/background_effects/data/all_slopes_nb.csv")
colnames(slopes_df) <- c("gene_id", "ore_slope", "sam_slope")
dim(slopes_df)

# Create a comparison distribution for each pathway, it will be the 2.5th and 97.5th percetiles of 1000 correlations between
# n random genes where n is the number of genes in each GO term group

datalist = list() # initiate a datalist 

iterations <- 1000 # set the number of iterations


for (n in 1: length(path_list)){
  pathway1 <- path_list[[n]]
  number_of_genes <- length(pathway1)
  #gene_vecs_subset <- slopes_df$gene_id[-which(slopes_df$gene_id %in% pathway1)]
  gene_vecs_subset <- slopes_df$gene_id
  iteration_vals <- matrix(nrow = iterations, ncol = 3)
  for (k in 1:iterations){
    iteration_vals[k,1] <- k
    iteration_vals[k,2] <- number_of_genes
    random_pathway <- sample(x = gene_vecs_subset, size = length(pathway1), replace = FALSE)
    ore_vec <- na.omit(slopes_df$ore_slope[match(random_pathway, slopes_df$gene_id)])
    sam_vec <- na.omit(slopes_df$sam_slope[match(random_pathway, slopes_df$gene_id)])
    iteration_vals[k,3] <- abs(cor(ore_vec, sam_vec))
  }
  
  datalist[[n]] <- iteration_vals
}

random_iterations = do.call(rbind, datalist) 
random_iterations <- as.data.frame(random_iterations)
colnames(random_iterations) <- c("iteration", "path_length", "correlation")

random_iterations$pathway <- rep(names(path_list), each = iterations) # add a column indicating which pathway the values come from 

# Extract 2.5% and 97.5% quantiles 
random_correlations <- do.call("rbind", tapply(random_iterations$correlation, random_iterations$pathway, 
                                           quantile, c(0.025, 0.975)))
random_correlations <- as.data.frame(random_correlations)
table(rownames(random_correlations))
random_correlations$pathway <- rownames(random_correlations)
rownames(random_correlations) <- NULL
colnames(random_correlations) <- c("bottom", "top", "pathwaynames")

### Plot with the new values!
# re-order to match the plot_df
random_correlations_1 <- random_correlations[match(plot_df$pathway, random_correlations$pathwaynames),]
plot_df_2 <- cbind.data.frame(plot_df, random_correlations_1)

# add in the perturbation genes row which was taken out by match

ggplot() + 
  geom_crossbar(data = plot_df_2, aes(x = xlabs, y = 1.5, ymin = bottom, ymax = top, fill = xlabs, colour = xlabs), alpha = 0.8) + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_classic() + 
  geom_crossbar(data = plot_df_2, aes(x = xlabs, y = abs_cor, ymin = abs_cor, ymax = abs_cor)) +
  ylab(expression(atop("Absolute value of correlation", paste("between SAM and ORE")))) + 
  xlab("Gene group") + 
  scale_fill_manual(values = c(viridis(8),"gray87")) +
  scale_colour_manual(values = c(viridis(8),"gray87")) +
  xlab("") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank())

ggsave("slope_vector_correlations_comparisons_no_subset.png", width = 11, height = 5, path = "/Users/amandan/Desktop/Dworkin/background_effects/outputs")

#Programming language - R

## Load packages and data files ------------------------------------------------

#Load packages
library(lme4)
library(dplyr)
library(writexl)
library(MASC)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
library(openxlsx)

#Load data
seurat_object <- readRDS("seurat_object.rds")

## Differential abundance of annotated clusters, using MASC---------------------

#Re-code phenotype so that controls are reference category
seurat_object$pheno2 = as.character(seurat_object$pheno)
seurat_object$pheno2[seurat_object$pheno2 == "control"] <- "0_control"
seurat_object$pheno2[seurat_object$pheno2 == "case"] <- "1_case"
seurat_object$pheno2 <- as.factor(seurat_object$pheno2)

#Create dataframe for masc command
masc_df <- data.frame(seurat_object$CELL.ID, seurat_object$pheno2, seurat_object$predicted.celltype.l2, seurat_object$predicted.celltype.l1, seurat_object$SingleR_labels, seurat_object$monaco_labels, stringsAsFactors = TRUE)

#Re-code specific characters to avoid errors later
masc_df$seurat_object.CELL.ID <- gsub("-", "_", masc_df$seurat_object.CELL.ID, fixed = TRUE)
masc_df$seurat_object.predicted.celltype.l2 <- gsub(" ", "_", masc_df$seurat_object.predicted.celltype.l2, fixed = TRUE)
masc_df$seurat_object.SingleR_labels <- gsub(" ", "_", masc_df$seurat_object.SingleR_labels, fixed = TRUE)
masc_df$seurat_object.SingleR_labels <- gsub("+", "", masc_df$seurat_object.SingleR_labels, fixed = TRUE)
masc_df$seurat_object.SingleR_labels <- gsub(",", "", masc_df$seurat_object.SingleR_labels, fixed = TRUE)
masc_df$seurat_object.monaco_labels <- gsub(" ", "_", masc_df$seurat_object.monaco_labels, fixed = TRUE)
masc_df$seurat_object.monaco_labels <- gsub("+", "", masc_df$seurat_object.monaco_labels, fixed = TRUE)
masc_df$seurat_object.monaco_labels <- gsub(",", "", masc_df$seurat_object.monaco_labels, fixed = TRUE)
masc_df$seurat_object.monaco_labels <- gsub("-", "_", masc_df$seurat_object.monaco_labels, fixed = TRUE)
masc_df$seurat_object.monaco_labels <- gsub("/", "_", masc_df$seurat_object.monaco_labels, fixed = TRUE)

#MASC command for main cluster annotations
masc_results_l2 <- MASC(data = masc_df, cluster = masc_df$seurat_object.predicted.celltype.l2, contrast = "seurat_object.pheno2", random_effects = "seurat_object.CELL.ID")

#MASC commands for Th cell subtypes (remove clusters with low cell count first)
masc_df_subset <- subset(x = masc_df, subset = seurat_object.monaco_labels != "Exhausted_B_cells")
masc_df_subset <- subset(x = masc_df_subset, subset = seurat_object.monaco_labels != "Intermediate_monocytes")
masc_df_subset <- subset(x = masc_df_subset, subset = seurat_object.monaco_labels != "Low_density_basophils")
masc_results_monaco <- MASC(data = masc_df_subset, cluster = masc_df_subset$seurat_object.monaco_labels, contrast = "seurat_object.pheno2", random_effects = "seurat_object.CELL.ID")
masc_results_monaco <- MASC(data = masc_df, cluster = masc_df$seurat_object.monaco_labels, contrast = "seurat_object.pheno2", random_effects = "seurat_object.CELL.ID")

#Document cell numbers by phenotype
cell_freq_l2 <- masc_df %>% group_by(seurat_object.predicted.celltype.l2, seurat_object.pheno2) %>% summarise(cell_count = n())
cell_freq_monaco <- masc_df %>% group_by(seurat_object.monaco_labels, seurat_object.pheno2) %>% summarise(cell_count = n())

#Save results outputs
write_xlsx(masc_results_l2, "masc_results_l2.xlsx")
write_xlsx(masc_results_monaco, "masc_results_monaco.xlsx")
write_xlsx(cell_freq_l2, "cell_freq_l2.xlsx")
write_xlsx(cell_freq_monaco, "cell_freq_monaco.xlsx")


## Differential abundance of cell neighbourhoods, using MiloR---------------------

#Load data
seurat_object <- readRDS("seurat_object.rds")
CD4 <- readRDS("seurat_object_CD4.rds")
CD8 <- readRDS("seurat_object_CD8.rds")
Mono <- readRDS("seurat_object_Mono.rds")
seurat_list <- list(CD4, CD8, Mono)

##Convert Seurat object to Single Cell Experiment object
sce <- as.SingleCellExperiment(seurat_object)
milo <- Milo(sce)

##Construct KNN graph
milo <- buildGraph(milo, k=30, d=30, reduced.dim = "PCAPREINTEGRATION")

##Define representative neighbourhoods
milo <- makeNhoods(milo, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "PCAPREINTEGRATION")

##Check how many cells form each neighbourhood (distribution peaking between 50 to 100 is good)
nhood_size_plot <- plotNhoodSizeHist(milo)

##Count number of cells from each sample in each neighbourhood
milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="CELL.ID")
head(nhoodCounts(milo))

##Differential abundance testing
milo_design <- data.frame(colData(milo))[,c("CELL.ID", "pheno")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$CELL.ID
milo_design
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCAPREINTEGRATION")
da_results <- testNhoods(milo, design = ~ pheno, design.df = milo_design, reduced.dim = "PCAPREINTEGRATION")
da_results %>%
  arrange(SpatialFDR) %>%
  head()

##Visualize neighbourhoods displaying differential abundance
milo <- buildNhoodGraph(milo)
plotUMAP(milo) + plotNhoodGraphDA(milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

##Volcano plot of differential abundance by phenotype
volcano_plot <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)

##Identify whether differential abundance particularly evident in certain cell types
da_results <- annotateNhoods(milo, da_results, coldata_col = "predicted.celltype.l2")
head(da_results)
ggplot(da_results, aes(predicted.celltype.l2_fraction)) + geom_histogram(bins=50)
da_results$predicted.celltype.l2 <- ifelse(da_results$predicted.celltype.l2_fraction < 0.7, "Mixed", da_results$predicted.celltype.l2)
plotDAbeeswarm(da_results, group.by = "predicted.celltype.l2")

##Identify if differential abundance evident in Th subsets
da_results_Th <- annotateNhoods(milo, da_results, coldata_col = "monaco_labels")
ggplot(da_results_Th, aes(monaco_labels_fraction)) + geom_histogram(bins=50)
da_results_Th$monaco_labels <- ifelse(da_results_Th$monaco_labels_fraction < 0.75, "Mixed", da_results_Th$monaco_labels)
plotDAbeeswarm(da_results_Th, group.by = "monaco_labels")


##Repeat entire code above but for monocyte and CD8/CD4 T cell subclusters
seurat_list <- readRDS("/mnt/jw01-aruk-home01/projects/single_cell/projects/APOLLO_Ali/analyses/primary/data/output/R/seurat_list.rds")
subcluster_cell_types <- c("CD4_T", "CD8_T", "Mono")
sce_list <- list()
milo_list <- list()
milo_design_list <- list()
da_results_list <- list()
nhood_hist_list <- list()
volcano_plot_list <- list()
celltype_histo_list <- list()
beeswarm_plot_list <- list()

for(i in subcluster_cell_types){
  DefaultAssay(object = seurat_list[[i]]) <- "SCT_pre_integration"
  sce_list[[i]] <- as.SingleCellExperiment(seurat_list[[i]])
  milo_list[[i]] <- Milo(sce_list[[i]])
  milo_list[[i]] <- buildGraph(milo_list[[i]], k=30, d=30, reduced.dim = "PCAPREINTEGRATION")
  milo_list[[i]] <- makeNhoods(milo_list[[i]], prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "PCAPREINTEGRATION")
  nhood_hist_list[[i]] <- plotNhoodSizeHist(milo_list[[i]])
  milo_list[[i]] <- countCells(milo_list[[i]], meta.data = data.frame(colData(milo_list[[i]])), sample="CELL.ID")
  #head(nhoodCounts(milo_list[[i]]))
  milo_design_list[[i]] <- data.frame(colData(milo_list[[i]]))[,c("CELL.ID", "pheno")]
  milo_design_list[[i]] <- distinct(milo_design_list[[i]])
  rownames(milo_design_list[[i]]) <- milo_design_list[[i]]$CELL.ID
  milo_list[[i]] <- calcNhoodDistance(milo_list[[i]], d=30, reduced.dim = "PCAPREINTEGRATION")
  da_results_list[[i]] <- testNhoods(milo_list[[i]], design = ~ pheno, design.df = milo_design_list[[i]], reduced.dim = "PCAPREINTEGRATION")
  da_results_list[[i]] %>%
    arrange(SpatialFDR) %>%
    head()
  milo_list[[i]] <- buildNhoodGraph(milo_list[[i]])
  volcano_plot_list[[i]] <- ggplot(da_results_list[[i]], aes(logFC, -log10(SpatialFDR))) +
    geom_point() +
    geom_hline(yintercept = 1)
  da_results_list[[i]] <- annotateNhoods(milo_list[[i]], da_results_list[[i]], coldata_col = "predicted.celltype.l2")
  celltype_histo_list[[i]] <- ggplot(da_results_list[[i]], aes(predicted.celltype.l2_fraction)) + geom_histogram(bins=50)
  da_results_list$i$predicted.celltype.l2 <- ifelse(da_results_list$i$predicted.celltype.l2_fraction < 0.7, "Mixed", da_results_list$i$predicted.celltype.l2)
  beeswarm_plot_list[[i]] <- plotDAbeeswarm(da_results_list[[i]], group.by = "predicted.celltype.l2")
  
##Save outputs
write.xlsx(da_results, "da_results.xlsx")
write.xlsx(da_results_list, "da_results_list.xlsx")
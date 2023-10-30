#Programming language - R
#Note: as described in the manuscript, initial data cleaning, debarcoding and batch normalisation was undertaken using FlowJo v10.9.

## Load packages and data files ------------------------------------------------
library(diffcyt)
library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(openCyto)
library(dplyr)
library(cowplot)
library(openxlsx)
library(ggplot2)

## Prepare data files ----------------------------------------------------------

#Read in csv files containing meta-data and panel data (created in Excel)
my_panel <- read.csv(file="panel.csv")
meta_data <- read.csv(file="meta_data.csv")

#Use flowCore to create a flowSet object from csv files exported from FlowJo (error if using .fcs files due to upper limit of some parameters being set to infinity)
data_path <- "/user/fcs_csv_directory"
data_files <- list.files("user/fcs_csv_directory")
data <- openCyto:::.read.flowSet.csv(paste(data_path, data_files, sep="/"))

#Use CATALYST to create single cell experiment object from flowset object
sce <- prepData(data, my_panel, meta_data)

#Split into stimulated and unstimulated
sce_unstim <- filterSCE(sce, condition %in% c("Case", "Control", "Healthy"))
sce_stim <- filterSCE(sce, condition %in% c("Case_stim", "Control_stim", "Healthy_stim"))

## Generate clusters -----------------------------------------------------------

sce_unstim <- filterSCE(sce_unstim, colSums(is.infinite(assay(sce_unstim, "exprs"))) == 0) #Preceding cluster with this command avoids an error
sce_stim <- filterSCE(sce_stim, colSums(is.infinite(assay(sce_stim, "exprs"))) == 0)
sce_unstim <- cluster(sce_unstim, features = "type", 
                      xdim = 12, ydim = 12, maxK = 60, 
                      verbose = FALSE, seed = 12643)

sce_stim <- cluster(sce_stim, features = "type", 
                    xdim = 12, ydim = 12, maxK = 60, 
                    verbose = FALSE, seed = 12643)

saveRDS(sce_unstim, file="sce_unstim.rds")
saveRDS(sce_stim, file="sce_stim.rds")

## Visualisation of data-------------------------------------------------------

#Plot of number of cells per barcoded sample

plotCounts(sce_unstim, 
           group_by = "sample_id", 
           color_by = "condition")

plotCounts(sce_stim, 
           group_by = "sample_id", 
           color_by = "condition")


#Heatmap of marker expression by cluster

cluster_heatmap_unstim <- plotExprHeatmap(sce_unstim, features = "type",
                                          by = "cluster_id", k = "meta60",
                                          scale = "first", q = 0.01, bars = FALSE, row_anno = FALSE, k_pal = CATALYST:::.cluster_cols)

cluster_heatmap_stim <- plotExprHeatmap(sce_stim, features = "type",
                                        by = "cluster_id", k = "meta60",
                                        scale = "first", q = 0.01, bars = FALSE, row_anno = FALSE, k_pal = CATALYST:::.cluster_cols)

setwd("/user/images")
png("60_cluster_heatmap_unstim.png", width=2000, height=3200, res=300)
cluster_heatmap_unstim
dev.off()

png("60_cluster_heatmap_stim.png", width=2000, height=3200, res=300)
cluster_heatmap_stim
dev.off()

#tSNE by cluster ID, including by each cluster and marker
set.seed(1601)
sce_unstim <- runDR(sce_unstim, dr = "TSNE", cells = 10000, features = "type")
saveRDS(sce_unstim, file="sce_unstim.rds")
sce_stim <- runDR(sce_stim, dr = "TSNE", cells = 10000, features = "type")
saveRDS(sce_stim, file="sce_stim.rds")

markers <- c("CD38", "CD56", "TCRgd", "CD294", "CD197", "CD14", "CD3", "CD161", "CD25", "CD57", "CD183", "CD185", "CD28", "CD123", "CD19", "CD4", "CD8a", "CD16", "CD45RA", "CD196", "CD11c", "CD45RO", "CD194", "CD27", "CD20", "CD66b", "HLA-DR", "IgD", "CD127")
cytokines <- c("TNFa", "IFNg", "IL-17A", "IL-10", "IFNa", "IL-4", "GM-CSF", "IL-13")
tSNE_unstim_markers_list <- list()
tSNE_unstim_cytokine_list <- list()
tSNE_unstim_pheno_cytokine_list <- list()
tSNE_stim_markers_list <- list()
tSNE_stim_cytokine_list <- list()
tSNE_stim_pheno_cytokine_list <- list()

tSNE_unstim_basic <- plotDR(sce_unstim, dr="TSNE") 						#Plot basic tSNE
tSNE_stim_basic <- plotDR(sce_stim, dr="TSNE")
tSNE_unstim_pheno <- plotDR(sce_unstim, dr="TSNE", facet_by = "condition") 			#Split by phenotype
tSNE_stim_pheno <- plotDR(sce_stim, dr="TSNE", facet_by = "condition")
tSNE_unstim_clusters_colours <- plotDR(sce_unstim, dr="TSNE", color_by = "meta60")		#Coloured by cluster ID
tSNE_stim_clusters_colours <- plotDR(sce_stim, dr="TSNE", color_by = "meta60")

setwd("/user/images/unstim_plots/cluster_tSNEs")
png("1.tSNE_unstim_clusters_colours.png", width=2000, height=2000, res=300)
tSNE_unstim_clusters_colours
dev.off()
setwd("/user/images/stim_plots/cluster_tSNEs")
png("1.tSNE_stim_clusters_colours.png", width=2000, height=2000, res=300)
tSNE_stim_clusters_colours
dev.off()

for(i in markers){	#Coloured by each cell surface marker
  tSNE_unstim_markers_list[[i]] <- plotDR(sce_unstim, dr="TSNE", color_by = i)
  tSNE_stim_markers_list[[i]] <- plotDR(sce_stim, dr="TSNE", color_by = i)
}
for(i in cytokines){	#Coloured by each cytokine marker, including split by phenotype
  tSNE_unstim_cytokine_list[[i]] <- plotDR(sce_unstim, dr="TSNE", color_by = i)
  tSNE_stim_cytokine_list[[i]] <- plotDR(sce_stim, dr="TSNE", color_by = i)
  tSNE_unstim_pheno_cytokine_list[[i]] <- plotDR(sce_unstim, dr="TSNE", color_by = i, facet_by = "condition")
  tSNE_stim_pheno_cytokine_list[[i]] <- plotDR(sce_stim, dr="TSNE", color_by = i, facet_by = "condition")
}

#Save tSNE plots
for(i in markers){
  setwd("E:/Ali/Documents/Work/SpR/PhD/Research Data/5.CyTOF/analysis/images/unstim_plots/surface_tSNEs")
  png(paste(i, "_unstim_tSNE", ".png", sep=""), width=2000, height=2000, res=300)
  print(tSNE_unstim_markers_list[[i]])
  dev.off()
  setwd("E:/Ali/Documents/Work/SpR/PhD/Research Data/5.CyTOF/analysis/images/stim_plots/surface_tSNEs")
  png(paste(i, "_stim_tSNE", ".png", sep=""), width=2000, height=2000, res=300)
  print(tSNE_stim_markers_list[[i]])
  dev.off()
}

for(i in cytokines){
  setwd("E:/Ali/Documents/Work/SpR/PhD/Research Data/5.CyTOF/analysis/images/unstim_plots/cytokine_tSNEs")
  png(paste(i, "_unstim_pheno_tSNE", ".png", sep=""), width=2000, height=2000, res=300)
  print(tSNE_unstim_pheno_cytokine_list[[i]])
  dev.off()
  setwd("E:/Ali/Documents/Work/SpR/PhD/Research Data/5.CyTOF/analysis/images/stim_plots/cytokine_tSNEs")
  png(paste(i, "_stim_pheno_tSNE", ".png", sep=""), width=2000, height=2000, res=300)
  print(tSNE_stim_pheno_cytokine_list[[i]])
  dev.off()
  setwd("E:/Ali/Documents/Work/SpR/PhD/Research Data/5.CyTOF/analysis")
}

## Clusters were manually annotated based on marker expression, as described in the manuscript.

## Manual merging of clusters---------------------------------------------------

#Read in Excel file with new and old cluster identities
merging_unstim1 <- read.xlsx("cluster_merging_unstim1.xlsx") #Main clusters
merging_unstim2 <- read.xlsx("cluster_merging_unstim2.xlsx") #Th subsets
merging_stim1 <- read.xlsx("cluster_merging_stim1.xlsx") #Main clusters
merging_stim2 <- read.xlsx("cluster_merging_stim2.xlsx") #Th subsets

#Merging command
sce_unstim <- mergeClusters(sce_unstim, k = "meta60", table = merging_unstim1, id = "merging1")
sce_unstim <- mergeClusters(sce_unstim, k = "meta60", table = merging_unstim2, id = "merging2")
sce_stim <- mergeClusters(sce_stim, k = "meta60", table = merging_stim1, id = "merging1")
sce_stim <- mergeClusters(sce_stim, k = "meta60", table = merging_stim2, id = "merging2")

#Remove excluded clusters
sce_stim <- filterSCE(sce_stim, k="meta60", cluster_id != 50)
sce_unstim_Th <- filterSCE(sce_unstim, k="merging2", cluster_id %in% c("Th1", "Th2", "Th17", "Tc1", "Tc2", "Tc17"))
sce_stim_Th <- filterSCE(sce_stim, k="merging2", cluster_id %in% c("Th17", "Th1_Th17"))

## Differential abundance and differential states using diffcyt-----------------

#Create design and contrast matrix
#Note: the design matrix includes group 1 (case non-stim) as an intercept. Therefore contrasts including this group code this group as 0.
design_unstim <- createDesignMatrix(ei(sce_unstim), cols_design = "condition")
design_stim <- createDesignMatrix(ei(sce_stim), cols_design = "condition")
design_unstim_Th <- createDesignMatrix(ei(sce_unstim_Th), cols_design = "condition")
design_stim_Th <- createDesignMatrix(ei(sce_stim_Th), cols_design = "condition")

contrast_unstim_1 <- createContrast(c(0, -1, 0))	#Case unstim vs control unstim - control is reference
contrast_unstim_2 <- createContrast(c(0, 0, -1))	#Case unstim vs healthy unstim - healthy is reference
contrast_unstim_3 <- createContrast(c(0, 1, -1))	#Control unstim vs healthy unstim - healthy is reference
contrast_stim_1 <- createContrast(c(0, -1, 0))		#Case stim vs control stim - control is reference
contrast_stim_2 <- createContrast(c(0, 0, -1))		#Case stim vs healthy stim - healthy is reference
contrast_stim_3 <- createContrast(c(0, 1, -1))		#Control stim vs healthy stim - healthy is reference

#Differential abundance - note: function automatically adjusts for number of cells inputted per donor. The "normalize" function additionally adjusts for compositional differences which can create false positives in non-differentially abundant populations. See paper: https://www.nature.com/articles/s42003-019-0415-5
#Note: transform set to FALSE - already arcsinh transformed during CATALYST's "prepData" step.

unstim_res_DA_1 <- diffcyt(sce_unstim,
                           analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                           design = design_unstim, contrast = contrast_unstim_1, verbose = FALSE, transform = FALSE)

unstim_res_DA_2 <- diffcyt(sce_unstim,
                           analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                           design = design_unstim, contrast = contrast_unstim_2, verbose = FALSE, transform = FALSE)

unstim_res_DA_3 <- diffcyt(sce_unstim,
                           analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                           design = design_unstim, contrast = contrast_unstim_3, verbose = FALSE, transform = FALSE)

stim_res_DA_1 <- diffcyt(sce_stim,
                         analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                         design = design_stim, contrast = contrast_stim_1, verbose = FALSE, transform = FALSE)

stim_res_DA_2 <- diffcyt(sce_stim,
                         analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                         design = design_stim, contrast = contrast_stim_2, verbose = FALSE, transform = FALSE)

stim_res_DA_3 <- diffcyt(sce_stim,
                         analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
                         design = design_stim, contrast = contrast_stim_3, verbose = FALSE, transform = FALSE)


#Repeat but for merging2 (Tc and Th subtypes)

unstim_Th_res_DA_1 <- diffcyt(sce_unstim_Th,
                              analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                              design = design_unstim_Th, contrast = contrast_unstim_1, verbose = FALSE, transform = FALSE)

unstim_Th_res_DA_2 <- diffcyt(sce_unstim_Th,
                              analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                              design = design_unstim, contrast = contrast_unstim_2, verbose = FALSE, transform = FALSE)

unstim_Th_res_DA_3 <- diffcyt(sce_unstim_Th,
                              analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                              design = design_unstim, contrast = contrast_unstim_3, verbose = FALSE, transform = FALSE)

stim_Th_res_DA_1 <- diffcyt(sce_stim_Th,
                            analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                            design = design_stim, contrast = contrast_stim_1, verbose = FALSE, transform = FALSE)

stim_Th_res_DA_2 <- diffcyt(sce_stim_Th,
                            analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                            design = design_stim, contrast = contrast_stim_2, verbose = FALSE, transform = FALSE)

stim_Th_res_DA_3 <- diffcyt(sce_stim_Th,
                            analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging2",
                            design = design_stim, contrast = contrast_stim_3, verbose = FALSE, transform = FALSE)


#Differential states
unstim_res_DS_1 <- diffcyt(sce_unstim, clustering_to_use = "merging1",
                           analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                           design = design_unstim, contrast = contrast_unstim_1, verbose = FALSE, transform = FALSE)

unstim_res_DS_2 <- diffcyt(sce_unstim, clustering_to_use = "merging1",
                           analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                           design = design_unstim, contrast = contrast_unstim_2, verbose = FALSE, transform = FALSE)

unstim_res_DS_3 <- diffcyt(sce_unstim, clustering_to_use = "merging1",
                           analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                           design = design_unstim, contrast = contrast_unstim_3, verbose = FALSE, transform = FALSE)

stim_res_DS_1 <- diffcyt(sce_stim, clustering_to_use = "merging1",
                         analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                         design = design_stim, contrast = contrast_stim_1, verbose = FALSE, transform = FALSE)

stim_res_DS_2 <- diffcyt(sce_stim, clustering_to_use = "merging1",
                         analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                         design = design_stim, contrast = contrast_stim_2, verbose = FALSE, transform = FALSE)

stim_res_DS_3 <- diffcyt(sce_stim, clustering_to_use = "merging1",
                         analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                         design = design_stim, contrast = contrast_stim_3, verbose = FALSE, transform = FALSE)


#Repeat but for merging2 (Tc and Th subtypes)
#Differential states
unstim_Th_res_DS_1 <- diffcyt(sce_unstim_Th, clustering_to_use = "merging2",
                              analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                              design = design_unstim, contrast = contrast_unstim_1, verbose = FALSE, transform = FALSE)
unstim_Th_tbl_DS_1 <- rowData(unstim_Th_res_DS_1$res)

unstim_Th_res_DS_2 <- diffcyt(sce_unstim_Th, clustering_to_use = "merging2",
                              analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                              design = design_unstim, contrast = contrast_unstim_2, verbose = FALSE, transform = FALSE)
unstim_Th_tbl_DS_2 <- rowData(unstim_Th_res_DS_2$res)

unstim_Th_res_DS_3 <- diffcyt(sce_unstim_Th, clustering_to_use = "merging2",
                              analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                              design = design_unstim, contrast = contrast_unstim_3, verbose = FALSE, transform = FALSE)
unstim_Th_tbl_DS_3 <- rowData(unstim_Th_res_DS_3$res)

stim_Th_res_DS_1 <- diffcyt(sce_stim_Th, clustering_to_use = "merging2",
                            analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                            design = design_stim, contrast = contrast_stim_1, verbose = FALSE, transform = FALSE)
stim_Th_tbl_DS_1 <- rowData(stim_Th_res_DS_1$res)

stim_Th_res_DS_2 <- diffcyt(sce_stim_Th, clustering_to_use = "merging2",
                            analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                            design = design_stim, contrast = contrast_stim_2, verbose = FALSE, transform = FALSE)
stim_Th_tbl_DS_2 <- rowData(stim_Th_res_DS_2$res)

stim_Th_res_DS_3 <- diffcyt(sce_stim_Th, clustering_to_use = "merging2",
                            analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                            design = design_stim, contrast = contrast_stim_3, verbose = FALSE, transform = FALSE)
stim_Th_tbl_DS_3 <- rowData(stim_Th_res_DS_3$res)

#Put all results into lists, save and export to Excel
DA_results_list <- list(unstim_res_DA_1, unstim_res_DA_2, unstim_res_DA_3, stim_res_DA_1, stim_res_DA_2, stim_res_DA_3, unstim_Th_res_DA_1, unstim_Th_res_DA_2, unstim_Th_res_DA_3, stim_Th_res_DA_1, stim_Th_res_DA_2, stim_Th_res_DA_3)
DS_results_list <- list(unstim_res_DS_1, unstim_res_DS_2, unstim_res_DS_3, stim_res_DS_1, stim_res_DS_2, stim_res_DS_3, unstim_Th_res_DS_1, unstim_Th_res_DS_2, unstim_Th_res_DS_3, stim_Th_res_DS_1, stim_Th_res_DS_2, stim_Th_res_DS_3)
results_names <- c("unstim_1", "unstim_2", "unstim_3", "stim_1", "stim_2", "stim_3", "unstim_Th_1", "unstim_Th_2", "unstim_Th_3", "stim_Th_1", "stim_Th_2", "stim_Th_3")
names(DA_results_list) <- results_names
names(DS_results_list) <- results_names

saveRDS(DA_results_list, file="DA_results_list.rds")
saveRDS(DS_results_list, file="DS_results_list.rds")

DA_tables_list <- list()
DS_tables_list <- list()
for(i in results_names){
  DA_tables_list[[i]] <- rowData(DA_results_list[[i]]$res)
  DS_tables_list[[i]] <- rowData(DS_results_list[[i]]$res)
}

sig_DA_results <- list()
sig_DS_results <- list()
all_DA_results <- list()
all_DS_results <- list()

for(i in results_names){
  sig_DA_results[[i]] <- data.frame(DA_tables_list[[i]]@listData) %>% filter(p_adj <= 0.1) %>% arrange(p_adj)
  sig_DS_results[[i]] <- data.frame(DS_tables_list[[i]]@listData) %>% filter(p_adj <= 0.1) %>% arrange(p_adj)
  all_DA_results[[i]] <- data.frame(DA_tables_list[[i]]@listData) %>% arrange(p_adj)
  all_DS_results[[i]] <- data.frame(DS_tables_list[[i]]@listData) %>% arrange(p_adj)
}


#Output all and significant results into Excel docs (multiple tabs in each)
setwd("/user/differential_results")
write.xlsx(sig_DA_results, "sig_DA_results.xlsx")
write.xlsx(sig_DS_results, "sig_DS_results.xlsx")
write.xlsx(all_DA_results, "all_DA_results.xlsx")
write.xlsx(all_DS_results, "all_DS_results.xlsx")

#Prep data for displaying results in heatmaps
sce_unstim$sample_id <- factor(sce_unstim$sample_id, levels = c("Patient1_non-stim", "Patient2_non-stim", "Patient3_non-stim", "Patient4_non-stim", "Patient5_non-stim", "Patient6_non-stim", "Patient7_non-stim", "Patient8_non-stim", "Patient9_non-stim", "Patient10_non-stim", "Patient11_non-stim", "Patient12_non-stim", "Patient13_non-stim", "Patient14_non-stim")) #Reorder factor levels, otherwise puts patients 10, 11 etc straight after patient 1

sce_stim$sample_id <- factor(sce_stim$sample_id, levels = c("Patient1_stim", "Patient2_stim", "Patient3_stim", "Patient4_stim", "Patient5_stim", "Patient6_stim", "Patient7_stim", "Patient8_stim", "Patient9_stim", "Patient10_stim", "Patient11_stim", "Patient12_stim", "Patient13_stim", "Patient14_stim"))  

sce_cases_controls <- filterSCE(sce_unstim, condition == c("Case", "Control"))
sce_cases_healthy <- filterSCE(sce_unstim, condition == c("Case", "Healthy"))
sce_controls_healthy <- filterSCE(sce_unstim, condition == c("Control", "Healthy"))
sce_stim_cases_controls <- filterSCE(sce_stim, condition == c("Case_stim", "Control_stim"))
sce_stim_cases_healthy <- filterSCE(sce_stim, condition == c("Case_stim", "Healthy_stim"))
sce_stim_controls_healthy <- filterSCE(sce_stim, condition == c("Control_stim", "Healthy_stim"))


##Modify catalyst package code to correctly display heatmap in descending logFC order.
#trace(plotDiffHeatmap, edit=T)
##Modified following lines:
#  if (sort_by != "none") {
#    o <- order(abs(y[[sort_by]]), decreasing = (sort_by == 
#      "lfc"))
##To following:
#  if (sort_by != "none") {
#    o <- order(y[[sort_by]], decreasing = TRUE)

#Heatmap commands (did not plot first two stim contrasts because no significant results)
case_control_heatmap <- plotDiffHeatmap(sce_cases_controls, DS_tables_list[["unstim_1"]], all = FALSE, fdr = 0.1, lfc = 0, col_anno = "condition", sort_by = "lfc", top_n=100)
case_healthy_heatmap <- plotDiffHeatmap(sce_cases_healthy, DS_tables_list[["unstim_2"]], all = FALSE, fdr = 0.1, lfc = 0, col_anno = "condition", sort_by = "lfc", top_n=100)
control_healthy_heatmap <- plotDiffHeatmap(sce_controls_healthy, DS_tables_list[["unstim_3"]], all = FALSE, fdr = 0.1, lfc = 0, col_anno = "condition", sort_by = "lfc", top_n=100)
stim_control_healthy_heatmap <- plotDiffHeatmap(sce_stim_controls_healthy, DS_tables_list[["stim_3"]], all = FALSE, fdr = 0.1, lfc = 0, col_anno = "condition", sort_by = "lfc", top_n=100)
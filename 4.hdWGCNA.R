#Programming language - R

## Load packages and data files ------------------------------------------------

library(Seurat)
library(WGCNA)
library(igraph)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qlcMatrix)
library(tidyverse)
library(ggrepel)
library(enrichR)
library(GeneOverlap)
library(sctransform)

seurat_object <- readRDS("seurat_object.rds")
CD4 <- readRDS("seurat_object_CD4.rds")
CD8 <- readRDS("seurat_object_CD8.rds")
Mono <- readRDS("seurat_object_Mono.rds")
seurat_list <- list(CD4, CD8, Mono)

# set random seed for reproducibility
set.seed(12345)



## Prep data sets --------------------------------------------------------------


## Remove "/" from the "Th1/Th17 cells" monaco cell type
new <- seurat_object$monaco_labels
new <- gsub("/", "_", new, fixed = TRUE)
seurat_object <- AddMetaData(seurat_object, new, col.name = "monaco_labels")

## Filter cell types (done separately for each method of cell labelling)

cells_multi <- c("B intermediate", "B memory", "B naive", "CD14 Mono", "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "cDC2", "gdT", "NK", "Treg")
cells_monaco <- c("Th2 cells", "Th17 cells", "Th1_Th17 cells")

Idents(seurat_object) <- "predicted.celltype.l2"
hdWGCNA_multi <- subset(x = seurat_object, idents = cells_multi)
Idents(seurat_object) <- "monaco_labels"
hdWGCNA_monaco <- subset(x = seurat_object, idents = cells_monaco)



## Analytical steps (multimodal reference annotation) ----------------------------------------------------------------------


## Set up Seurat object for WGCNA (select all genes, previously outputted from SCT_pre_integration assay)


power_plots_multi <- list()
dendrograms_multi <- list()
TOM_list_multi <- list()
kme_plots_multi <- list()
modules_multi <- list()
powers_multi <- list("B intermediate" = 5, "B memory" = 10, "B naive" = 10, "CD14 Mono" = 5, "CD16 Mono" = 6, "CD4 CTL" = 10, "CD4 Naive" = 8, "CD4 TCM" = 6, "CD4 TEM" = 16, "CD8 Naive" = 14, "CD8 TCM" = 12, "CD8 TEM" = 12, "cDC2" = 8, "gdT" = 12, "NK" = 10, "Treg" = 6)

for(i in cells_multi){
  hdWGCNA_multi <- SetupForWGCNA(
    hdWGCNA_multi,
    features = GetWGCNAGenes(hdWGCNA_multi, "SCT_pre_integration"),
    wgcna_name = i)
  
  ## Construct metacells
  
  hdWGCNA_multi <- MetacellsByGroups(
    seurat_obj = hdWGCNA_multi,
    group.by = c("predicted.celltype.l2", "CELL.ID"),
    k = 20,
    max_shared=10,
    min_cells = 50,
    ident.group = 'predicted.celltype.l2',
    slot = 'scale.data',
    assay = 'SCT_pre_integration',
    reduction = "pcapreintegration",
    wgcna_name = i
  )
  
  ## Set up expression matrix for network analysis on selected cell types (note, this needs to be repeated on individual cell types unless reason to combine them).
  
  hdWGCNA_multi  <- SetDatExpr(
    hdWGCNA_multi,
    group_name = i,
    group.by='predicted.celltype.l2',
    wgcna_name = i
  )
  
  hdWGCNA_multi  <- TestSoftPowers(
    hdWGCNA_multi,
    group.by = "predicted.celltype.l2",
    group_name = i
  )
  
  ## Select soft power threshold (identified as soft power threshold using below commands, then saved in powers_multi list as above)
  
  power_plots_multi[[i]] <- PlotSoftPowers(hdWGCNA_multi, wgcna_name = i)
  
  ##View plot
  #wrap_plots(power_plots_multi[[i]], ncol=2)
  
  ## Construct co-expression network (need separate soft power for each cell type, identified after inspecting above plots)
  
  hdWGCNA_multi <- ConstructNetwork(
    hdWGCNA_multi,
    setDatExpr=FALSE,
    wgcna_name = i,
    tom_name = i,
    soft_power = powers_multi[i],
    overwrite_tom = TRUE
  )
  
  dendrograms_multi[[i]] <- PlotDendrogram(hdWGCNA_multi, main=i, wgcna_name = i)
  
  ## Inspect topological overlap matrix (TOM) which can be used for custom downstream analyses later
  
  TOM_list_multi[[i]] <- GetTOM(hdWGCNA_multi, wgcna_name = i)
  
  
  ## Compute module eigengenes (first PC of gene expression in each module) and connectivity (how connected each gene is to the hub/main gene in each module)
  hdWGCNA_multi <- ModuleEigengenes(
    hdWGCNA_multi,
    wgcna_name = i
  )
  
  ## Compute model connectivity
  hdWGCNA_multi <- ModuleConnectivity(
    hdWGCNA_multi,
    group.by = 'predicted.celltype.l2', group_name = i,
    wgcna_name = i
  )
  
  ## kme plots
  kme_plots_multi[[i]] <- PlotKMEs(hdWGCNA_multi, ncol=5, wgcna_name = i)
  
  ## Get module assignment table
  modules_multi[[i]] <- GetModules(hdWGCNA_multi, wgcna_name = i)
  
  ## Gene scoring of top 25 genes in each set
  hdWGCNA_multi <- ModuleExprScore(
    hdWGCNA_multi,
    n_genes = 25,
    method='Seurat',
    wgcna_name = i
  )
}


## Analytical steps (Th subsets) ----------------------------------------------------------------------

## Set up Seurat object for WGCNA (select all genes, previously outputted from SCT_pre_integration assay)

power_plots_monaco <- list()
dendrograms_monaco <- list()
TOM_list_monaco <- list()
kme_plots_monaco <- list()
modules_monaco <- list()
powers_monaco <- list("Th2 cells" = 9, "Th17 cells" = 6, "Th1_Th17 cells" = 16)

for(i in cells_monaco){
  hdWGCNA_monaco <- SetupForWGCNA(
    hdWGCNA_monaco,
    features = GetWGCNAGenes(hdWGCNA_monaco, "SCT_pre_integration"),
    wgcna_name = i)
  
  ## Construct metacells
  
  hdWGCNA_monaco <- MetacellsByGroups(
    seurat_obj = hdWGCNA_monaco,
    group.by = c("monaco_labels", "CELL.ID"),
    k = 25,
    max_shared=10,
    min_cells = 50,
    ident.group = 'monaco_labels',
    slot = 'scale.data',
    assay = 'SCT_pre_integration',
    reduction = "pcapreintegration",
    wgcna_name = i
  )
  
  ## Set up expression matrix for network analysis on selected cell types (note, this needs to be repeated on individual cell types unless reason to combine them).
  
  hdWGCNA_monaco  <- SetDatExpr(
    hdWGCNA_monaco,
    group_name = i,
    group.by='monaco_labels',
    wgcna_name = i
  )
  
  hdWGCNA_monaco  <- TestSoftPowers(
    hdWGCNA_monaco,
    group.by = "monaco_labels",
    group_name = i
  )
  
  ## Select soft power threshold (identified different soft powers for each cell type [added to powers_monaco list above])
  
  power_plots_monaco[[i]] <- PlotSoftPowers(hdWGCNA_monaco, wgcna_name = i)
  
  ##View plot
  #wrap_plots(power_plots_monaco[[i]], ncol=2)
  
  
  ## Construct co-expression network (need separate soft power for each cell type, identified after inspecting above plots)
  
  hdWGCNA_monaco <- ConstructNetwork(
    hdWGCNA_monaco,
    setDatExpr=FALSE,
    wgcna_name = i,
    tom_name = i,
    soft_power = powers_monaco[i],
    overwrite_tom = TRUE
  )
  
  dendrograms_monaco[[i]] <- PlotDendrogram(hdWGCNA_monaco, main=i, wgcna_name = i)
  
  ## Inspect topological overlap matrix (TOM) which can be used for custom downstream analyses later
  
  TOM_list_monaco[[i]] <- GetTOM(hdWGCNA_monaco, wgcna_name = i)
  
  
  ## Compute module eigengenes (first PC of gene expression in each module) and connectivity (how connected each gene is to the hub/main gene in each module)
  hdWGCNA_monaco <- ModuleEigengenes(
    hdWGCNA_monaco,
    wgcna_name = i
  )
  
  ## Compute model connectivity
  hdWGCNA_monaco <- ModuleConnectivity(
    hdWGCNA_monaco,
    group.by = 'monaco_labels', group_name = i,
    wgcna_name = i
  )
  
  ## kme plots
  kme_plots_monaco[[i]] <- PlotKMEs(hdWGCNA_monaco, ncol=5, wgcna_name = i)
  
  ## Get module assignment table
  modules_monaco[[i]] <- GetModules(hdWGCNA_monaco, wgcna_name = i)
  
  ## Gene scoring of top 25 genes in each set
  hdWGCNA_monaco <- ModuleExprScore(
    hdWGCNA_monaco,
    n_genes = 25,
    method='Seurat',
    wgcna_name = i
  )
}




## Visualisation of results (multimodal)-----------------------------------------------------------

## Set up necessary lists
correlogram_list_multi <- list()
hMEs_list_multi <- list()
MEs_list_multi <- list()
mods_list_multi <- list()
module_feature_plots_multi <- list()
module_network_plots_multi <- list()
module_dotplots_list_multi <- list()

## Create output directory for network plots
dir.create("/user/ModuleNetworks")

for(i in cells_multi){
  
  ## Get module eigengenes (harmonized and non-harmonized); add harmonized to Seurat meta-data (needed for dot plots below) after renaming, to ensure that two cell types don't have the same module name.
  hMEs_list_multi[[i]] <- GetMEs(hdWGCNA_multi, wgcna_name = i)
  MEs_list_multi <- GetMEs(hdWGCNA_multi, wgcna_name = i, harmonized=FALSE)
  colnames(hMEs_list_multi[[i]]) <- paste(i, colnames(hMEs_list_multi[[i]]), sep="_")
  mods_list_multi[[i]] <- colnames(hMEs_list_multi[[i]]); mods_list_multi[[i]] <- mods_list_multi[[i]][mods_list_multi[[i]] != 'grey']
  hdWGCNA_multi@meta.data <- cbind(hdWGCNA_multi@meta.data, hMEs_list_multi[[i]])
  
  
  ## Correlation between modules
  correlogram_list_multi[[i]] <- ModuleCorrelogram(hdWGCNA_multi, wgcna_name = i)
  
  ## Cluster plot by cell cell modules (module feature plot)
  
  module_feature_plots_multi[[i]] <- ModuleFeaturePlot(
    hdWGCNA_multi,
    features='hMEs',
    order=TRUE,
    wgcna_name = i,
    reduction="umappreintegration"
  )
  
  #wrap_plots(module_feature_plots_multi[[i]], ncol=6)
  
  ## Module network plot
  ModuleNetworkPlot(hdWGCNA_multi, wgcna_name = i, outdir = paste("/user/ModuleNetworks", i, sep = "/"))
  
}

## Dot plot of cell types vs each other cell type's modules
for(i in cells_multi){
  hdWGCNA_multi <- hdWGCNA_multi ; 
  hdWGCNA_multi@meta.data <- cbind(hdWGCNA_multi@meta.data, hMEs_list_multi[[i]]) ; 
  module_dotplots_list_multi[[i]] <- DotPlot(hdWGCNA_multi, features=mods_list_multi[[i]], group.by = "predicted.celltype.l2") ; 
  module_dotplots_list_multi[[i]] <- module_dotplots_list_multi[[i]] +   coord_flip() + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue')
}

## Save outputs
saveRDS(correlogram_list_multi, file="correlogram_list_multi.rds")
#(opted not to save this one, because took over 100GB space) saveRDS(module_feature_plots_multi, file="module_feature_plots_multi.rds")
saveRDS(module_network_plots_multi, file="module_network_plots_multi.rds")
#(opted not to save this one, because took over 100GB space) saveRDS(module_dotplots_list_multi, file="module_dotplots_list_multi.rds")


## Differential module eigengene (DME) analysis ----------------------------------------------------------------------

cases_list_multi <- list()
controls_list_multi <- list()
DMEs_list_multi <- list()
DME_volcano_plots_multi <- list()

for(i in cells_multi){
  ## Extract cases and controls for each cell type
  cases_list_multi[[i]] <- hdWGCNA_multi@meta.data %>% subset(predicted.celltype.l2 == i & pheno == "case") %>% rownames
  controls_list_multi[[i]] <- hdWGCNA_multi@meta.data %>% subset(predicted.celltype.l2 == i & pheno == "control") %>% rownames
  
  ## Run differential expression command
  DMEs_list_multi[[i]] <- FindDMEs(
    hdWGCNA_multi,
    barcodes1 = cases_list_multi[[i]],
    barcodes2 = controls_list_multi[[i]],
    wgcna_name = i,
    test.use='wilcox',
    add_missing = TRUE
  )
  
  ## Visualize with volcano plot
  DME_volcano_plots_multi[[i]] <- PlotDMEsVolcano(hdWGCNA_multi, DMEs_list_multi[[i]], wgcna_name = i)
}



## Enrichment analysis ----------------------------------------------------------------------

dbs <- c('GO_Biological_Process_2021','KEGG_2021_Human','MSigDB_Hallmark_2020')

#Create output lists for results, and results directory
enrich_list_multi <- list()
enrich_dotplots_multi <- list()
dir.create("/user/enrichr_plots")

for(i in cells_multi){
  # perform enrichment tests
  hdWGCNA_multi <- RunEnrichr(
    hdWGCNA_multi,
    dbs=dbs,
    max_genes = 100,
    wgcna_name = i
  )
  
  
  enrich_list_multi[[i]] <- GetEnrichrTable(hdWGCNA_multi, wgcna_name = i)
  
  #Bar plots of enrichment
  EnrichrBarPlot(
    hdWGCNA_multi,
    outdir = paste("/user/enrichr_plots", i, sep = "/"),
    n_terms = 10,
    plot_size = c(5,7), # width, height of the output .pdfs
    logscale=TRUE,
    wgcna_name = i
  )
  
  #Enrichr dot plots
  enrich_dotplots_multi[[i]] <- EnrichrDotPlot(
    hdWGCNA_multi,
    mods = "all",
    database = "GO_Biological_Process_2021",
    n_terms=1,
    wgcna_name = i
  )
}


## Visualisation of results (monaco)-----------------------------------------------------------

## Set up necessary lists
correlogram_list_monaco <- list()
hMEs_list_monaco <- list()
MEs_list_monaco <- list()
mods_list_monaco <- list()
module_feature_plots_monaco <- list()
module_network_plots_monaco <- list()
module_dotplots_list_monaco <- list()

## Create output directory for network plots
dir.create("/user/ModuleNetworks")

for(i in cells_monaco){
  
  ## Get module eigengenes (harmonized and non-harmonized); add harmonized to Seurat meta-data (needed for dot plots below)
  hMEs_list_monaco[[i]] <- GetMEs(hdWGCNA_monaco, wgcna_name = i)
  MEs_list_monaco <- GetMEs(hdWGCNA_monaco, wgcna_name = i, harmonized=FALSE)
  colnames(hMEs_list_monaco[[i]]) <- paste(i, colnames(hMEs_list_monaco[[i]]), sep="_")
  mods_list_monaco[[i]] <- colnames(hMEs_list_monaco[[i]]); mods_list_monaco[[i]] <- mods_list_monaco[[i]][mods_list_monaco[[i]] != 'grey']
  hdWGCNA_monaco@meta.data <- cbind(hdWGCNA_monaco@meta.data, hMEs_list_monaco[[i]])
  
  
  ## Correlation between modules
  correlogram_list_monaco[[i]] <- ModuleCorrelogram(hdWGCNA_monaco, wgcna_name = i)
  
  ## Cluster plot by cell cell modules (module feature plot)
  
  module_feature_plots_monaco[[i]] <- ModuleFeaturePlot(
    hdWGCNA_monaco,
    features='hMEs',
    order=TRUE,
    wgcna_name = i,
    reduction="umappreintegration"
  )
  
  #wrap_plots(module_feature_plots_monaco[[i]], ncol=6)
  
  ## Module network plot
  ModuleNetworkPlot(hdWGCNA_monaco, wgcna_name = i, outdir = paste("/user/ModuleNetworks", i, sep = "/"))
  
}

## Dot plot of cell types vs each other cell type's modules
for(i in cells_monaco){
  hdWGCNA_monaco <- hdWGCNA_monaco ; 
  hdWGCNA_monaco@meta.data <- cbind(hdWGCNA_monaco@meta.data, hMEs_list_monaco[[i]]) ; 
  module_dotplots_list_monaco[[i]] <- DotPlot(hdWGCNA_monaco, features=mods_list_monaco[[i]], group.by = "monaco_labels") ; 
  module_dotplots_list_monaco[[i]] <- module_dotplots_list_monaco[[i]] +   coord_flip() + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue')
}


## Differential module eigengene (DME) analysis ----------------------------------------------------------------------

cases_list_monaco <- list()
controls_list_monaco <- list()
DMEs_list_monaco <- list()
DME_volcano_plots_monaco <- list()

for(i in cells_monaco){
  ## Extract cases and controls for each cell type
  cases_list_monaco[[i]] <- hdWGCNA_monaco@meta.data %>% subset(monaco_labels == i & pheno == "case") %>% rownames
  controls_list_monaco[[i]] <- hdWGCNA_monaco@meta.data %>% subset(monaco_labels == i & pheno == "control") %>% rownames
  
  ## Run differential expression command
  DMEs_list_monaco[[i]] <- FindDMEs(
    hdWGCNA_monaco,
    barcodes1 = cases_list_monaco[[i]],
    barcodes2 = controls_list_monaco[[i]],
    wgcna_name = i,
    test.use='wilcox',
    add_missing = TRUE
  )
  
  ## Visualize with volcano plot
  DME_volcano_plots_monaco[[i]] <- PlotDMEsVolcano(hdWGCNA_monaco, DMEs_list_monaco[[i]], wgcna_name = i)
}


saveRDS(mtcor_list_monaco, file="mtcor_list_monaco.rds")

## Enrichment analysis ----------------------------------------------------------------------

dbs <- c('GO_Biological_Process_2021','KEGG_2021_Human','MSigDB_Hallmark_2020')

#Create output lists for results, and results directory
enrich_list_monaco <- list()
enrich_dotplots_monaco <- list()
dir.create("/user/enrichr_plots")

for(i in cells_monaco){
  # perform enrichment tests
  hdWGCNA_monaco <- RunEnrichr(
    hdWGCNA_monaco,
    dbs=dbs,
    max_genes = 100,
    wgcna_name = i
  )
  
  
  enrich_list_monaco[[i]] <- GetEnrichrTable(hdWGCNA_monaco, wgcna_name = i)
  
  #Bar plots of enrichment
  EnrichrBarPlot(
    hdWGCNA_monaco,
    outdir = paste("/user/enrichr_plots", i, sep = "/"),
    n_terms = 10,
    plot_size = c(5,7), # width, height of the output .pdfs
    logscale=TRUE,
    wgcna_name = i
  )
  
  #Enrichr dot plots
  enrich_dotplots_monaco[[i]] <- EnrichrDotPlot(
    hdWGCNA_monaco,
    mods = "all",
    database = "GO_Biological_Process_2021",
    n_terms=1,
    wgcna_name = i
  )
}


## Combined outputs ------------------------------------------------------------

##Extract top 10 hub genes for each cell type, combine into one frame and export to Excel
hub_df_multi <- list()
hub_df_monaco <- list()
for(i in cells_multi){
  hub_df_multi[[i]] <- GetHubGenes(hdWGCNA_multi, wgcna_name = i, n_hubs = 10)
}
for(i in cells_monaco){
  hub_df_monaco[[i]] <- GetHubGenes(hdWGCNA_monaco, wgcna_name = i, n_hubs = 10)
}

hub_df_combined <-c(hub_df_multi, hub_df_monaco)
hub_df_combined <- dplyr::bind_rows(hub_df_combined, .id = "index")
write.xlsx(hub_df_combined, "hub_df_combined.xlsx")

##Enrichment analysis results per module for each cell type
enrich_list_combined <- c(enrich_list_multi, enrich_list_monaco)
enrich_list_combined <- dplyr::bind_rows(enrich_list_combined, .id = "index")
write.xlsx(enrich_list_combined, "enrich_list_combined.xlsx")


##Pull differential module eigengenes (all results then only significant results)
DMEs_combined <- c(DMEs_list_multi, DMEs_list_monaco)
DMEs_combined <- dplyr::bind_rows(DMEs_combined, .id = "index")

#Add in top 3 significant gene sets per module
enrich_list_combined_significant <- enrich_list_combined %>% dplyr::filter(Adjusted.P.value <= 0.05)
enrich_list_combined_significant <- enrich_list_combined_significant %>% arrange(index, module, Adjusted.P.value)
enrich_val1 <- enrich_list_combined_significant %>% group_by(index,module) %>% filter(row_number()==1) 
enrich_val2 <- enrich_list_combined_significant %>% group_by(index,module) %>% filter(row_number()==2)
enrich_val3 <- enrich_list_combined_significant %>% group_by(index,module) %>% filter(row_number()==3)
DMEs_combined <- left_join(DMEs_combined, enrich_val1[,c("index", "module", "Term")], by = c("index", "module")) %>% rename(Term1 = Term)
DMEs_combined <- left_join(DMEs_combined, enrich_val2[,c("index", "module", "Term")], by = c("index", "module")) %>% rename(Term2 = Term)
DMEs_combined <- left_join(DMEs_combined, enrich_val3[,c("index", "module", "Term")], by = c("index", "module")) %>% rename(Term3 = Term)

DMEs_combined_significant <- DMEs_combined %>% dplyr::filter(p_val_adj <= 0.05)
write.xlsx(DMEs_combined_significant, "DMEs_combined_significant.xlsx")
write.xlsx(DMEs_combined, "DMEs_combined.xlsx")

#Save volcano plots
ggsave(DME_volcano_plots_monaco[["Th17 cells"]], file="DME_Th17_volcano.png", dpi=300, height=4, width=4)
ggsave(DME_volcano_plots_multi[["CD14 Mono"]], file="DME_CD14Mono_volcano.png", dpi=300, height=4, width=4)
ggsave(DME_volcano_plots_monaco[["Th1_Th17 cells"]], file="DME_Th1_Th17_volcano.png", dpi=300, height=4, width=4)
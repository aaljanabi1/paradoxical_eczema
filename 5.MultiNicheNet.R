#Programming language - R

## Load packages and data files ------------------------------------------------
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)
library(RColorBrewer)
library(openxlsx)

seurat_object <- readRDS("seurat_object.rds")

##Load NicheNet ligand-receptor network and ligand-target matrix
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
options(timeout = max(1000, getOption("timeout"))) #Need to ensure URL doesn't time out for read of large matrix file
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()

##Convert Seurat object to SCE object
sce <- as.SingleCellExperiment(seurat_object, assay = "RNA")

##Convert to most recent version of official gene symbols
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

##Define metadata columns containing group, sample and cell type IDs, and ensure syntactically valid
sample_id    <- "CELL.ID"
group_id     <- "pheno"
celltype_id  <- "celltype"
covariates   <- NA
batches      <- NA
SummarizedExperiment::colData(sce)$CELL.ID = SummarizedExperiment::colData(sce)$CELL.ID %>% make.names()
SummarizedExperiment::colData(sce)$pheno = SummarizedExperiment::colData(sce)$pheno %>% make.names()
SummarizedExperiment::colData(sce)$celltype = SummarizedExperiment::colData(sce)$celltype %>% make.names()

##Define sender and receiver cell types (all cell types for both in this case)
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

## Extract expression information and link ligand expression in senders to receptor expression in receivers---------------------------------------------

##Set minimum number of cells per cluster
min_cells = 10

##Calculate abundance/expression for each cell type/sample/group combination
abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

##View abundance plots to see which cell types are excluded from which donors
abundance_expression_info$abund_plot_sample

## Perform differential expression between receiver and sender cell types---------------------------------------------------

##Define contrasts of interest
contrasts_oi <- c("'case-control','control-case'")
contrast_tbl = tibble(contrast = c("case-control","control-case"), group = c("case","control"))

##Perform DE analysis
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

##Check DE results
DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()
celltype_de = DE_info$celltype_de$de_output_tidy

##Combine DE information for ligand-senders and receptor-receivers
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head(20)

## Predict NicheNet ligand activities and ligand-target links based on DE data---------------------------------------------------

##Define thresholds
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05
p_val_adj = FALSE
top_n_target = 250

##Set number of cores
verbose = TRUE
cores_system = 8
n.cores = 8

##Run the NicheNet ligand activity analysis
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))


##Check DE genes used for activity analysis
ligand_activities_targets_DEgenes$de_genes_df %>% head(20)

##Check the output of the activity analysis
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

## Prioritise sender-ligand:receiver-receptor pairs---------------------------------------------------------------------------------

##Define weights
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

##Make grouping dataframe
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

##Run the prioritisation
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

##Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

## Add information on prior knowledge and expression correlation between LR and target expression-------------------------------------

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

## Save all outputs-------------------------------------------------------------------------------------------------------------------

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

saveRDS(multinichenet_output, "multinichenet_output.rds")

## Visualisation of results-------------------------------------------------------------------------------------------------------------

## Heatmap of ligand activities in each receiver-group combination (i.e. which ligands are most active in which cell types)
ligands_oi = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% inner_join(contrast_tbl) %>% 
  group_by(group, receiver) %>% distinct(ligand, receiver, group, activity) %>% 
  top_n(5, activity) %>% pull(ligand) %>% unique()

plot_oi = make_ligand_activity_plots(multinichenet_output$prioritization_tables, ligands_oi, contrast_tbl, widths = NULL)
plot_oi
ggsave(plot_oi, file="ligand_heatmap.png", dpi=300, height=10, width=21)


##Top 500 ligand receptor interactions in cases (search IFN, TNF, IL7 and others) and then controls
group_oi = "case"
prioritized_tbl_case = get_top_n_lr_pairs(prioritization_tables, 500, groups_oi = group_oi, receivers_oi = receivers_oi)
group_oi = "control"
prioritized_tbl_control = get_top_n_lr_pairs(prioritization_tables, 500, groups_oi = group_oi, receivers_oi = receivers_oi)
group_oi = "case"
write.xlsx(prioritized_tbl_case, "prioritized_tbl_case.xlsx")
write.xlsx(prioritized_tbl_control, "prioritized_tbl_control.xlsx")


##Circos plot of top 50 predictions
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = FALSE)
prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()
nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
colors_sender = mycolors %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors %>% magrittr::set_names(senders_receivers)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

png(file="circos_case.png", res=300, width=1900, height=1900)
circos_list$case
dev.off()
png(file="circos_control.png", res=300, width=1900, height=1900)
circos_list$control
dev.off()
png(file="circos_legend.png", res=300, width=1900, height=1900)
circos_list$legend
dev.off()
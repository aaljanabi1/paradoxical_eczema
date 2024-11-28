#Programming language - R

## Load packages and data files --------------------------------------------

#Load packages
library(dplyr)
library(patchwork)
library(Seurat)
library(tidyverse)
library(sctransform)
library(scds)
library(tidyr)
library(SingleCellExperiment)
library(SeuratDisk)
library(cowplot)
library(tidyseurat)
library(openxlsx)
library(xlsx)
library(lme4)
library(celldex)
library(SingleR)
library(presto)

#Set working directory:
setwd("/user/directory")

#Load Cellranger filtered feature-barcode matrices (one for each batch) and initialise Seurat objects
batch1.data <- Read10X(data.dir = "/user/cellranger/batch1/outs/filtered_feature_bc_matrix/")
batch2.data <- Read10X(data.dir = "/user/cellranger/batch2/outs/filtered_feature_bc_matrix/")
batch1 <- CreateSeuratObject(counts = batch1.data, project = "batch1.apollo", min.cells = 3, min.features = 200)
batch2 <- CreateSeuratObject(counts = batch2.data, project = "batch2.apollo", min.cells = 3, min.features = 200)

#SCDS to identify same-sample heterotypic doublets
batch1.sce <- as.SingleCellExperiment(batch1)
batch2.sce <- as.SingleCellExperiment(batch2)
batch1.sce = bcds(batch1.sce, estNdbl = TRUE)
batch2.sce = bcds(batch2.sce, estNdbl = TRUE)
batch1.sce = cxds(batch1.sce, estNdbl = TRUE)
batch2.sce = cxds(batch2.sce, estNdbl = TRUE)
batch1.sce = cxds_bcds_hybrid(batch1.sce, estNdbl = TRUE)
batch2.sce = cxds_bcds_hybrid(batch2.sce, estNdbl = TRUE)
doublets_batch1 <- as.data.frame(cbind(rownames(colData(batch1.sce)), colData(batch1.sce)$hybrid_score, colData(batch1.sce)$hybrid_call))
doublets_batch2 <- as.data.frame(cbind(rownames(colData(batch2.sce)), colData(batch2.sce)$hybrid_score, colData(batch2.sce)$hybrid_call))
colnames(doublets_batch1) <- c("BARCODE","scds_score","scds_DropletType")
doublets_batch1$scds_DropletType <- gsub("FALSE","singlet",doublets_batch1$scds_DropletType)
doublets_batch1$scds_DropletType <- gsub("TRUE","doublet",doublets_batch1$scds_DropletType)
colnames(doublets_batch2) <- c("BARCODE","scds_score","scds_DropletType")
doublets_batch2$scds_DropletType <- gsub("FALSE","singlet",doublets_batch2$scds_DropletType)
doublets_batch2$scds_DropletType <- gsub("TRUE","doublet",doublets_batch2$scds_DropletType)
batch1_scds_doublet_count <- table(unlist(doublets_batch1))
batch1_scds_doublet_count["doublet"]
batch2_scds_doublet_count <- table(unlist(doublets_batch2))
batch2_scds_doublet_count["doublet"]

#Merge SCDS output with demuxlet output and input as meta data into Seurat object
batch1.best <- read.csv("/user/demuxlet/batch1.best", sep = "\t")
batch2.best <- read.csv("/user/demuxlet/batch2.best", sep = "\t")
batch1.best <- left_join(batch1.best, doublets_batch1, by = "BARCODE")
batch2.best <- left_join(batch2.best, doublets_batch2, by = "BARCODE")
batch1.best$BARCODES = batch1.best$BARCODE
batch2.best$BARCODES = batch2.best$BARCODE
batch1.best <- batch1.best %>% remove_rownames %>% column_to_rownames(var="BARCODE")
batch2.best <- batch2.best %>% remove_rownames %>% column_to_rownames(var="BARCODE")
batch1.best$DROPLET.FINAL = "SNG"
batch1.best$DROPLET.FINAL[batch1.best$DROPLET.TYPE=="DBL"] <- "DBL"
batch1.best$DROPLET.FINAL[batch1.best$DROPLET.TYPE=="AMB"] <- "AMB"
batch1.best$DROPLET.FINAL[batch1.best$scds_DropletType=="doublet"] <- "DBL"
batch2.best$DROPLET.FINAL = "SNG"
batch2.best$DROPLET.FINAL[batch2.best$DROPLET.TYPE=="DBL"] <- "DBL"
batch2.best$DROPLET.FINAL[batch2.best$DROPLET.TYPE=="AMB"] <- "AMB"
batch2.best$DROPLET.FINAL[batch2.best$scds_DropletType=="doublet"] <- "DBL"
batch1.best$CELL.ID = batch1.best$SNG.BEST.GUESS
batch1.best$CELL.ID[batch1.best$DROPLET.FINAL=="DBL"] = "DBL"
batch1.best$CELL.ID[batch1.best$DROPLET.FINAL=="AMB"] = "AMB"
batch2.best$CELL.ID = batch2.best$SNG.BEST.GUESS
batch2.best$CELL.ID[batch2.best$DROPLET.FINAL=="DBL"] = "DBL"
batch2.best$CELL.ID[batch2.best$DROPLET.FINAL=="AMB"] = "AMB"
batch1 <- CreateSeuratObject(counts = batch1.data, project = "batch1.apollo", min.cells = 3, min.features = 200, meta.data = batch1.best)
batch2 <- CreateSeuratObject(counts = batch2.data, project = "batch2.apollo", min.cells = 3, min.features = 200, meta.data = batch2.best)



## QC ----------------------------------------------------------------------

#Visualise features
batch1[["percent.mt"]] <- PercentageFeatureSet(batch1, pattern = "^MT-")
batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT-")
Idents(batch1) <- "orig.ident"
Idents(batch2) <- "orig.ident"
batch1_violinplot1 <- VlnPlot(batch1, features = "nFeature_RNA", split.by = "CELL.ID", cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
batch1_violinplot2 <- VlnPlot(batch1, features = "nCount_RNA", split.by = "CELL.ID", y.max = 50000, cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
batch1_violinplot3 <- VlnPlot(batch1, features = "percent.mt", split.by = "CELL.ID", y.max = 15, cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
batch2_violinplot1 <- VlnPlot(batch2, features = "nFeature_RNA", split.by = "CELL.ID", cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
batch2_violinplot2 <- VlnPlot(batch2, features = "nCount_RNA", split.by = "CELL.ID", y.max = 50000, cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
batch2_violinplot3 <- VlnPlot(batch2, features = "percent.mt", split.by = "CELL.ID", y.max = 15, cols = c("blue", "red", "purple", "yellow", "orange", "brown", "grey", "green"))
plot1 <- FeatureScatter(batch1, feature1 = "nCount_RNA", "percent.mt")
plot2 <- FeatureScatter(batch1, feature1 = "nCount_RNA", "nFeature_RNA")
plot3 <- FeatureScatter(batch2, feature1 = "nCount_RNA", "percent.mt")
plot4 <- FeatureScatter(batch2, feature1 = "nCount_RNA", "nFeature_RNA")

#Identify doublet or ambiguous clusters in batch1
batch1 <- SCTransform(batch1, vars.to.regress = "percent.mt", verbose = FALSE)
batch1 <- RunPCA(batch1, verbose = FALSE)
batch1 <- RunUMAP(object = batch1, dims = 1:20, verbose = FALSE)
batch1 <- FindNeighbors(object = batch1, dims = 1:20, verbose = FALSE)
batch1 <- FindClusters(object = batch1, verbose = FALSE)
batch1_raw_clusters <- DimPlot(batch1, reduction = "umap", label=TRUE)
batch1_raw_clusters_celltype <- DimPlot(batch1, reduction = "umap", split.by = "DROPLET.FINAL", label=TRUE)
batch1_raw_clusters_participants <- DimPlot(batch1, reduction = "umap", split.by = "CELL.ID", label=TRUE)
batch1_raw_clusters_participants2 <- DimPlot(batch1, reduction = "umap", group.by = "CELL.ID")

#Identify doublet or ambiguous clusters in batch2
batch2 <- SCTransform(batch2, vars.to.regress = "percent.mt", verbose = FALSE)
batch2 <- RunPCA(batch2, verbose = FALSE)
batch2 <- RunUMAP(object = batch2, dims = 1:20, verbose = FALSE)
batch2 <- FindNeighbors(object = batch2, dims = 1:20, verbose = FALSE)
batch2 <- FindClusters(object = batch2, verbose = FALSE)
batch2_raw_clusters <- DimPlot(batch2, reduction = "umap", label=TRUE)
batch2_raw_clusters_celltype <- DimPlot(batch2, reduction = "umap", split.by = "DROPLET.FINAL", label=TRUE)
batch2_raw_clusters_participants <- DimPlot(batch2, reduction = "umap", split.by = "CELL.ID", label=TRUE)
batch2_raw_clusters_participants2 <- DimPlot(batch2, reduction = "umap", group.by = "CELL.ID")

#Pull cluster identities and apply to an untransformed Seurat object
clusters_raw_batch1 <- batch1$seurat_clusters 
clusters_raw_batch1 <- as.data.frame(clusters_raw_batch1)
clusters_raw_batch1 <- tibble::rownames_to_column(clusters_raw_batch1, "BARCODES")
clusters_raw_batch2 <- batch2$seurat_clusters
clusters_raw_batch2 <- as.data.frame(clusters_raw_batch2) 
clusters_raw_batch2 <- tibble::rownames_to_column(clusters_raw_batch2, "BARCODES")
batch1.best <- left_join(batch1.best, clusters_raw_batch1, by="BARCODES")
batch2.best <- left_join(batch2.best, clusters_raw_batch2, by="BARCODES")
batch1.best <- batch1.best %>% dplyr::rename(clusters_raw = clusters_raw_batch1)
batch2.best <- batch2.best %>% dplyr::rename(clusters_raw = clusters_raw_batch2)
batch1.best$BARCODE = batch1.best$BARCODES
batch2.best$BARCODE = batch2.best$BARCODES
batch1.best <- batch1.best %>% remove_rownames %>% column_to_rownames(var="BARCODE")
batch2.best <- batch2.best %>% remove_rownames %>% column_to_rownames(var="BARCODE")
batch1_filtered <- CreateSeuratObject(counts = batch1.data, project = "batch1.apollo", min.cells = 3, min.features = 200, meta.data = batch1.best)
batch2_filtered <- CreateSeuratObject(counts = batch2.data, project = "batch2.apollo", min.cells = 3, min.features = 200, meta.data = batch2.best)

#Remove clusters that are predominantly doublet or ambiguous cells
batch1_filtered[["percent.mt"]] <- PercentageFeatureSet(batch1_filtered, pattern = "^MT-")
batch2_filtered[["percent.mt"]] <- PercentageFeatureSet(batch2_filtered, pattern = "^MT-")
batch1_filtered <- subset(x = batch1_filtered, subset = clusters_raw == 7,invert=TRUE)
batch1_filtered <- subset(x = batch1_filtered, subset = clusters_raw == 15,invert=TRUE)
batch1_filtered <- subset(x = batch1_filtered, subset = clusters_raw == 16,invert=TRUE)
batch1_filtered <- subset(x = batch1_filtered, subset = clusters_raw == 19,invert=TRUE)
batch2_filtered <- subset(x = batch2_filtered, subset = clusters_raw == 14, invert=TRUE)
batch2_filtered <- subset(x = batch2_filtered, subset = clusters_raw == 16, invert=TRUE)
batch2_filtered <- subset(x = batch2_filtered, subset = clusters_raw == 18, invert=TRUE)
batch2_filtered <- subset(x = batch2_filtered, subset = clusters_raw == 20, invert=TRUE)
batch2_filtered <- subset(x = batch2_filtered, subset = clusters_raw == 23, invert=TRUE)

#Filter doublets and ambiguous from demuxlet output (identifiable doublets, i.e from different samples)
batch1_filtered <- subset(batch1_filtered, subset = DROPLET.FINAL == "DBL", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = DROPLET.FINAL == "AMB", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 7 & CELL.ID == "Case1", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 9 & CELL.ID == "Case2", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 10 & CELL.ID == "Case1", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 10 & CELL.ID == "Control1", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 10 & CELL.ID == "Control2", invert=TRUE)
batch1_filtered <- subset(batch1_filtered, subset = percent.mt > 9 & CELL.ID == "Control3", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = DROPLET.FINAL == "DBL", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = DROPLET.FINAL == "AMB", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 7 & CELL.ID == "Case1", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 9 & CELL.ID == "Case2", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 10 & CELL.ID == "Case3", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 10 & CELL.ID == "Control1", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 10 & CELL.ID == "Control2", invert=TRUE)
batch2_filtered <- subset(batch2_filtered, subset = percent.mt > 9 & CELL.ID == "Control3", invert=TRUE)

#Merge Seurat objects
seurat_object <- merge(batch1_filtered, y = batch2_filtered, add.cell.ids = c("batch1", "batch2"), project = "AA_combined")

##  Split and reintegrate data (needed for clustering)--------------------------

#Split object by sample ID
sample.list <- SplitObject(seurat_object, split.by = "CELL.ID")
case1 <- sample.list[["Case1"]]
case2 <- sample.list[["Case2"]]
case3 <- sample.list[["Case3"]]
control1 <- sample.list[["Control1"]]
control2 <- sample.list[["Control2"]]
control3 <- sample.list[["Control3"]]

#Transform data (separately for each sample) using sctransform version 2
case1 <- SCTransform(case1, vst.flavor = "v2", verbose = FALSE)
case2 <- SCTransform(case2, vst.flavor = "v2", verbose = FALSE)
case3 <- SCTransform(case3, vst.flavor = "v2", verbose = FALSE)
control1 <- SCTransform(control1, vst.flavor = "v2", verbose = FALSE)
control2 <- SCTransform(control2, vst.flavor = "v2", verbose = FALSE)
control3 <- SCTransform(control3, vst.flavor = "v2", verbose = FALSE)
sample.list <- list(case1 = case1, case2 = case2, case3 = case3, control1 = control1, control2 = control2, control3 = control3)

#Integrate data by sample ID
features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors (object.list = sample.list, normalization.method = "SCT", anchor.features = features)
seurat_object <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:30)
seurat_object <- FindNeighbors(object = seurat_object, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat_object <- FindClusters(object = seurat_object, verbose = FALSE)


## Clustering method 1 - Multimodal reference mapping--------------------------

#Referencing commands
reference <- LoadH5Seurat("/user_directory/pbmc_multimodal.h5seurat")
mapping.anchors <- FindTransferAnchors(reference = reference, query = seurat_object, normalization.method = "SCT", reference.reduction = "spca", dims = 1:50)
seurat_object <- MapQuery(anchorset = mapping.anchors, query = seurat_object, reference = reference, refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", predicted_ADT = "ADT"), reference.reduction = "spca", reduction.model = "wnn.umap")
markers_list <- list()



## Clustering method 2 for Th subsets - SingleR--------------------------

monaco <- celldex::MonacoImmuneData(ensembl=FALSE, cell.ont="nonna")
sce_object <- as.SingleCellExperiment(seurat_object, assay = "RNA")
set.seed(2000)
singler_monaco <- SingleR(sce_object, ref=monaco, labels=monaco$label.fine, assay.type.test=1, aggr.ref=TRUE)
seurat_object <- AddMetaData(seurat_object, singler_monaco$labels, col.name = "monaco_labels")



##Sub cluster T cells and monocytes ------------------------------------

#Repeat splitting, transformation and integration for only these cell types.
subcluster_cell_types <- c("CD4_T", "CD8_T", "Mono")
DefaultAssay(object = seurat_object) <- "RNA"
seurat_object$predicted.celltype.l1 <- gsub(" ", "_", seurat_object$predicted.celltype.l1, fixed=TRUE)
Idents(seurat_object) <- 'predicted.celltype.l1'
sample.list <- list()
seurat_list <- list()
case1_list <- list()
case2_list <- list()
case3_list <- list()
control1_list <- list()
control2_list <- list()
control3_list <- list()
features_list <- list()
anchors_list <- list()

for(i in subcluster_cell_types){
  seurat_list[[i]] <- subset(seurat_object, idents = i)
  sample.list[[i]] <- SplitObject(seurat_list[[i]], split.by = "CELL.ID")
  case1_list[[i]] <- sample.list[[i]][["case1"]]
  case2_list[[i]] <- sample.list[[i]][["case2"]]
  case3_list[[i]] <- sample.list[[i]][["case3"]]
  control1_list[[i]] <- sample.list[[i]][["control1"]]
  control2_list[[i]] <- sample.list[[i]][["control2"]]
  control3_list[[i]] <- sample.list[[i]][["control3"]]
  case1_list[[i]] <- SCTransform(case1_list[[i]], vst.flavor = "v2", verbose = FALSE)
  case2_list[[i]] <- SCTransform(case2_list[[i]], vst.flavor = "v2", verbose = FALSE)
  case3_list[[i]] <- SCTransform(case3_list[[i]], vst.flavor = "v2", verbose = FALSE)
  control1_list[[i]] <- SCTransform(control1_list[[i]], vst.flavor = "v2", verbose = FALSE)
  control2_list[[i]] <- SCTransform(control2_list[[i]], vst.flavor = "v2", verbose = FALSE)
  control3_list[[i]] <- SCTransform(control3_list[[i]], vst.flavor = "v2", verbose = FALSE)
  sample.list[[i]] <- list(case1 = case1_list[[i]], case2 = case2_list[[i]], case3 = case3_list[[i]], control1 = control1_list[[i]], control2 = control2_list[[i]], control3 = control3_list[[i]])
  features_list[[i]] <- SelectIntegrationFeatures(object.list = sample.list[[i]], nfeatures = 3000)
  sample.list[[i]] <- PrepSCTIntegration(object.list = sample.list[[i]], anchor.features = features_list[[i]])
  anchors_list[[i]] <- FindIntegrationAnchors(object.list = sample.list[[i]], normalization.method = "SCT", anchor.features = features_list[[i]])
  seurat_list[[i]] <- IntegrateData(anchorset = anchors_list[[i]], normalization.method = "SCT")
  seurat_list[[i]] <- RunPCA(seurat_list[[i]], verbose = FALSE, reduction.name = paste(i, "pca", sep ="_"))
  seurat_list[[i]] <- RunUMAP(seurat_list[[i]], reduction = paste(i, "pca", sep ="_"), dims = 1:30, reduction.name = paste(i, "umap", sep ="_"))
  seurat_list[[i]] <- FindNeighbors(object = seurat_list[[i]], reduction = paste(i, "pca", sep ="_"), dims = 1:30, verbose = FALSE)
  seurat_list[[i]] <- FindClusters(object = seurat_list[[i]], verbose = FALSE)
}

#Manually assigned sub-clusters based on expression of canonical markers (visualised in violin plots).

#Allocate new cluster identities to Seurat object
new.CD4.ids <- c("CD4 Naive 1", "CD4 Naive 2", "CD4 TCM 1", "CD4 TCM 2", "CD4 Naive 3", "CD4 TCM 3", "Treg-FOXP3hi", "CD4 TEM 1", "CD4 Naive 4", "CD4 Naive 5", "CD4 TEM 2", "Treg-FOXP3lo", "CD4 TCM 4", "CD4 Naive 6", "CD4 Naive 7")
names(new.CD4.ids) <- levels(seurat_list[["CD4_T"]])seurat_list[["CD4_T"]] <- RenameIdents(seurat_list[["CD4_T"]], new.CD4.ids)
new.CD8.ids <- c("CD8 Naive 1", "CD8 Naive 2", "CD8 TEM 1", "CD8 TEM 2", "CD8 TEM 3", "CD8 TCM", "CD8 TEM 4", "CD8 Naive 3", "CD8 TEM 5", "CD8 Naive 4")
names(new.CD8.ids) <- levels(seurat_list[["CD8_T"]])
seurat_list[["CD8_T"]] <- RenameIdents(seurat_list[["CD8_T"]], new.CD8.ids)
new.Mono.ids <- c("CD14 Mono 1", "CD14 Mono 2", "CD14 Mono 3", "CD16 Mono 1", "CD14 Mono 4", "CD14 Mono 5", "CD16 Mono 2", "CD14 Mono 6", "CD14 Mono 7", "CD14 Mono 8", "CD14 Mono 9", "CD16 Mono 3", "CD14 Mono 10")
names(new.Mono.ids) <- levels(seurat_list[["Mono"]])
seurat_list[["Mono"]] <- RenameIdents(seurat_list[["Mono"]], new.Mono.ids)


##Identify and save cluster markers for all clusters and sub-clusters

#Broad clusters
Idents(seurat_object) <- 'predicted.celltype.l1'
broad_cluster_markers_all <- FindAllMarkers(seurat_object, assay = "SCT", only.pos = TRUE)
broad_cluster_markers <- broad_cluster_markers_all %>% group_by(cluster) %>% dplyr::filter(pct.1 > 0.25) %>% dplyr::slice (1:6) %>% ungroup()
broad_cluster_markers <- broad_cluster_markers[order(as.character(broad_cluster_markers$cluster)),c(6,7,2,3,4,1,5)] 

#Th subtypes
Idents(seurat_object) <- 'monaco_labels'
Th_clusters <- c("Th2 cells", "Th17 cells", "Th1/Th17 cells")
Th_markers_all <- FindAllMarkers(subset(seurat_object, assay = "SCT", subset = monaco_labels %in% Th_clusters), only.pos = TRUE)
Th_markers <- Th_markers_all %>% group_by(cluster) %>% dplyr::filter(pct.1 > 0.25) %>% dplyr::slice (1:6) %>% ungroup()
Th_markers <- Th_markers[order(as.character(Th_markers$cluster)),c(6,7,2,3,4,1,5)] 

#Narrow clusters
Idents(seurat_object) <- 'predicted.celltype.l2'
narrow_cluster_markers_all <- FindAllMarkers(seurat_object, assay = "SCT", only.pos = TRUE)
narrow_cluster_markers <- narrow_cluster_markers_all %>% group_by(cluster) %>% dplyr::filter(pct.1 > 0.25) %>% dplyr::slice (1:6) %>% ungroup()
Idents(seurat_object) <- factor(seurat_object@active.ident, sort(decreasing = TRUE, levels(seurat_object@active.ident))) #Re-order clusters alphabetically
narrow_cluster_markers <- narrow_cluster_markers[order(as.character(narrow_cluster_markers$cluster)),c(6,7,2,3,4,1,5)] 

#Sub-clusters
subcluster_markers_all <- list()
subcluster_markers <- list()
mono_levels <- c("CD16 Mono 3", "CD16 Mono 2", "CD16 Mono 1", "CD14 Mono 10", "CD14 Mono 9", "CD14 Mono 8", "CD14 Mono 7", "CD14 Mono 6", "CD14 Mono 5", "CD14 Mono 4", "CD14 Mono 3", "CD14 Mono 2", "CD14 Mono 1")

for(i in subcluster_cell_types){
   Idents(seurat_list[[i]]) <- 'manual_clusters'
   Idents(seurat_list[[i]]) <- factor(seurat_list[[i]]@active.ident, sort(decreasing = TRUE, levels(seurat_list[[i]]@active.ident))) #Re-order clusters alphabetically
   Idents(seurat_list[["Mono"]]) <- factor(seurat_list[["Mono"]]@active.ident, levels = mono_levels)
   subcluster_markers_all[[i]] <- FindAllMarkers(seurat_list[[i]], assay = "SCT", only.pos = TRUE)
   subcluster_markers[[i]] <- subcluster_markers_all[[i]] %>% group_by(cluster) %>% dplyr::filter(pct.1 > 0.25, pct.2 < 0.25) %>% dplyr::slice (1:6) %>% ungroup()
   subcluster_markers[[i]] <- subcluster_markers[[i]][order(subcluster_markers[[i]]$cluster,decreasing=TRUE),]
}
mono_levels <- c("CD14 Mono 1", "CD14 Mono 2", "CD14 Mono 3", "CD14 Mono 4", "CD14 Mono 5", "CD14 Mono 6", "CD14 Mono 7", "CD14 Mono 8", "CD14 Mono 9", "CD14 Mono 10", "CD16 Mono 1", "CD16 Mono 2", "CD16 Mono 3")
subcluster_markers[["Mono"]] <- subcluster_markers[["Mono"]] %>% mutate(cluster=factor(cluster)) %>% mutate(cluster=fct_relevel(cluster, mono_levels)) %>%
  arrange(cluster) #Re-arranging mono factor levels so "CD14 Mono 10" comes after 9

subcluster_markers_joined <- rbind(subcluster_markers[["CD4_T"]], subcluster_markers[["CD8_T"]]) %>% rbind(subcluster_markers[["Mono"]])
subcluster_markers_joined <- subcluster_markers_joined[,c(6,7,2,3,4,1,5)]

#Output markers onto spreadsheet
write.xlsx(broad_cluster_markers, "top_markers_new.xlsx", sheetName="broad")
write.xlsx(main_cluster_markers, "top_markers_new.xlsx", sheetName="narrow", append=TRUE)
write.xlsx(subcluster_markers_joined, "top_markers_new.xlsx", sheetName="subclusters",append=TRUE)
write.xlsx(Th_markers, "top_markers_new.xlsx", sheetName="Th",append=TRUE)


##Create cluster plots ------------------------------------

#Create phenotype variable
seurat_object <- seurat_object %>% tidyseurat::mutate(pheno = if_else(`CELL.ID` %in% c("Case1", "Case2", "Case3"), "case", "control"))
saveRDS(seurat_object, file="seurat_object.rds")

#Broad cluster categories
broad_plot = DimPlot(seurat_object, reduction = "umappreintegration", group.by = "predicted.celltype.l1", split.by = "pheno")
broad_plot <- broad_plot + ggtitle("") + ylab("UMAP-2") + xlab("UMAP-1") + guides(color=guide_legend(override.aes = list(size=2), byrow=TRUE, keyheight=0.2)) + theme(strip.text.x=element_text(size=6), axis.text.x = element_text(size=6), axis.title.x=element_text(size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.title = element_text(size=5), legend.text = element_text(size=5))

#Narrow cluster categories
detailed_plot = DimPlot(seurat_object, reduction = "umappreintegration", group.by = "predicted.celltype.l2")
detailed_plot <- detailed_plot + ggtitle("") + ylab("UMAP-2") + xlab("UMAP-1") + guides(color=guide_legend(override.aes = list(size=2), byrow=TRUE, keyheight=0.2, ncol=2)) + theme(strip.text.x=element_text(size=6), axis.text.x = element_text(size=6), axis.title.x=element_text(size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.title = element_text(size=5), legend.text = element_text(size=5), legend.position="bottom")

#CD4 T cell, CD8 T cell and monocyte sub-clusters
cd4_subclusters <- DimPlot(seurat_list[["CD4_T"]] , reduction = "umappreintegration", group.by = "manual_clusters", label = FALSE)
cd8_subclusters <- DimPlot(seurat_list[["CD8_T"]] , reduction = "umappreintegration", group.by = "manual_clusters", label = FALSE)
mono_subclusters <- DimPlot(seurat_list[["Mono"]] , reduction = "umappreintegration", group.by = "manual_clusters", label = FALSE)
cd4_subclusters <- cd4_subclusters + ggtitle("CD4 T cells") + ylab("UMAP-2") + xlab("UMAP-1") + guides(color=guide_legend(override.aes = list(size=1.5), byrow=TRUE, keyheight=0.1)) + theme(title=element_text(size=6), axis.text.x = element_text(size=6), axis.title.x=element_text(size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.title = element_text(size=5), legend.text = element_text(size=5))
cd8_subclusters <- cd8_subclusters + ggtitle("CD8 T cells") + ylab("UMAP-2") + xlab("UMAP-1") + guides(color=guide_legend(override.aes = list(size=1.5), byrow=TRUE, keyheight=0.2)) + theme(title=element_text(size=6), axis.text.x = element_text(size=6), axis.title.x=element_text(size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.title = element_text(size=5), legend.text = element_text(size=5))
mono_subclusters <- mono_subclusters + ggtitle("Monocytes") + ylab("UMAP-2") + xlab("UMAP-1") + guides(color=guide_legend(override.aes = list(size=1.5), byrow=TRUE, keyheight=0.2)) + theme(title=element_text(size=6), axis.text.x = element_text(size=6), axis.title.x=element_text(size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.title = element_text(size=5), legend.text = element_text(size=5))


##Differential expression analysis---------------------------

#Ensure corrected counts are set properly
DefaultAssay(object = seurat_object) <- "integrated"
seurat_object <- PrepSCTFindMarkers(seurat_object)

#Differential expression data for each cell type, excluding those with very low numbers
Idents(seurat_object) <- 'predicted.celltype.l2'
seurat_object$celltype <- Idents(seurat_object)
seurat_object$celltype.pheno <- paste(Idents(seurat_object), seurat_object$pheno, sep = "_")
Idents(seurat_object) <- "celltype.pheno"
cell_types_filtered <- c("ASDC", "B intermediate", "B memory", "B naive", "CD14 Mono", "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "cDC2", "dnT", "gdT", "HSPC", "ILC", "MAIT", "NK", "NK_CD56bright", "pDC", "Treg")
cell_diff_exp <- list()
for(i in cell_types_filtered) {
  cell_diff_exp[[i]] <- FindMarkers(seurat_object, assay = "SCT", ident.1 = paste(i,"_case", sep = ""), ident.2 = paste(i,"_control", sep = ""), verbose = FALSE)
}

#Save each list to tab in an excel document
cell_diff_exp <- lapply(cell_diff_exp, function (x) rownames_to_column (x, var="gene"))
write.xlsx(cell_diff_exp, "cell_diff_exp.xlsx")

#Differential expression analysis within Th subtypes using SingleR labels-------------

DefaultAssay(object = seurat_object) <- "SCT"
Idents(seurat_object) <- 'monaco_labels'
seurat_object$monaco_celltype <- Idents(seurat_object)
seurat_object$monaco_celltype.pheno <- paste(Idents(seurat_object), seurat_object$pheno, sep = "_")
Idents(seurat_object) <- "monaco_celltype.pheno"
Th_cell_types <- c("Th1/Th17 cells", "Th2 cells", "Th17 cells")
cell_diff_exp_Th <- list()
for(i in Th_cell_types) {
  cell_diff_exp_Th[[i]] <- FindMarkers(seurat_object, assay = "SCT", ident.1 = paste(i,"_case", sep = ""), ident.2 = paste(i,"_control", sep = ""), verbose = FALSE)
}

cell_diff_exp_Th <- lapply(cell_diff_Th, function (x) rownames_to_column (x, var="gene"))
write.xlsx(cell_diff_exp_Th, "cell_diff_exp_Th.xlsx")

##Differential expression analysis within subclusters-------------
cell_types_filtered_subclusters <- list()
for(i in subcluster_cell_types){
  DefaultAssay(object = seurat_list[[i]]) <- "SCT"
  Idents(seurat_list[[i]]) <- 'manual_clusters'
  seurat_list[[i]]$seurat_celltype <- Idents(seurat_list[[i]])
  seurat_list[[i]]$seurat_celltype.pheno <- paste(Idents(seurat_list[[i]]), seurat_list[[i]]$pheno, sep = "_")
  Idents(seurat_list[[i]]) <- "seurat_celltype.pheno"
  cell_types_filtered_subclusters[[i]] <- unique(seurat_list[[i]]$manual_clusters)
  seurat_list[[i]] <- PrepSCTFindMarkers(seurat_list[[i]])
}

cell_diff_exp_CD4 <- list()
cell_diff_exp_CD8 <- list()
cell_diff_exp_Mono <- list()
for(i in cell_types_filtered_subclusters[["CD4_T"]]) {
  cell_diff_exp_CD4[[i]] <- FindMarkers(seurat_list[["CD4_T"]], assay = "SCT", ident.1 = paste(i,"_case", sep = ""), ident.2 = paste(i,"_control", sep = ""), verbose = FALSE)
}

#Remove CD8 Naive 4 as no control cells in this category
cell_types_filtered_subclusters[["CD8_T"]] <- cell_types_filtered_subclusters[["CD8_T"]][-1]

for(i in cell_types_filtered_subclusters[["CD8_T"]]) {
  cell_diff_exp_CD8[[i]] <- FindMarkers(seurat_list[["CD8_T"]], assay = "SCT", ident.1 = paste(i,"_case", sep = ""), ident.2 = paste(i,"_control", sep = ""), verbose = FALSE)
}
for(i in cell_types_filtered_subclusters[["Mono"]]) {
  cell_diff_exp_Mono[[i]] <- FindMarkers(seurat_list[["Mono"]], assay = "SCT", ident.1 = paste(i,"_case", sep = ""), ident.2 = paste(i,"_control", sep = ""), verbose = FALSE)
}

cell_diff_exp_CD4 <- lapply(cell_diff_exp_CD4, function (x) rownames_to_column (x, var="gene"))
cell_diff_exp_CD8 <- lapply(cell_diff_exp_CD8, function (x) rownames_to_column (x, var="gene"))
cell_diff_exp_Mono <- lapply(cell_diff_exp_Mono, function (x) rownames_to_column (x, var="gene"))
DEGs_CD4 <- dplyr::bind_rows(cell_diff_exp_CD4, .id = "index")
DEGs_CD8 <- dplyr::bind_rows(cell_diff_exp_CD8, .id = "index")
DEGs_Mono <- dplyr::bind_rows(cell_diff_exp_Mono, .id = "index")
write.xlsx(DEGs_CD4, "DEGs_CD4.xlsx")
write.xlsx(DEGs_CD8, "DEGs_CD8.xlsx")
write.xlsx(DEGs_Mono, "DEGs_Mono.xlsx")

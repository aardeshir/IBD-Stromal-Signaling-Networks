#!/usr/bin/env Rscript
# 01-preprocessing.R
# Preprocessing and data preparation for the Kinchen IBD stromal cell analysis

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set random seed for reproducibility
set.seed(42)

# Create directories for results and figures if they don't exist
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

# Load Kinchen dataset
message("Loading Kinchen dataset...")
kinchen_data_path <- "data/kinchen_ibd_stromal_seurat.rds"
kinchen_data <- readRDS(kinchen_data_path)

# Basic dataset information
message("Dataset information:")
message("Number of cells: ", ncol(kinchen_data))
message("Number of genes: ", nrow(kinchen_data))

# Apply quality control filters
kinchen_data <- subset(kinchen_data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# Run standard Seurat workflow
kinchen_data <- NormalizeData(kinchen_data)
kinchen_data <- FindVariableFeatures(kinchen_data, selection.method = "vst", nfeatures = 2000)
kinchen_data <- ScaleData(kinchen_data)
kinchen_data <- RunPCA(kinchen_data, features = VariableFeatures(object = kinchen_data))
kinchen_data <- RunUMAP(kinchen_data, dims = 1:30)
kinchen_data <- FindNeighbors(kinchen_data, dims = 1:30)
kinchen_data <- FindClusters(kinchen_data, resolution = 0.8)

# Create UMAP visualization
p1 <- DimPlot(kinchen_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
      ggtitle("Kinchen IBD Stromal Clusters")
p2 <- DimPlot(kinchen_data, reduction = "umap", group.by = "condition") + 
      ggtitle("Kinchen IBD Stromal Conditions")

umaps <- p1 + p2
ggsave("figures/kinchen_umap_overview.png", umaps, width = 12, height = 6, dpi = 300)

# Save processed dataset
saveRDS(kinchen_data, "results/kinchen_data_processed.rds")

message("Preprocessing complete. Processed data saved to results/kinchen_data_processed.rds") 
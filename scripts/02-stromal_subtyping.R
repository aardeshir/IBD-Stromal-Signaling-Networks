#!/usr/bin/env Rscript
# 02-stromal_subtyping.R
# Identification of stromal cell subtypes in the Kinchen IBD dataset

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Define stromal cell markers for subtyping
stromal_markers <- list(
  general = c("COL1A1", "VIM", "DCN", "LUM"),
  inflammatory = c("IL11", "TNFSF14", "IL33", "IL6", "CCL2"),
  myofibroblasts = c("ACTA2", "TAGLN", "MYH11"),
  wnt_secreting = c("WNT2B", "WNT5A", "RSPO3"),
  mesothelial = c("MSLN", "UPK3B", "KRT19")
)

# Function to calculate module scores for each stromal subtype
calculate_module_scores <- function(seurat_obj, markers_list) {
  message("Calculating module scores for stromal subtypes...")
  
  # Check for existing module scores and remove them if present
  meta_cols <- colnames(seurat_obj@meta.data)
  existing_scores <- meta_cols[grepl("^Stromal_", meta_cols)]
  if (length(existing_scores) > 0) {
    message("Removing existing stromal scores...")
    seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% existing_scores]
  }
  
  # Process each stromal subtype
  for (subtype in names(markers_list)) {
    marker_genes <- markers_list[[subtype]]
    genes_present <- marker_genes[marker_genes %in% rownames(seurat_obj)]
    
    if (length(genes_present) == 0) {
      warning(paste0("No genes found for subtype ", subtype, ". Skipping."))
      next
    }
    
    # Calculate module score
    score_name <- paste0("Stromal_", subtype)
    message(paste0("Calculating scores for ", subtype, " using ", length(genes_present), " genes"))
    
    # Handle different Seurat versions gracefully
    tryCatch({
      seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(genes_present),
        name = score_name,
        ctrl = min(100, nrow(seurat_obj) - length(genes_present)),
        seed = 12345
      )
      # Rename the score column (Seurat adds "1" to the end of the name)
      colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == paste0(score_name, "1")] <- score_name
    }, error = function(e) {
      message("Error in AddModuleScore, using manual calculation method...")
      
      # Manual calculation as fallback
      expr_data <- FetchData(seurat_obj, vars = genes_present)
      score <- rowMeans(expr_data)
      seurat_obj@meta.data[[score_name]] <<- score
    })
  }
  
  return(seurat_obj)
}

# Function to classify cells into stromal subtypes
classify_stromal_subtypes <- function(seurat_obj) {
  message("Classifying cells into stromal subtypes...")
  
  # Get score columns
  score_cols <- colnames(seurat_obj@meta.data)[grepl("^Stromal_", colnames(seurat_obj@meta.data))]
  
  if (length(score_cols) == 0) {
    stop("No stromal score columns found. Run calculate_module_scores first.")
  }
  
  # Extract scores
  scores <- seurat_obj@meta.data[, score_cols]
  
  # Determine the subtype with the highest score for each cell
  max_scores <- apply(scores, 1, which.max)
  subtypes <- colnames(scores)[max_scores]
  
  # Extract subtype name (remove "Stromal_" prefix)
  subtypes <- gsub("Stromal_", "", subtypes)
  
  # Add to metadata
  seurat_obj$stromal_subtype <- subtypes
  
  # Count cells by subtype
  subtype_counts <- table(seurat_obj$stromal_subtype)
  message("Stromal subtype counts:")
  print(subtype_counts)
  
  return(seurat_obj)
}

# Load preprocessed Kinchen dataset
message("Loading Kinchen dataset...")
kinchen_data <- readRDS("results/kinchen_data_processed.rds")

# Calculate stromal module scores
kinchen_data <- calculate_module_scores(kinchen_data, stromal_markers)

# Classify cells into stromal subtypes
kinchen_data <- classify_stromal_subtypes(kinchen_data)

# Visualize stromal subtypes
p1 <- DimPlot(kinchen_data, 
             reduction = "umap", 
             group.by = "stromal_subtype", 
             cols = c("general" = "grey", 
                      "inflammatory" = "yellow", 
                      "myofibroblasts" = "blue", 
                      "wnt_secreting" = "darkblue", 
                      "mesothelial" = "lightblue"),
             label = TRUE) + 
      ggtitle("Stromal Subtypes")

p2 <- DimPlot(kinchen_data, 
             reduction = "umap", 
             group.by = "condition") + 
      ggtitle("Condition")

p3 <- FeaturePlot(kinchen_data, 
                 features = "BMP2", 
                 min.cutoff = "q10", 
                 max.cutoff = "q90") + 
      ggtitle("BMP2 Expression")

# Combine plots and save
combined_plot <- (p1 + p2) / p3
ggsave("figures/kinchen_stromal_umap.png", combined_plot, width = 12, height = 10, dpi = 300)

# Create violin plot of BMP2 expression by stromal subtype
p4 <- VlnPlot(kinchen_data, 
             features = "BMP2", 
             group.by = "stromal_subtype", 
             pt.size = 0) + 
      ggtitle("BMP2 Expression by Stromal Subtype") +
      theme(legend.position = "none")

p5 <- VlnPlot(kinchen_data, 
             features = "BMP2", 
             group.by = "condition", 
             pt.size = 0) + 
      ggtitle("BMP2 Expression by Condition") +
      theme(legend.position = "none")

# Combine plots and save
bmp2_violin_plot <- p4 + p5
ggsave("figures/kinchen_bmp2_expression.png", bmp2_violin_plot, width = 10, height = 6, dpi = 300)

# Save updated Seurat object with stromal subtype annotations
saveRDS(kinchen_data, "results/kinchen_data_with_subtypes.rds")

message("Stromal subtyping complete. Updated data saved to results/kinchen_data_with_subtypes.rds") 
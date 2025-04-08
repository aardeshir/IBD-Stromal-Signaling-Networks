#!/usr/bin/env Rscript
# 03-pathway_correlation_analysis.R
# Analysis of correlations between BMP2 and key signaling pathways in IBD

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(reshape2)
library(viridis)

# Define pathway genes to analyze
pathway_genes <- list(
  tgfb = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7"),
  wnt = c("WNT2B", "WNT5A", "RSPO3", "LRP5", "LRP6", "FZD1", "FZD2", "CTNNB1", "LEF1", "TCF7"),
  hedgehog = c("SHH", "IHH", "DHH", "PTCH1", "PTCH2", "SMO", "GLI1", "GLI2", "GLI3"),
  notch = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", "HES1"),
  bmp = c("BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR1B", "BMPR2", "SMAD1", "SMAD5", "SMAD9"),
  ecm = c("COL1A1", "COL3A1", "COL6A1", "COL6A3", "FN1", "LAMA4", "LAMB1", "LAMC1"),
  growth_factors = c("EGF", "PDGFA", "PDGFB", "PDGFC", "VEGFA", "IGF1", "FGF2", "HGF"),
  inflammation = c("IL1B", "IL6", "TNF", "CXCL8", "CCL2", "TGFB1", "IL10", "IL17A", "IL22")
)

# Function to calculate correlations between BMP2 and pathway genes
calculate_pathway_correlations <- function(seurat_obj, pathway_list, output_dir = "results") {
  message("Calculating correlations between BMP2 and pathway genes...")
  
  # Check if BMP2 is in the dataset
  if (!"BMP2" %in% rownames(seurat_obj)) {
    stop("BMP2 gene not found in the dataset")
  }
  
  # Get BMP2 expression
  bmp2_expr <- FetchData(seurat_obj, vars = "BMP2")
  
  # Initialize results list
  pathway_results <- list()
  all_correlations <- data.frame(
    gene = character(),
    pathway = character(),
    correlation = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each pathway
  for (pathway_name in names(pathway_list)) {
    message(paste0("Processing ", pathway_name, " pathway..."))
    
    pathway_genes <- pathway_list[[pathway_name]]
    genes_in_data <- pathway_genes[pathway_genes %in% rownames(seurat_obj)]
    
    if (length(genes_in_data) == 0) {
      warning(paste0("No genes found for pathway ", pathway_name, ". Skipping."))
      next
    }
    
    # Get expression data for pathway genes
    pathway_expr <- FetchData(seurat_obj, vars = genes_in_data)
    
    # Calculate correlations
    correlations <- data.frame(
      gene = character(),
      correlation = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (gene in genes_in_data) {
      # Skip BMP2 correlating with itself
      if (gene == "BMP2") next
      
      # Calculate Spearman correlation
      cor_test <- cor.test(bmp2_expr$BMP2, pathway_expr[[gene]], method = "spearman")
      
      # Add to results
      correlations <- rbind(correlations, data.frame(
        gene = gene,
        correlation = cor_test$estimate,
        p_value = cor_test$p.value,
        stringsAsFactors = FALSE
      ))
    }
    
    # Sort by absolute correlation
    correlations <- correlations[order(abs(correlations$correlation), decreasing = TRUE), ]
    
    # Add pathway info
    correlations$pathway <- pathway_name
    
    # Store results
    pathway_results[[pathway_name]] <- correlations
    all_correlations <- rbind(all_correlations, correlations)
    
    # Create bar plot for this pathway
    if (nrow(correlations) > 0) {
      # Color by correlation sign
      colors <- ifelse(correlations$correlation > 0, "firebrick", "steelblue")
      
      # Create plot
      p <- ggplot(correlations, aes(x = reorder(gene, correlation), y = correlation, fill = correlation > 0)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("steelblue", "firebrick")) +
        theme_minimal() +
        labs(
          title = paste0("BMP2 correlation with ", pathway_name, " pathway genes"),
          x = "Gene",
          y = "Spearman correlation"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      
      # Save plot
      ggsave(
        paste0(output_dir, "/kinchen_", pathway_name, "_correlation.png"),
        p,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  }
  
  # Adjust p-values for multiple testing
  all_correlations$adjusted_p_value <- p.adjust(all_correlations$p_value, method = "BH")
  
  # Add significance indicators
  all_correlations$significance <- ifelse(all_correlations$adjusted_p_value < 0.001, "***",
    ifelse(all_correlations$adjusted_p_value < 0.01, "**",
      ifelse(all_correlations$adjusted_p_value < 0.05, "*", "ns")
    )
  )
  
  # Save overall results
  write.csv(all_correlations, paste0(output_dir, "/bmp2_pathway_correlations.csv"), row.names = FALSE)
  
  # Create combined pathway correlation plot
  # Select top N correlations by absolute value from each pathway
  top_n_per_pathway <- 5
  top_correlations <- data.frame()
  
  for (pathway_name in names(pathway_list)) {
    if (pathway_name %in% unique(all_correlations$pathway)) {
      pathway_data <- all_correlations[all_correlations$pathway == pathway_name, ]
      if (nrow(pathway_data) > 0) {
        # Sort by absolute correlation and take top N
        pathway_data <- pathway_data[order(abs(pathway_data$correlation), decreasing = TRUE), ]
        top_pathway <- pathway_data[1:min(top_n_per_pathway, nrow(pathway_data)), ]
        top_correlations <- rbind(top_correlations, top_pathway)
      }
    }
  }
  
  # Create combined plot
  if (nrow(top_correlations) > 0) {
    # Create faceted plot
    p <- ggplot(top_correlations, aes(x = reorder(gene, correlation), y = correlation, fill = correlation > 0)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("steelblue", "firebrick")) +
      facet_wrap(~pathway, scales = "free_y") +
      theme_minimal() +
      labs(
        title = "BMP2 correlation with key signaling pathways",
        x = "Gene",
        y = "Spearman correlation"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold")
      )
    
    # Save plot
    ggsave(
      paste0(output_dir, "/kinchen_pathway_correlations.png"),
      p,
      width = 12,
      height = 10,
      dpi = 300
    )
  }
  
  return(all_correlations)
}

# Function to create correlation heatmap
create_correlation_heatmap <- function(correlation_data, output_dir = "results") {
  message("Creating correlation heatmap...")
  
  # Prepare data for heatmap
  corr_matrix <- dcast(correlation_data, gene ~ pathway, value.var = "correlation")
  rownames(corr_matrix) <- corr_matrix$gene
  corr_matrix$gene <- NULL
  
  # Replace NA with 0
  corr_matrix[is.na(corr_matrix)] <- 0
  
  # Set up colors
  colors <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
  
  # Create heatmap plot
  heatmap_data <- as.matrix(corr_matrix)
  
  # Create heatmap using ggplot2
  melted_data <- melt(heatmap_data)
  colnames(melted_data) <- c("Gene", "Pathway", "Correlation")
  
  p <- ggplot(melted_data, aes(x = Pathway, y = Gene, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "steelblue",
      mid = "white",
      high = "firebrick",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    theme_minimal() +
    labs(
      title = "BMP2 correlation with pathway genes",
      x = "Pathway",
      y = "Gene"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  # Save plot
  ggsave(
    paste0(output_dir, "/kinchen_correlation_heatmap.png"),
    p,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  # Create BMP2-specific correlation heatmap
  bmp2_correlations <- correlation_data[correlation_data$gene == "BMP2", ]
  
  if (nrow(bmp2_correlations) > 0) {
    p <- ggplot(bmp2_correlations, aes(x = pathway, y = "BMP2", fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "steelblue",
        mid = "white",
        high = "firebrick",
        midpoint = 0,
        limits = c(-1, 1)
      ) +
      theme_minimal() +
      labs(
        title = "BMP2 correlation with pathways",
        x = "Pathway",
        y = ""
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    # Save plot
    ggsave(
      paste0(output_dir, "/kinchen_bmp2_correlation_heatmap.png"),
      p,
      width = 8,
      height = 4,
      dpi = 300
    )
  }
}

# Load Kinchen dataset with stromal subtypes
message("Loading Kinchen dataset with stromal subtypes...")
kinchen_data <- readRDS("results/kinchen_data_with_subtypes.rds")

# Calculate pathway correlations
correlation_results <- calculate_pathway_correlations(kinchen_data, pathway_genes, "figures")

# Create correlation heatmap
create_correlation_heatmap(correlation_results, "figures")

message("Pathway correlation analysis complete. Results saved to figures/ and results/ directories.") 
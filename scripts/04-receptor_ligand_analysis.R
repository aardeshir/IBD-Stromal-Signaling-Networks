#!/usr/bin/env Rscript
# 04-receptor_ligand_analysis.R
# Analysis of BMP2 receptor expression and potential ligand-receptor interactions

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(viridis)

# Define BMP2 receptors
bmp2_receptors <- c("BMPR1A", "BMPR1B", "BMPR2", "ACVR1", "ACVR2A", "ACVR2B")

# Function to analyze ligand-receptor interactions
analyze_ligand_receptor <- function(seurat_obj, ligand = "BMP2", receptors = bmp2_receptors, 
                                   group_var = "stromal_subtype", output_dir = "figures") {
  message(paste0("Analyzing ", ligand, " interactions with receptors..."))
  
  # Check if directory exists, create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Verify ligand and receptors exist in the dataset
  genes_to_check <- c(ligand, receptors)
  genes_present <- genes_to_check[genes_to_check %in% rownames(seurat_obj)]
  missing_genes <- setdiff(genes_to_check, genes_present)
  
  if (length(missing_genes) > 0) {
    message("The following genes are not found in the dataset: ", paste(missing_genes, collapse = ", "))
  }
  
  if (!ligand %in% genes_present) {
    stop(paste0(ligand, " not found in dataset. Cannot perform receptor analysis."))
  }
  
  receptors_present <- intersect(receptors, genes_present)
  if (length(receptors_present) == 0) {
    stop("No receptors found in dataset. Cannot perform receptor analysis.")
  }
  
  message(paste0("Found ", length(receptors_present), " receptors in the dataset: ", 
                 paste(receptors_present, collapse = ", ")))
  
  # Check if group_var exists
  if (!group_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Group variable '", group_var, "' not found in metadata"))
  }
  
  # Get unique groups
  groups <- unique(seurat_obj@meta.data[[group_var]])
  groups <- groups[!is.na(groups)]
  
  # Create violin plots for ligand and receptor expression by group
  genes_to_plot <- c(ligand, receptors_present)
  
  # Plot only cells with valid grouping
  cells_with_group <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data[[group_var]])]
  plot_obj <- subset(seurat_obj, cells = cells_with_group)
  
  # Create violin plots
  message("Creating violin plots for receptor expression...")
  violin_plots <- VlnPlot(plot_obj, 
                         features = genes_to_plot, 
                         group.by = group_var, 
                         pt.size = 0,
                         combine = FALSE)
  
  # Customize plots
  for (i in seq_along(violin_plots)) {
    violin_plots[[i]] <- violin_plots[[i]] + 
      theme(legend.position = "none",
            plot.title = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Combine and save violin plots
  if (length(violin_plots) > 0) {
    ncols <- min(3, length(violin_plots))
    combined_plot <- patchwork::wrap_plots(violin_plots, ncol = ncols)
    
    # Save combined plot
    output_file <- paste0(output_dir, "/", tolower(ligand), "_receptor_expression.png")
    ggsave(output_file, combined_plot, width = 3 * ncols, height = 2.5 * ceiling(length(violin_plots) / ncols), dpi = 300)
    message("Saved receptor expression plot to: ", output_file)
  }
  
  # Calculate average expression by cell type
  message("Calculating average expression by ", group_var, "...")
  avg_expr <- AverageExpression(plot_obj, 
                               features = genes_to_plot, 
                               group.by = group_var,
                               assays = DefaultAssay(plot_obj))
  
  # Convert to data frame for easier handling
  avg_expr_df <- as.data.frame(avg_expr[[DefaultAssay(plot_obj)]])
  
  # Create heatmap of average expression
  if (length(avg_expr_df) > 0) {
    avg_expr_df <- as.data.frame(t(avg_expr_df))
    
    # Melt for ggplot
    melted_data <- reshape2::melt(avg_expr_df)
    colnames(melted_data) <- c("Group", "Gene", "Expression")
    
    # Create heatmap
    p <- ggplot(melted_data, aes(x = Gene, y = Group, fill = Expression)) +
      geom_tile() +
      scale_fill_viridis(option = "plasma") +
      theme_minimal() +
      labs(
        title = paste0("Average expression of ", ligand, " and receptors by ", group_var),
        x = "Gene",
        y = group_var
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    # Save heatmap
    output_file <- paste0(output_dir, "/", tolower(ligand), "_receptor_heatmap.png")
    ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    message("Saved receptor expression heatmap to: ", output_file)
  }
  
  # Calculate interaction scores between ligand-expressing and receptor-expressing cell types
  message("Calculating interaction scores...")
  interaction_scores <- matrix(0, nrow = length(groups), ncol = length(groups))
  rownames(interaction_scores) <- groups
  colnames(interaction_scores) <- groups
  
  # Get ligand expression by group
  ligand_expr <- avg_expr_df[[ligand]]
  
  # For each receptor, calculate interaction scores
  receptor_scores <- list()
  for (receptor in receptors_present) {
    # Get receptor expression by group
    receptor_expr <- avg_expr_df[[receptor]]
    
    # Calculate interaction score for each source-target pair
    score_matrix <- matrix(0, nrow = length(groups), ncol = length(groups))
    rownames(score_matrix) <- groups
    colnames(score_matrix) <- groups
    
    for (source_group in groups) {
      for (target_group in groups) {
        # Interaction score = ligand expression in source * receptor expression in target
        score_matrix[source_group, target_group] <- 
          ligand_expr[source_group] * receptor_expr[target_group]
      }
    }
    
    # Store scores for this receptor
    receptor_scores[[receptor]] <- score_matrix
    
    # Add to total interaction scores
    interaction_scores <- interaction_scores + score_matrix
  }
  
  # Create interaction heatmap for total interactions
  melted_scores <- reshape2::melt(interaction_scores)
  colnames(melted_scores) <- c("Source", "Target", "Interaction")
  
  p <- ggplot(melted_scores, aes(x = Target, y = Source, fill = Interaction)) +
    geom_tile() +
    scale_fill_viridis() +
    theme_minimal() +
    labs(
      title = paste0(ligand, "-Receptor Interaction Potential"),
      x = paste0(group_var, " (Receptor-expressing)"),
      y = paste0(group_var, " (", ligand, "-expressing)")
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  # Save interaction heatmap
  output_file <- paste0(output_dir, "/", tolower(ligand), "_interaction_potential.png")
  ggsave(output_file, p, width = 8, height = 7, dpi = 300)
  message("Saved interaction heatmap to: ", output_file)
  
  # Create receptor-specific interaction heatmaps
  if (length(receptors_present) > 1) {
    receptor_plots <- list()
    
    for (receptor in receptors_present) {
      melted_receptor <- reshape2::melt(receptor_scores[[receptor]])
      colnames(melted_receptor) <- c("Source", "Target", "Interaction")
      
      receptor_plots[[receptor]] <- ggplot(melted_receptor, aes(x = Target, y = Source, fill = Interaction)) +
        geom_tile() +
        scale_fill_viridis() +
        theme_minimal() +
        labs(
          title = paste0(ligand, "-", receptor, " Interaction"),
          x = paste0(group_var, " (", receptor, "-expressing)"),
          y = paste0(group_var, " (", ligand, "-expressing)")
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10)
        )
    }
    
    # Combine receptor plots
    ncols <- min(3, length(receptor_plots))
    combined_receptor_plot <- patchwork::wrap_plots(receptor_plots, ncol = ncols)
    
    # Save combined receptor interaction plots
    output_file <- paste0(output_dir, "/", tolower(ligand), "_receptor_interactions.png")
    ggsave(output_file, combined_receptor_plot, 
          width = 4 * ncols, 
          height = 3.5 * ceiling(length(receptor_plots) / ncols), 
          dpi = 300)
    message("Saved receptor-specific interaction plots to: ", output_file)
  }
  
  # Return results
  return(list(
    receptors_present = receptors_present,
    avg_expression = avg_expr_df,
    interaction_scores = interaction_scores,
    receptor_scores = receptor_scores
  ))
}

# Load Kinchen dataset with stromal subtypes
message("Loading Kinchen dataset...")
kinchen_data <- readRDS("results/kinchen_data_with_subtypes.rds")

# Run receptor-ligand analysis
message("Starting receptor-ligand interaction analysis for BMP2 signaling...")
result <- analyze_ligand_receptor(
  kinchen_data,
  ligand = "BMP2",
  receptors = bmp2_receptors,
  group_var = "stromal_subtype",
  output_dir = "figures"
)

# Also analyze by condition
if ("condition" %in% colnames(kinchen_data@meta.data)) {
  message("Analyzing BMP2 receptor interactions by condition...")
  condition_result <- analyze_ligand_receptor(
    kinchen_data,
    ligand = "BMP2",
    receptors = bmp2_receptors,
    group_var = "condition",
    output_dir = "figures"
  )
}

message("Receptor-ligand analysis complete. Results saved to figures/ directory.") 
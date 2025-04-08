#!/usr/bin/env Rscript
# 05-statistical_tests.R
# Statistical tests for BMP2 expression and signaling in IBD

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(broom)

# Function to perform differential expression tests for BMP2
perform_bmp2_diff_expression_tests <- function(seurat_obj) {
  message("Performing statistical tests on BMP2 expression...")
  
  # Extract metadata
  meta_data <- seurat_obj@meta.data
  
  # 1. Test BMP2 expression difference between conditions (Healthy vs UC)
  if("condition" %in% colnames(meta_data)) {
    # Get BMP2 expression
    bmp2_expr <- FetchData(seurat_obj, vars = "BMP2")
    
    # Combine with condition
    test_data <- cbind(bmp2_expr, condition = meta_data$condition)
    
    # Wilcoxon test
    wilcox_result <- wilcox.test(BMP2 ~ condition, data = test_data)
    
    # Format results
    condition_test_result <- data.frame(
      test = "Wilcoxon rank-sum",
      comparison = "BMP2 expression: Healthy vs UC",
      statistic = wilcox_result$statistic,
      p_value = wilcox_result$p.value,
      significance = ifelse(wilcox_result$p.value < 0.001, "***", 
                     ifelse(wilcox_result$p.value < 0.01, "**",
                     ifelse(wilcox_result$p.value < 0.05, "*", "ns")))
    )
    
    print(condition_test_result)
  }
  
  # 2. Test BMP2 expression differences across stromal subtypes
  if("stromal_subtype" %in% colnames(meta_data)) {
    # Get BMP2 expression
    bmp2_expr <- FetchData(seurat_obj, vars = "BMP2")
    
    # Combine with stromal subtype
    test_data <- cbind(bmp2_expr, stromal_subtype = meta_data$stromal_subtype)
    
    # Filter out NAs
    test_data <- test_data[!is.na(test_data$stromal_subtype),]
    
    # Kruskal-Wallis test
    kruskal_result <- kruskal.test(BMP2 ~ stromal_subtype, data = test_data)
    
    # Format results
    subtype_test_result <- data.frame(
      test = "Kruskal-Wallis",
      comparison = "BMP2 expression across stromal subtypes",
      statistic = kruskal_result$statistic,
      p_value = kruskal_result$p.value,
      significance = ifelse(kruskal_result$p.value < 0.001, "***", 
                     ifelse(kruskal_result$p.value < 0.01, "**",
                     ifelse(kruskal_result$p.value < 0.05, "*", "ns")))
    )
    
    print(subtype_test_result)
    
    # Pairwise Wilcoxon tests for subtypes
    if(length(unique(test_data$stromal_subtype)) > 1) {
      pairwise_result <- pairwise.wilcox.test(
        test_data$BMP2, test_data$stromal_subtype,
        p.adjust.method = "BH"
      )
      
      # Print pairwise results
      cat("\nPairwise comparisons between stromal subtypes (p-values with BH correction):\n")
      print(pairwise_result$p.value)
    }
  }
  
  # Return list of test results
  return(list(
    condition_test = if(exists("condition_test_result")) condition_test_result else NULL,
    subtype_test = if(exists("subtype_test_result")) subtype_test_result else NULL
  ))
}

# Function to test significance of BMP2 correlations with pathway genes
test_bmp2_pathway_correlations <- function(seurat_obj, pathway_genes, min_cells_expressing = 10) {
  message("Testing significance of BMP2 correlations with pathway genes...")
  
  # Get expression data for BMP2 and pathway genes
  genes_to_fetch <- c("BMP2", pathway_genes)
  
  # Check which genes are in the dataset
  genes_in_data <- genes_to_fetch[genes_to_fetch %in% rownames(seurat_obj)]
  missing_genes <- setdiff(genes_to_fetch, genes_in_data)
  
  if(length(missing_genes) > 0) {
    message(paste0("Genes not found in dataset: ", paste(missing_genes, collapse=", ")))
  }
  
  if(!"BMP2" %in% genes_in_data) {
    stop("BMP2 not found in dataset")
  }
  
  # Fetch expression data
  expr_data <- FetchData(seurat_obj, vars = genes_in_data)
  
  # Filter genes with too few cells expressing
  genes_to_analyze <- c()
  for(gene in genes_in_data) {
    if(sum(expr_data[[gene]] > 0) >= min_cells_expressing) {
      genes_to_analyze <- c(genes_to_analyze, gene)
    } else {
      message(paste0(gene, " is expressed in fewer than ", min_cells_expressing, " cells and will be excluded"))
    }
  }
  
  # If BMP2 is not expressed in enough cells, stop
  if(!"BMP2" %in% genes_to_analyze) {
    stop("BMP2 is not expressed in enough cells")
  }
  
  # Calculate correlations and p-values
  correlation_results <- data.frame(
    gene = character(),
    correlation = numeric(),
    p_value = numeric(),
    adjusted_p_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  for(gene in setdiff(genes_to_analyze, "BMP2")) {
    cor_test <- cor.test(expr_data$BMP2, expr_data[[gene]], method = "spearman")
    
    correlation_results <- rbind(correlation_results, data.frame(
      gene = gene,
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      adjusted_p_value = NA,  # Will be filled later
      significance = NA       # Will be filled later
    ))
  }
  
  # Adjust p-values for multiple testing
  if(nrow(correlation_results) > 0) {
    correlation_results$adjusted_p_value <- p.adjust(correlation_results$p_value, method = "BH")
    
    # Add significance indicators
    correlation_results$significance <- ifelse(correlation_results$adjusted_p_value < 0.001, "***", 
                                      ifelse(correlation_results$adjusted_p_value < 0.01, "**",
                                      ifelse(correlation_results$adjusted_p_value < 0.05, "*", "ns")))
    
    # Sort by absolute correlation
    correlation_results <- correlation_results[order(abs(correlation_results$correlation), decreasing = TRUE),]
  }
  
  return(correlation_results)
}

# Function to test proportion of BMP2+ cells between conditions or subtypes
test_bmp2_positive_proportions <- function(seurat_obj, group_var, bmp2_threshold = 0) {
  message(paste0("Testing differences in BMP2+ cell proportions by ", group_var, "..."))
  
  # Check if group variable exists
  if(!group_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Group variable ", group_var, " not found in metadata"))
  }
  
  # Get BMP2 expression
  bmp2_expr <- FetchData(seurat_obj, vars = "BMP2")
  
  # Get grouping variable
  groups <- seurat_obj@meta.data[[group_var]]
  
  # Combine data
  test_data <- data.frame(
    BMP2 = bmp2_expr$BMP2,
    group = groups
  )
  
  # Remove NAs
  test_data <- test_data[!is.na(test_data$group),]
  
  # Mark cells as BMP2+ or BMP2-
  test_data$BMP2_status <- ifelse(test_data$BMP2 > bmp2_threshold, "BMP2+", "BMP2-")
  
  # Create contingency table
  contingency_table <- table(test_data$group, test_data$BMP2_status)
  print(contingency_table)
  
  # Chi-square test if all expected frequencies >= 5
  expected <- chisq.test(contingency_table)$expected
  
  if(all(expected >= 5)) {
    test_result <- chisq.test(contingency_table)
    test_name <- "Chi-square"
  } else {
    # Fisher's exact test for small expected frequencies
    test_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
    test_name <- "Fisher's exact"
  }
  
  # Format results
  proportion_test_result <- data.frame(
    test = test_name,
    comparison = paste0("Proportion of BMP2+ cells by ", group_var),
    statistic = ifelse(test_name == "Chi-square", test_result$statistic, NA),
    p_value = test_result$p.value,
    significance = ifelse(test_result$p.value < 0.001, "***", 
                   ifelse(test_result$p.value < 0.01, "**",
                   ifelse(test_result$p.value < 0.05, "*", "ns")))
  )
  
  # Calculate percentages
  percentages <- prop.table(contingency_table, margin = 1) * 100
  percentages <- as.data.frame(percentages)
  colnames(percentages) <- c("Group", "BMP2_Status", "Percentage")
  
  # Filter for BMP2+ cells
  bmp2_pos_percentages <- percentages[percentages$BMP2_Status == "BMP2+",]
  
  return(list(
    test_result = proportion_test_result,
    contingency_table = contingency_table,
    percentages = bmp2_pos_percentages
  ))
}

# Function to run all statistical tests
run_all_statistical_tests <- function(seurat_obj) {
  message("Running all statistical tests for BMP2 analysis...")
  
  # Create results directory if it doesn't exist
  if(!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
  }
  
  # 1. Differential expression tests
  diff_expr_results <- perform_bmp2_diff_expression_tests(seurat_obj)
  saveRDS(diff_expr_results, "results/bmp2_diff_expression_tests.rds")
  
  # 2. Pathway correlation tests
  # Define pathway genes
  tgfb_genes <- c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7")
  wnt_genes <- c("WNT2B", "WNT5A", "RSPO3", "LRP5", "LRP6", "FZD1", "FZD2", "CTNNB1", "LEF1", "TCF7")
  bmp_genes <- c("BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR1B", "BMPR2", "SMAD1", "SMAD5", "SMAD9")
  ecm_genes <- c("COL1A1", "COL3A1", "COL6A1", "COL6A3", "FN1", "LAMA4", "LAMB1", "LAMC1")
  growth_factors <- c("EGF", "PDGFA", "PDGFB", "PDGFC", "VEGFA", "IGF1", "FGF2", "HGF")
  
  # Combine all genes
  all_pathway_genes <- unique(c(tgfb_genes, wnt_genes, bmp_genes, ecm_genes, growth_factors))
  
  # Test correlations
  correlation_results <- test_bmp2_pathway_correlations(seurat_obj, all_pathway_genes)
  saveRDS(correlation_results, "results/bmp2_correlation_tests.rds")
  
  # Also save as CSV for easier viewing
  write.csv(correlation_results, "results/bmp2_correlation_tests.csv", row.names = FALSE)
  message("Saved correlation test results to: results/bmp2_correlation_tests.csv")
  
  # 3. Test BMP2+ proportions by condition
  if("condition" %in% colnames(seurat_obj@meta.data)) {
    proportion_by_condition <- test_bmp2_positive_proportions(seurat_obj, "condition")
    saveRDS(proportion_by_condition, "results/bmp2_positive_proportion_by_condition.rds")
  }
  
  # 4. Test BMP2+ proportions by stromal subtype
  if("stromal_subtype" %in% colnames(seurat_obj@meta.data)) {
    proportion_by_subtype <- test_bmp2_positive_proportions(seurat_obj, "stromal_subtype")
    saveRDS(proportion_by_subtype, "results/bmp2_positive_proportion_by_subtype.rds")
  }
  
  message("All statistical tests completed and results saved to results/ directory.")
}

# Load Kinchen dataset
message("Loading Kinchen dataset...")
kinchen_data <- readRDS("results/kinchen_data_with_subtypes.rds")

# Run all statistical tests
message("Running statistical analysis on Kinchen dataset...")
run_all_statistical_tests(kinchen_data)

message("Statistical analysis complete. Results saved to results/ directory.") 
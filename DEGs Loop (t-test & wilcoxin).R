library(stats)
library(ggplot2)
library(ggrepel)

process_samples <- function(normal_file, disease_file, output_prefix) {
  # Load the normal data
  env <- new.env()
  load(normal_file, envir = env)
  normal_objects <- mget(ls(envir = env), envir = env)
  Normal <- normal_objects[[1]]
  
  # Load the disease data
  Disease <- readRDS(disease_file)
  
  # Extract the gene symbols
  colnames(Disease) <- paste("Disease", seq_along(Disease), sep = "_")
  gene_symbols <- rownames(Disease)
  
  # Ensure only genes present in both datasets are used
  Normal <- Normal[rownames(Normal) %in% gene_symbols, ]
  Disease <- Disease[rownames(Disease) %in% gene_symbols, ]
  
  # Initialize results storage
  results <- data.frame(Gene = character(),
                        Normal_Shapiro_p = numeric(),
                        Disease_Shapiro_p = numeric(),
                        Test = character(),
                        Statistic = numeric(),
                        p_value = numeric(),
                        Adjusted_p_value = numeric(),
                        Log_FC = numeric(),
                        stringsAsFactors = FALSE)
  
  # Perform Shapiro-Wilk test, appropriate test (t-test or Wilcoxon test), and log fold change for each gene
  for (gene in gene_symbols) {
    normal_values <- as.numeric(Normal[gene, ])
    disease_values <- as.numeric(Disease[gene, ])
    
    # Ensure there are no non-finite values
    if (any(!is.finite(normal_values)) | any(!is.finite(disease_values))) next
    
    # Shapiro-Wilk test for normality
    normal_shapiro_test <- shapiro.test(normal_values)
    disease_shapiro_test <- shapiro.test(disease_values)
    
    # Determine which test to use based on normality
    if (normal_shapiro_test$p.value >= 0.05 & disease_shapiro_test$p.value >= 0.05) {
      # Perform t-test if data is normally distributed
      test <- t.test(normal_values, disease_values, paired = FALSE)
      test_type <- "t-test"
    } else {
      # Perform Wilcoxon test if data is not normally distributed
      test <- wilcox.test(normal_values, disease_values, paired = FALSE)
      test_type <- "Wilcoxon test"
    }
    
    # Calculate log fold change, adding a small constant to avoid division by zero
    log_fc <- log2((mean(disease_values) + 1) / (mean(normal_values) + 1))
    
    # Store the results
    results <- rbind(results, data.frame(Gene = gene,
                                         Normal_Shapiro_p = normal_shapiro_test$p.value,
                                         Disease_Shapiro_p = disease_shapiro_test$p.value,
                                         Test = test_type,
                                         Statistic = as.numeric(test$statistic),
                                         p_value = test$p.value,
                                         Adjusted_p_value = NA,  # Placeholder for now
                                         Log_FC = log_fc,
                                         stringsAsFactors = FALSE))
  }
  
  # Perform multiple testing correction using Benjamini-Hochberg method
  results$Adjusted_p_value <- p.adjust(results$p_value, method = "BH")
  
  # Save the results
  write.csv(results, paste0(output_prefix, "_test_results_with_correction.csv"), row.names = FALSE)
  
  # Filter significant genes
  significant_genes <- subset(results, Adjusted_p_value < 0.05 & abs(Log_FC) > 1)
  
  # Save the significant results
  write.csv(significant_genes, paste0(output_prefix, "_significant_genes.csv"), row.names = FALSE)
  
  # Print results
  print(results)
  
  # Count upregulated and downregulated genes
  upregulated_genes <- sum(significant_genes$Log_FC > 1)
  downregulated_genes <- sum(significant_genes$Log_FC < -1)
  
  # Create a volcano plot
  volcano_plot <- ggplot(results, aes(x = Log_FC, y = -log10(Adjusted_p_value), label = Gene)) +
    geom_point(aes(color = ifelse(Log_FC > 1 & Adjusted_p_value < 0.05, "Upregulated",
                                  ifelse(Log_FC < -1 & Adjusted_p_value < 0.05, "Downregulated", "Not Significant"))), alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = "Volcano Plot", 
         x = "Log2 Fold Change", 
         y = "-log10 Adjusted P-value",
         subtitle = paste("Upregulated genes:", upregulated_genes, 
                          "Downregulated genes:", downregulated_genes)) +
    theme(legend.position = "none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(results, Adjusted_p_value < 0.05 & abs(Log_FC) > 1), 
                    aes(label = Gene), size = 3)
  
  # Save the volcano plot
  ggsave(paste0(output_prefix, "_volcano_plot.png"), plot = volcano_plot, bg = "white")
  
  
  # Display the volcano plot
  print(volcano_plot)
}

# Define your samples and file names
samples <- list(
  list(normal_file = "LumA_Normal.rda", disease_file = "LumA_Predicted.rda", output_prefix = "LumA"),
  list(normal_file = "LumB_Normal.rda", disease_file = "LumB_Predicted.rda", output_prefix = "LumB")
)

# Process each sample
for (sample in samples) {
  process_samples(sample$normal_file, sample$disease_file, sample$output_prefix)
}

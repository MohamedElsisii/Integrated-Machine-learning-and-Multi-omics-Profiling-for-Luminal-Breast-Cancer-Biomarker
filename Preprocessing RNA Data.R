#1] Import libraries
library(edgeR)
library(biomaRt)
library(dplyr)
library(e1071)
library(limma)
############################################
#2] Functions
# Main function to preprocess the RNA-seq count matrix
preprocess_count_matrix <- function(file_path, dataset = "hsapiens_gene_ensembl") {
  # Load the RNA-Seq count matrix
  count_matrix <- load_count_matrix(file_path)
  
  # Map Ensembl IDs to gene symbols and aggregate duplicate symbols
  count_matrix <- map_and_aggregate(count_matrix, dataset)
  
  # Filter lowly-expressed genes and handle any potential zeros
  filtered_log_count_matrix <- filter_lowly_expressed_genes_and_handle_zeros(count_matrix)
  
  # Handle zeros if there is any
  if (sum(filtered_log_count_matrix == 0) != 0){
    cat("Processing warning: Zeros detected in the dataset after processing.")
  }
  
  # Adjust for outliers
  final_count_matrix <- t(apply(filtered_log_count_matrix, 1, adjust_outliers))
  
  return(final_count_matrix)
}

# Function to load count matrix from CSV or TSV
load_count_matrix <- function(file_path) {
  file_extension <- tools::file_ext(file_path)
  
  if (file_extension == "csv") {
    count_matrix <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (file_extension == "tsv" || file_extension == "txt") {
    count_matrix <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported file type. Please provide a CSV or TSV file.")
  }
  
  return(count_matrix)
}

# Function to map Ensembl IDs to gene symbols and aggregate duplicates
map_and_aggregate <- function(count_matrix, dataset) {
  ensembl_ids <- sub("\\..*$", "", count_matrix[, 1]) 
  ensembl <- biomaRt::useMart("ensembl", dataset = dataset)
  
  gene_symbols <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                                 filters = 'ensembl_gene_id', 
                                 values = ensembl_ids, 
                                 mart = ensembl)
  
  gene_symbols$external_gene_name[gene_symbols$external_gene_name == ""] <- NA
  id_to_symbol_map <- setNames(gene_symbols$external_gene_name, gene_symbols$ensembl_gene_id)
  count_matrix[, 1] <- id_to_symbol_map[ensembl_ids]
  count_matrix <- count_matrix[!is.na(count_matrix[, 1]), ]
  
  # Aggregate duplicate gene symbols
  colnames(count_matrix)[1] <- "Gene"
  if (anyDuplicated(count_matrix[, "Gene"])) {
    count_matrix <- count_matrix %>%
      group_by(Gene) %>%
      summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')
  }
  
  # Convert to a matrix for further analysis
  count_matrix <- as.data.frame(count_matrix)
  rownames(count_matrix) <- count_matrix$Gene
  count_matrix <- as.matrix(count_matrix[, -1])
  
  return(count_matrix)
}

# Function to filter lowly-expressed genes
filter_lowly_expressed_genes_and_handle_zeros <- function(count_matrix){
  # Normalize the data using the TMM method from edgeR
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  log_count_matrix <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Calculate median expression level across genes and median absolute deviation (MAD) for each gene
  median_expression <- apply(log_count_matrix, 1, median)
  mad_expression <- apply(log_count_matrix, 1, mad)
  
  # Dynamically set mad_multiplier based on skewness
  data_skewness <- skewness(median_expression)
  mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))
  
  # Combine median and MAD for a robust measure of gene expression variability
  threshold_expression <- median_expression + (mad_multiplier * mad_expression)
  
  # Filter genes with low expression variability
  keep_genes <- mad_expression > quantile(mad_expression, 0.25)  # Keep genes with highest variability
  filtered_log_count_matrix <- log_count_matrix[keep_genes, ]
  
  # RNA-seq data typically should not have missing values post-normalization
  if(any(is.na(filtered_log_count_matrix))) {
    print ("hi")
    # Dynamic thresholding based on the distribution of missing values
    missing_per_gene <- colSums(is.na(filtered_log_count_matrix)) / nrow(filtered_log_count_matrix)
    missing_per_sample <- rowSums(is.na(filtered_log_count_matrix)) / ncol(filtered_log_count_matrix)
    gene_missing_threshold <- mean(missing_per_gene) + sd(missing_per_gene)
    sample_missing_threshold <- mean(missing_per_sample) + sd(missing_per_sample)
    
    # Filter based on missing data thresholds
    filtered_log_count_matrix <- filtered_log_count_matrix[rowSums(is.na(filtered_log_count_matrix)) / ncol(filtered_log_count_matrix) < sample_missing_threshold, ]
    filtered_log_count_matrix <- filtered_log_count_matrix[, colSums(is.na(filtered_log_count_matrix)) / nrow(filtered_log_count_matrix) < gene_missing_threshold]
  }
  
  return(filtered_log_count_matrix)
}

# Function to adjust for outliers
adjust_outliers <- function(filtered_log_count_matrix) {
  calculate_outlier_bounds <- function(filtered_log_count_matrix, k) {
    Q1 <- quantile(filtered_log_count_matrix, 0.25, na.rm = TRUE)
    Q3 <- quantile(filtered_log_count_matrix, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(c(lower_bound, upper_bound))
  }
  
  # Initial outlier detection with standard k=1.5
  k <- 1.5
  bounds <- calculate_outlier_bounds(filtered_log_count_matrix, k)
  outlier_proportion <- mean(filtered_log_count_matrix < bounds[1] | filtered_log_count_matrix > bounds[2], na.rm = TRUE)
  
  # Adjust k based on the outlier proportion
  target_proportion <- 0.05 # Targeting 5% outliers
  while (outlier_proportion > target_proportion && k < 10) {
    k <- k + 0.5
    bounds <- calculate_outlier_bounds(filtered_log_count_matrix, k)
    outlier_proportion <- mean(filtered_log_count_matrix < bounds[1] | filtered_log_count_matrix > bounds[2], na.rm = TRUE)
  }
  
  # Adjust the data points considered as outliers
  adjusted_matrix <- ifelse(filtered_log_count_matrix < bounds[1], bounds[1], ifelse(filtered_log_count_matrix > bounds[2], bounds[2], filtered_log_count_matrix))
  
  return(adjusted_matrix)
}
###############################################
#3] Loop on downloaded files and then save them as .rda

# List all files in the directory
files_in_directory <- list.files()

# Filter files with names matching the pattern "cancer name_RNA"
# Adjust the pattern if needed
pattern <- "_RNA\\.csv$" # This assumes your files end with "_RNA.csv"
files_to_process <- grep(pattern, files_in_directory, value = TRUE)

# Loop over the filtered files, process, and save them
for (file_name in files_to_process) {
  # Construct the name of the RDA file
  rda_file_name <- sub("\\.csv$", ".rda", file_name) # Replace .csv with .rda
  
  # Check if the RDA file already exists
  if (file.exists(rda_file_name)) {
    cat("Skipping", file_name, "- RDA file already exists.\n")
    next
  }
  
  # Process the count matrix for each file
  processed_count_matrix <- preprocess_count_matrix(file_name)
  
  # Save the processed count matrix as an RDA file
  saveRDS(processed_count_matrix, rda_file_name)
}
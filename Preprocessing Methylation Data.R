# Load necessary libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(impute)
library(dplyr)
library(e1071)
library(preprocessCore)
library(data.table)
library(BiocParallel)
library(lumi)

# Functions
preprocess_methylation_matrix <- function(file_path) {
  beta_matrix <- load_beta_value_matrix(file_path)
  imputed_filtered_beta_matrix <- filter_and_impute(beta_matrix)
  normalized_beta_matrix <- normalize(imputed_filtered_beta_matrix)
  final_beta_matrix <- t(apply(normalized_beta_matrix, 1, adjust_outliers))
  return(final_beta_matrix)
}

load_beta_value_matrix <- function(file_path) {
  file_extension <- tools::file_ext(file_path)
  
  # Use fread for efficient data reading
  beta_matrix <- if (file_extension == "csv") {
    fread(file_path, data.table = FALSE, check.names = FALSE)
  } else if (file_extension == "tsv" || file_extension == "txt") {
    fread(file_path, sep = "\t", data.table = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported file type. Please provide a CSV or TSV file.")
  }
  
  colnames(beta_matrix)[1] <- "CpG_Site"
  
  if (anyDuplicated(beta_matrix[, "CpG_Site"])) {
    beta_matrix <- beta_matrix %>%
      group_by(CpG_Site) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')
  }
  
  beta_matrix <- as.data.frame(beta_matrix)
  row.names(beta_matrix) <- beta_matrix$CpG_Site
  beta_matrix <- as.matrix(beta_matrix[, -1])
  
  return(beta_matrix)
}

filter_and_impute <- function(beta_matrix) {
  missing_per_site <- colSums(is.na(beta_matrix)) / nrow(beta_matrix)
  missing_per_sample <- rowSums(is.na(beta_matrix)) / ncol(beta_matrix)
  
  site_missing_threshold <- mean(missing_per_site) + sd(missing_per_site)
  sample_missing_threshold <- mean(missing_per_sample) + sd(missing_per_sample)
  
  filtered_beta_matrix <- beta_matrix[rowSums(is.na(beta_matrix)) / ncol(beta_matrix) < sample_missing_threshold, ]
  filtered_beta_matrix <- filtered_beta_matrix[, colSums(is.na(filtered_beta_matrix)) / nrow(filtered_beta_matrix) < site_missing_threshold]
  
  median_methylation <- apply(filtered_beta_matrix, 1, median, na.rm = TRUE)
  mad_methylation <- apply(filtered_beta_matrix, 1, function(x) mad(x, constant = 1, na.rm = TRUE))
  
  data_skewness <- skewness(log1p(rowMedians(filtered_beta_matrix, na.rm = TRUE)))
  mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))
  
  threshold_methylation <- median_methylation + (mad_multiplier * mad_methylation)
  
  overall_non_zero_proportion <- mean(rowSums(filtered_beta_matrix > 0, na.rm = TRUE)) / ncol(filtered_beta_matrix)
  samples_threshold <- 0.05 + (0.15 * (1 - overall_non_zero_proportion))
  
  keep_sites <- rowSums(filtered_beta_matrix >= threshold_methylation, na.rm = TRUE) >= (samples_threshold * ncol(filtered_beta_matrix))
  filtered_beta_matrix <- filtered_beta_matrix[keep_sites, ]
  
  k <- min(10, ncol(filtered_beta_matrix) - 1)
  imputed_filtered_beta_matrix <- impute.knn(as.matrix(filtered_beta_matrix), k = k)$data
  
  return(imputed_filtered_beta_matrix)
}

normalize <- function(imputed_filtered_beta_matrix){
  normalized_beta_matrix <- as.data.frame(normalize.quantiles(as.matrix(imputed_filtered_beta_matrix)))
  
  rownames(normalized_beta_matrix) <- rownames(imputed_filtered_beta_matrix)
  colnames(normalized_beta_matrix) <- colnames(imputed_filtered_beta_matrix)
  
  return(as.matrix(normalized_beta_matrix))
}

adjust_outliers <- function(normalized_beta_matrix) {
  calculate_outlier_bounds <- function(normalized_beta_matrix, k) {
    Q1 <- quantile(normalized_beta_matrix, 0.25, na.rm = TRUE)
    Q3 <- quantile(normalized_beta_matrix, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(c(lower_bound, upper_bound))
  }
  
  k <- 1.5
  bounds <- calculate_outlier_bounds(normalized_beta_matrix, k)
  outlier_proportion <- mean(normalized_beta_matrix < bounds[1] | normalized_beta_matrix > bounds[2], na.rm = TRUE)
  
  target_proportion <- 0.05
  while (outlier_proportion > target_proportion && k < 10) {
    k <- k + 0.5
    bounds <- calculate_outlier_bounds(normalized_beta_matrix, k)
    outlier_proportion <- mean(normalized_beta_matrix < bounds[1] | normalized_beta_matrix > bounds[2], na.rm = TRUE)
  }
  
  adjusted_matrix <- ifelse(normalized_beta_matrix < bounds[1], bounds[1], ifelse(normalized_beta_matrix > bounds[2], bounds[2], normalized_beta_matrix))
  
  return(adjusted_matrix)
}

# Main processing loop
setwd("")
files_in_directory <- list.files()
pattern <- "_Methylation\\.csv$"
files_to_process <- grep(pattern, files_in_directory, value = TRUE)

for (file_name in files_to_process) {
  cancer_name <- gsub(pattern, "", file_name)
  rda_file_name <- paste0(cancer_name, "_Methylation.rda")
  
  # Check if the .rda file already exists
  if (file.exists(rda_file_name)) {
    cat("Skipping", rda_file_name, "- file already exists.\n")
    next
  }
  
  # Process in chunks to handle large data
  try({
    processed_beta_matrix <- preprocess_methylation_matrix(file_name)
    # Convert beta values to M values
    m_matrix <- beta2m(processed_beta_matrix)
    m_df <- as.data.frame(m_matrix)
    
    saveRDS(m_df, rda_file_name)
    cat("Processed and saved", rda_file_name, "\n")
  }, silent = TRUE)
  
  if (exists("processed_beta_matrix")) {
    rm(processed_beta_matrix)
  }
  gc()
}

#######################################################

#Check Row Numbers the resulted Methylome Data (RDA files)
# Load necessary libraries
library(data.table)

# Function to check the number of rows in an RDA file
check_rows_in_rda <- function(rda_file) {
  # Load the RDA file
  data <- readRDS(rda_file)
  
  # Get the number of rows
  num_rows <- nrow(data)
  
  return(num_rows)
}

# List all .rda files in the directory with the specific pattern
rda_files <- list.files(pattern = "_Methylation\\.rda$")

# Create a data frame to store the results
results <- data.frame(CancerType = character(), File = character(), NumRows = integer(), stringsAsFactors = FALSE)

# Loop through each .rda file and check the number of rows
for (rda_file in rda_files) {
  num_rows <- check_rows_in_rda(rda_file)
  cancer_type <- gsub("_Methylation\\.rda$", "", rda_file)
  results <- rbind(results, data.frame(CancerType = cancer_type, File = rda_file, NumRows = num_rows))
}

# Print the results
print(results)




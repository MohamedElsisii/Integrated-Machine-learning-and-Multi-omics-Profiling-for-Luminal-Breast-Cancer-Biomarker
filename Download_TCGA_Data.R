# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
############################################
# Install the libraries from Bioconductor
BiocManager::install(c("TCGAbiolinks", "maftools", "SummarizedExperiment", "sesameData", "sesame","tidyverse"))
############################################
# Import libraries
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(SummarizedExperiment)
library(sesameData)
library(sesame)
############################################
# Set working directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
############################################
# Set functions to download and label the TCGA data
clean_and_rename_columns <- function(df, common_samples) {
  # Extract the simplified column names (up to -01, inclusive)
  simplified_colnames <- sapply(colnames(df), function(name) {
    sub("^(TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-01).*", "\\1", name)
  })
  
  # Create a logical vector indicating whether each simplified name is in common_samples
  in_common_samples <- simplified_colnames %in% common_samples
  
  # Filter the dataframe columns based on the match to common_samples
  df <- df[, in_common_samples]
  
  # Rename columns to their simplified form, ensuring they exactly match common_samples
  colnames(df) <- simplified_colnames[in_common_samples]
  
  # Remove duplicate columns by keeping only the first occurrence of each unique identifier
  df <- df[, !duplicated(colnames(df))]
  
  return(df)
}
############################################
# Define a vector of cancer types to loop over
cancer_types <- c("BLCA","COAD","KIRP","LGG","LIHC","LUAD","MESO","PAAD","READ","SARC","STAD","UCEC","PRAD", "THCA", "SKCM","BRCA")
############################################
# Loop through each cancer type and download its RNA and Methylation data common samples
for(cancer_type in cancer_types) {
  common_samples <- matchedMetExp(paste0("TCGA-", cancer_type))
  
  Methylation <- GDCquery(project = paste0("TCGA-", cancer_type),
                          data.category = "DNA Methylation",
                          platform = "Illumina Human Methylation 450",
                          data.type = "Methylation Beta Value",
                          barcode = common_samples)
  
  mRNA <- GDCquery(project = paste0("TCGA-", cancer_type),
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   barcode = common_samples)
  
  GDCdownload(Methylation)
  GDCdownload(mRNA)
  
  data_rnaseq <- GDCprepare(mRNA)
  data_methylation <- GDCprepare(Methylation)
  
  geneExp <- SummarizedExperiment::assay(data_rnaseq)
  geneExp = as.data.frame(geneExp)
  geneExp <- clean_and_rename_columns(geneExp, common_samples)
  
  MethExp <- SummarizedExperiment::assay(data_methylation)
  MethExp = as.data.frame(MethExp)
  MethExp <- clean_and_rename_columns(MethExp, common_samples)
  
  write.csv(geneExp, file = paste0("D:/Data/", cancer_type, "_RNA.csv"))
  write.csv(MethExp, file = paste0("D:/Data/", cancer_type, "_Methylation.csv"))
}
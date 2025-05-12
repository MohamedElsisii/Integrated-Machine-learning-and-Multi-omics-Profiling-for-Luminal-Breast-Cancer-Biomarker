Integrated Machine-learning and Multi-omics Profiling for Luminal Breast Cancer Biomarker Identification Pipeline
![Integrated machine-learning and multi-omics profiling for luminal breast cancer biomarker identification_page-0001](https://github.com/user-attachments/assets/a4ad48d6-bec8-4247-8cf9-96ea91a4c134)

## Execution Order & Description

### 1. Download_TCGA_Data Script
- **Purpose**: Automates the download of RNA and methylation data from TCGA for multiple cancer types.
- **Main Functions**: 
  - Queries TCGA
  - Downloads data
  - Retains common samples
  - Renames columns for consistency
- **Output**: `<Cancer_Type>_RNA.csv`, `<Cancer_Type>_Methylation.csv`

### 2. Preprocessing RNA Data Script
- **Purpose**: Processes RNA-seq count matrices, including normalization, gene mapping, and filtering.
- **Main Functions**: 
  - Loads data
  - Maps Ensembl IDs
  - Filters low-expression genes
  - Normalizes data
  - Adjusts outliers
- **Output**: Processed RNA count matrix (`.rda` format)

### 3. Preprocessing Methylation Data Script
- **Purpose**: Processes DNA methylation beta value matrices, including imputation, normalization, and outlier handling.
- **Main Functions**: 
  - Loads beta values
  - Imputes missing data
  - Applies quantile normalization
  - Adjusts outliers
- **Output**: Processed Methylation Matrix (`.rda` format)

### 4. BRCA Samples Subtyping Script
- **Purpose**: Classifies BRCA samples into molecular subtypes (LumA, LumB, Basal, Her2, Normal) using the PAM50 classifier.
- **Main Functions**: 
  - Maps gene names to probe IDs
  - Prepares expression data
  - Performs PAM50 subtyping
- **Output**: Subtype sample lists (`.rda` format)

### 5. Mapping Script
- **Purpose**: Generates a mapping between genes and CpG sites based on genomic coordinates.
- **Main Functions**: 
  - Imports reference data
  - Processes genomic locations
  - Creates gene-to-CpG site mapping
- **Output**: Gene-to-CpG site mapping (`.rda` format)

### 6. Model Script
- **Purpose**: Trains predictive models using Lasso regression to predict gene expression from methylation data.
- **Main Functions**: 
  - Loads data
  - Applies Lasso regression
  - Stores model outputs
- **Output**: Regression model file (`.rda` format)

### 7. Prediction Script
- **Purpose**: Uses trained Lasso models to predict gene expression levels from new methylation data.
- **Main Functions**: 
  - Loads trained models
  - Applies them to new methylation data
  - Generates predicted expression values
- **Output**: Predicted expression data (`.rda` format)

### 8. DEGs Loop (t-test & Wilcoxon) Script
- **Purpose**: Identifies differentially expressed genes (DEGs) between normal and cancer samples using statistical tests.
- **Main Functions**: 
  - Performs normality testing
  - t-tests or Wilcoxon tests
  - Multiple testing correction
  - Generates volcano plots
- **Output**: 
  - DEG results (`.csv`)
  - Significant genes (`.csv`)
  - Volcano plots (`.png`)

## Usage Instructions
- Ensure all required R packages are installed before running the scripts.
- Follow the execution order to maintain data integrity and consistency.
- Modify cancer types, parameters, and paths as needed.
- Use the generated `.rda` and `.csv` files for downstream analysis and visualization.

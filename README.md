# 🧬 Integrated Machine-learning and Multi-omics Profiling for Luminal Breast Cancer Biomarker Identification Pipeline
![Integrated machine-learning and multi-omics profiling for luminal breast cancer biomarker identification_page-0001](https://github.com/user-attachments/assets/a4ad48d6-bec8-4247-8cf9-96ea91a4c134)

## 🚀 Execution Order & Description  

### 1️⃣ **`Download_TCGA_Data` Script**  
🔹 **Purpose**: Automates the download of RNA and methylation data from TCGA for multiple cancer types.  
🔹 **Main Functions**:  
  - 📥 Queries TCGA  
  - ⬇️ Downloads data  
  - 🔄 Retains common samples  
  - 🏷️ Renames columns for consistency  
🔹 **Output**: `📂 <Cancer_Type>_RNA.csv`, `📂 <Cancer_Type>_Methylation.csv`  

### 2️⃣ **`Preprocessing RNA Data` Script**  
🔹 **Purpose**: Processes RNA-seq count matrices (normalization, gene mapping, filtering).  
🔹 **Main Functions**:  
  - 🧹 Loads & cleans data  
  - 🧬 Maps Ensembl IDs  
  - 🔍 Filters low-expression genes  
  - ⚖️ Normalizes data  
  - 📊 Adjusts outliers  
🔹 **Output**: `📂 Processed RNA count matrix (.rda)`  

### 3️⃣ **`Preprocessing Methylation Data` Script**  
🔹 **Purpose**: Processes DNA methylation beta values (imputation, normalization, outlier handling).  
🔹 **Main Functions**:  
  - 🧪 Loads beta values  
  - 🔄 Imputes missing data  
  - 📏 Quantile normalization  
  - 🎯 Adjusts outliers  
🔹 **Output**: `📂 Processed Methylation Matrix (.rda)`  

### 4️⃣ **`BRCA Samples Subtyping` Script**  
🔹 **Purpose**: Classifies BRCA samples into molecular subtypes (**LumA, LumB, Basal, Her2, Normal**) using PAM50.  
🔹 **Main Functions**:  
  - 🧬 Maps gene names → probe IDs  
  - 📊 Prepares expression data  
  - 🔍 Performs PAM50 subtyping  
🔹 **Output**: `📂 Subtype sample lists (.rda)`  

### 5️⃣ **`Mapping` Script**  
🔹 **Purpose**: Creates gene-to-CpG site mappings based on genomic coordinates.  
🔹 **Main Functions**:  
  - 🗺️ Imports reference data  
  - 🧬 Processes genomic locations  
  - 🔗 Generates gene-CpG mappings  
🔹 **Output**: `📂 Gene-to-CpG mapping (.rda)`  

### 6️⃣ **`Model` Script**  
🔹 **Purpose**: Trains **Lasso regression** models to predict gene expression from methylation data.  
🔹 **Main Functions**:  
  - 🤖 Loads data  
  - 📉 Applies Lasso regression  
  - 💾 Stores model outputs  
🔹 **Output**: `📂 Trained model (.rda)`  

### 7️⃣ **`Prediction` Script**  
🔹 **Purpose**: Predicts gene expression from new methylation data using trained models.  
🔹 **Main Functions**:  
  - 🔮 Loads trained models  
  - 🎯 Applies to new data  
  - 📈 Generates predicted expression  
🔹 **Output**: `📂 Predicted expression (.rda)`  

### 8️⃣ **`DEGs Loop (t-test & Wilcoxon)` Script**  
🔹 **Purpose**: Identifies **Differentially Expressed Genes (DEGs)** between normal/cancer samples.  
🔹 **Main Functions**:  
  - 📊 Normality testing  
  - ⚖️ t-test / Wilcoxon test  
  - 🔢 Multiple testing correction  
  - 🌋 Volcano plots  
🔹 **Output**:  
  - `📂 DEG results (.csv)`  
  - `📂 Significant genes (.csv)`  
  - `📂 Volcano plots (.png)`  

---

## 🛠️ **Usage Instructions**  
✅ **Before running**:  
   - Install all required **R packages** (`📦 tidyverse`, `📦 glmnet`, etc.)  
   - Check file paths and parameters.  

🔁 **Execution flow**:  
   **1 → 2 → 3 → 4 → 5 → 6 → 7 → 8** (Follow order for consistency!)  

⚙️ **Customization**:  
   - Modify `cancer_types`, `parameters`, and `output paths` as needed.  

📤 **Outputs**:  
   - Use `.rda` / `.csv` files for **downstream analysis** & visualization.  

library(glmnet)
###########################################################
#1]Set directory to the working location
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
###########################################################
#2]Read expression, methylation, gene_to_cpg_site_mapping files
# Expression and Methylation data must be the same samples#
gene_to_cpg_site_mapping =readRDS("gene_to_cpg_site_mapping.rda")
expression <- read.csv("expression.csv", header = TRUE, row.names = 1)
methylation <- read.csv("methylation.csv", header = TRUE, row.names = 1)
###########################################################
#3]Model training function
train_predict_lasso <- function(gene, cg_sites, expression_data, methylation_data) {
  # Subset methylation data for relevant CpG sites and transpose for glmnet
  X_train <- t(methylation_data[cg_sites, , drop = FALSE])
  
  # Get corresponding expression data for the gene
  Y_train <- expression_data[gene, ]
  
  # Ensure Y_train is a numeric vector
  Y_train <- as.numeric(Y_train)
  
  # Train Lasso model with cross-validation
  model <- cv.glmnet(x = as.matrix(X_train), y = Y_train, alpha = 1, nfolds = 10)
  
  # Predict gene expression
  predictions <- predict(model, newx = as.matrix(X_train), s = "lambda.min")
  
  # Calculate R-squared
  R_squared <- cor(Y_train, predictions)^2
  
  # Return a named list including the gene, model, predictions, and R-squared
  return(list(Gene = gene,
              Fit_Parameters = model, 
              Predicted_Expression = predictions, 
              R_squared = R_squared))
}

# Wrapper function to apply the Lasso model to each gene and return a data frame
apply_lasso_to_genes <- function(expression_data, methylation_data, gene_to_cpg_site_mapping) {
  results_list <- lapply(1:nrow(gene_to_cpg_site_mapping), function(i) {
    gene <- gene_to_cpg_site_mapping$gene[i]
    cg_sites <- gene_to_cpg_site_mapping$cg_sites[[i]]
    result <- train_predict_lasso(gene, cg_sites, expression_data, methylation_data)
    result
  })
  
  # Convert the list of results to a data frame
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(Gene = x$Gene,
               Fit_Parameters = I(list(x$Fit_Parameters)), 
               Predicted_Expression = I(list(x$Predicted_Expression)), 
               R_squared = x$R_squared)
  }))
  colnames(results_df) <- c("Gene", "Fit_Parameters", "Predicted_Expression", "R_squared")
  
  colnames(results_df) <- c("Gene", "Fit_Parameters", "Predicted_Expression", "R_squared")
  return(results_df)
}
#############################################################
#4]Call the function
regression.output <- apply_lasso_to_genes(expression_data = expression, methylation_data = methylation, gene_to_cpg_site_mapping = gene_to_cpg_site_mapping)
#############################################################
#5]Save the model
saveRDS(regression.output, "regression.output.rda")


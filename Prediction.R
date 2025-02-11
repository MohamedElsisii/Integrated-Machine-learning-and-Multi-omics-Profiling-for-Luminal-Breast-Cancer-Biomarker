library(glmnet)
library(dplyr)
library(parallel)
#################################################################
#1]Set directory to the working location
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
#################################################################
#2]Read gene_to_cpg_site_mapping, methylation data for the prediction of the expression data
gene_to_cpg_site_mapping <- readRDS("gene_to_cpg_site_mapping.rda")
BRCA_Methylation <- read.csv("BRCA_Methylation.csv", header = TRUE, row.names = 1)
#################################################################
#3]Read regression model and convert the gene column into row names
regression.output=readRDS("regression.output.rda")
row.names(regression.output) <- regression.output[, 1]; regression.output <- regression.output[, -1]
#################################################################
#4] Prediction function
predict_expression <- function(methylation, gene_to_cpg_site_mapping, regression.output) {
  predicted_expression_list <- list()
  
  for (gene in rownames(regression.output)) {
    cg_sites <- gene_to_cpg_site_mapping$cg_sites[gene_to_cpg_site_mapping$gene == gene][[1]]
    
    if (length(cg_sites) == 0) {
      cat("No CpG sites found for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    fit <- regression.output$Fit_Parameters[[which(rownames(regression.output) == gene)]]

    if (is.null(fit)) {
      cat("Invalid or NULL model for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    if (!all(cg_sites %in% rownames(methylation))) {
      cat("Not all CpG sites found in methylation data for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    X <- t(methylation[cg_sites, , drop = FALSE])
    gene_prediction <- predict(fit, newx = X, s = "lambda.min")
    predicted_expression_list[[gene]] <- as.vector(gene_prediction)
  }
  
  predicted_expression <- do.call(cbind, predicted_expression_list)
  rownames(predicted_expression) <- colnames(methylation)
  colnames(predicted_expression) <- names(predicted_expression_list)
  
  return(predicted_expression)
}
#####################################################
#5] Call function
predicted_expression <- predict_expression(methylation=BRCA_Methylation, gene_to_cpg_site_mapping, regression.output=regression.output)
#####################################################
#6] Transpose the matrix and then change it into a dataframe 
predicted_expression=t(predicted_expression)
predicted_expression=as.data.frame(predicted_expression)
#####################################################
#7] Save the predicted dataframe
saveRDS(predicted_expression,"Predicted.rda")


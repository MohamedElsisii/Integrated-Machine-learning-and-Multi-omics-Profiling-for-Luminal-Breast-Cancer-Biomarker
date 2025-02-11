######################## Libraries #########################
library(limma)
library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(stringr)
library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(caret)
library(biomaRt)  
library(dplyr)
library(tibble)
#######################################################
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
BRCA <- read.csv("BRCA_RNA.csv")
colnames(BRCA)[1]<- "Genes"
####################################################################
data(pam50)
data("pam50.robust")
#################################### Probe ID ########################################
ensembl_dataset <- "hsapiens_gene_ensembl"  
ensembl_host <- "https://www.ensembl.org"

# Create a biomaRt object
mart <- useMart(biomart = "ensembl", host = ensembl_host, dataset = ensembl_dataset)

# Extract gene names from row names of exprset
gene_names <- BRCA[, 1]

# Use biomaRt to convert gene names to probe IDs
probe_ids <- getBM(attributes = c("external_gene_name", "affy_hg_u133_plus_2"),
                   filters = "external_gene_name",
                   values = gene_names,
                   mart = mart)

#Remove empty rows and choose only one probe and changing the column names
probe_ids <- probe_ids[probe_ids$affy_hg_u133_plus_2 != "", ]
probe_ids <- probe_ids %>% distinct(external_gene_name, affy_hg_u133_plus_2, .keep_all = TRUE)
colnames(probe_ids) <- c("Gene_Name", "Probe_id")

# Match the values of the Genes column with the Gene_Name column in probe_ids
match_vec <- match(BRCA$Genes, probe_ids$Gene_Name)

# Extract matched probe_ids using match_vec
matched_probe_ids <- probe_ids$Probe_id[match_vec]

# Replace Genes column in BRCA with matched_probe_ids
BRCA$Genes <- matched_probe_ids

#Remove NA rows
BRCA <- BRCA[complete.cases(BRCA), ]

#Remove duplicates
BRCA <- BRCA %>% distinct(Genes, .keep_all = TRUE)

#Convert first column into row names
row.names(BRCA) <- BRCA$Genes
BRCA <- BRCA[, -1, drop = FALSE]
######################################Entrez id########################################################

probe_ids <- rownames(BRCA)

# Use biomaRt to convert probe IDs to Entrez Gene IDs
gene_ids <- getBM(attributes = c("affy_hg_u133_plus_2", "entrezgene_id"),
                  filters = "affy_hg_u133_plus_2",
                  values = probe_ids,
                  mart = mart)

annot_df <- merge(data.frame(probe = probe_ids), gene_ids, by.x = "probe", by.y = "affy_hg_u133_plus_2", all.x = TRUE)
# Create the desired data frame
annot_df <- data.frame(probe = annot_df$probe, EntrezGene.ID = annot_df$entrezgene_id, row.names = NULL)

##############################################Remove Null values###############################################
complete_rows_annot <- complete.cases(annot_df)
annot_df <- annot_df[complete_rows_annot, ]
common_values <- intersect(annot_df$probe, rownames(BRCA))
New <- BRCA[common_values, ]
############################################Subtyping###################################################
TRANS <- t(New)
PAM50Preds <- molecular.subtyping(sbt.model = "pam50", data = TRANS, annot=annot_df, do.mapping=TRUE)
table(PAM50Preds$subtype)
LumA<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumA")]
LumB<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumB")]
Basal<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Basal")]
Her2<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Her2")]
Normal<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "Normal")]

saveRDS(LumA, "LumA_Samples.rda")
saveRDS(LumB, "LumB_Samples.rda")

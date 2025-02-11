library(rtracklayer)
library(dplyr)
library(parallel)
library(pracma)
####################################################
#Set directory to the working location
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
####################################################
HG19.REF.GENE.PATH <- "hg19.refGene.gtf"
GENE.COORDINATES.PATH <- "gene.coordinates.rda"
HUMAN.METHYLATION.450.MANIFEST.PATH <-"HumanMethylation450_15017482_v1-2.csv"
GENE.TO.CG.SITE.MAPPING.PATH.FORMAT <- "gene_to_cpg_site_mapping.rda"
####################################################
#1]Get genes coordinates 
####################################################
get.full.ref.gene.df <- function() 
{
  hg19.ref.gene <- rtracklayer::import(HG19.REF.GENE.PATH)
  hg19.ref.gene.df <- as.data.frame(hg19.ref.gene)
  return(hg19.ref.gene.df)
}
get.gene.coordinates.df <- function(force.update = FALSE) 
{
  if (file.exists(GENE.COORDINATES.PATH) & !force.update) {
    return(readRDS(file = GENE.COORDINATES.PATH))
  }
  
  hg19.ref.gene.df <- get.full.ref.gene.df()
  
  first.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chr1")
  last.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chrY")
  
  gene.coordinates.df <- hg19.ref.gene.df %>% 
    filter(type == "transcript" &
             as.integer(seqnames) >= first.chr.index & as.integer(seqnames) <= last.chr.index) %>%
    select(gene_id, chr = seqnames, strand, start, end) %>% 
    distinct()
  
  saveRDS(gene.coordinates.df, file = GENE.COORDINATES.PATH)
  return(gene.coordinates.df)
}
####################################################
#Call function
gene_coords_df <- get.gene.coordinates.df()
####################################################
#2]Create gene to cpg site mapping file
####################################################
is.valid.value <- function(value)
{
  return(!is.na(value) & value != "")
}
utils.log <- function(log.str) {
  
  print(paste(format(Sys.time(), "[%X]"), log.str))
}
create.gene.to.cg.site.mapping <- function(
    gene.coordinates.df, upstream.margin.in.bases, downstream.margin.in.bases, force.update = FALSE)
{
  mapping.file.path <- 
    sprintf(GENE.TO.CG.SITE.MAPPING.PATH.FORMAT, upstream.margin.in.bases, downstream.margin.in.bases)
  
  if (file.exists(mapping.file.path) & !force.update) {
    return(readRDS(file = mapping.file.path))
  }
  
  methylation.manifest <- read.csv(file = HUMAN.METHYLATION.450.MANIFEST.PATH)
  
  methylation.manifest.filtered <- methylation.manifest %>%
    filter(is.valid.value(CHR) & is.valid.value(MAPINFO) & Genome_Build == 37)
  
  gene.to.cg.site.mapping <- data.frame(gene = unique(gene.coordinates.df$gene_id))
  gene.to.cg.site.mapping$cg_sites <- lapply(1:nrow(gene.to.cg.site.mapping), function(x) c())
  
  cg.sites.per.gene <- mclapply(1:nrow(gene.to.cg.site.mapping), function(gene_idx) {
    gene <- gene.to.cg.site.mapping$gene[gene_idx]
    gene.coordinates.filtered <- gene.coordinates.df %>% filter(gene_id == gene)
    
    #convert chromosome to string comparable with the one in methylation.manifest.filtered$CHR
    gene.coordinates.filtered$chr <- sapply(gene.coordinates.filtered$chr, function(chr.str) 
      substr(as.character(chr.str), start = nchar("chr") + 1, stop = nchar(as.character(chr.str))))
    
    gene.cg.sites <- c()
    
    for (i in 1:nrow(gene.coordinates.filtered)) {
      chr <- gene.coordinates.filtered$chr[i]
      
      if (gene.coordinates.filtered$strand[i] == "+") {
        
        start <- gene.coordinates.filtered$start[i] - upstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + downstream.margin.in.bases
        
      } else {
        
        start <- gene.coordinates.filtered$start[i] - downstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + upstream.margin.in.bases
      }
      
      cg.sites.in.range <- (methylation.manifest.filtered %>% 
                              filter(tolower(CHR) == tolower(chr) & MAPINFO >= start & MAPINFO <= end) %>%
                              select(Name))$Name
      
      gene.cg.sites <- c(gene.cg.sites, cg.sites.in.range)
    }
    
    if (gene_idx %% 1000 == 1) {
      utils.log(sprintf("Done mapping gene number %d", gene_idx - 1))
    }
    
    return(gene.cg.sites)
  })
  
  gene.to.cg.site.mapping$cg_sites <- lapply(cg.sites.per.gene, unique)
  
  saveRDS(gene.to.cg.site.mapping, file = mapping.file.path)
  return(gene.to.cg.site.mapping)
}
####################################################
# Define upstream and downstream margins "10Kb"
upstream_margin <- 10000
downstream_margin <- 10000
####################################################
# Call the function to create gene-to-CpG site mapping
gene_to_cpg_site_mapping <- create.gene.to.cg.site.mapping(gene.coordinates.df =  gene_coords_df,upstream.margin.in.bases =  upstream_margin,downstream.margin.in.bases =  downstream_margin)
####################################################
#3]Filter the gene_to_cpg_site_mapping file from NA, genes with less than 2 CpG sites and with expression and methylation features
####################################################
# Extract features "gene names" from the 'expression' training data frame
expression <- read.csv("expression.csv", header = TRUE, row.names = 1)
expression_features <- rownames(expression)
# Extract features "cg_sites" from the 'methylation' training data frame
methylation <- read.csv("methylation.csv", header = TRUE, row.names = 1)
methylation_features <- rownames(methylation)
#Filtration
gene_to_cpg_site_mapping <- gene_to_cpg_site_mapping %>% 
  filter(gene %in% expression_features)

gene_to_cpg_site_mapping$cg_sites <- lapply(gene_to_cpg_site_mapping$cg_sites, function(cpg_sites_for_gene) 
  intersect(cpg_sites_for_gene, methylation_features))

cpg_sites_lengths <- sapply(gene_to_cpg_site_mapping$cg_sites, length)

gene_to_cpg_site_mapping <- gene_to_cpg_site_mapping[which(cpg_sites_lengths > 1), ]
####################################################
#Save final file
saveRDS(gene_to_cpg_site_mapping, "gene_to_cpg_site_mapping.rda")

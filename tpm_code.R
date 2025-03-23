# Function Information
# ------------------
# Function Name: CountsToTPM.R
# Description: Given a count data matrix, returns TPM normalized data
## Input: 
### counts matrix: a m x n matrix with m genes and n samples
### gene_id_col: Name of the column where the geneIds are
### datset: User defined organism dataset (e.g. hsapiens_gene_ensembl)
### attribute: biomaR_attribute. Should be of the same type of the counts matrix GeneID
## 
# Author: José Basílio
# Date: 2023-08-17

CountsToTPM <- function(attribute, counts, gene_id_col, dataset) {
  if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
  library(pacman)
  pacman::p_load(tidyverse, biomaRt)
  
  # Function to retrieve gene length data
  retrieve_gene_length <- function(attribute, dataset) {
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
    
    attributes <- c(attribute, "start_position", "end_position")
    
    gene_data <- getBM(attributes = attributes, mart = ensembl)
    
    gene_dataframe <- as.data.frame(gene_data) %>% 
      dplyr::mutate(gene_length = end_position - start_position) %>% 
      tidyr::drop_na(!!sym(attribute)) %>% 
      dplyr::distinct(!!sym(attribute), .keep_all = TRUE) %>% 
      dplyr::select(!!sym(attribute), gene_length) %>%
      dplyr::arrange(!!sym(attribute))
    
    return(gene_dataframe)
  }
  
  # Function to calculate TPM matrix
  func_tpm <- function(counts, gene_lengths_df, gene_id_col) {
    common_gene_ids <- intersect(counts[, gene_id_col], gene_lengths_df[[attribute]])
    
    counts_common_genes <- counts[counts[, gene_id_col] %in% common_gene_ids, ]
    lengths_common_genes <- gene_lengths_df[gene_lengths_df[[attribute]] %in% common_gene_ids, ]
    
    gene_lengths <- lengths_common_genes$gene_length
    x <- counts_common_genes[, -which(names(counts_common_genes) == gene_id_col)] / gene_lengths
    tpm_matrix <- t(t(x) * 1e6 / colSums(x, na.rm = TRUE))
    
    rownames(tpm_matrix) <- common_gene_ids
    colnames(tpm_matrix) <- colnames(counts_common_genes)[-which(names(counts_common_genes) == gene_id_col)]
    
    return(tpm_matrix)
  }
  
  # Calling the retrieve_gene_length function
  gene_lengths_df <- retrieve_gene_length(attribute, dataset)
  
  # Calling the func_tpm function with retrieved gene length data
  tpm_matrix <- as.data.frame(func_tpm(counts, gene_lengths_df, gene_id_col))
  
  return(tpm_matrix)
}

countData<- read.xlsx("C://Users/Anna/Documents/CVID_project/countData.xlsx")

# Example Usage:

# gene_counts dataframe
## counts_df <- data.frame(
##  gene_symbol = c("TSPAN6", "TNMD", "DPM1", "SCYL3"),
##  sample1 = c(100, 50, 200, 300),
##  sample2 = c(150, 30, 250, 280),
##  sample3 = c(180, 40, 220, 320))

# User-defined dataset, attribute and column for gene IDs
## dataset <- "hsapiens_gene_ensembl"
## attribute <- "external_gene_name"
## gene_id_col <- "gene_symbol"

# Call the combined_gene_analysis function
tpm_result <- CountsToTPM(attribute = "external_gene_name", counts = countData, gene_id_col = "row.names", dataset = "hsapiens_gene_ensembl")

write.xlsx(tpm_result,"C:/Users/Anna/Documents/CVID_project/tpm_countTable.xlsx", rowNames=TRUE)
# Now you have the TPM matrix in the 'tpm_result' variable
## print(tpm_result)
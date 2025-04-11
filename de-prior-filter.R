## Testing out different quartiles for filtering the DE prior list

library(dplyr)

find_common_genes <- function(de_prior, gene_dimnames, dataset_name, genes_info_table) {
  # Identify common genes
  common_genes <- de_prior[de_prior$Gene_Name %in% gene_dimnames, ]
  
  # Identify genes only in the assay
  genes_only_in_assay <- setdiff(gene_dimnames, de_prior$Gene_Name)
  
  # Update the genes_info_table
  genes_info_table <- rbind(genes_info_table, 
                            data.frame(dataset = dataset_name, 
                                       common_genes_count = nrow(common_genes),
                                       only_assay_count = length(genes_only_in_assay)))
  
  return(genes_info_table)
}

genes_info_table <- data.frame()

# Filter out genes with DE prior rank greater than or equal to 0.95

de_prior_95 <- de_prior %>%
  filter(DE_Prior_Rank >= 0.95)

results_filtered_95 <- find_common_genes(de_prior_95, gene_dimnames_blca, "BLCA", genes_info_table)
results_filtered_95 <- find_common_genes(de_prior_95, gene_dimnames_brca, "BRCA", results_filtered_95)
results_filtered_95 <- find_common_genes(de_prior_95, gene_dimnames_thca, "THCA", results_filtered_95)
results_filtered_95 <- find_common_genes(de_prior_95, gene_dimnames_kipan, "KIPAN", results_filtered_95)
results_filtered_95 <- find_common_genes(de_prior_95, gene_dimnames_GSE71669$hgnc_symbol, "GSE Bladder Cancer", results_filtered_95)

# Filter out genes with DE prior rank greater than or equal to 0.9

de_prior_90 <- de_prior %>%
  filter(DE_Prior_Rank >= 0.90)

results_filtered_90 <- find_common_genes(de_prior_90, gene_dimnames_blca, "BLCA", genes_info_table)
results_filtered_90 <- find_common_genes(de_prior_90, gene_dimnames_brca, "BRCA", results_filtered_90)
results_filtered_90 <- find_common_genes(de_prior_90, gene_dimnames_thca, "THCA", results_filtered_90)
results_filtered_90 <- find_common_genes(de_prior_90, gene_dimnames_kipan, "KIPAN", results_filtered_90)
results_filtered_90 <- find_common_genes(de_prior_90, gene_dimnames_GSE71669$hgnc_symbol, "GSE Bladder Cancer", results_filtered_90)

# Filter out genes with DE prior rank greater than or equal to 0.8

de_prior_80 <- de_prior %>%
  filter(DE_Prior_Rank >= 0.80)

results_filtered_80 <- find_common_genes(de_prior_80, gene_dimnames_blca, "BLCA", genes_info_table)
results_filtered_80 <- find_common_genes(de_prior_80, gene_dimnames_brca, "BRCA", results_filtered_80)
results_filtered_80 <- find_common_genes(de_prior_80, gene_dimnames_thca, "THCA", results_filtered_80)
results_filtered_80 <- find_common_genes(de_prior_80, gene_dimnames_kipan, "KIPAN", results_filtered_80)
results_filtered_80 <- find_common_genes(de_prior_80, gene_dimnames_GSE71669$hgnc_symbol, "GSE Bladder Cancer", results_filtered_80)

# Filter out genes with DE prior rank greater than or equal to 0.75

de_prior_75 <- de_prior %>%
  filter(DE_Prior_Rank >= 0.75)

results_filtered_75 <- find_common_genes(de_prior_75, gene_dimnames_blca, "BLCA", genes_info_table)
results_filtered_75 <- find_common_genes(de_prior_75, gene_dimnames_brca, "BRCA", results_filtered_75)
results_filtered_75 <- find_common_genes(de_prior_75, gene_dimnames_thca, "THCA", results_filtered_75)
results_filtered_75 <- find_common_genes(de_prior_75, gene_dimnames_kipan, "KIPAN", results_filtered_75)
results_filtered_75 <- find_common_genes(de_prior_75, gene_dimnames_GSE71669$hgnc_symbol, "GSE Bladder Cancer", results_filtered_75)

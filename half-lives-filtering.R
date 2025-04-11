half_lives <- readr::read_delim("/scratch/st-singha53-1/nair1602/multiomics_project/method_2_mapping_output.txt", delim = "\t")

df_hl <-
  half_lives %>%
  mutate(Protein_Half_Life = log(as.numeric(`Protein half-life average [h]`)),
         mRNA_Half_Life = log(as.numeric(`mRNA half-life average [h]`))) %>%
  select(Human.gene.name, `Protein Names`, `Protein Descriptions`,
         `Protein half-life average [h]`, Protein_Half_Life, 
         `mRNA half-life average [h]`, mRNA_Half_Life)
  

df_hl %>%
  filter(!is.na(mRNA_Half_Life))

hist(df_hl$mRNA_Half_Life) 

# Define the quantiles for 5%, 10%, 20%, and 25%
quantiles <- c(0.025, 0.05, 0.1, 0.125, 0.2, 0.25)

# Define the function to generate datasets with both the top and bottom quantile shortest-living and longest-living proteins
get_extreme_living_proteins <- function(df_hl, quantiles) {
  
  # Initialize an empty list to store the resulting datasets
  result_list <- list()
  
  # Loop through the requested quantiles
  for (quant in quantiles) {
    
    # Calculate the quantile value corresponding to the top x% (e.g., 0.05 = top 5%)
    quantile_value_top <- quantile(df_hl$Protein_Half_Life, 1 - quant)
    quantile_value_bottom <- quantile(df_hl$Protein_Half_Life, quant)
    
    # Filter the proteins with half-lives greater than or equal to the quantile value (for top proteins)
    top_proteins <- df_hl[df_hl$Protein_Half_Life >= quantile_value_top, ]
    
    # Filter the proteins with half-lives less than or equal to the quantile value (for bottom proteins)
    bottom_proteins <- df_hl[df_hl$Protein_Half_Life <= quantile_value_bottom, ]
    
    # Add the result to the list with a named entry for top and bottom proteins
    result_list[[paste0(quant * 100, "%_Longest")]] <- top_proteins
    result_list[[paste0(quant * 100, "%_Shortest")]] <- bottom_proteins
  }
  
  return(result_list)
}

# Generate datasets for both top and bottom quantiles 5%, 10%, 20%, and 25%
extreme_living_proteins_datasets <- get_extreme_living_proteins(df_hl, quantiles)


# Print the resulting datasets for top and bottom quantiles
extreme_living_proteins_datasets


# Define the function to generate datasets with both the top and bottom quantile shortest-living and longest-living proteins
get_extreme_living_mrna <- function(df_hl, quantiles) {
  
  df_hl_filtered <- df_hl %>%
    filter(!is.na(`mRNA half-life average [h]`))
  
  # Initialize an empty list to store the resulting datasets
  result_list <- list()
  
  # Loop through the requested quantiles
  for (quant in quantiles) {
    
    # Calculate the quantile value corresponding to the top x% (e.g., 0.05 = top 5%)
    quantile_value_top <- quantile(df_hl_filtered$mRNA_Half_Life, 1 - quant)
    quantile_value_bottom <- quantile(df_hl_filtered$mRNA_Half_Life, quant)
    
    # Filter the mrna with half-lives greater than or equal to the quantile value (for top mrna)
    top_mrna <- df_hl_filtered[df_hl_filtered$mRNA_Half_Life >= quantile_value_top, ]
    
    # Filter the mrna with half-lives less than or equal to the quantile value (for bottom mrna)
    bottom_mrna <- df_hl_filtered[df_hl_filtered$mRNA_Half_Life <= quantile_value_bottom, ]
    
    # Add the result to the list with a named entry for top and bottom mrna
    result_list[[paste0(quant * 100, "%_Longest")]] <- top_mrna
    result_list[[paste0(quant * 100, "%_Shortest")]] <- bottom_mrna
  }
  
  return(result_list)
}

# Generate datasets for both top and bottom quantiles 5%, 10%, 20%, and 25%
extreme_living_mrna_datasets <- get_extreme_living_mrna(df_hl, quantiles)

# Print the resulting datasets for top and bottom quantiles
extreme_living_mrna_datasets


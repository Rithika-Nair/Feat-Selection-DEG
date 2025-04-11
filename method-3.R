# Load necessary libraries
library(glmnet)
library(caret)
library(dplyr)
library(pROC)
library(tibble)
library(ggplot2)
library(ggrain)

# Creating target df
samples <- rownames(colData(mae_tcga_blca))
response <- factor(colData(mae_tcga_blca)$response, levels = levels(colData(mae_tcga_blca)$response)) 
target <- data.frame(samples = samples, response = response)

# Prepare data 
x <- as.matrix(blca_assay_df)  # Genes as features
y <- target  # Group labels (binary classification)

# Split the data into training and testing sets
set.seed(123)  
trainIndex <- sample(1:nrow(x), size = 0.8 * nrow(x))
x_train <- x[trainIndex, ]
x_test <- x[-trainIndex, ]

y_train <- y %>% filter(samples %in% rownames(x_train))
y_train$response <- factor(y_train$response, levels = c(1, 0))
y_test <- y %>% filter(samples %in% rownames(x_test))
y_test$response <- factor(y_test$response, levels = c(1, 0))

# Selecting only genes found in mrna half lives list
df_hl_filtered <- df_hl %>%
  filter(!is.na(`mRNA half-life average [h]`))

selected_genes <- intersect(colnames(x), df_hl_filtered$Human.gene.name)
x_train_100_percent <- x_train[, selected_genes]
x_test_100_percent <- x_test[, selected_genes]

# List of possible alpha values
alpha_values <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

# mRNA half lives
quantile_mrnas <- extreme_living_mrna_datasets

######################
# Initialize an empty data frame
results_df <- data.frame()

for (alpha_val in alpha_values) {
  cat("\nRunning Elastic Net for alpha =", alpha_val, "\n") 
  
  ###CONTROL CONDITION: Full Gene Set (Top 100% of genes from mrnas list)
  # Train Elastic Net Model
  cv_model_full <- cv.glmnet(x_train_100_percent, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_full <- which.max(cv_model_full$cvm)
  best_lambda_full <- cv_model_full$lambda[best_lambda_index_full]
  best_auc_mean_full <- cv_model_full$cvm[best_lambda_index_full]
  best_auc_lo_full <- cv_model_full$cvlo[best_lambda_index_full]
  best_auc_hi_full <- cv_model_full$cvup[best_lambda_index_full]
  
  cat("Best Lambda (Full):", best_lambda_full, "\n")
  cat("AUC Mean (Full):", best_auc_mean_full, "\n")
  cat("AUC CI (Full):", best_auc_lo_full, "-", best_auc_hi_full, "\n")
  
  ### CONDITION A: Feature Shortest (Shortest Living n% Genes)
  shortest_genes <- intersect(colnames(x), quantile_mrnas[["25%_Shortest"]]$Human.gene.name)
  x_train_shortest <- x_train[, shortest_genes]
  x_test_shortest <- x_test[, shortest_genes]
  
  # Train Elastic Net Model
  cv_model_shortest <- cv.glmnet(x_train_shortest, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_shortest <- which.max(cv_model_shortest$cvm)
  best_lambda_shortest <- cv_model_shortest$lambda[best_lambda_index_shortest]
  best_auc_mean_shortest <- cv_model_shortest$cvm[best_lambda_index_shortest]
  best_auc_lo_shortest <- cv_model_shortest$cvlo[best_lambda_index_shortest]
  best_auc_hi_shortest <- cv_model_shortest$cvup[best_lambda_index_shortest]
  
  cat("Best Lambda (Shortest):", best_lambda_shortest, "\n")
  cat("AUC Mean (Shortest):", best_auc_mean_shortest, "\n")
  cat("AUC CI (Shortest):", best_auc_lo_shortest, "-", best_auc_hi_shortest, "\n")
  
  ### **CONDITION B: Feature Longest (Longest Living n% Genes)**
  longest_genes <- intersect(colnames(x), quantile_mrnas[["25%_Longest"]]$Human.gene.name)
  x_train_longest <- x_train[, longest_genes]
  x_test_longest <- x_test[, longest_genes]
  
  # Train Elastic Net Model
  cv_model_longest <- cv.glmnet(x_train_longest, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_longest <- which.max(cv_model_longest$cvm)
  best_lambda_longest <- cv_model_longest$lambda[best_lambda_index_longest]
  best_auc_mean_longest <- cv_model_longest$cvm[best_lambda_index_longest]
  best_auc_lo_longest <- cv_model_longest$cvlo[best_lambda_index_longest]
  best_auc_hi_longest <- cv_model_longest$cvup[best_lambda_index_longest]
  
  cat("Best Lambda (Longest):", best_lambda_longest, "\n")
  cat("AUC Mean (Longest):", best_auc_mean_longest, "\n")
  cat("AUC CI (Longest):", best_auc_lo_longest, "-", best_auc_hi_longest, "\n")
  
  ### **CONDITION C: Feature Longest and Shortest (Select Longest n/2% and Shortest n/2% Genes)**
  no_middle_genes <- intersect(colnames(x), c(quantile_mrnas[["12.5%_Longest"]]$Human.gene.name, 
                                              quantile_mrnas[["12.5%_Shortest"]]$Human.gene.name))
  x_train_no_middle<- x_train[, no_middle_genes]
  x_test_no_middle <- x_test[, no_middle_genes]
  
  # Train Elastic Net Model
  cv_model_no_middle <- cv.glmnet(x_train_no_middle, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_no_middle <- which.max(cv_model_no_middle$cvm)
  best_lambda_no_middle <- cv_model_no_middle$lambda[best_lambda_index_no_middle]
  best_auc_mean_no_middle <- cv_model_no_middle$cvm[best_lambda_index_no_middle]
  best_auc_lo_no_middle <- cv_model_no_middle$cvlo[best_lambda_index_no_middle]
  best_auc_hi_no_middle<- cv_model_no_middle$cvup[best_lambda_index_no_middle]
  
  cat("Best Lambda (Shortest-and-Longest):", best_lambda_no_middle, "\n")
  cat("AUC Mean (Shortest-and-Longest):", best_auc_mean_no_middle, "\n")
  cat("AUC CI (Shortest-and-Longest):", best_auc_lo_no_middle, "-", best_auc_hi_no_middle, "\n")
  
  # Append results to the dataframe
  results_df <- rbind(results_df, data.frame(
    Alpha = rep(alpha_val, 4),  # Repeat alpha for each condition
    Condition = c("Full", "Shortest", "Longest", "Shortest-and-Longest"),
    Lambda = c(best_lambda_full, best_lambda_shortest, best_lambda_longest, best_lambda_no_middle),
    AUC_Mean = c(best_auc_mean_full, best_auc_mean_shortest, best_auc_mean_longest, best_auc_mean_no_middle),
    AUC_LowerCI = c(best_auc_lo_full, best_auc_lo_shortest, best_auc_lo_longest, best_auc_lo_no_middle),
    AUC_UpperCI = c(best_auc_hi_full, best_auc_hi_shortest, best_auc_hi_longest, best_auc_hi_no_middle)
  ))
}

results_df_25 <- results_df

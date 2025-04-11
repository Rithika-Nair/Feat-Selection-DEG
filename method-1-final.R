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

# Selecting only genes found in de prior list
selected_genes <- intersect(colnames(x), de_prior$Gene_Name)
x_train_100_percent <- x_train[, selected_genes]
x_test_100_percent <- x_test[, selected_genes]

# List of possible alpha values
alpha_values <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

#################################
# Initialize an empty data frame
results_df <- data.frame()

### Loop through the elastic net code for different alpha values

for (alpha_val in alpha_values) {
  cat("\nRunning Elastic Net for alpha =", alpha_val, "\n")
  
  ###CONTROL CONDITION: Full Gene Set (Top 100% of genes from Pavilidis list)
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
  
  ### CONDITION A: Feature Removal (Removing Top n% Genes)
  removed_genes <- setdiff(colnames(x), de_prior_75$Gene_Name)
  x_train_removed <- x_train[, removed_genes]
  x_test_removed <- x_test[, removed_genes]
  
  # Train Elastic Net Model
  cv_model_removed <- cv.glmnet(x_train_removed, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_removed <- which.max(cv_model_removed$cvm)
  best_lambda_removed <- cv_model_removed$lambda[best_lambda_index_removed]
  best_auc_mean_removed <- cv_model_removed$cvm[best_lambda_index_removed]
  best_auc_lo_removed <- cv_model_removed$cvlo[best_lambda_index_removed]
  best_auc_hi_removed <- cv_model_removed$cvup[best_lambda_index_removed]
  
  cat("Best Lambda (Removed):", best_lambda_removed, "\n")
  cat("AUC Mean (Removed):", best_auc_mean_removed, "\n")
  cat("AUC CI (Removed):", best_auc_lo_removed, "-", best_auc_hi_removed, "\n")
  
  ### **CONDITION B: Feature Selection (Selecting Top n% Genes)**
  selected_genes <- intersect(colnames(x), de_prior_75$Gene_Name)
  x_train_selected <- x_train[, selected_genes]
  x_test_selected <- x_test[, selected_genes]
  
  # Train Elastic Net Model
  cv_model_selected <- cv.glmnet(x_train_selected, as.character(y_train$response), alpha = alpha_val, family = 'binomial', type.measure = "auc", nfolds = 5)
  
  # Extract best lambda and corresponding AUC
  best_lambda_index_selected <- which.max(cv_model_selected$cvm)
  best_lambda_selected <- cv_model_selected$lambda[best_lambda_index_selected]
  best_auc_mean_selected <- cv_model_selected$cvm[best_lambda_index_selected]
  best_auc_lo_selected <- cv_model_selected$cvlo[best_lambda_index_selected]
  best_auc_hi_selected <- cv_model_selected$cvup[best_lambda_index_selected]
  
  cat("Best Lambda (Selected):", best_lambda_selected, "\n")
  cat("AUC Mean (Selected):", best_auc_mean_selected, "\n")
  cat("AUC CI (Selected):", best_auc_lo_selected, "-", best_auc_hi_selected, "\n")
  
  # Append results to the dataframe
  results_df <- rbind(results_df, data.frame(
    Alpha = rep(alpha_val, 3),  # Repeat alpha for each condition
    Condition = c("Full", "Removed", "Selected"),
    Lambda = c(best_lambda_full, best_lambda_removed, best_lambda_selected),
    AUC_Mean = c(best_auc_mean_full, best_auc_mean_removed, best_auc_mean_selected),
    AUC_LowerCI = c(best_auc_lo_full, best_auc_lo_removed, best_auc_lo_selected),
    AUC_UpperCI = c(best_auc_hi_full, best_auc_hi_removed, best_auc_hi_selected)
  ))
}


result_df_25 <- results_df




####################





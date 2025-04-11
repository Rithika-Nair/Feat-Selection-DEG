# results_df_5$Quantile <- "5"
# results_df_10$Quantile <- "10"
# results_df_20$Quantile <- "20"
# results_df_25$Quantile <- "25"
# 
# combined_result_df_blca <- bind_rows(results_df_5, results_df_10, results_df_20, results_df_25)
# 
# head(combined_result_df_blca)
# 
# write.csv(combined_result_df_blca, file = "combined_result_method_3_blca.csv")

combined_result_df_blca <- read.csv("~/plots-final-3/combined_result_method_3_blca.csv")

combined_result_df_kipan <-
  combined_result_df_kipan %>%
  filter(Condition != "Shortest-and-Longest")
# Ensure Alpha and Quantile are factors to control their ordering
combined_result_df_blca$Alpha <- factor(combined_result_df_blca$Alpha, levels = unique(combined_result_df_blca$Alpha))
combined_result_df_blca$Quantile <- factor(combined_result_df_blca$Quantile, levels = c("5", "10", "20", "25"))

# Custom labeller for Quantile
quantile_labeller <- c("5" = "Top 5%", "10" = "Top 10%", "20" = "Top 20%", "25" = "Top 25%")

# Create the plot
plot_1_blca <-
  ggplot(combined_result_df_blca, aes(x = Condition, y = AUC_Mean, color = as.factor(Condition), group = Alpha)) +
  # Error bars (shaded region for confidence intervals)
  geom_errorbar(aes(ymin = AUC_LowerCI, ymax = AUC_UpperCI), 
                width = 0.2, position = position_dodge(0.4), size = 0.5) +  # Standard error bars
  # Points for mean AUC
  geom_point(position = position_dodge(0.4), size = 2) +
  labs(title = "Elastic Net Cross Validation AUC ",
       subtitle = "For Different Quantiles and Alpha values on BLCA Dataset",
       x = "Feature Selection Condition", y = "Mean AUC (Cross-Validation)") +
  # facet_wrap(~ Alpha, scales = "fixed", 
  #            labeller = labeller(Alpha = function(x) paste("alpha =", x))) +  # Custom facet labels for alpha
  theme_bw() +
  theme(legend.position = "none", 
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, lineheight = 1)) +  # Increase facet header size
  scale_color_manual(values = c("Full" = "#1710D4", 
                                "Shortest" = "#D10D0D", 
                                "Longest" = "#079228",
                                "Shortest-and-Longest" = "#8e44ad")) +
  facet_grid(~ Quantile, scales = "fixed", 
             labeller = labeller(Quantile = quantile_labeller))

ggsave("blca_plot_1m3.png", width = 14, height = 8, dpi = 300)

# Rank by cvm (mean AUC) for each condition and quantile, and select the best model (highest AUC)
best_results_df <- combined_result_df_blca %>%
  group_by(Condition, Quantile) %>%
  arrange(desc(AUC_Mean)) %>%
  slice(1) %>%  # Select the top row (best AUC) for each condition and quantile
  ungroup()

plot_2_blca <-
  ggplot(best_results_df, aes(x = Condition, y = AUC_Mean, color = Condition)) +
  # Error bars (shaded region for confidence intervals)
  geom_errorbar(aes(ymin = AUC_LowerCI, ymax = AUC_UpperCI), 
                width = 0.2, position = position_dodge(0.4), size = 0.5) +  # Standard error bars
  # Points for mean AUC
  geom_point(position = position_dodge(0.4), size = 2) +
  labs(title = "Best AUC Model Ranked by Mean AUC",
       subtitle = "For Different Quantiles on BLCA Dataset",
       x = "Feature Selection Condition", y = "Mean AUC (Cross-Validation)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, lineheight = 1)) +
  scale_color_manual(values = c("Full" = "#1710D4", 
                                "Shortest" = "#D10D0D", 
                                "Longest" = "#079228",
                                "Shortest-and-Longest" = "#8e44ad")) +
  facet_grid(~ Quantile, scales = "fixed", labeller = labeller(Quantile = quantile_labeller)) # Add facet for quantile

ggsave("blca_plot_2m3.png", width = 12, height = 6, dpi = 300)

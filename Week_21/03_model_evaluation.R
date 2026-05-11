#' 03_Model_Evaluation.R
#' 
#' This script evaluates the trained machine learning models and generates 
#' comparative figures and tables.
#' It covers the third phase of the lab practical:
#' Model Evaluation and Performance Metrics.
#'
#' Steps:
#' 1. Aggregating Performance Metrics
#' 2. Visualizing Performance Comparison (AUC, prAUC)
#' 3. Plotting ROC and Precision-Recall Curves
#' 4. Exporting Summary Tables

# Load configuration and dependencies
script_file <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_file[grep("^--file=", script_file)])
if (length(script_file) == 1) {
  base_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
} else {
  base_dir <- getwd()
}
source(file.path(base_dir, "config.R"))
library(mikropml)
library(dplyr)
library(ggplot2)
library(tidyr)

# 0. Setup Logging and Seed
log_file <- file.path(logs_dir, "03_model_evaluation.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting Model Evaluation Phase...\n")

# 1. Load Trained Models
input_file <- file.path(results_dir, "model_results.rds")
if (!file.exists(input_file)) {
  stop("Error: model_results.rds not found. Run 02_model_training.R first.")
}
model_list <- readRDS(input_file)

# 2. Combine Performance Data
# We extract the 'performance' dataframe from each model result.
# Novice Tip: Sometimes R packages return "NA" as a string instead of a missing value.
# We ensure all metrics are numeric before combining.
fix_perf_types <- function(df) {
  df %>%
    mutate(across(where(is.character) & !c(method), as.numeric))
}

perf_glmnet <- fix_perf_types(model_list$glmnet$performance)
perf_rf <- fix_perf_types(model_list$rf$performance)

combined_perf <- bind_rows(perf_glmnet, perf_rf)

cat("\nStep 1: Comparative Performance Table\n")
# Note: In mikropml, 'Mean_F1' might be named 'F1' depending on the version.
# Let's check for F1 or Mean_F1 dynamically.
f1_col <- if ("Mean_F1" %in% names(combined_perf)) "Mean_F1" else "F1"
print(combined_perf %>% select(method, AUC, prAUC, Accuracy, !!sym(f1_col)))

# Save table to results
write.csv(combined_perf, file.path(results_dir, "performance_comparison.csv"), row.names = FALSE)

# 3. Visualizing Comparison
cat("\nStep 2: Generating performance comparison plots...\n")

perf_long <- combined_perf %>%
  select(method, AUC, prAUC, Accuracy, !!sym(f1_col)) %>%
  pivot_longer(-method, names_to = "metric", values_to = "value")

p_comp <- ggplot(perf_long, aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  labs(title = "Model Performance Comparison", 
       subtitle = "Higher values indicate better performance",
       y = "Score", x = "ML Method") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path(figures_dir, "model_comparison.png"), p_comp, width = 8, height = 6)

# 4. ROC and Precision-Recall Curves
cat("\nStep 3: Calculating ROC and PRC curves...\n")

# Helper function to get sens/spec for plotting
get_curves_data <- function(result, method_name) {
  calc_model_sensspec(result$trained_model, result$test_data, "dx") %>%
    mutate(method = method_name)
}

roc_prc_glmnet <- get_curves_data(model_list$glmnet, "glmnet")
roc_prc_rf <- get_curves_data(model_list$rf, "rf")
combined_curves <- bind_rows(roc_prc_glmnet, roc_prc_rf)

# Plot ROC Curve
# Novice Tip: ROC curve shows the trade-off between Sensitivity and Specificity.
# A curve closer to the top-left corner is better.
p_roc <- combined_curves %>%
  ggplot(aes(x = specificity, y = sensitivity, color = method)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey") +
  scale_x_reverse() +
  theme_bw() +
  coord_fixed() +
  labs(title = "Receiver Operating Characteristic (ROC) Curve",
       x = "Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)")

ggsave(file.path(figures_dir, "roc_curve.png"), p_roc, width = 6, height = 6)

# Plot Precision-Recall Curve
# Novice Tip: PRC is often better for imbalanced biological data.
p_prc <- combined_curves %>%
  rename(recall = sensitivity) %>%
  ggplot(aes(x = recall, y = precision, color = method)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  coord_fixed() +
  labs(title = "Precision-Recall Curve (PRC)",
       x = "Recall (Sensitivity)",
       y = "Precision")

ggsave(file.path(figures_dir, "prc_curve.png"), p_prc, width = 6, height = 6)

cat("\nSuccess: Evaluation tables and figures saved to results/ and figures/\n")
cat("Model Evaluation Phase Complete.\n")

sink()

#' 04_Hyperparameter_Tuning.R
#' 
#' This script demonstrates the impact of hyperparameter tuning on model performance.
#' It covers the fourth phase of the lab practical:
#' Hyperparameter Tuning.
#'
#' Steps:
#' 1. Defining Custom Hyperparameter Grids
#' 2. Re-training models with custom grids
#' 3. Visualizing tuning performance (Grid Search results)
#' 4. Comparing "Default" vs "Tuned" performance

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
library(future)

# 0. Setup Logging, Seed, and Parallelization
log_file <- file.path(logs_dir, "04_hyperparameter_tuning.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting Hyperparameter Tuning Phase...\n")
future::plan(future::multicore, workers = 4)

# 1. Load Prepared Data and Previous Results
input_data <- readRDS(file.path(results_dir, "prepared_data.rds"))
previous_results <- readRDS(file.path(results_dir, "model_results.rds"))

# Keep learned preprocessing inside caret::train() to avoid estimating scaling,
# variance filters, or correlation filters from held-out test samples.
model_preprocess <- c("nzv", "corr", "center", "scale")

# 2. Define Custom Grids
# Novice Tip: Hyperparameters are the "knobs" we turn. 
# For GLMNET, 'lambda' controls the strength of the penalty.
# For RF, 'mtry' is the number of features sampled at each split.
cat("\nStep 1: Defining custom hyperparameter grids...\n")

# Custom lambda for GLMNET (Search a more granular range)
custom_hp_glmnet <- list(
  alpha = 0, # Keep L2 (Ridge)
  lambda = c(0.0001, 0.001, 0.01, 0.1, 1)
)

# Custom mtry for RF (Normally defaults to sqrt(features))
# Our dataset has 10 features, so default mtry is ~3. Let's test 2, 5, 10.
custom_hp_rf <- list(
  mtry = c(2, 5, 10)
)

# 3. Re-train with Custom Tuning
cat("\nStep 2: Re-training models with custom tuning grids...\n")

tuned_glmnet <- run_ml(
  dataset = input_data,
  method = "glmnet",
  outcome_colname = "dx",
  hyperparameters = custom_hp_glmnet,
  seed = seed,
  preProcess = model_preprocess,
  cv_times = 5 # Higher than before to see tuning stability
)

tuned_rf <- run_ml(
  dataset = input_data,
  method = "rf",
  outcome_colname = "dx",
  hyperparameters = custom_hp_rf,
  seed = seed,
  preProcess = model_preprocess,
  cv_times = 5
)

# 4. Visualize Tuning Progress
cat("\nStep 3: Generating tuning performance plots...\n")

# Use mikropml's built-in plotting for tuning
p_tune_glmnet <- plot_hp_performance(tuned_glmnet$trained_model$results, lambda, AUC) +
  labs(title = "GLMNET Tuning: AUC vs Lambda")
ggsave(file.path(figures_dir, "tuning_glmnet.png"), p_tune_glmnet)

p_tune_rf <- plot_hp_performance(tuned_rf$trained_model$results, mtry, AUC) +
  labs(title = "RF Tuning: AUC vs mtry")
ggsave(file.path(figures_dir, "tuning_rf.png"), p_tune_rf)

# 5. Compare Default vs Tuned
cat("\nStep 4: Comparing Default vs Tuned models...\n")

# Helper to fix types as seen in Script 03
fix_perf_types <- function(df) {
  df %>% mutate(across(where(is.character) & !c(method), as.numeric))
}

comparison_table <- bind_rows(
  fix_perf_types(previous_results$glmnet$performance) %>% mutate(version = "Default"),
  fix_perf_types(tuned_glmnet$performance) %>% mutate(version = "Tuned"),
  fix_perf_types(previous_results$rf$performance) %>% mutate(version = "Default"),
  fix_perf_types(tuned_rf$performance) %>% mutate(version = "Tuned")
)

print(comparison_table %>% select(method, version, AUC, prAUC))

# Save results
saveRDS(list(glmnet = tuned_glmnet, rf = tuned_rf), 
        file.path(results_dir, "tuned_model_results.rds"))
write.csv(comparison_table, 
          file.path(results_dir, "tuning_comparison_summary.csv"), row.names = FALSE)

cat("\nSuccess: Tuning figures and comparison saved.\n")
cat("Hyperparameter Tuning Phase Complete.\n")

sink()

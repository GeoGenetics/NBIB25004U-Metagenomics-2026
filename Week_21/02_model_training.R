#' 02_Model_Training.R
#' 
#' This script performs machine learning model training and experimental design.
#' It covers the second phase of the lab practical:
#' Model Training and Experimental Design.
#'
#' Steps:
#' 1. Data Partitioning (Internal vs. External Validation)
#' 2. Leakage-safe preprocessing inside model training
#' 3. k-Fold Cross-Validation (Hyperparameter Tuning)
#' 4. Model Execution (Baseline: Logistic Regression)
#' 5. Model Execution (Non-linear: Random Forest)
#' 6. Parallelization for Speed

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
library(future)

# 0. Setup Logging, Seed, and Parallelization
log_file <- file.path(logs_dir, "02_model_training.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting Model Training Phase...\n")

# Use multiple cores for training speed.
# Novice Tip: Training models can be computationally expensive. 
# Parallelization allows R to train multiple versions of the model simultaneously.
n_cores <- 4 # Adjust based on environment; students may use more.
future::plan(future::multicore, workers = n_cores)
cat("Parallel processing enabled with", n_cores, "cores.\n")

# 1. Load Prepared Data
input_file <- file.path(results_dir, "prepared_data.rds")
if (!file.exists(input_file)) {
  stop("Error: prepared_data.rds not found. Run 01_data_preparation.R first.")
}
dat <- readRDS(input_file)

# 2. Leakage-safe learned preprocessing
# caret estimates these operations from training data/resamples only, then applies
# the fitted transformation during prediction.
model_preprocess <- c("nzv", "corr", "center", "scale")
cat("Learned preprocessing inside caret::train():", 
    paste(model_preprocess, collapse = ", "), "\n")

# 3. Model Training Strategy: Logistic Regression (glmnet)
# Novice Tip: Logistic regression with L2 regularization (Ridge) is a great baseline.
# It is fast, interpretable, and less prone to overfitting than complex models.
cat("\nStep 3: Training Logistic Regression (glmnet)...\n")
cat("Strategy: 80/20 Train-Test split, 5-fold Cross-Validation (repeated 2x for this tutorial)\n")

# run_ml() handles partitioning and CV automatically.
results_glmnet <- run_ml(
  dataset = dat,
  method = "glmnet",
  outcome_colname = "dx",
  seed = seed,
  preProcess = model_preprocess,
  cv_times = 2 # Reduced for the tutorial run; should be 100 in real research
)

# 4. Model Training Strategy: Random Forest (rf)
# Novice Tip: Random Forest is a non-linear ensemble method. It's robust 
# but can take much longer to train.
cat("\nStep 4: Training Random Forest (rf)...\n")
results_rf <- run_ml(
  dataset = dat,
  method = "rf",
  outcome_colname = "dx",
  seed = seed,
  preProcess = model_preprocess,
  cv_times = 2 # Reduced for the tutorial run
)

# 4. Inspecting Initial Results
cat("\n--- Training Results Summary ---\n")
cat("GLMNET Performance (Internal CV AUC):", results_glmnet$performance$cv_metric_AUC, "\n")
cat("GLMNET Performance (Test Set AUC):", results_glmnet$performance$AUC, "\n")
cat("\nRF Performance (Internal CV AUC):", results_rf$performance$cv_metric_AUC, "\n")
cat("RF Performance (Test Set AUC):", results_rf$performance$AUC, "\n")

# Novice Tip: If CV AUC is 1.0 but Test AUC is 0.5, the model is OVERFITTING.
# It has memorized the training data but cannot generalize to new samples.

# 5. Save Model Results
output_file <- file.path(results_dir, "model_results.rds")
model_list <- list(
  glmnet = results_glmnet,
  rf = results_rf
)
saveRDS(model_list, file = output_file)

cat("\nSuccess: Trained models saved to results/model_results.rds\n")
cat("Model Training Phase Complete.\n")

sink()

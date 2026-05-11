#' 01_Data_Preparation.R
#' 
#' This script performs initial data cleaning for microbial metagenomics 
#' machine learning analysis. It covers the first phase of the lab practical:
#' Data Preparation.
#'
#' Steps:
#' 1. Cleaning Names and Columns
#' 2. Standardizing Feature Data Types
#' 3. Saving a prepared dataset for leakage-safe model training

# Load configuration and dependencies
script_file <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_file[grep("^--file=", script_file)])
if (length(script_file) == 1) {
  base_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
} else {
  base_dir <- getwd()
}
source(file.path(base_dir, "config.R"))
library(dplyr)
library(janitor)

# Ensure output directories exist before logging or saving artifacts.
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# 0. Setup Logging and Seed
log_file <- file.path(logs_dir, "01_data_preparation.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting Data Preparation Phase...\n")
cat("Using random seed:", seed, "\n")

# 1. Cleaning Names and Columns
# Novice Tip: Special characters or spaces in column names can break R functions.
# We use janitor::clean_names() to ensure all headers are snake_case and ASCII-safe.
cat("Step 1: Cleaning column names...\n")
cleaned_data <- raw_data %>%
  clean_names()

if (!"dx" %in% names(cleaned_data)) {
  stop("Error: expected outcome column `dx` was not found after cleaning names.")
}

# 2. Standardizing Feature Data Types
# Novice Tip: Models require numeric inputs. We ensure all OTU columns are numeric.
# In this example dataset they already are, but we apply this as a best practice.
# We leave dx as character because mikropml handles outcome conversion internally.
cat("Step 2: Ensuring numeric data types for features...\n")
cleaned_data <- cleaned_data %>%
  mutate(across(-dx, as.numeric))

# Learned preprocessing such as scaling, near-zero variance filtering, and
# correlation filtering is handled inside model training to avoid data leakage.
cat("Step 3: Deferring learned preprocessing to model training...\n")

# Summary of changes
cat("\n--- Preprocessing Summary ---\n")
cat("Sample count:", nrow(cleaned_data), "\n")
cat("Feature count:", ncol(cleaned_data) - 1, "\n")
cat("Outcome counts:\n")
print(table(cleaned_data$dx))

# Save intermediate output
output_file <- file.path(results_dir, "prepared_data.rds")
saveRDS(cleaned_data, file = output_file)

cat("\nSuccess: Prepared data saved to results/prepared_data.rds\n")
cat("Data Preparation Phase Complete.\n")

sink()

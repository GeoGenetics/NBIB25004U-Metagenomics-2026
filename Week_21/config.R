# Project Configuration

# Global seed for reproducibility
seed <- 42

if (!exists("base_dir", inherits = TRUE)) {
  base_dir <- getwd()
}

# Input data
# In a real scenario, this might be a path to a CSV. 
# For this lab, we use the built-in mikropml dataset.
library(mikropml)
data("otu_mini_bin")
raw_data <- otu_mini_bin

# Output paths
results_dir <- file.path(base_dir, "results")
logs_dir <- file.path(base_dir, "logs")
figures_dir <- file.path(base_dir, "figures")

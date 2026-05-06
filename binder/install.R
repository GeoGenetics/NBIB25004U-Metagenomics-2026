# Install CRAN packages
install.packages(c(
  "tidyverse",
  "ggplot2",
  "dplyr",
  "tidyr",
  "janitor",
  "future",
  "caret",
  "mikropml",
  "vegan",
  "reshape2",
  "pheatmap",
  "RColorBrewer",
  "knitr",
  "rmarkdown"
), repos = "https://cloud.r-project.org")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "phyloseq",
  "microbiome",
  "DESeq2",
  "apeglm"
), ask = FALSE)

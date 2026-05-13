# Install CRAN packages
cran_repo <- "https://packagemanager.posit.co/cran/__linux__/jammy/2024-01-16"
options(repos = c(CRAN = cran_repo))

cran_packages <- c(
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
  "ggh4x",
  "rmarkdown",
  "broom.mixed",
  "lme4",
  "lmerTest"
)

install.packages(cran_packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "phyloseq",
  "microbiome",
  "DESeq2",
  "apeglm"
)

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

required_packages <- c(cran_packages, bioc_packages)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages failed to install: ",
    paste(missing_packages, collapse = ", ")
  )
}

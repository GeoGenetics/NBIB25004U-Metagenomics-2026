options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "tidyverse",
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
  "lmerTest",
  "igraph"
)

install.packages(cran_packages, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "microbiome",
  "DESeq2",
  "apeglm"
)

BiocManager::install(
  bioc_packages,
  ask = FALSE,
  update = FALSE
)

required_packages <- c(cran_packages, bioc_packages, "BiocManager")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages failed to install: ",
    paste(missing_packages, collapse = ", ")
  )
}
# NBIB25004U Metagenomics — R Studio Environment

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GeoGenetics/NBIB25004U-Metagenomics-2026/r-studio?urlpath=rstudio)

This branch provides a cloud-based RStudio environment for the R-based practical sessions of the NBIB25004U Metagenomics course. It is pre-configured with all packages needed to work with phyloseq objects and microbiome data.

---

## Launching the RStudio environment

1. Click the **launch binder** badge above (or the button below).
2. Wait for the environment to build — this can take **5–10 minutes** on the first launch. Subsequent launches are faster if the image is cached.
3. Once loaded, an RStudio interface will open directly in your browser. No installation is required on your computer.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GeoGenetics/NBIB25004U-Metagenomics-2026/r-studio?urlpath=rstudio)

> **Note:** Binder sessions are temporary. Any changes you make will be lost when the session ends. Download your modified `.Rmd` files or rendered outputs before closing the browser tab.

---

## Working with the Rmd templates

Three R Markdown templates are provided in this repository:

| File | Topic |
|------|-------|
| [01_template.Rmd](01_template.Rmd) | Template 1 |
| [02_template.Rmd](02_template.Rmd) | Template 2 |
| [03_template.Rmd](03_template.Rmd) | Template 3 |

Open a template in RStudio, fill in the code chunks as instructed, and knit the document to produce an HTML report.

---

## Pre-installed R packages

The environment is built from R 4.3 and includes the following packages:

### Bioconductor
| Package | Purpose |
|---------|---------|
| [phyloseq](https://joey711.github.io/phyloseq/) | Core data structure for microbiome data (OTU tables, taxonomy, sample metadata, phylogenetic trees) |
| [microbiome](https://microbiome.github.io/tutorials/) | Additional utilities for phyloseq objects |
| [DESeq2](https://bioconductor.org/packages/DESeq2) | Differential abundance analysis |
| [apeglm](https://bioconductor.org/packages/apeglm) | Log-fold change shrinkage for DESeq2 |

### CRAN
| Package | Purpose |
|---------|---------|
| tidyverse | Data wrangling and plotting (dplyr, ggplot2, tidyr, ...) |
| vegan | Multivariate ecology statistics (ordinations, diversity) |
| reshape2 | Data reshaping |
| pheatmap | Heatmap visualisation |
| RColorBrewer | Colour palettes |
| knitr / rmarkdown | Reproducible reports |

---

## Troubleshooting

- **Build takes too long / times out:** Binder has a 10-minute build limit. Try again; the image may already be cached.
- **Session disconnected:** Re-open the Binder link. Work is not saved automatically — always download files you want to keep.
- **Package not found:** Open an issue in this repository and the missing package will be added to `binder/install.R`.

---

## Course materials

Full course instructions are available in the Wiki:

👉 [NBIB25004U-Metagenomics-2026 Wiki](https://github.com/GeoGenetics/NBIB25004U-Metagenomics-2026/wiki)

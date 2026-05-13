# 🧬 Differential Abundance Analysis Workflow

## 📌 Overview

In this practical you will analyse differential abundance patterns in microbiome data using an already assembled **phyloseq** object. The session builds directly on the preceding theory session: raw read counts must first be normalised, microbiome abundance values violate simple test assumptions, thousands of taxon-level tests require FDR correction, and relative abundances are compositional. The practical starts with basic data inspection and exploratory visualisation, then compares four increasingly sophisticated approaches for identifying taxa that differ between **infant** and **adult** gut microbiomes.

The workflow has eight conceptual phases:

1. **[Load and inspect the phyloseq object](#-1-load-and-inspect-the-data)**
2. **[Understand read counts and choose an appropriate normalisation](#-2-visualise-community-composition)**
3. **[Visualise community composition](#-2-visualise-community-composition) and [individual taxa](#-3-explore-individual-taxa)**
4. **[Diagnose the statistical properties of microbiome abundance data](#-4-diagnose-the-data-before-testing)**
5. **[Test age-associated taxa with a Wilcoxon rank-sum test](#-5-wilcoxon-rank-sum-test)**
6. **[Fit linear models](#-6-age-difference-using-log-transformed-linear-models) with and [without health-status adjustment](#-7-adjust-the-linear-model-for-health-status)**
7. **[Control false discoveries](#-8-control-false-discoveries-with-fdr) and [account for compositionality](#-9-compositional-analysis-with-clr-transformation)**
8. **[Use constrained ordination to identify community-level drivers](#-11-constrained-ordination-with-rda)**

---

## 🚀 Launch the Analysis Environment

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GeoGenetics/NBIB25004U-Metagenomics-2026/r-studio?urlpath=rstudio)

You can run this practical in either a local RStudio installation or in the prepared **Binder RStudio environment** for the course.

| Option | What you need | When to use it |
| --- | --- | --- |
| **Local RStudio** | Clone this repository and install all required R dependencies | Best if you want a persistent local copy and already have R/RStudio configured |
| **Binder RStudio** | No local setup; click the Binder badge above | Best for the practical session because all dependencies, data, and code are already available |

---

### 🔹 Step 1 — Choose where to run the practical

**Option A: Run locally in RStudio**

Clone the repository to your computer, open the project folder in RStudio, and make sure all packages used by the practical are installed. The required packages are listed in the setup chunk of:

```text
week21_differential_abundance/differential_abundance.Rmd
```

**Option B: Run in Binder**

Click the **launch binder** badge above. The first launch can take **5-10 minutes** while the environment builds; later launches are usually faster if Binder has cached the image.

---

### 🔹 Step 2 — Open the practical Rmd

When RStudio opens in the browser, use the **Files** pane to open:

```text
week21_differential_abundance/differential_abundance.Rmd
```

You can run the analysis in either of two ways:

- copy-paste code chunks from the Rmd into the R console,
- or use the Rmd directly and run chunks inside the document with the green chunk-run buttons.

If you are not using the green chunk-run buttons in RStudio, you can also copy-paste the `r` code chunks directly from this wiki page into the R console. Use the main visible chunks for running the analysis; the collapsed **Annotated code** chunks are there to help you understand what each line does.

After working through the analysis, you can click **Knit** to render the full practical report.

---

### 🔹 Step 3 — Save anything you want to keep

Binder sessions are temporary. Any edits, rendered HTML files, or notes created inside the session will disappear when the session shuts down.

> ⚠️ **Before closing Binder**
> Download your modified `.Rmd`, rendered `.html`, or any notes/results you want to keep.

---

> 🧬 **In-class vs after-class**
> - **During the class** you will run `week21_differential_abundance/differential_abundance.Rmd` chunk by chunk in RStudio. The main objective is to understand what each method assumes, what each result table contains, and how the plots should be interpreted.
> - **After the class**, re-run the practical by changing the example taxa selected in the exploratory plots. Compare whether the same biological conclusions hold for different taxa and across different statistical methods.
> - The formal differential abundance tests focus on **gut samples** and compare `age = infant` against `age = adult`. Other body sites are used at the beginning to understand the broader community structure.

> ⚙️ **Working directory layout used by the Rmd**
> ```
> NBIB25004U-Metagenomics-2026/
> ├── data/
> │   ├── phyloseq_augmented.rds          # Main input object for this practical
> │   ├── mag_table.tsv                   # Abundance table used in earlier sessions
> │   ├── sample_data.tsv                 # Sample metadata
> │   └── tax_table.tsv                   # Taxonomic annotation
> │
> └── week21_differential_abundance/
>     └── differential_abundance.Rmd      # Practical script
> ```

> 💡 **Important**
> Run the chunks in order. Several later analyses depend on objects created earlier, especially `dat_rel`, `wilcox_results`, `glm_age_results`, `glm_age_adjusted_status_results`, `dat_clr`, `clr_age_adjusted_status_results`, and `comparison_results`.

---

## 🧪 1. Load and Inspect the Data

The practical starts from a saved **phyloseq** object. A phyloseq object bundles together the abundance table, taxonomy table, sample metadata, and optionally a phylogenetic tree. This keeps the different layers of a microbiome dataset synchronised during filtering, transformation, and plotting.

### 🔹 Step 1 — Open the practical in RStudio

Open:

```text
week21_differential_abundance/differential_abundance.Rmd
```

The repository Binder/RStudio environment already contains the packages needed for the R-based practical sessions.

---

### 🔹 Step 2 — Load the required R packages

The Rmd begins by loading the packages used for microbiome data handling, data wrangling, plotting, model fitting, and ordination:

```r
library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggh4x)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(microbiome)
```

<details>
<summary>Annotated code</summary>

```r
# Load phyloseq data structures and microbiome-specific helper functions.
library(phyloseq)

# Load ecological analysis tools used for ordination and community statistics.
library(vegan)

# Load core data-wrangling verbs such as filter(), mutate(), and summarise().
library(dplyr)

# Load tidy data-frame helpers, including row-name and tibble conversion tools.
library(tibble)

# Load functional programming helpers used for repeated model fitting and mapping.
library(purrr)

# Load string manipulation helpers for cleaning and matching taxonomic labels.
library(stringr)

# Load the main plotting system used throughout the practical.
library(ggplot2)

# Load plot-composition tools for combining multiple ggplot panels.
library(patchwork)

# Load the broader tidyverse toolkit used for wrangling, plotting, and importing data.
library(tidyverse)

# Load extra ggplot scale and facet tools used to customise figures.
library(ggh4x)

# Load mixed-effect model fitting functions.
library(lme4)

# Load p-values and test summaries for mixed-effect models.
library(lmerTest)

# Load tidy summaries for mixed-model outputs.
library(broom.mixed)

# Load microbiome-specific transformations and utility functions.
library(microbiome)
```

</details>

The most important packages for this practical are `phyloseq`, `tidyverse`, `vegan`, `ggh4x`, and the model-summary tools used through `broom::tidy()`.

---

### 🔹 Step 3 — Load the phyloseq object

```r
ps <- readRDS("data/phyloseq_augmented.rds")

sample_metadata <- sample_data(ps) %>%
  as("data.frame") %>%
  rownames_to_column("sample_id") %>%
  mutate(
    body_site_original = as.character(body_site),
    body_site = if_else(str_detect(body_site_original, "^gut_"), "gut", body_site_original),
    age = case_when(
      str_detect(body_site_original, "infant") ~ "infant",
      str_detect(body_site_original, "adult") ~ "adult",
      TRUE ~ NA_character_
    ),
    status = if_else(as.logical(simulated), "sick", "healthy"),
    age = factor(age, levels = c("infant", "adult")),
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  column_to_rownames("sample_id")

sample_data(ps) <- sample_data(sample_metadata)

ps
sample_data(ps)
rank_names(ps)
sample_variables(ps)
```

<details>
<summary>Annotated code</summary>

```r
# Read the saved phyloseq object from the course data directory.
ps <- readRDS("data/phyloseq_augmented.rds")

# Convert the compact teaching metadata into the variables used below.
sample_metadata <- sample_data(ps) %>%
  as("data.frame") %>%
  rownames_to_column("sample_id") %>%
  mutate(
    # Preserve the original body-site labels from the input object.
    body_site_original = as.character(body_site),
    # Collapse gut_adult and gut_infant into the common body_site value gut.
    body_site = if_else(str_detect(body_site_original, "^gut_"), "gut", body_site_original),
    # Derive the age comparison used in the gut differential-abundance tests.
    age = case_when(
      str_detect(body_site_original, "infant") ~ "infant",
      str_detect(body_site_original, "adult") ~ "adult",
      TRUE ~ NA_character_
    ),
    # Derive the teaching health-status variable used for covariate adjustment.
    status = if_else(as.logical(simulated), "sick", "healthy"),
    # Set reference levels used by the model formulas later in the practical.
    age = factor(age, levels = c("infant", "adult")),
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  column_to_rownames("sample_id")

# Replace the original sample metadata with the augmented table.
sample_data(ps) <- sample_data(sample_metadata)

# Print a compact summary of the phyloseq object, including sample and taxon counts.
ps

# Inspect the sample metadata table attached to the phyloseq object.
sample_data(ps)

# List the available taxonomic ranks, such as kingdom, phylum, genus, and species.
rank_names(ps)

# List the metadata variables that can be used for filtering, plotting, or modelling.
sample_variables(ps)
```

</details>

This first inspection tells you:

- how many samples are present,
- how many taxa are present,
- which taxonomic ranks are available,
- which sample metadata variables can be used in downstream analyses.

The key metadata variables for this practical are:

| Variable | Role in the practical |
| --- | --- |
| `body_site_original` | Original input label, such as `gut_adult`, `gut_infant`, or `vaginal` |
| `body_site` | Used to compare broad community composition and then focus on gut samples |
| `age` | Main variable of interest: adult vs infant |
| `status` | Potential confounder: healthy vs sick |

The input object stores the teaching groups compactly. The practical therefore derives `age` from the original gut labels, collapses `gut_adult` and `gut_infant` into `body_site = gut`, and creates a teaching `status` variable from the `simulated` flag so that later models can demonstrate covariate adjustment.

> 💡 **Why start with inspection?**
> Many downstream errors in microbiome analysis come from mismatched sample names, missing metadata, or unexpected factor levels. Inspecting the object before modelling confirms that the data structure matches the biological question.

---

### 🔹 Step 4 — Remember what read counts measure

The abundance layer in a microbiome object is a count matrix: taxa are rows, samples are columns, and the values are integer read counts. A read count does **not** directly measure the number of microbial cells or genomes in a sample. It measures how many sequencing reads happened to be assigned to a taxon after extraction, library preparation, sequencing, and classification.

Two samples with identical biological composition can have very different total read counts because of technical factors such as pooling, PCR efficiency, multiplexing, and sequencing depth. This is why raw counts should not be compared sample-to-sample without normalisation.

> 💡 **Theory link**
> Library size is technical, not biology. The count of taxon *i* in sample *j* only becomes interpretable relative to the other counts in sample *j*.

---

## 📊 2. Visualise Community Composition

Before testing individual taxa, we first look at the overall community composition. This is a sanity check: if gut, oral, skin, or vaginal samples are compositionally distinct, that should be visible before formal modelling.

### 🔹 Step 1 — Choose the normalisation for the question

The theory session introduced three common normalisation strategies:

| Normalisation | What it does | Useful for | Main caution |
| --- | --- | --- | --- |
| Relative abundance | Divides each sample by its total count | Intuitive visualisation and first-pass summaries | Creates compositional data constrained to sum to one |
| Rarefaction | Subsamples every sample to the same depth | Some alpha-diversity estimators | Throws away reads and is not used here for differential abundance |
| CLR | Compares each taxon to the sample geometric mean | Compositionally aware modelling and ordination | Requires a pseudo-count for zeros and changes the scale |

In this practical, relative abundance is used for visualisation and some first-pass tests. CLR is used later as the compositionally aware reference analysis. Rarefaction is deliberately not used for differential abundance because discarding reads can increase both false positives and false negatives.

---

### 🔹 Step 2 — Convert counts to relative abundance

```r
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_rel <- as.data.frame(t(otu_table(ps_rel)))
meta_rel <- as.data.frame(sample_data(ps_rel))
tax_rel <- as.data.frame(tax_table(ps_rel))
```

<details>
<summary>Annotated code</summary>

```r
# Convert raw counts to relative abundances within each sample.
# For each sample, every count is divided by that sample's total count.
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Extract the relative-abundance table and transpose it so samples are rows.
otu_rel <- as.data.frame(t(otu_table(ps_rel)))

# Extract sample metadata as a regular data frame for joining and plotting.
meta_rel <- as.data.frame(sample_data(ps_rel))

# Extract the taxonomy table as a regular data frame for adding taxonomic labels.
tax_rel <- as.data.frame(tax_table(ps_rel))
```

</details>

Raw read counts are affected by sequencing depth. Relative abundance places every sample on the same 0-1 scale by dividing each taxon count by the total count in that sample.

> ⚠️ **Relative abundance is useful, but not perfect**
> Relative abundances are easy to visualise, but they are **compositional**: all taxa in a sample must sum to one. This means an apparent decrease in one taxon can be caused by an increase in another. The practical returns to this issue in the CLR section.

---

### 🔹 Step 3 — Aggregate taxa to genus level

```r
ps_genus <- tax_glom(ps_rel, taxrank = "genus")

plot_bar(ps_genus, fill = "genus") +
  facet_nested(~ body_site + age + status, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank())
```

<details>
<summary>Annotated code</summary>

```r
# Collapse taxa that share the same genus into a single genus-level abundance.
ps_genus <- tax_glom(ps_rel, taxrank = "genus")

# Draw a stacked barplot where each colour represents a genus.
plot_bar(ps_genus, fill = "genus") +
  # Split the plot by body site, age group, and health status.
  facet_nested(~ body_site + age + status, scales = "free_x") +
  # Use a simple high-contrast ggplot theme.
  theme_bw() +
  # Hide sample labels on the x-axis because there are too many to read clearly.
  theme(axis.text.x = element_blank())
```

</details>

Aggregating to genus level reduces visual clutter while preserving a biologically interpretable taxonomic resolution. The stacked barplot lets you compare major community shifts across body sites, age groups, and health status.

> 🧠 Discussion
> Look for genera that dominate specific body sites. Do gut samples look compositionally different from other body sites? Are infant gut samples visually distinct from adult gut samples before any statistical test is run?

---

### 🔹 Step 4 — Reshape the data into long format

The Rmd creates a tidy long-format table called `dat_rel`:

```r
dat_rel <- otu_rel %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "taxon", values_to = "rel_abund") %>%
  left_join(meta_rel %>% as("data.frame") %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  left_join(tax_rel %>% rownames_to_column("taxon"), by = "taxon")
```

<details>
<summary>Annotated code</summary>

```r
# Start from the relative-abundance matrix, where rows are samples and columns are taxa.
dat_rel <- otu_rel %>%
  # Move sample IDs from row names into an explicit column.
  rownames_to_column("sample_id") %>%
  # Convert the wide matrix into long format: one row per sample-taxon pair.
  pivot_longer(-sample_id, names_to = "taxon", values_to = "rel_abund") %>%
  # Add sample metadata such as body_site, age, and status.
  left_join(meta_rel %>% as("data.frame") %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  # Add taxonomic annotations for each taxon.
  left_join(tax_rel %>% rownames_to_column("taxon"), by = "taxon")
```

</details>

Each row now represents one **sample-taxon combination**, with the taxon's relative abundance, taxonomy, and sample metadata all in the same table.

> 💡 **Why long format?**
> Long format is ideal for `dplyr` and `ggplot2`: it makes it easy to filter one taxon, group by taxon, split by age, and join taxonomic labels to statistical results.

---

## 🔍 3. Explore Individual Taxa

Formal tests are easier to interpret if you first look at individual taxon abundance profiles. The Rmd selects one taxon by index and plots it across body sites and then across age groups within the gut.

### 🔹 Step 1 — Pick one taxon

```r
one_taxon <- dat_rel$taxon[3]
# one_taxon <- dat_rel$taxon[5]
# one_taxon <- dat_rel$taxon[20]
```

<details>
<summary>Annotated code</summary>

```r
# Select the third taxon in dat_rel as the example taxon for exploratory plots.
one_taxon <- dat_rel$taxon[3]

# Alternative example taxa: remove the # from one line to explore a different taxon.
# one_taxon <- dat_rel$taxon[5]
# one_taxon <- dat_rel$taxon[20]
```

</details>

Changing the index lets you explore how different taxa behave. Some taxa will be widespread; others will be nearly absent outside a specific niche.

---

### 🔹 Step 2 — Compare the taxon across body sites

```r
dat_rel %>%
  filter(taxon == one_taxon) %>%
  ggplot(aes(x = body_site, y = rel_abund, fill = body_site)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    theme_bw()
```

<details>
<summary>Annotated code</summary>

```r
# Start from the long relative-abundance table.
dat_rel %>%
  # Keep only the taxon selected above.
  filter(taxon == one_taxon) %>%
  # Plot relative abundance by body site, with matching fill colours.
  ggplot(aes(x = body_site, y = rel_abund, fill = body_site)) +
    # Draw boxplots but suppress outlier dots because raw points are added below.
    geom_boxplot(outlier.shape = NA) +
    # Add individual samples with slight horizontal jitter to reduce overlap.
    geom_jitter(width = 0.15) +
    # Use a simple high-contrast theme.
    theme_bw()
```

</details>

The boxplot summarises the distribution in each body site, while the jittered points show the raw sample-level values.

> 💡 **Why show points on top of boxplots?**
> Microbiome datasets often have small sample sizes and many zeros. The raw points reveal whether a pattern is shared by many samples or driven by one or two high-abundance outliers.

---

### 🔹 Step 3 — Compare the same taxon between infant and adult gut samples

```r
dat_rel %>%
  filter(taxon == one_taxon, body_site == "gut") %>%
  ggplot(aes(x = age, y = rel_abund, fill = body_site)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    theme_bw()
```

<details>
<summary>Annotated code</summary>

```r
# Start from the long relative-abundance table.
dat_rel %>%
  # Keep the selected taxon and restrict the comparison to gut samples.
  filter(taxon == one_taxon, body_site == "gut") %>%
  # Plot relative abundance by age group within the gut samples.
  ggplot(aes(x = age, y = rel_abund, fill = body_site)) +
    # Show the group-level distribution while hiding automatic outlier dots.
    geom_boxplot(outlier.shape = NA) +
    # Overlay each sample so sparse or outlier-driven patterns are visible.
    geom_jitter(width = 0.15) +
    # Use a simple high-contrast theme.
    theme_bw()
```

</details>

This plot gives a first visual impression of whether a taxon is age-associated in the gut. It does not replace statistical testing, but it helps you understand what the later model estimates mean.

---

## 🧪 4. Diagnose the Data Before Testing

Microbiome abundance data rarely satisfy the assumptions of simple parametric tests. The practical demonstrates three common features:

1. **Zero-inflation**: many taxa are absent from many samples.
2. **Right-skewed distributions**: most non-zero values are small, with a few large values.
3. **Departure from normality**: taxon abundances do not follow a Gaussian distribution.

> 🧠 **Why not just run a t-test?**
> A t-test compares the difference between two group means to the spread within the groups. It is calibrated for roughly continuous, symmetric distributions with comparable variance. Microbiome relative abundances usually have a spike at zero, a long right tail, and strong non-normality, so the t-test's assumptions are often visibly wrong before any formal test is run.

### 🔹 Step 1 — Check zero-inflation

```r
dat_rel %>%
  filter(body_site == "gut") %>%
  group_by(taxon) %>%
  summarise(prop_zero = mean(rel_abund == 0), .groups = "drop") %>%
  ggplot(aes(x = prop_zero)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
    theme_bw()
```

<details>
<summary>Annotated code</summary>

```r
# Start from the long relative-abundance table.
dat_rel %>%
  # Focus on gut samples because the differential abundance tests compare gut samples.
  filter(body_site == "gut") %>%
  # Calculate zero frequency separately for each taxon.
  group_by(taxon) %>%
  # mean(rel_abund == 0) gives the proportion of samples where the taxon is absent.
  summarise(prop_zero = mean(rel_abund == 0), .groups = "drop") %>%
  # Plot the distribution of zero proportions across taxa.
  ggplot(aes(x = prop_zero)) +
    # Use bins of width 0.05 so each bin represents a 5 percentage-point interval.
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
    # Use a simple high-contrast theme.
    theme_bw()
```

</details>

This histogram shows the proportion of gut samples where each taxon has zero relative abundance.

> ⚠️ **Why zero-inflation matters**
> A taxon present in only a few samples cannot be well represented by a symmetric continuous distribution. This is one reason why microbiome differential abundance is more difficult than applying a simple `t.test()` to every taxon.

---

### 🔹 Step 2 — Inspect skewness among non-zero values

```r
dat_rel %>%
  filter(body_site == "gut", taxon %in% example_taxa, rel_abund > 0) %>%
  mutate(label = paste(genus, species)) %>%
  ggplot(aes(x = rel_abund)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    facet_wrap(~ label, scales = "free") +
    theme_bw()
```

<details>
<summary>Annotated code</summary>

```r
# Start from the long relative-abundance table.
dat_rel %>%
  # Keep gut samples, selected example taxa, and non-zero abundance values.
  filter(body_site == "gut", taxon %in% example_taxa, rel_abund > 0) %>%
  # Combine genus and species into a readable facet label.
  mutate(label = paste(genus, species)) %>%
  # Plot the distribution of non-zero relative abundances.
  ggplot(aes(x = rel_abund)) +
    # Draw one histogram per taxon.
    geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    # Give each taxon its own panel and scale because abundances can differ strongly.
    facet_wrap(~ label, scales = "free") +
    # Use a simple high-contrast theme.
    theme_bw()
```

</details>

Even when zeros are removed, most taxa remain strongly right-skewed. This means that means and variances can be heavily influenced by a few samples.

---

### 🔹 Step 3 — Use Q-Q plots to assess normality

```r
dat_rel %>%
  filter(body_site == "gut", taxon %in% example_taxa) %>%
  mutate(label = paste(genus, species)) %>%
  ggplot(aes(sample = rel_abund)) +
    stat_qq(size = 0.8, alpha = 0.6) +
    stat_qq_line(color = "red", linetype = "dashed") +
    facet_wrap(~ label, scales = "free") +
    theme_bw()
```

<details>
<summary>Annotated code</summary>

```r
# Start from the long relative-abundance table.
dat_rel %>%
  # Keep gut samples and the selected example taxa.
  filter(body_site == "gut", taxon %in% example_taxa) %>%
  # Combine genus and species into a readable facet label.
  mutate(label = paste(genus, species)) %>%
  # Build Q-Q plots by comparing observed abundances to normal-theory quantiles.
  ggplot(aes(sample = rel_abund)) +
    # Plot the observed quantiles.
    stat_qq(size = 0.8, alpha = 0.6) +
    # Add the reference line expected under approximate normality.
    stat_qq_line(color = "red", linetype = "dashed") +
    # Show each taxon separately with its own scale.
    facet_wrap(~ label, scales = "free") +
    # Use a simple high-contrast theme.
    theme_bw()
```

</details>

If points followed the dashed line, a normal model would be plausible. Strong curvature or J-shaped patterns indicate non-normality.

> 🧠 Discussion
> Which assumption fails most clearly: absence of zeros, symmetry, or normality? How might that affect p-values from a standard parametric test?

---

## 🧬 5. Wilcoxon Rank-Sum Test

The **Wilcoxon rank-sum test** is a non-parametric alternative to a two-sample t-test. Instead of modelling the raw abundance values directly, it ranks all observations and asks whether one group tends to have higher ranks than the other.

This makes it more robust to skewness, outliers, and non-normality.

> 🧠 **What the Wilcoxon test is doing**
> The test pools both groups, ranks all observations from smallest to largest, sums the ranks for one group, and asks whether that rank sum is more extreme than expected under the null hypothesis that group labels are exchangeable. With small samples the exact null distribution can be enumerated; with realistic sample sizes, `wilcox.test(..., exact = FALSE)` uses a normal approximation to the rank-sum statistic.

> ⚠️ **Effect size note**
> The Wilcoxon test itself gives a p-value, not a biological fold change. The practical therefore computes `log2_fc` separately from group means so that significance and enrichment magnitude can be interpreted together.

### 🔹 Step 1 — Run one Wilcoxon test per taxon

```r
wilcox_results <- dat_rel %>%
  filter(body_site == "gut") %>%
  group_by(taxon, genus, species) %>%
  summarise(
    p_value = wilcox.test(rel_abund ~ age, exact = FALSE)$p.value,
    mean_adult = mean(rel_abund[age == "adult"]),
    mean_infant = mean(rel_abund[age == "infant"]),
    log2_fc = log2((mean_adult + 1e-6) / (mean_infant + 1e-6)),
    .groups = "drop"
  ) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))
```

<details>
<summary>Annotated code</summary>

```r
# Create a result table with one Wilcoxon test per taxon.
wilcox_results <- dat_rel %>%
  # Restrict the analysis to gut samples.
  filter(body_site == "gut") %>%
  # Run the same summary and test separately for every taxon.
  group_by(taxon, genus, species) %>%
  summarise(
    # Test whether relative abundance ranks differ between age groups.
    p_value = wilcox.test(rel_abund ~ age, exact = FALSE)$p.value,
    # Calculate the mean relative abundance in adult gut samples.
    mean_adult = mean(rel_abund[age == "adult"]),
    # Calculate the mean relative abundance in infant gut samples.
    mean_infant = mean(rel_abund[age == "infant"]),
    # Convert the adult/infant abundance ratio to a log2 fold change.
    log2_fc = log2((mean_adult + 1e-6) / (mean_infant + 1e-6)),
    # Drop the grouping structure after summarising.
    .groups = "drop"
  ) %>%
  # Correct raw p-values for multiple testing using Benjamini-Hochberg FDR.
  mutate(q_value = p.adjust(p_value, method = "BH"))
```

</details>

The result table contains:

| Column | Meaning |
| --- | --- |
| `p_value` | Raw Wilcoxon p-value for adult vs infant |
| `q_value` | Benjamini-Hochberg adjusted p-value |
| `mean_adult` | Mean relative abundance in adult gut samples |
| `mean_infant` | Mean relative abundance in infant gut samples |
| `log2_fc` | Log2 fold change: adult / infant |

> 💡 **Interpreting `log2_fc`**
> Positive values indicate adult enrichment. Negative values indicate infant enrichment. The pseudo-count `1e-6` prevents division by zero.

---

### 🔹 Step 2 — Inspect the strongest infant- and adult-enriched taxa

```r
wilcox_infant_enriched <- wilcox_results %>%
  select(taxon, genus, species, mean_adult, mean_infant, log2_fc, p_value, q_value) %>%
  arrange(log2_fc) %>%
  head(10)

wilcox_adult_enriched <- wilcox_results %>%
  select(taxon, genus, species, mean_adult, mean_infant, log2_fc, p_value, q_value) %>%
  arrange(-log2_fc) %>%
  head(10)
```

<details>
<summary>Annotated code</summary>

```r
# Find the taxa with the strongest negative log2 fold changes.
wilcox_infant_enriched <- wilcox_results %>%
  # Keep the identifiers, effect size, and significance columns needed for interpretation.
  select(taxon, genus, species, mean_adult, mean_infant, log2_fc, p_value, q_value) %>%
  # Sort from most infant-enriched to least infant-enriched.
  arrange(log2_fc) %>%
  # Keep the top 10 rows.
  head(10)

# Find the taxa with the strongest positive log2 fold changes.
wilcox_adult_enriched <- wilcox_results %>%
  # Keep the same columns so the two tables can be compared directly.
  select(taxon, genus, species, mean_adult, mean_infant, log2_fc, p_value, q_value) %>%
  # Sort from most adult-enriched to least adult-enriched.
  arrange(-log2_fc) %>%
  # Keep the top 10 rows.
  head(10)
```

</details>

The top infant-enriched taxa are those with the most negative fold changes. The top adult-enriched taxa are those with the most positive fold changes.

> 🧠 Discussion
> Do the top infant-associated taxa match what you expect biologically? In human gut studies, infant-associated taxa often include early colonisers and milk-associated bacteria, while adult-associated taxa often include strict anaerobes adapted to a more complex diet.

---

### 🔹 Step 3 — Use volcano plots and ranked bar charts

The volcano plot combines effect size and significance:

- x-axis: `log2_fc`
- y-axis: `-log10(q_value)`
- colour: adult-enriched, infant-enriched, or not significant

Taxa in the upper-left are significant and infant-enriched. Taxa in the upper-right are significant and adult-enriched.

> ⚠️ **Limitation of the Wilcoxon test**
> The Wilcoxon test treats every taxon independently, does not model additional covariates such as health status, and does not address the compositional structure of relative abundance data. The next sections address these limitations.

---

## 📈 6. Age Difference Using Log-Transformed Linear Models

The next method fits a Gaussian model to log-transformed relative abundance. This is more parametric than the Wilcoxon test, but it has an important advantage: model formulas can be extended to include additional covariates.

### 🔹 Step 1 — Understand logs and pseudo-counts

Relative abundances span orders of magnitude. Taking logs compresses the right tail: a 10-fold change becomes the same-sized step on the log scale whether it occurs at high or low abundance.

The problem is that microbiome data contain many zeros, and `log(0)` is undefined. The standard workaround is to add a small pseudo-count before logging:

```text
log10(relative abundance + ε)
```

| Pseudo-count choice | Consequence |
| --- | --- |
| Very small `ε` | Zeros sit far below non-zero values and may dominate sparse taxa |
| Larger `ε` | Low real abundances and zeros become harder to distinguish |
| `1e-6` in this Rmd | A common compromise for relative-abundance-scale data |

> ⚠️ **Sanity check**
> Pseudo-counts are not biologically measured values. If a conclusion depends strongly on the exact value of `ε`, treat it as fragile.

---

### 🔹 Step 2 — Fit a model for age only

```r
glm_age_results <- dat_rel %>%
  filter(body_site == "gut") %>%
  mutate(age = factor(age, levels = c("infant", "adult"))) %>%
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    fit <- glm(
      log10(rel_abund + 1e-6) ~ age,
      data = .x,
      family = gaussian()
    )

    broom::tidy(fit)
  }) %>%
  ungroup() %>%
  filter(term == "ageadult") %>%
  mutate(
    q_value = p.adjust(p.value, method = "BH"),
    fold_change_adult_vs_infant = 10^estimate
  )
```

<details>
<summary>Annotated code</summary>

```r
# Create a result table with one age-only linear model per taxon.
glm_age_results <- dat_rel %>%
  # Restrict the analysis to gut samples.
  filter(body_site == "gut") %>%
  # Set infant as the reference level, so the ageadult coefficient means adult minus infant.
  mutate(age = factor(age, levels = c("infant", "adult"))) %>%
  # Fit the same model separately for every taxon.
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    # Model log10-transformed relative abundance as a function of age.
    fit <- glm(
      log10(rel_abund + 1e-6) ~ age,
      data = .x,
      family = gaussian()
    )

    # Convert the model summary into a tidy data frame.
    broom::tidy(fit)
  }) %>%
  # Remove grouping so later operations act on the full result table.
  ungroup() %>%
  # Keep only the adult-vs-infant age coefficient.
  filter(term == "ageadult") %>%
  mutate(
    # Correct p-values across taxa using Benjamini-Hochberg FDR.
    q_value = p.adjust(p.value, method = "BH"),
    # Back-transform the log10 estimate into an adult/infant fold change.
    fold_change_adult_vs_infant = 10^estimate
  )
```

</details>

The model estimates the difference between adults and infants on the `log10(relative abundance)` scale.

| Quantity | Interpretation |
| --- | --- |
| `estimate > 0` | Higher abundance in adults |
| `estimate < 0` | Higher abundance in infants |
| `10^estimate` | Fold change on the original abundance scale |
| `q_value < 0.05` | Significant after FDR correction |

---

### 🔹 Step 3 — Compare model results to Wilcoxon results

The GLM and Wilcoxon approaches should often agree on the strongest taxa, but they are not testing exactly the same thing:

| Method | What it uses | Main strength | Main limitation |
| --- | --- | --- | --- |
| Wilcoxon | Rank order of abundances | Robust to skewness and outliers | Cannot include covariates |
| GLM | Log-transformed abundance values | Provides model estimates and standard errors | Still sensitive to model assumptions |

> 🧠 Discussion
> Which taxa are found by both methods? Which taxa appear in one method only? Method-specific hits deserve closer inspection because they may be driven by zeros, outliers, or distributional assumptions.

---

## ⚖️ 7. Adjust the Linear Model for Health Status

In observational data, variables can be correlated with each other. Here, `status` may be a confounder if healthy and sick samples are not evenly distributed between infants and adults.

The adjusted model estimates the age effect **while holding health status constant**.

### 🔹 Step 1 — Fit the adjusted model

```r
glm_age_adjusted_status_results <- dat_rel %>%
  filter(body_site == "gut") %>%
  mutate(
    age = factor(age, levels = c("infant", "adult")),
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    fit <- glm(
      log10(rel_abund + 1e-6) ~ age + status,
      data = .x,
      family = gaussian()
    )

    broom::tidy(fit)
  }) %>%
  ungroup() %>%
  filter(term == "ageadult") %>%
  mutate(
    q_value = p.adjust(p.value, method = "BH"),
    fold_change_adult_vs_infant_adjusted_for_status = 10^estimate
  )
```

<details>
<summary>Annotated code</summary>

```r
# Create a result table with one age-plus-status model per taxon.
glm_age_adjusted_status_results <- dat_rel %>%
  # Restrict the analysis to gut samples.
  filter(body_site == "gut") %>%
  mutate(
    # Set infant as the reference level for the age comparison.
    age = factor(age, levels = c("infant", "adult")),
    # Set healthy as the reference level for the health-status adjustment.
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  # Fit the same adjusted model separately for every taxon.
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    # Model log10-transformed relative abundance using age and health status.
    fit <- glm(
      log10(rel_abund + 1e-6) ~ age + status,
      data = .x,
      family = gaussian()
    )

    # Convert the model summary into a tidy data frame.
    broom::tidy(fit)
  }) %>%
  # Remove grouping so filtering and p-value adjustment use the full result table.
  ungroup() %>%
  # Keep only the age effect, because status is included here as an adjustment variable.
  filter(term == "ageadult") %>%
  mutate(
    # Correct age-effect p-values across taxa using Benjamini-Hochberg FDR.
    q_value = p.adjust(p.value, method = "BH"),
    # Back-transform the log10 estimate into an adjusted adult/infant fold change.
    fold_change_adult_vs_infant_adjusted_for_status = 10^estimate
  )
```

</details>

Only the `ageadult` term is retained because the practical focuses on the age effect. The model also estimates a `statussick` coefficient internally, but that coefficient is not the target here.

---

### 🔹 Step 2 — Interpret changes after adjustment

Compare the adjusted and unadjusted GLM results:

- taxa with stable estimates are robust to health-status adjustment,
- taxa whose effects shrink may be partly confounded by health status,
- taxa that change direction require careful interpretation,
- taxa that lose significance after adjustment should not be described as clean age effects.

> 💡 **Key concept: confounding**
> A confounder is a variable associated with both the predictor and the outcome. Adjusting for `status` helps separate age-associated microbial patterns from disease-associated microbial patterns.

---

## 🧾 8. Control False Discoveries with FDR

Differential abundance analysis usually tests many taxa. If 1,000 independent null hypotheses are tested at raw `p < 0.05`, about 50 false positives are expected by chance alone. This is why every result table in the practical reports a BH-adjusted `q_value`.

### 🔹 Step 1 — Distinguish FWER and FDR

| Error-control strategy | What it controls | Practical consequence |
| --- | --- | --- |
| FWER / Bonferroni | Probability of making any false positive | Very strict; often removes most microbiome signal |
| FDR / Benjamini-Hochberg | Expected proportion of false positives among discoveries | More appropriate for high-dimensional screening |

The practical uses **Benjamini-Hochberg FDR correction**:

```r
mutate(q_value = p.adjust(p_value, method = "BH"))
```

<details>
<summary>Annotated code</summary>

```r
# Add a q_value column by applying Benjamini-Hochberg correction to raw p-values.
mutate(q_value = p.adjust(p_value, method = "BH"))
```

</details>

or equivalently for model outputs:

```r
mutate(q_value = p.adjust(p.value, method = "BH"))
```

<details>
<summary>Annotated code</summary>

```r
# Model summaries from broom use p.value, so this corrects that column instead.
mutate(q_value = p.adjust(p.value, method = "BH"))
```

</details>

> 💡 **How to report significance**
> Use `q_value`, not the raw p-value, when deciding which taxa are significant. A threshold such as `q_value < 0.05` means that roughly 5% of the declared discoveries are expected to be false positives under the assumptions of the correction.

---

### 🔹 Step 2 — Interpret q-values together with effect sizes

A small q-value means the association is statistically reliable under a method's assumptions. It does not guarantee that the effect is biologically large or robust to modelling choices. This is why the practical always pairs q-values with effect-size plots such as volcano plots, ranked bar charts, and cross-method comparisons.

> ⚠️ **Practical rule**
> Report the direction, effect size, and q-value together. A taxon with `q < 0.05` but tiny effect size is less biologically compelling than a taxon with a large, consistent effect across methods.

---

## 🧪 9. Compositional Analysis with CLR Transformation

Relative abundances are **compositional** because all taxa in each sample sum to one. This creates dependence between taxa: if one taxon increases in relative abundance, at least one other taxon must decrease in relative terms, even if its absolute abundance did not change.

The **centred log-ratio (CLR)** transformation addresses this by expressing each taxon relative to the geometric mean of the whole community in that sample:

```text
CLR(x_i) = log(x_i / geometric_mean(x))
```

> 🧠 **The closure trap**
> Imagine three taxa where only taxon A truly increases in absolute abundance. After converting to relative abundance, taxa B and C can appear to decrease even if their absolute amounts did not change. Someone else's increase can look like your decrease. CLR is designed to reduce this artefact by analysing log-ratios rather than closed proportions.

> 💡 **Why the geometric mean?**
> The geometric mean is multiplicative, which is natural for compositions where ratios matter more than raw differences. On the log scale, `log(geometric_mean(x))` is the average of the log values. CLR therefore takes logs and subtracts the sample's log-scale centre.

| Mean | How to think about it | Why it matters here |
| --- | --- | --- |
| Arithmetic mean | Adds values and divides by count | Pulled upward by very abundant taxa |
| Median | Middle sorted value | Robust, but not tied to log-ratio geometry |
| Geometric mean | Multiplies values and takes the nth root | Natural centre for ratios and CLR |

### 🔹 Step 1 — Transform the phyloseq object

```r
ps_clr <- transform_sample_counts(ps, function(x) {
  x <- x + 1e-9
  log(x) - mean(log(x))
})

otu_clr <- as.data.frame(t(otu_table(ps_clr)))
meta_clr <- as.data.frame(sample_data(ps_clr))
tax_clr <- as.data.frame(tax_table(ps_clr))
```

<details>
<summary>Annotated code</summary>

```r
# Apply the centred log-ratio transformation to each sample in the phyloseq object.
ps_clr <- transform_sample_counts(ps, function(x) {
  # Add a tiny pseudo-count so zeros can be used in logarithms.
  x <- x + 1e-9
  # Subtract the sample's mean log abundance from each taxon's log abundance.
  log(x) - mean(log(x))
})

# Extract the CLR-transformed abundance table and transpose it so samples are rows.
otu_clr <- as.data.frame(t(otu_table(ps_clr)))

# Extract sample metadata for later joins and model design matrices.
meta_clr <- as.data.frame(sample_data(ps_clr))

# Extract taxonomic annotations for labelling CLR results.
tax_clr <- as.data.frame(tax_table(ps_clr))
```

</details>

CLR values are no longer constrained to the 0-1 range:

- positive CLR values mean a taxon is above the sample's geometric mean,
- negative CLR values mean a taxon is below the sample's geometric mean,
- values are comparable across taxa in log-ratio space.

> ⚠️ **Pseudo-counts and zeros**
> CLR uses logarithms, so zeros must be handled before transformation. The practical adds a very small pseudo-count (`1e-9`) before taking logs.

---

### 🔹 Step 2 — Reshape CLR values into long format

```r
dat_clr <- otu_clr %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "taxon", values_to = "clr") %>%
  left_join(meta_clr %>% as("data.frame") %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  left_join(tax_clr %>% rownames_to_column("taxon"), by = "taxon")
```

<details>
<summary>Annotated code</summary>

```r
# Start from the CLR abundance matrix, where rows are samples and columns are taxa.
dat_clr <- otu_clr %>%
  # Move sample IDs from row names into an explicit column.
  rownames_to_column("sample_id") %>%
  # Convert the wide CLR matrix into long format.
  pivot_longer(-sample_id, names_to = "taxon", values_to = "clr") %>%
  # Add sample metadata such as body_site, age, and status.
  left_join(meta_clr %>% as("data.frame") %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  # Add taxonomic labels for each taxon.
  left_join(tax_clr %>% rownames_to_column("taxon"), by = "taxon")
```

</details>

This mirrors the earlier `dat_rel` object, but the response variable is now `clr` instead of `rel_abund`.

---

### 🔹 Step 3 — Fit a CLR-based linear model

```r
clr_age_adjusted_status_results <- dat_clr %>%
  filter(body_site == "gut") %>%
  mutate(
    age = factor(age, levels = c("infant", "adult")),
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    fit <- lm(
      clr ~ age + status,
      data = .x
    )

    broom::tidy(fit)
  }) %>%
  ungroup() %>%
  filter(term == "ageadult") %>%
  mutate(q_value = p.adjust(p.value, method = "BH"))
```

<details>
<summary>Annotated code</summary>

```r
# Create a result table with one CLR-based adjusted linear model per taxon.
clr_age_adjusted_status_results <- dat_clr %>%
  # Restrict the analysis to gut samples.
  filter(body_site == "gut") %>%
  mutate(
    # Set infant as the reference age level.
    age = factor(age, levels = c("infant", "adult")),
    # Set healthy as the reference status level.
    status = factor(status, levels = c("healthy", "sick"))
  ) %>%
  # Fit the same CLR model separately for every taxon.
  group_by(taxon, genus, species) %>%
  group_modify(~ {
    # Model CLR abundance as a function of age and health status.
    fit <- lm(
      clr ~ age + status,
      data = .x
    )

    # Convert the model summary into a tidy data frame.
    broom::tidy(fit)
  }) %>%
  # Remove grouping before filtering and multiple-testing correction.
  ungroup() %>%
  # Keep the adult-vs-infant age coefficient.
  filter(term == "ageadult") %>%
  # Correct the age-effect p-values across taxa using Benjamini-Hochberg FDR.
  mutate(q_value = p.adjust(p.value, method = "BH"))
```

</details>

The `ageadult` estimate is now a difference in CLR values between adults and infants, adjusted for health status.

> 💡 **Why this is the most compositionally appropriate model in the practical**
> The CLR model no longer treats raw relative abundances as independent absolute measurements. Instead, it models each taxon relative to the whole-community geometric mean. Distances in CLR space correspond to Aitchison distance, the standard geometry for compositional data.

---

## 🔁 10. Compare Differential Abundance Methods

No single method gives the full truth. The practical therefore compares all four approaches:

1. Wilcoxon rank-sum test
2. GLM on log-transformed relative abundance
3. GLM adjusted for health status
4. CLR-based linear model adjusted for health status

### 🔹 Step 1 — Merge all result tables

```r
comparison_results <- wilcox_results %>%
  select(taxon, genus, species,
         wilcox_log2fc = log2_fc,
         wilcox_q = q_value) %>%
  left_join(
    glm_age_results %>%
      select(taxon, glm_estimate = estimate, glm_q = q_value),
    by = "taxon"
  ) %>%
  left_join(
    glm_age_adjusted_status_results %>%
      select(taxon, glm_adj_estimate = estimate, glm_adj_q = q_value),
    by = "taxon"
  ) %>%
  left_join(
    clr_age_adjusted_status_results %>%
      select(taxon, clr_estimate = estimate, clr_q = q_value),
    by = "taxon"
  )
```

<details>
<summary>Annotated code</summary>

```r
# Start the comparison table from the Wilcoxon results.
comparison_results <- wilcox_results %>%
  # Keep taxon labels plus Wilcoxon effect sizes and q-values.
  select(taxon, genus, species,
         wilcox_log2fc = log2_fc,
         wilcox_q = q_value) %>%
  left_join(
    # Add the age-only GLM estimate and q-value.
    glm_age_results %>%
      select(taxon, glm_estimate = estimate, glm_q = q_value),
    by = "taxon"
  ) %>%
  left_join(
    # Add the status-adjusted GLM estimate and q-value.
    glm_age_adjusted_status_results %>%
      select(taxon, glm_adj_estimate = estimate, glm_adj_q = q_value),
    by = "taxon"
  ) %>%
  left_join(
    # Add the status-adjusted CLR model estimate and q-value.
    clr_age_adjusted_status_results %>%
      select(taxon, clr_estimate = estimate, clr_q = q_value),
    by = "taxon"
  )
```

</details>

This creates one comparison table where each taxon has effect sizes and q-values from every method.

---

### 🔹 Step 2 — Count significant taxa per method

```r
sig_summary <- tibble(
  method = c("Wilcoxon", "GLM", "GLM adjusted", "CLR"),
  n_significant = c(
    sum(wilcox_results$q_value < 0.05, na.rm = TRUE),
    sum(glm_age_results$q_value < 0.05, na.rm = TRUE),
    sum(glm_age_adjusted_status_results$q_value < 0.05, na.rm = TRUE),
    sum(clr_age_adjusted_status_results$q_value < 0.05, na.rm = TRUE)
  )
)

sig_summary
```

<details>
<summary>Annotated code</summary>

```r
# Build a small summary table with one row per method.
sig_summary <- tibble(
  # Name the four differential-abundance approaches compared in the practical.
  method = c("Wilcoxon", "GLM", "GLM adjusted", "CLR"),
  # Count how many taxa pass the q-value threshold in each method.
  n_significant = c(
    sum(wilcox_results$q_value < 0.05, na.rm = TRUE),
    sum(glm_age_results$q_value < 0.05, na.rm = TRUE),
    sum(glm_age_adjusted_status_results$q_value < 0.05, na.rm = TRUE),
    sum(clr_age_adjusted_status_results$q_value < 0.05, na.rm = TRUE)
  )
)

# Print the summary table.
sig_summary
```

</details>

Methods that detect many taxa are not automatically better. They may be more sensitive, but they may also be more prone to false positives. Conservative methods may miss real biology but provide stronger confidence in the taxa they do detect.

---

### 🔹 Step 3 — Identify taxa significant in all four methods

```r
all_methods_sig <- comparison_results %>%
  filter(wilcox_q < 0.05, glm_q < 0.05, glm_adj_q < 0.05, clr_q < 0.05) %>%
  select(genus, species, wilcox_q, glm_q, glm_adj_q, clr_q) %>%
  arrange(wilcox_q)

all_methods_sig
```

<details>
<summary>Annotated code</summary>

```r
# Start from the merged method-comparison table.
all_methods_sig <- comparison_results %>%
  # Keep only taxa that are significant after FDR correction in all four methods.
  filter(wilcox_q < 0.05, glm_q < 0.05, glm_adj_q < 0.05, clr_q < 0.05) %>%
  # Keep readable taxon labels and the q-values from each method.
  select(genus, species, wilcox_q, glm_q, glm_adj_q, clr_q) %>%
  # Sort by the Wilcoxon q-value as one simple way to order the table.
  arrange(wilcox_q)

# Print the taxa supported by all four methods.
all_methods_sig
```

</details>

Taxa significant across all four methods are the strongest candidates for robust age-associated microbes in this dataset.

> 🧠 Discussion
> Do the taxa significant in all methods make biological sense? Are they enriched in the expected age group? Do they also appear among the strongest RDA loadings later?

---

### 🔹 Step 4 — Compare effect sizes and q-values visually

The Rmd uses several plots to compare methods:

| Plot | Question it answers |
| --- | --- |
| GLM vs adjusted GLM scatter plot | Does health-status adjustment change the age effect? |
| Adjusted GLM vs CLR scatter plot | Does compositional correction change the inferred effect? |
| Significant-method count bar chart | How many methods support each taxon? |
| Directional heatmap | Are taxa consistently adult- or infant-enriched across methods? |
| Wilcoxon vs CLR scatter plot | Do the most different methods agree on effect direction? |

> 💡 **Most robust interpretation strategy**
> Prioritise taxa with consistent direction, low q-values, and agreement across multiple methods. Treat taxa detected by only one method as hypotheses that require follow-up, not as final conclusions.

---

## 🧭 11. Constrained Ordination with RDA

The univariate tests above analyse one taxon at a time. Microbiomes are communities, so the practical ends with a multivariate method: **redundancy analysis (RDA)**.

RDA asks whether the whole community composition is structured by explanatory variables such as `age` and `status`.

> 💡 **Alpha vs beta/community-level analyses**
> Alpha diversity asks how rich or even one sample is and returns one value per sample. Beta diversity asks how different samples are from each other and works with distances or ordination. RDA belongs to the second family: it tests whether whole-community structure aligns with metadata.

### 🔹 Step 1 — Build the CLR community matrix

```r
gut_clr_mat <- dat_clr %>%
  filter(body_site == "gut") %>%
  select(sample_id, taxon, clr) %>%
  pivot_wider(names_from = taxon, values_from = clr) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

gut_meta <- meta_clr[rownames(gut_clr_mat), ] %>% as("data.frame")
```

<details>
<summary>Annotated code</summary>

```r
# Build the community matrix used as the RDA response.
gut_clr_mat <- dat_clr %>%
  # RDA is focused on gut samples in this practical.
  filter(body_site == "gut") %>%
  # Keep only the sample identifier, taxon identifier, and CLR abundance.
  select(sample_id, taxon, clr) %>%
  # Convert from long format to a matrix-like wide format: samples x taxa.
  pivot_wider(names_from = taxon, values_from = clr) %>%
  # Move sample IDs into row names because vegan ordination functions expect row names.
  column_to_rownames("sample_id") %>%
  # Convert the data frame to a numeric matrix for RDA.
  as.matrix()

# Reorder the metadata so its rows match the sample order in the community matrix.
gut_meta <- meta_clr[rownames(gut_clr_mat), ] %>% as("data.frame")
```

</details>

RDA needs:

- a response matrix with samples as rows and taxa as columns,
- a metadata table with the same samples in the same order.

---

### 🔹 Step 2 — Fit the constrained ordination

```r
rda_age <- rda(gut_clr_mat ~ age + status, data = gut_meta)
```

<details>
<summary>Annotated code</summary>

```r
# Fit an RDA model that explains CLR community composition using age and status.
rda_age <- rda(gut_clr_mat ~ age + status, data = gut_meta)
```

</details>

This model constrains the ordination axes to combinations of `age` and `status`. In other words, the plot is not just describing total community variation; it is specifically showing the part of the variation explained by the predictors in the formula.

> 🧠 **PCA vs RDA**
> PCA finds the directions of greatest overall variance. RDA finds the directions of community variance that are best explained by the predictors in the model formula.

---

### 🔹 Step 3 — Test the model with permutation ANOVA

```r
anova(rda_age, permutations = 999)
anova(rda_age, by = "term", permutations = 999)
```

<details>
<summary>Annotated code</summary>

```r
# Test whether the full constrained RDA model explains more variation than expected by chance.
anova(rda_age, permutations = 999)

# Test each explanatory variable separately to compare age and status effects.
anova(rda_age, by = "term", permutations = 999)
```

</details>

The first test asks whether the overall constrained model explains more community variation than expected by chance. The second test evaluates predictors term by term, allowing you to compare the contribution of `age` and `status`.

> 💡 **Why permutations?**
> Ordination tests often use permutations because they do not rely on simple parametric assumptions. Sample labels are shuffled repeatedly to build a null distribution for the explained variance.

---

### 🔹 Step 4 — Extract taxa driving the primary RDA axis

```r
taxon_scores <- scores(rda_age, display = "species")[, 1]

rda_taxa <- tibble(
  taxon = names(taxon_scores),
  rda1_loading = taxon_scores
) %>%
  arrange(desc(abs(rda1_loading)))
```

<details>
<summary>Annotated code</summary>

```r
# Extract taxon scores from the first RDA axis.
taxon_scores <- scores(rda_age, display = "species")[, 1]

# Turn the named score vector into a tidy table.
rda_taxa <- tibble(
  # Keep the taxon identifier from the score names.
  taxon = names(taxon_scores),
  # Store each taxon's loading on RDA1.
  rda1_loading = taxon_scores
) %>%
  # Sort taxa by the size of their loading, regardless of direction.
  arrange(desc(abs(rda1_loading)))
```

</details>

Taxa with large absolute RDA1 loadings contribute most strongly to the primary constrained axis. Positive and negative loadings point toward opposite ends of the community gradient.

---

### 🔹 Step 5 — Interpret the RDA biplot

The Rmd builds a biplot with:

- sample points coloured by `age`,
- sample shapes indicating `status`,
- black arrows for explanatory variables,
- grey arrows for the highest-loading taxa.

In the biplot:

| Visual element | Interpretation |
| --- | --- |
| Distance between sample points | Difference in community composition |
| Direction of the `age` arrow | Main axis of age-associated community change |
| Direction of the `status` arrow | Main axis of health-status-associated community change |
| Taxon arrows | Taxa contributing strongly to the constrained axes |
| Clustering by colour | Evidence that age structures the gut microbiome |

> 🧠 Discussion
> Do the top RDA taxa overlap with the taxa significant in all four univariate methods? Agreement between the RDA loadings and the differential abundance tests gives stronger evidence that those taxa are genuine drivers of the infant-adult gut microbiome difference.

---

## 📊 12. How to Read the Practical Results

The practical deliberately shows several methods because differential abundance analysis is sensitive to assumptions. Use the following hierarchy when interpreting results:

| Evidence level | Interpretation |
| --- | --- |
| Significant in all four methods, same direction | Strong candidate age-associated taxon |
| Significant in adjusted GLM and CLR, same direction | Good candidate, especially if biologically plausible |
| Significant before but not after status adjustment | Possible confounding by health status |
| Significant only in Wilcoxon | May reflect rank differences, zero patterns, or outliers |
| Significant only in CLR | May reflect compositional effects hidden in relative abundance |
| Opposite directions across methods | Do not interpret without deeper inspection |

> ⚠️ **Avoid over-interpreting p-values**
> A low q-value indicates statistical evidence under a specific model. It does not by itself prove biological importance. Always consider effect size, direction, distribution shape, and consistency across methods.

---

## 🧠 Discussion Questions

Use these questions while going through the Rmd:

1. Why is it useful to visualise relative abundance before running statistical tests?
2. Which feature of the data is most problematic for a standard t-test: zeros, skewness, non-normality, or compositionality?
3. Do the Wilcoxon and GLM methods identify the same top infant- and adult-enriched taxa?
4. Which taxa change most after adjusting for health status?
5. Does the CLR model change the direction or significance of any major taxa?
6. Which taxa are significant across all four methods?
7. Do the strongest RDA loadings support the same biological story as the univariate tests?
8. Which results would you trust most, and why?

---

## 📤 Suggested Outputs to Keep

By the end of the practical, keep a record of:

- the rendered HTML version of `differential_abundance.Rmd`,
- the `sig_summary` table,
- the `all_methods_sig` table,
- the top infant- and adult-enriched taxa from each method,
- the method-comparison heatmap,
- the RDA biplot and top RDA loading bar chart,
- a short written interpretation of which taxa are most robustly associated with infant vs adult gut microbiomes.

> 📂 **Suggested interpretation note**
> Write 5-10 lines summarising the biological conclusion. Focus on taxa that are consistent across methods, explain whether health-status adjustment changed the result, and state whether the CLR/RDA analyses support the same conclusion.

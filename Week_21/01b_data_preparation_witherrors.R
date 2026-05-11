#' 01b_Data_Preparation_WithErrors.R
#'
#' This optional teaching script intentionally demonstrates incorrect data
#' preparation patterns. It should not be used to feed the main ML pipeline.
#'
#' The goal is to help students recognize characteristic signatures of:
#' 1. Fitting learned preprocessing on held-out data
#' 2. Accidentally creating outcome/label leakage features
#' 3. Coercing the outcome column during feature conversion
#' 4. Breaking sample-label alignment
#'
#' All outputs are prefixed with `witherrors_` so they do not overwrite the
#' leakage-safe pipeline artifacts.

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
library(ggplot2)
library(tidyr)
library(caret)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "01b_data_preparation_witherrors.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting intentional data preparation error demo...\n")
cat("These artifacts are for instruction only and should not feed the main pipeline.\n\n")

# Start from the same structurally clean data used by the corrected pipeline.
prepared_data <- raw_data %>%
  clean_names() %>%
  mutate(across(-dx, as.numeric))

if (!"dx" %in% names(prepared_data)) {
  stop("Error: expected outcome column `dx` was not found after cleaning names.")
}

feature_cols <- setdiff(names(prepared_data), "dx")
outcome_binary <- if_else(prepared_data$dx == "cancer", 1, 0)

# Use an explicit split for diagnostics so students can compare train-only
# preprocessing against incorrect full-data and test-data preprocessing.
train_idx <- caret::createDataPartition(prepared_data$dx, p = 0.8, list = FALSE) %>%
  as.vector()
test_idx <- setdiff(seq_len(nrow(prepared_data)), train_idx)
train_raw <- prepared_data[train_idx, ]
test_raw <- prepared_data[test_idx, ]

scale_using <- function(df, centers, scales) {
  safe_scales <- scales
  safe_scales[is.na(safe_scales) | safe_scales == 0] <- 1
  scaled <- sweep(as.matrix(df[, names(centers), drop = FALSE]), 2, centers, "-")
  scaled <- sweep(scaled, 2, safe_scales, "/")
  as.data.frame(scaled)
}

summarize_scaled_test <- function(scaled_df, method, interpretation) {
  tibble(
    method = method,
    feature = names(scaled_df),
    test_mean = vapply(scaled_df, mean, numeric(1), na.rm = TRUE),
    test_sd = vapply(scaled_df, sd, numeric(1), na.rm = TRUE),
    interpretation = interpretation
  ) %>%
    mutate(
      abs_test_mean = abs(test_mean),
      abs_test_sd_minus_1 = abs(test_sd - 1)
    )
}

feature_cor_scan <- function(df, label) {
  y <- if_else(df$dx == "cancer", 1, 0)
  df %>%
    select(-dx) %>%
    summarise(across(everything(), ~ abs(cor(.x, y, use = "complete.obs")))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "abs_cor_with_dx") %>%
    mutate(dataset = label) %>%
    arrange(dataset, desc(abs_cor_with_dx))
}

cat("Case 1: Learned preprocessing fit on the wrong data boundary...\n")

train_centers <- vapply(train_raw[feature_cols], mean, numeric(1), na.rm = TRUE)
train_scales <- vapply(train_raw[feature_cols], sd, numeric(1), na.rm = TRUE)
full_centers <- vapply(prepared_data[feature_cols], mean, numeric(1), na.rm = TRUE)
full_scales <- vapply(prepared_data[feature_cols], sd, numeric(1), na.rm = TRUE)
test_centers <- vapply(test_raw[feature_cols], mean, numeric(1), na.rm = TRUE)
test_scales <- vapply(test_raw[feature_cols], sd, numeric(1), na.rm = TRUE)

correct_test_scaled <- scale_using(test_raw, train_centers, train_scales)
bad_full_data_scaled <- scale_using(test_raw, full_centers, full_scales)
bad_test_fit_scaled <- scale_using(test_raw, test_centers, test_scales)

scaling_diagnostics <- bind_rows(
  summarize_scaled_test(
    correct_test_scaled,
    "correct_train_fit",
    "Expected: held-out test data is transformed with training parameters."
  ),
  summarize_scaled_test(
    bad_full_data_scaled,
    "error_full_dataset_fit",
    "Leakage: preprocessing parameters used held-out rows."
  ),
  summarize_scaled_test(
    bad_test_fit_scaled,
    "error_test_fit_separately",
    "Leakage signature: test means become exactly 0 and SDs exactly 1."
  )
)

scaling_param_drift <- tibble(
  feature = feature_cols,
  train_mean = train_centers,
  full_data_mean = full_centers,
  train_sd = train_scales,
  full_data_sd = full_scales
) %>%
  mutate(
    mean_shift_from_using_full_data = full_data_mean - train_mean,
    sd_shift_from_using_full_data = full_data_sd - train_sd
  )

write.csv(
  scaling_diagnostics,
  file.path(results_dir, "witherrors_scaling_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  scaling_param_drift,
  file.path(results_dir, "witherrors_preprocessing_parameter_drift.csv"),
  row.names = FALSE
)

p_scaling <- scaling_diagnostics %>%
  select(method, feature, abs_test_mean, abs_test_sd_minus_1) %>%
  pivot_longer(
    cols = c(abs_test_mean, abs_test_sd_minus_1),
    names_to = "diagnostic",
    values_to = "value"
  ) %>%
  ggplot(aes(x = reorder(feature, value), y = value, fill = method)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~diagnostic, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Preprocessing Boundary Error Signatures",
    subtitle = "A test set with exactly zero means and unit SDs was fit separately.",
    x = "Feature",
    y = "Absolute diagnostic value",
    fill = "Preprocessing method"
  )

ggsave(
  file.path(figures_dir, "witherrors_preprocessing_signatures.png"),
  p_scaling,
  width = 10,
  height = 6
)

cat("Case 2: Outcome leakage feature accidentally added...\n")

bad_leakage_data <- prepared_data %>%
  mutate(
    dx_leak_numeric = if_else(dx == "cancer", 1, 0),
    dx_leak_almost_numeric = dx_leak_numeric + rnorm(n(), mean = 0, sd = 0.001)
  )

leakage_scan <- feature_cor_scan(bad_leakage_data, "outcome_leakage_error") %>%
  mutate(
    suspicious = abs_cor_with_dx >= 0.95,
    interpretation = if_else(
      suspicious,
      "Near-perfect association with outcome; inspect for label leakage.",
      "No near-perfect outcome association."
    )
  )

write.csv(
  leakage_scan,
  file.path(results_dir, "witherrors_outcome_leakage_scan.csv"),
  row.names = FALSE
)

p_leakage <- leakage_scan %>%
  ggplot(aes(x = reorder(feature, abs_cor_with_dx), y = abs_cor_with_dx, fill = suspicious)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  labs(
    title = "Outcome Leakage Signature",
    subtitle = "Features near correlation 1.0 with the label are usually not biology.",
    x = "Feature",
    y = "Absolute correlation with dx",
    fill = "Suspicious"
  )

ggsave(
  file.path(figures_dir, "witherrors_outcome_leakage_signature.png"),
  p_leakage,
  width = 8,
  height = 6
)

cat("Case 3: Outcome column accidentally coerced during numeric feature conversion...\n")

bad_outcome_conversion <- suppressWarnings(
  prepared_data %>%
    mutate(across(everything(), as.numeric))
)

outcome_conversion_diagnostics <- tibble(
  version = c("correct_prepared_data", "error_across_everything_as_numeric"),
  dx_missing_count = c(
    sum(is.na(prepared_data$dx)),
    sum(is.na(bad_outcome_conversion$dx))
  ),
  dx_unique_non_missing = c(
    n_distinct(prepared_data$dx, na.rm = TRUE),
    n_distinct(bad_outcome_conversion$dx, na.rm = TRUE)
  ),
  interpretation = c(
    "Outcome labels are intact.",
    "Outcome labels became missing; across(everything(), as.numeric) touched dx."
  )
)

write.csv(
  outcome_conversion_diagnostics,
  file.path(results_dir, "witherrors_outcome_conversion_diagnostics.csv"),
  row.names = FALSE
)

p_conversion <- ggplot(
  outcome_conversion_diagnostics,
  aes(x = version, y = dx_missing_count, fill = version)
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Outcome Conversion Error Signature",
    subtitle = "If dx becomes missing during preparation, feature conversion touched the label.",
    x = "Dataset version",
    y = "Missing dx labels"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(figures_dir, "witherrors_outcome_conversion_signature.png"),
  p_conversion,
  width = 8,
  height = 4
)

cat("Case 4: Sample-label alignment broken by shuffling labels...\n")

bad_label_alignment_data <- prepared_data %>%
  mutate(dx = sample(dx))

alignment_scan <- bind_rows(
  feature_cor_scan(prepared_data, "correct_labels"),
  feature_cor_scan(bad_label_alignment_data, "shuffled_label_error")
) %>%
  group_by(dataset) %>%
  mutate(rank = row_number()) %>%
  ungroup()

alignment_summary <- alignment_scan %>%
  group_by(dataset) %>%
  summarise(
    max_abs_cor_with_dx = max(abs_cor_with_dx, na.rm = TRUE),
    mean_abs_cor_with_dx = mean(abs_cor_with_dx, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = if_else(
      dataset == "shuffled_label_error",
      "Signal collapses because labels no longer match samples.",
      "Baseline signal scan before intentional label shuffling."
    )
  )

write.csv(
  alignment_scan,
  file.path(results_dir, "witherrors_label_alignment_scan.csv"),
  row.names = FALSE
)
write.csv(
  alignment_summary,
  file.path(results_dir, "witherrors_label_alignment_summary.csv"),
  row.names = FALSE
)

p_alignment <- alignment_scan %>%
  filter(rank <= 10) %>%
  ggplot(aes(x = reorder(feature, abs_cor_with_dx), y = abs_cor_with_dx, fill = dataset)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~dataset) +
  theme_bw() +
  labs(
    title = "Sample-Label Alignment Error Signature",
    subtitle = "Shuffling labels erases real feature-label structure.",
    x = "Top features within each scan",
    y = "Absolute correlation with dx",
    fill = "Dataset"
  )

ggsave(
  file.path(figures_dir, "witherrors_label_alignment_signature.png"),
  p_alignment,
  width = 10,
  height = 6
)

bad_examples <- list(
  notes = c(
    "These datasets are intentionally wrong.",
    "Use them only to demonstrate diagnostics; do not feed them to scripts 02-05."
  ),
  split_indices = list(train_idx = train_idx, test_idx = test_idx),
  bad_full_data_scaled_test = bad_full_data_scaled,
  bad_test_fit_scaled = bad_test_fit_scaled,
  bad_leakage_data = bad_leakage_data,
  bad_outcome_conversion = bad_outcome_conversion,
  bad_label_alignment_data = bad_label_alignment_data
)

saveRDS(bad_examples, file.path(results_dir, "witherrors_bad_examples.rds"))

cat("\n--- Error Demo Summary ---\n")
cat("Scaling diagnostics saved to results/witherrors_scaling_diagnostics.csv\n")
cat("Outcome leakage scan saved to results/witherrors_outcome_leakage_scan.csv\n")
cat("Outcome conversion diagnostics saved to results/witherrors_outcome_conversion_diagnostics.csv\n")
cat("Label alignment scan saved to results/witherrors_label_alignment_scan.csv\n")
cat("Figures saved to figures/witherrors_*.png\n")
cat("Intentional data preparation error demo complete.\n")

sink()

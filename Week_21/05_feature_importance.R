#' 05_Feature_Importance.R
#' 
#' This script identifies which features (OTUs) are most influential for model predictions.
#' It covers the fifth phase of the lab practical:
#' Biological Interpretation and Feature Importance.
#'
#' Steps:
#' 1. Manual Permutation Feature Importance (Custom Implementation)
#' 2. Visualizing Top Features (Decrease in Performance)
#' 3. Linear Model Interpretation (Directionality via Coefficients)

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
library(caret)

# 0. Setup Logging and Seed
log_file <- file.path(logs_dir, "05_feature_importance.log")
sink(log_file, append = FALSE, split = TRUE)
set.seed(seed)

cat("Starting Feature Importance Phase (Manual Implementation)...\n")

# 1. Load Tuned Models
input_file <- file.path(results_dir, "tuned_model_results.rds")
tuned_models <- readRDS(input_file)

# 2. Manual Permutation Importance Function
# Novice Tip: We create our own importance calculator to see exactly how it works!
# We shuffle one column at a time and see how much the AUC drops.
get_retained_predictors <- function(trained_model, test_data, outcome_col) {
  all_predictors <- setdiff(colnames(test_data), outcome_col)
  
  if (!is.null(trained_model$trainingData)) {
    retained_predictors <- setdiff(colnames(trained_model$trainingData), ".outcome")
    return(intersect(all_predictors, retained_predictors))
  }
  
  all_predictors
}

calc_manual_importance <- function(model_obj, test_data, outcome_col, n_perms = 5) {
  
  # Calculate Baseline AUC
  actual_probs <- predict(model_obj$trained_model, test_data, type = "prob")
  baseline_perf <- twoClassSummary(data.frame(obs = factor(test_data[[outcome_col]]), 
                                              pred = predict(model_obj$trained_model, test_data),
                                              actual_probs), 
                                   lev = levels(factor(test_data[[outcome_col]])))
  baseline_auc <- baseline_perf["ROC"]
  
  features <- get_retained_predictors(model_obj$trained_model, test_data, outcome_col)
  imp_results <- data.frame()
  
  for (f in features) {
    perm_aucs <- numeric(n_perms)
    for (i in 1:n_perms) {
      perm_data <- test_data
      perm_data[[f]] <- sample(perm_data[[f]]) # SHUFFLE
      
      perm_probs <- predict(model_obj$trained_model, perm_data, type = "prob")
      perm_perf <- twoClassSummary(data.frame(obs = factor(test_data[[outcome_col]]), 
                                              pred = predict(model_obj$trained_model, perm_data),
                                              perm_probs), 
                                   lev = levels(factor(test_data[[outcome_col]])))
      perm_aucs[i] <- perm_perf["ROC"]
    }
    avg_perm_auc <- mean(perm_aucs)
    imp_results <- rbind(imp_results, data.frame(feat = f, perf_metric_diff = baseline_auc - avg_perm_auc))
  }
  return(imp_results)
}

cat("\nStep 1: Calculating permutation importance (Manual shuffle)...\n")
importance_rf <- calc_manual_importance(tuned_models$rf, tuned_models$rf$test_data, "dx")
importance_glmnet <- calc_manual_importance(tuned_models$glmnet, tuned_models$glmnet$test_data, "dx")

# 3. Visualize Importance
cat("\nStep 2: Visualizing top features...\n")
combined_importance <- bind_rows(
  importance_rf %>% mutate(method = "rf"),
  importance_glmnet %>% mutate(method = "glmnet")
)

p_imp <- combined_importance %>%
  group_by(method) %>%
  slice_max(order_by = perf_metric_diff, n = 10) %>%
  ggplot(aes(x = reorder(feat, perf_metric_diff), y = perf_metric_diff, fill = method)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~method, scales = "free") +
  theme_bw() +
  labs(title = "Top Important Features (Manual Calculation)",
       subtitle = "Decrease in AUC when feature is shuffled",
       x = "Feature (OTU)", y = "Decrease in AUC")

ggsave(file.path(figures_dir, "feature_importance.png"), p_imp, width = 10, height = 6)

# 4. Directionality (GLMNET Coefficients)
cat("\nStep 3: Extracting directionality from GLMNET...\n")
final_model <- tuned_models$glmnet$trained_model$finalModel
best_lambda <- tuned_models$glmnet$trained_model$bestTune$lambda
coefs <- as.matrix(coef(final_model, s = best_lambda))
coef_df <- data.frame(feat = rownames(coefs), weight = coefs[,1]) %>%
  filter(feat != "(Intercept)", weight != 0)

p_coef <- coef_df %>%
  ggplot(aes(x = reorder(feat, weight), y = weight, fill = weight > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), 
                    labels = c("TRUE" = "Associated with Normal", "FALSE" = "Associated with Cancer")) +
  theme_bw() +
  labs(title = "GLMNET Feature Weights (Directionality)",
       x = "Feature (OTU)", y = "Model Coefficient (Weight)",
       fill = "Direction")

ggsave(file.path(figures_dir, "feature_weights_direction.png"), p_coef, width = 8, height = 6)

saveRDS(combined_importance, file.path(results_dir, "feature_importance.rds"))
cat("\nSuccess: Manual importance figures saved.\n")

sink()

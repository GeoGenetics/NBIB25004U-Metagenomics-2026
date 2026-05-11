# Lab Practical: Machine Learning for Microbial Metagenomics
**Objective**: Learn to build, tune, and interpret machine learning models for microbiome-based classification using the `mikropml` package.

This practical is designed to live inside `Week_21/` in a shared course repo.

---

## Introduction
In this practical, we will predict host health (Diagnosis: `dx`) based on the relative abundance of Operational Taxonomic Units (OTUs). We will follow a rigorous pipeline designed to avoid overfitting and data leakage.

### Dataset
We are using a reduced version of the colorectal cancer dataset from Topçuoğlu et al. (2020).
- **Outcome (`dx`)**: `cancer` or `normal`.
- **Features**: Abundance of 10 specific OTUs.

---

## Phase 1: Data Preparation
**Script**: `01_data_preparation.R`
**Optional Diagnostic Script**: `01b_data_preparation_witherrors.R`

Before any ML begins, the data must be structurally clean. 
- **Cleaning Names**: Ensure R can read headers without error.
- **Type Checks**: Ensure OTU features are numeric and the outcome labels are present.
- **Leakage Boundary**: Do not scale, filter low-variance features, or collapse correlated features here. Those learned transformations are estimated inside model training.

**Key Concept**: *Data Leakage*. Why do we estimate scaling, variance filtering, and correlation filtering from the training data only? (Hint: To avoid the model "seeing" the test data too early).

Run `01b_data_preparation_witherrors.R` when you want to inspect characteristic failure signatures, including test-set preprocessing leakage, outcome leakage features, broken outcome conversion, and shuffled sample labels. It writes only `witherrors_*` outputs and should not feed the main pipeline.

---

## Phase 2: Model Training & Experimental Design
**Script**: `02_model_training.R`

We train two types of models:
1. **GLMNET**: Logistic Regression (Baseline). Fast and interpretable.
2. **RF**: Random Forest (Non-linear). Robust ensemble of trees.

**Experimental Setup**:
- **80/20 Split**: 80% for training, 20% for testing.
- **k-Fold Cross-Validation**: The model iterates through 5 "folds" of training data to find the best settings internally.
- **Preprocessing Inside Training**: `caret::train()` estimates near-zero variance filtering, correlation filtering, centering, and scaling within the training/resampling workflow, then applies those fitted transformations to held-out samples.

---

## Phase 3: Model Evaluation
**Script**: `03_model_evaluation.R`

How do we know if our model is good? We look at the **AUC (Area Under the Curve)**.
- **AUC = 0.5**: No better than a coin flip.
- **AUC = 1.0**: Perfect prediction.

**Expected Results**:
You should see that **Random Forest** generally outperforms Logistic Regression on this biological dataset.
- *GLMNET AUC*: ~0.63
- *RF AUC*: ~0.70

---

## Phase 4: Hyperparameter Tuning
**Script**: `04_hyperparameter_tuning.R`

Hyperparameters are the "knobs" of the model. 
- **RF**: We tune `mtry` (how many features are checked at each split).
- **GLMNET**: We tune `lambda` (the strength of the penalty).

**Results**:
Notice how the RF performance increases when we search for the optimal `mtry`.
- *Tuned RF AUC*: ~0.72 (vs 0.70 default).

---

## Phase 5: Biological Interpretation
**Script**: `05_feature_importance.R`

ML isn't just about prediction; it's about discovery.
- **Permutation Importance**: We shuffle an OTU's values. If the AUC drops significantly, that OTU was important for the prediction.
- **Directionality**: Using GLMNET coefficients, we can see if an OTU is associated with `cancer` (Positive weight) or `normal` (Negative weight).

**Challenge**: Look at `feature_importance.png`. Which OTU is the most important for the Random Forest model?

---

## How to Run
Change into the `Week_21` folder, then run the scripts in order:
```bash
cd Week_21
Rscript 01_data_preparation.R
# Optional diagnostic lab, not part of the main pipeline:
Rscript 01b_data_preparation_witherrors.R
Rscript 02_model_training.R
Rscript 03_model_evaluation.R
Rscript 04_hyperparameter_tuning.R
Rscript 05_feature_importance.R
```

Check the `results/` folder for data and `figures/` for plots!

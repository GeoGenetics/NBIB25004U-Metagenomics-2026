# Lab Practical Plan: Machine Learning for Bioinformatics

This practical is intended to live in `Week_21/` within the shared course repository.

## 1. Data Preparation and Feature Engineering
*This stage separates global data hygiene from learned transformations, so the model can interpret the data without leaking information from held-out samples.*

1. **Cleaning Names and Columns**: Remove special characters, spaces, and non-ASCII characters from headers. This prevents syntax errors in downstream R/Python functions that treat column names as variables.
2. **Standardizing Data Types**: Ensure numeric values are not stored as strings (characters). Machine learning models require quantitative inputs; a single "N/A" in a column can force the entire vector to become a character type, breaking the model.
3. **Assessing Data Completeness (Imputation)**:
    * **Median/Mean Imputation**: Replacing missing values with the average of the column.
    * **Zero/NA Replacement**: Deciding if an "NA" in bioinformatics (e.g., proteomics or transcriptomics) represents a true zero (below detection limit) or a missing observation.
4. **Redundant and Correlated Feature Removal**: Removing features that are highly correlated ($r > 0.90$). If two genes always vary together, keeping both provides no new info and can confuse model weights (Multicollinearity). This should be estimated inside the training/resampling workflow, not on the full dataset before splitting.
5. **Target Variable Curation**: Remove any observations (rows) that are missing the "Outcome" (Labels). A model cannot learn from a sample if it doesn't know the truth it is supposed to predict.
6. **Feature Variance Filtering**:
    * **Zero Variance**: Removing columns that are the same value for every sample.
    * **Near-Zero Variance**: Removing features with very few unique values or extremely skewed distributions, as they offer no predictive power. Like correlation filtering, this is a learned training step.
7. **Data Normalization and Scaling**:
    * **Z-score Scaling**: Centering data so the mean is 0 and standard deviation is 1.
    * **Min-Max Scaling**: Squishing data between 0 and 1. This prevents features with large raw numbers (e.g., 5000 units) from dominating features with small numbers (e.g., 0.05 units). Scaling parameters must come from the training data only.
8. **Categorical and Binary Data Handling**:
    * **One-Hot Encoding**: Converting "Wildtype/Mutant" into numeric 0 and 1.
    * **Dummy Variables**: Ensuring multi-level categories (e.g., Tissue Type: Lung, Heart, Liver) are represented numerically.

## 2. Model Training and Experimental Design
*This stage defines how we teach the model while ensuring we don't "overfit" (memorize the data).*

1. **Data Partitioning Strategy**:
    * **Training/Validation Split**: Holding out a portion (e.g., 20%) of data that the model never sees during training to test its real-world performance.
    * **k-Fold Cross-Validation**: Dividing data into $k$ groups (folds). The model trains on $k-1$ folds and tests on the remaining one, repeating this until every fold has been the "test" set. This ensures the results aren't a fluke of a lucky split.
    * **Exemplar Data/Downsampling**: Balancing classes if you have 1000 healthy samples but only 10 sick ones, ensuring the model doesn't just learn to guess "healthy" every time.
2. **Model Execution**:
    * **Algorithm Selection**: Choosing between Random Forest, XGBoost, or Support Vector Machines.
    * **Parallelization**: Using multiple CPU cores (via `mclapply` in R or `n_jobs` in Python) to speed up training on large genomic datasets.
3. **Inspecting Model Outputs**:
    * **Convergence**: Checking if the model actually finished learning.
    * **Probability Scores**: Looking at the raw confidence (e.g., "This sample is 85% likely to be Mutant") rather than just the final 0/1 label.

## 3. Model Evaluation and Performance Metrics
*Evaluating if the model's predictions are reliable and useful for biological discovery.*

1. **Classification Performance**:
    * **ROC/AUC**: Area Under the Receiver Operating Characteristic curve. A score of 1.0 is perfect; 0.5 is as good as a coin flip. It measures performance across all classification thresholds.
    * **prAUC (Precision-Recall AUC)**: Often better for bioinformatics when classes are imbalanced (e.g., a rare disease).
2. **Clinical and Biological Utility**:
    * **Sensitivity (Recall)**: The ability to find all the "sick" samples. Crucial for medical screening.
    * **Specificity**: The ability to correctly identify "healthy" samples. Crucial for avoiding false alarms.
    * **F1-Score**: The harmonic mean of Precision and Recall. High F1 means the model is balanced and not biased toward one error type.
    * **Confusion Matrix**: A table showing True Positives, False Positives, True Negatives, and False Negatives.

## 4. Hyperparameter Tuning
*The "fine-tuning" phase to squeeze the best performance out of an algorithm.*

1. **Grid Search vs. Random Search**: Automatically testing different "knobs" on the model (e.g., the depth of a decision tree or the learning rate) to find the optimal combination.
2. **Overfitting Detection**: Monitoring if performance is great on training data but poor on validation data—a sign the model is memorizing noise rather than learning biology.

## 5. Biological Interpretation and Feature Importance
*Moving from a "black box" to understanding the actual biology learned by the model.*

1. **Feature Importance (Variable Importance)**: Ranking which genes, metabolites, or clinical markers were most influential in making the prediction.
2. **Biological Context**:
    * **Pathway Analysis**: Checking if the "important" genes belong to a specific biological pathway (e.g., inflammation).
    * **Directionality**: Understanding if a feature's increase or decrease leads to the predicted outcome.
3. **Latent Feature Review**: Checking if the model picked up on "Batch Effects" (e.g., samples from Lab A vs. Lab B) rather than actual biology.

# remove any environment variables to prevent any code reading 
rm(list = ls())

# Load required libraries
library(corrplot) # For plotting a correlation  
library(reticulate) # For Python connection 
library(ggplot2) # For plotting colour graph + boxplot + chart 
library(moments)  # For skewness calculation



# use reticulate to connect to local python environment, in this case, connect to conda environment  
reticulate::use_python("USER/python.exe", required = TRUE)
# check if reticulate is connect to local environment and version of the python libraries 
reticulate::py_config()

# Check and set a working directory in Rstudio
setwd("USER/Work_dir")
# Import dataset
data <- read.csv("plagiodontia_database.csv")

# Extract numeric data
numeric_data <- data[, sapply(data, is.numeric)]

# Remove constant variables
constant_vars <- sapply(numeric_data, function(x) {
  sd(x, na.rm = TRUE) < .Machine$double.eps
})

if(any(constant_vars)) {
  cat("Removing constant variables:", names(constant_vars)[constant_vars], "\n")
  numeric_data <- numeric_data[, !constant_vars]
}

# Get a numerical data 
data_num <- numeric_data

# Clone GAIN repository if not exists
if(!dir.exists("GAIN-master")) {
  system("git clone https://github.com/jsyoon0823/GAIN.git")
}


# import GAIN from directory 
gain <- import_from_path("gain", path = "USER/GAIN-master")


# Metrics Calculation Function

calculate_rsquared <- function(observed, predicted) { # create function with 2 variables: observed (same data) and predicted (new data)
  valid_idx <- complete.cases(observed, predicted) # create a vector where only return TRUE if observed and predicted is not missing
  # keep only valid idx variables for both
  observed <- observed[valid_idx]
  predicted <- predicted[valid_idx]
  
  # If fewer than 2 valid observations remain, returns a list with all metrics as NA
  if(length(observed) < 2) return(list(R2 = NA, Adj_R2 = NA, MSE = NA, MAE = NA, NRMSE = NA))
  
  # Calculate metrics
  ss_res <- sum((observed - predicted)^2) # sum of square residues,  measures how much the predictions deviate from the observed values
  ss_tot <- sum((observed - mean(observed))^2) # sum of total square, measures the total variance in the observed data around its mean
  r2 <- 1 - (ss_res / ss_tot) # R-square 
  
  # Adjusted R-square 
  n <- length(observed)
  adj_r2 <- 1 - (1 - r2) * ((n - 1) / (n - 2))
  
  # Error metrics
  mse <- mean((observed - predicted)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(observed - predicted))
  mape <- mean(abs((observed - predicted) / observed)) * 100
  
  # Range-normalized NRMSE
  data_range <- max(observed, na.rm = TRUE) - min(observed, na.rm = TRUE)
  nrmse_range <- rmse / data_range
  
  # Mean-normalized NRMSE 
  nrmse_mean <- rmse / mean(observed, na.rm = TRUE)
  
  # SD-normalized NRMSE
  nrmse_sd <- rmse / sd(observed, na.rm = TRUE)
  
  # Correlation
  pearson_r <- cor(observed, predicted, method = "pearson")
  spearman_rho <- cor(observed, predicted, method = "spearman")
  
  return(list(
    R2 = r2,
    Adj_R2 = adj_r2,
    RMSE = rmse,
    NRMSE_range = nrmse_range,  
    NRMSE_mean = nrmse_mean,    
    NRMSE_sd = nrmse_sd,      
    MAE = mae,
    MSE = mse,
    MAPE = mape,
    Pearson_r = pearson_r,
    Spearman_rho = spearman_rho
  ))
}  






# setting GAIN parameters 

n_samples <- nrow(data_num) # number of samples 
cat("Number of samples:", n_samples, "\n")

# Find valid batch size
valid_batch_sizes <- c(8, 12, 14, 21, 24, 28, 42, 56, 84, 168)
valid_batch_sizes <- valid_batch_sizes[valid_batch_sizes <= n_samples]
optimal_batch <- valid_batch_sizes[which.min(abs(valid_batch_sizes - 32))]
cat("Using batch size:", optimal_batch, "\n")


# GAIN parameters
gain_parameters <- dict(
  batch_size = as.integer(optimal_batch),
  hint_rate = 0.3,        # Lower hint rate preserves distribution
  alpha = 20,             # Balanced generator/discriminator
  iterations = 5000L
)

# scaling function 

# using z-score for scaling 
# using mean_value and sd_value further to prevent data leakage 

scale_zscore <- function(x, mean_val = NULL, sd_val = NULL) { # if mean_value and sd_value are not provided, using data to calculate
  
  # setting condition where both values are not provided 
  if(is.null(mean_val)) {
    # calculate mean and SD, also ignore NA value 
    mean_val <- apply(x, 2, mean, na.rm = TRUE)
    sd_val <- apply(x, 2, sd, na.rm = TRUE)
  }
  
  # handle zero variance by setting any smaller SD to 1 to avoid division by zero
  sd_val[sd_val < .Machine$double.eps] <- 1
  
  # perform Z-score 
  scaled <- sweep(sweep(x, 2, mean_val, "-"), 2, sd_val, "/")
  
  # return values 
  return(list(
    scaled = scaled,
    mean_val = mean_val,
    sd_val = sd_val
  ))
}





# perform artificial masking test data using MAR protocol 

# set random seed to ensure reproducibility
set.seed(123)

# create MAR mask using same dimension with data_num
# initially false(not marked) and will marked true for values to be artificially removed
test_mask <- matrix(FALSE, nrow = nrow(data_num), ncol = ncol(data_num))

# Loop funtion for MAR masking 
for(j in 1:ncol(data_num)) { # for each column j, decideds which rows will be masked
  predictor_col <- 1 # declare first column
  # skip columns that are 
  # 1. same as column j (can't use the same column to predict its own missingness)
  # 2. constant column (prevent zero variances(no predictive power))
  while(predictor_col <= ncol(data_num) && 
        (predictor_col == j || sd(data_num[, predictor_col], na.rm = TRUE) == 0)) {
    predictor_col <- predictor_col + 1 # predict only last column
  }
  
  
  if(predictor_col <= ncol(data_num)) { # if valid predictor found
    # extract the column and preform z-scaling (standardize)
    predictor <- data_num[, predictor_col]
    predictor_scaled <- (predictor - mean(predictor, na.rm = TRUE)) / sd(predictor, na.rm = TRUE)
    # calculate probability of missingness of each row 
    # use plogis to convert scale between 0 and 1
    # minimum missingness = 0.05 + 0.25 * 0 = 0.05 = 5% 
    # maximum missingness = 0.05 + 0.25 * 1 = 0.30 = 30%
    # at mean = 0.05 + 0.25 * 0.5 = 0.175 = 17.5% 
    prob_missing <- 0.05 + 0.25 * plogis(predictor_scaled)
    # setting the missingness to depend on predictor value (MAR)
    # stochastic pattern
    test_mask[, j] <- runif(nrow(data_num)) < prob_missing # if random number < probability, mark as TRUE = masked 
  } else { # if not a suitable predictor
    test_mask[, j] <- runif(nrow(data_num)) < 0.1 # use MCAR with 10% constant probability
  }
}

# identifies columns that already contain NA in original data 
original_missing <- is.na(data_num)
# only masked if already masked with loop function and not contain NA in original dataset
test_mask <- test_mask & !original_missing
# safety check 
test_mask[is.na(test_mask)] <- FALSE

# count total TRUE values in MAR masking step
total_masked <- sum(test_mask, na.rm = TRUE)
cat("Number of values artificially masked:", total_masked, "\n")

# Create MAR test dataset
data_test <- data_num
data_test[test_mask] <- NA


# Performing GAIN imputation on test dataset 

# Scale test data to incresing stability
scaling_result <- scale_zscore(data_test) # standardise with z-scaling function 
scaled_test <- scaling_result$scaled # using only scaled data, into vector
scaled_test_matrix <- as.matrix(scaled_test) # convert as matrix for GAIN

# Create missingness mask for GAIN: 0 = missing, 1 = observed
# Let GAIN know which value are real and which to impute
mask_matrix <- matrix(1, nrow = nrow(scaled_test_matrix), ncol = ncol(scaled_test_matrix))
mask_matrix[is.na(scaled_test_matrix)] <- 0 # set NA to 0

# Convert to Python format
scaled_test_py <- r_to_py(scaled_test_matrix) # use reticulate to converts R objects to python object
mask_py <- r_to_py(mask_matrix)



# Run GAIN
cat("\nRunning GAIN imputation on masked test data\n")
imputed_scaled <- tryCatch({ # use trycatch to prevent any errors 
  # call gain from python on standardised data with NAs with parameters and missingness indicator matrix 
  gain$gain(scaled_test_py, gain_parameters, mask = mask_py) 
}, error = function(e) { # fallback function 
  cat("Trying without mask\n") # retry without masking parameter
  gain$gain(scaled_test_py, gain_parameters)
})


# Convert back to original scale (Reverse Standardisation)
imputed_test <- imputed_scaled * scaling_result$sd_val + scaling_result$mean_val

# Force non-negative globally (domain-specific constraints)
# Apply soft constraints (winsorization)
for(j in 1:ncol(imputed_test)) { # Loop through each columns in imputed test dataset
  obs_vals <- data_num[, j][!is.na(data_num[, j])] # Extracts originally observed values
  
  # Calculate bounds using IQR (InterQuartile Range)
  if(length(obs_vals) > 0) { # Only calculate if there are observed values
    # Q1 = 25% and Q3 = 75%
    Q1 <- quantile(obs_vals, 0.25, na.rm = TRUE)
    Q3 <- quantile(obs_vals, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    
    # Use Tukey's fences = value outset bound is an extreme outlier
    lower_bound <- Q1 - 3 * IQR_val
    upper_bound <- Q3 + 3 * IQR_val
    
    # Apply winsorization = caps extreme values 
    # Replace value > upper bound with upper value, same goes with lower bound
    imputed_test[, j] <- pmax(pmin(imputed_test[, j], upper_bound), lower_bound)
    
    # Check if all observed values are non-negative 
    if(all(obs_vals >= 0, na.rm = TRUE)) {
      imputed_test[, j] <- pmax(imputed_test[, j], 0)
    }
  }
}

# Overall GAIn performance

true_values <- data_num[test_mask]
pred_values <- imputed_test[test_mask]

overall_metrics <- calculate_rsquared(true_values, pred_values)


cat(sprintf("R²:    %.4f\n", overall_metrics$R2))
cat(sprintf("Adjusted R²:                           %.4f\n", overall_metrics$Adj_R2))
cat(sprintf("Pearson correlation:               %.4f\n", overall_metrics$Pearson_r))
cat(sprintf("Spearman rank correlation:         %.4f\n", overall_metrics$Spearman_rho))
cat(sprintf("RMSE:                                  %.4f\n", overall_metrics$RMSE))
cat(sprintf("NRMSE:              %.4f\n", overall_metrics$NRMSE_range))  # Add this
cat(sprintf("MAE:                                   %.4f\n", overall_metrics$MAE))
cat(sprintf("MAPE:                                  %.2f%%\n", overall_metrics$MAPE))


# Re-Imputation on dataset to create a fully populated version

# Scale full dataset (same method as test dataset)
scaling_full <- scale_zscore(data_num) # Standardise on original data
scaled_full <- scaling_full$scaled
scaled_full_matrix <- as.matrix(scaled_full)

# Create mask for full dataset
mask_full <- matrix(1, nrow = nrow(scaled_full_matrix), ncol = ncol(scaled_full_matrix))
mask_full[is.na(scaled_full_matrix)] <- 0

# Run GAIN on full dataset
cat("\nRunning GAIN imputation on full data\n")
imputed_full_scaled <- tryCatch({
  gain$gain(scaled_full_matrix, gain_parameters, mask = mask_full)
}, error = function(e) {
  gain$gain(scaled_full_matrix, gain_parameters)
})

# Convert back to original scale (reverse standardisation)
imputed_full <- imputed_full_scaled * scaling_full$sd_val + scaling_full$mean_val

# Apply soft constraints
for(j in 1:ncol(imputed_full)) {
  obs_vals <- data_num[, j][!is.na(data_num[, j])]
  
  if(length(obs_vals) > 0) {
    
    Q1 <- quantile(obs_vals, 0.25, na.rm = TRUE)
    Q3 <- quantile(obs_vals, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    
    lower_bound <- Q1 - 3 * IQR_val
    upper_bound <- Q3 + 3 * IQR_val
    
    # Winsorization 
    imputed_full[, j] <- pmax(pmin(imputed_full[, j], upper_bound), lower_bound)
    

  }
}

# Add random noise to imputed values 
missing_idx <- is.na(data_num) 

for(j in 1:ncol(imputed_full)) { 
  sd_j <- sd(data_num[, j], na.rm = TRUE) # Calculate the sd of observed values only (shows the natural variability of the variables)
  
  # Add random noise to reflect imputation uncertainty and prevent identical imputation
  imputed_full[missing_idx[, j], j] <- 
    imputed_full[missing_idx[, j], j] +
    rnorm(sum(missing_idx[, j]), 0, 0.05 * sd_j)
}

# Final check if model collapse of not (no large/constant value)
apply(imputed_full, 2, function(x) max(table(round(x, 4))))






# Save output 

# Create final dataset
final_data <- data
numeric_cols <- sapply(data, is.numeric)
final_data[, numeric_cols] <- imputed_full

# Save files
write.csv(final_data, "GAIN_imputed_full_complete_dataset.csv", row.names = FALSE)
write.csv(imputed_full, "GAIN_imputed_full_numeric_dataset.csv", row.names = FALSE)

# Save detailed metrics
metrics_summary <- data.frame(
  Metric = c("R2", "Adj_R2", "Pearson_r", "Spearman_rho", "RMSE", "MAE", "MAPE"),
  Value = c(overall_metrics$R2, overall_metrics$Adj_R2, overall_metrics$Pearson_r,
            overall_metrics$Spearman_rho, overall_metrics$RMSE, overall_metrics$MAE,
            overall_metrics$MAPE)
)

write.csv(metrics_summary, "GAIN_full_metrics_summary.csv", row.names = FALSE)
write.csv(per_var_metrics, "GAIN_full_per_variable_metrics.csv", row.names = FALSE)

# shows summaries


cat(sprintf("Overall R²: %.4f (%.2f%% variance explained)\n", 
            overall_metrics$R2, overall_metrics$R2 * 100))
cat(sprintf("Overall correlation: r = %.4f\n", overall_metrics$Pearson_r))
cat(sprintf("RMSE: %.4f\n", overall_metrics$RMSE))
cat("\nNegative values check:\n")
cat("Test data:", sum(imputed_test < 0, na.rm = TRUE), "\n")
cat("Full data:", sum(imputed_full < 0, na.rm = TRUE), "\n")
cat("\nVariances check:\n")
dup_ratio <- apply(imputed_full, 2, function(x) {
  max(table(round(x, 3))) / length(x)
})
summary(dup_ratio) # Max value should not exceed 0.7 ++ 




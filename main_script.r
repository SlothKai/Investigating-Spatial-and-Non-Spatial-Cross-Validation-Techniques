# Load the necessary libraries
library(ranger)
library(caret)
library(readr)
library(dplyr)
library(spatialsample)
library(sf)

source("general_functions.r")

source("Cross_Validation_Technique/Random_KFoldCV.r")
source("Cross_Validation_Technique/BootstrapCV.r")
source("Cross_Validation_Technique/Importance_Weighted_CV.r")

source("Cross_Validation_Technique/Spatial_KFoldCV.r")
source("Cross_Validation_Technique/BlockedCV.r")
source("Cross_Validation_Technique/BufferedCV.r")


results_K <- data.frame(
  Range = integer(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

results_Boot <- data.frame(
  Range = integer(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

results_IWCV <- data.frame(
  Range = integer(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

results_SpatialK <- data.frame(
  Range = integer(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

results_BlockedCV <- data.frame(
  Method = character(),
  Block_Size = numeric(),
  Range = numeric(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

results_BlockedBuffCV <- data.frame(
  Method = character(),
  Block_Size = numeric(),
  Range = numeric(),
  Buffer = numeric(),
  Dataset = integer(),
  Fold = integer(),
  RMSE_Test = numeric(),
  MAE_Test = numeric(),
  R_Squared_Test = numeric(),
  MAPE_Test = numeric(),
  BIAS_Test = numeric(),
  RMSE_Validation = numeric(),
  MAE_Validation = numeric(),
  R_Squared_Validation = numeric(),
  MAPE_Validation = numeric(),
  BIAS_Validation = numeric()
)

# Define range
range_array <- c(1, 6, 12)
# Define scenarios
scenarios <- c("sd", "sdcs", "si", "sics")
# Define number of datasets
number_of_Dataset <- 50
# Coordinates columns (in this case 'lat' and 'lon')
coords_cols <- c("lon", "lat")
# Define the folds to evaluate, Random K fold and Spatial K fold
folds_array <- c(2, 5, 10)
# Number of resamples. For bootstrap
numbers_array <- c(50, 75, 100)


# Create a tuning grid with default values for regression
cat("Creating the tuning grid...\n")

tune_grid <- expand.grid(
  mtry = NA, # Using all variables
  splitrule = "variance", # Default for regression
  min.node.size = 1 # Default minimum node size
)

for (scene in scenarios) {
  for (range in range_array) {
    cat("Range: ", range, "...\n")
    validation_set <- paste("simulated_data/validation_", scene, range, ".csv", sep = "")
    validation_data <- invisible(read.csv(validation_set))
    validation_X <- validation_data[, -which(names(validation_data) == "z")]
    validation_y <- validation_data$z

    # For block_cv
    if (range == 1) {
      block_size <- 2
    } else {
      block_size <- range
    }

    for (dataset in 1:number_of_Dataset) {
      cat(" Dataset: ", dataset, "...\n")
      spatial_file <- paste("simulated_data/spatial/", scene, "/complete_", range, "_", dataset, ".csv", sep = "")
      non_spatial_train <- paste("simulated_data/non_spatial/", scene, "/train_", range, "_", dataset, ".csv", sep = "")
      non_spatial_test <- paste("simulated_data/non_spatial/", scene, "/test_", range, "_", dataset, ".csv", sep = "")


      # Load the train and test data
      train_data <- invisible(read.csv(non_spatial_train))
      test_data <- invisible(read.csv(non_spatial_test))
      spatial_data <- invisible(read.csv(spatial_file))

      # Preparing Train and Test for modelling.
      train_X <- train_data[, -which(names(train_data) == "z")]
      train_y <- train_data$z

      test_X <- test_data[, -which(names(test_data) == "z")]
      test_y <- test_data$z

      spatial_X <- spatial_data[, -which(names(spatial_data) == "z")]
      spatial_y <- spatial_data$z
      # Initialize features for model training
      features_to_include <- coords_cols # Start with 'lat' and 'lon'

      for (feature in setdiff(names(train_data), c("z", coords_cols))) {
        features_to_include <- c(features_to_include, feature)
        cat("   Features: ", features_to_include, "...\n")
        # Prepare training and test sets with the current features
        train_X <- train_data[, features_to_include]
        test_X <- test_data[, features_to_include]
        spatial_X <- spatial_data[, features_to_include]

        # Set mtry based on the number of features
        mtry_value <- length(features_to_include)
        tune_grid$mtry[1] <- mtry_value
        # Run the evaluation with Random K-Fold CV
        temp_results_K <- evaluate_with_k_folds(train_X, train_y, test_X, test_y, validation_X, validation_y, tune_grid, range, dataset, folds_array, scene)

        # Run the evaluation with Bootstrap CV
        temp_results_Boot <- evaluate_with_Bootstrap(train_X, train_y, test_X, test_y, validation_X, validation_y, tune_grid, range, dataset, numbers_array, scene)

        # Run the evaluation with Importance Weighted CV
        temp_results_IWCV <- evaluate_with_ImportanceWeightedCV(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, folds_array, scene)

        # Run the evaluation with Spatial K-Fold CV
        temp_results_SpatialK <- evaluate_with_spatial_k_folds(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, folds_array, coords_cols, scene)

        # Run the evaluation with Blocked CV
        temp_results_BlockedCV <- evaluate_with_blocked_cv(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, folds_array, block_size, coords_cols, scene)

        # Run the evaluation with Buffered CV
        temp_results_BlockedBuffCV <- evaluate_with_blockedBuff_cv(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, folds_array, block_size, coords_cols, scene)

        # Record results
        # Add scenario and features to the results, ensuring completeness of information
        temp_results_K$Scenario <- paste(scene, collapse = ", ")
        temp_results_K$Features <- paste(features_to_include, collapse = ", ")
        temp_results_Boot$Scenario <- paste(scene, collapse = ", ")
        temp_results_Boot$Features <- paste(features_to_include, collapse = ", ")
        temp_results_IWCV$Scenario <- paste(scene, collapse = ", ")
        temp_results_IWCV$Features <- paste(features_to_include, collapse = ", ")
        temp_results_SpatialK$Scenario <- paste(scene, collapse = ", ")
        temp_results_SpatialK$Features <- paste(features_to_include, collapse = ", ")
        temp_results_BlockedCV$Scenario <- paste(scene, collapse = ", ")
        temp_results_BlockedCV$Features <- paste(features_to_include, collapse = ", ")
        temp_results_BlockedBuffCV$Scenario <- paste(scene, collapse = ", ")
        temp_results_BlockedBuffCV$Features <- paste(features_to_include, collapse = ", ")

        # Append the results to the main data frame of each specific technique
        results_K <- rbind(results_K, temp_results_K)
        results_Boot <- rbind(results_Boot, temp_results_Boot)
        results_IWCV <- rbind(results_IWCV, temp_results_IWCV)
        results_SpatialK <- rbind(results_SpatialK, temp_results_SpatialK)
        results_BlockedCV <- rbind(results_BlockedCV, temp_results_BlockedCV)
        results_BlockedBuffCV <- rbind(results_BlockedBuffCV, temp_results_BlockedBuffCV)
      }
    }
  }
}

# Write results to CSV
RandomKFoldCV_CSV <- paste("RandomKFoldCV.csv")
BootstrapCSV <- paste("BootstrapCV_Results.csv")
ImportanceWeightedCSV <- paste("ImportanceWeightedCV_results.csv")
SpatialK_CSV <- paste("SptialKFoldCV.csv")
BlockCV <- paste("BlockedCV.csv")
BlockBuffCV <- paste("BlockedBuffCV.csv")

write_csv(results_K, RandomKFoldCV_CSV)
write_csv(results_Boot, BootstrapCSV)
write.csv(results_IWCV, ImportanceWeightedCSV)
write.csv(results_SpatialK, SpatialK_CSV)
write.csv(results_BlockedCV, BlockCV)
write.csv(results_BlockedBuffCV, BlockBuffCV)

cat("Results written")

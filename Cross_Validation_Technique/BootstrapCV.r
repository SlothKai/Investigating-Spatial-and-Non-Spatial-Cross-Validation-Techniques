# Load the necessary libraries
library(ranger)
library(caret)
library(readr)
library(dplyr)
library(doParallel)
library(parallel)


# Source functions
source("general_functions.r")

# Main function to evaluate models with varying folds
evaluate_with_Bootstrap <- function(train_X, train_y, test_X, test_y, validation_X, validation_y, tune_grid, range, dataset, fold_array, scene) {
  results <- data.frame(
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

  # Detect the number of CPU cores and set up a parallel backend
  num_cores <- parallel::detectCores()
  cat("Setting up parallel processing with", num_cores, "cores...\n")
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Ensure that the cluster is stopped on exit, even if an error occurs
  on.exit(stopCluster(cl))

  for (fold in fold_array) {
    cat("   Evaluating with", fold, "resamples... \n")

    # Set up cross-validation
    train_control <- trainControl(method = "boot", number = fold, summaryFunction = custom_metrics)

    # Train the random forest model
    cat("     Training model...\n")
    model_rf <- train(train_X, train_y,
      method = "ranger",
      trControl = train_control,
      metric = "RMSE",
      tuneGrid = tune_grid
    )

    # Make predictions on the test set
    cat("     Making predictions on the test set...\n")
    predictions_test <- predict(model_rf, newdata = test_X)

    # Make predictions on the validation set
    cat("     Making predictions on the validation set...\n")
    predictions_validation <- predict(model_rf, newdata = validation_X)

    # Calculate metrics on the test set
    rmse_test <- rmse(test_y, predictions_test)
    mae_test <- mae(test_y, predictions_test)
    r2_test <- r_squared(test_y, predictions_test)
    mape_test <- mape(test_y, predictions_test)
    bias_test <- bias(test_y, predictions_test)

    # Calculate metrics on the validation set
    rmse_validation <- rmse(validation_y, predictions_validation)
    mae_validation <- mae(validation_y, predictions_validation)
    r2_validation <- r_squared(validation_y, predictions_validation)
    mape_validation <- mape(validation_y, predictions_validation)
    bias_validation <- bias(validation_y, predictions_validation)

    # Results DF
    results <- rbind(results, data.frame(
      Range = range,
      Dataset = dataset,
      Fold = fold,
      RMSE_Test = rmse_test,
      MAE_Test = mae_test,
      R_Squared_Test = r2_test,
      MAPE_Test = mape_test,
      BIAS_Test = bias_test,
      RMSE_Validation = rmse_validation,
      MAE_Validation = mae_validation,
      R_Squared_Validation = r2_validation,
      MAPE_Validation = mape_validation,
      BIAS_Validation = bias_validation,
      stringsAsFactors = FALSE
    ))
  }

  # Stop the parallel backend
  stopCluster(cl)
  cat("BOOTSTRAP Parallel processing completed and cluster stopped", range, dataset, scene)

  return(results)
}

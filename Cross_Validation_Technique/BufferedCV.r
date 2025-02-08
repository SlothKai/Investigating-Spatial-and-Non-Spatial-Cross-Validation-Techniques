# Load necessary libraries
library(ranger)
library(caret)
library(readr)
library(dplyr)
library(spatialsample) # for spatial k-fold
library(sf) # for handling spatial data (simple features)
library(ggplot2)
library(parallel)
library(doParallel)

# Source functions
source("general_functions.r")

block_method <- c("random", "snake", "continuous")

# Main function to evaluate models with spatial k-fold cross validation
evaluate_with_blockedBuff_cv <- function(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, fold_array, block_size, coords_cols, scene) {
    # Convert train and test data to sf objects
    cat("Converting data to sf objects...\n")
    train_X_sf <- st_as_sf(spatial_X, coords = coords_cols, crs = 3857)

    results <- data.frame(
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
        BIAS_Validation = numeric(),
        # Track the number of training points, reason for this is that I've encountered cases where
        # there were no more training points left due to the buffer zone. Tracking the points helps to
        # understand the reason for the error, as well as the number of training points used in the model.
        Analysis_Size = integer(), # Track the number of training points
        Assessment_Size = integer(), # Track the number of test points
        Buffer_Size = integer() # Track the number of buffer points
    )

    # Detect the number of CPU cores and set up a parallel backend
    num_cores <- parallel::detectCores()
    cat("Setting up parallel processing with", num_cores, "cores...\n")
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # Ensure that the cluster is stopped on exit, even if an error occurs
    on.exit(stopCluster(cl))

    for (fold in fold_array) {
        cat("Evaluating with ", fold, "fold...\n")

        for (method in block_method) {
            cat("method: ", method, "...\n")

            # Try-catch for the spatial fold creation
            spatial_folds <- tryCatch(
                {
                    spatial_block_cv(train_X_sf, cellsize = 12.5, v = fold, radius = NULL, buffer = range, method = method)
                },
                error = function(e) {
                    cat("Error during spatial fold creation:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            # Check if the spatial_folds object is NULL (due to error)
            if (is.null(spatial_folds)) {
                next
            }

            # Try-catch for the extraction of data points (analysis, assessment, buffer)
            fold_data <- tryCatch(
                {
                    lapply(spatial_folds$splits, function(split) {
                        analysis_size <- nrow(analysis(split)) # Number of analysis (training) points
                        assessment_size <- nrow(assessment(split)) # Number of assessment (testing) points
                        total_points <- nrow(spatial_X) # Total number of data points
                        buffer_size <- total_points - analysis_size - assessment_size # Calculate buffer size
                        list(analysis_size = analysis_size, assessment_size = assessment_size, buffer_size = buffer_size)
                    })
                },
                error = function(e) {
                    cat("Error during extraction of rsplit data points:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            # Extract training and test indices from spatial_folds
            train_indices <- tryCatch(
                {
                    lapply(spatial_folds$splits, function(split) {
                        as.integer(row.names(analysis(split))) # Correctly extract row indices for training
                    })
                },
                error = function(e) {
                    cat("Error during training index extraction:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            test_indices <- tryCatch(
                {
                    lapply(spatial_folds$splits, function(split) {
                        as.integer(row.names(assessment(split))) # Correctly extract row indices for testing
                    })
                },
                error = function(e) {
                    cat("Error during test index extraction:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            # Train control setup
            train_control <- trainControl(
                method = "cv",
                index = train_indices,
                indexOut = test_indices,
                summaryFunction = custom_metrics
            )

            # Try-catch for model training
            model_rf <- tryCatch(
                {
                    train(
                        x = spatial_X,
                        y = spatial_y,
                        method = "ranger",
                        trControl = train_control,
                        metric = "RMSE",
                        tuneGrid = tune_grid
                    )
                },
                error = function(e) {
                    cat("Error during model training:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            # Try-catch for validation prediction
            predictions_validation <- tryCatch(
                {
                    predict(model_rf, newdata = validation_X)
                },
                error = function(e) {
                    cat("Error during validation prediction:", conditionMessage(e), "\n")
                    next # Skip this method if there is an error
                }
            )

            # Calculate metrics on the validation set
            rmse_validation <- rmse(validation_y, predictions_validation)
            mae_validation <- mae(validation_y, predictions_validation)
            r2_validation <- r_squared(validation_y, predictions_validation)
            mape_validation <- mape(validation_y, predictions_validation)
            bias_validation <- bias(validation_y, predictions_validation)

            # Extract test metrics
            rmse_test <- mean(model_rf$results$RMSE) # Mean RMSE from cross-validation
            mae_test <- mean(model_rf$results$MAE) # Mean MAE from cross-validation
            r2_test <- mean(model_rf$results$R2) # Mean R² from cross-validation
            mape_test <- mean(model_rf$results$MAPE) # Assuming you've implemented MAPE in the model
            bias_test <- mean(model_rf$results$Bias) # Assuming you've implemented Bias in the model

            # Add results for each fold
            for (i in seq_along(fold_data)) {
                results <- rbind(results, data.frame(
                    Method = method,
                    Block_Size = block_size,
                    Range = range,
                    Buffer = range,
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
                    Analysis_Size = fold_data[[i]]$analysis_size,
                    Assessment_Size = fold_data[[i]]$assessment_size,
                    Buffer_Size = fold_data[[i]]$buffer_size,
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    cat("BLOCKED BUFFERED CV Parallel processing completed and cluster stopped.", range, dataset, scene, "\n")

    return(results)
}

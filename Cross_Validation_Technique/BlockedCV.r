# Load necessary libraries
library(ranger)
library(caret)
library(readr)
library(dplyr)
library(spatialsample) # for spatial k-fold
library(sf) # for handling spatial data (simple features)
library(ggplot2)

# Source functions
source("general_functions.r")

block_method <- c("random", "snake", "continuous")

# Main function to evaluate models with spatial k-fold cross validation
evaluate_with_blocked_cv <- function(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, fold_array, block_size, coords_cols, scene) {
    # Convert train and test data to sf objects
    cat("Converting data to sf objects...\n")
    train_X_sf <- st_as_sf(spatial_X, coords = coords_cols, crs = 3857)

    results <- data.frame(
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
            print(paste("method: ", method, "...\n"))
            # Set up spatial cross-validation using the sf object
            spatial_folds <- spatial_block_cv(train_X_sf, cellsize = block_size, v = fold, method = method) # Create spatial folds
            # Extract training and test indices from spatial_folds
            train_indices <- lapply(spatial_folds$splits, function(split) {
                as.integer(row.names(analysis(split))) # Correctly extract row indices for training
            })
            test_indices <- lapply(spatial_folds$splits, function(split) {
                as.integer(row.names(assessment(split))) # Correctly extract row indices for testing
            })

            train_control <- trainControl(
                method = "cv",
                index = train_indices,
                indexOut = test_indices,
                summaryFunction = custom_metrics
            )

            # Train the random forest model
            cat("  Training model...\n")
            model_rf <- train(
                x = spatial_X,
                y = spatial_y,
                method = "ranger",
                trControl = train_control,
                metric = "RMSE",
                tuneGrid = tune_grid
            )

            # Extract the metrics from the trained model
            rmse_test <- mean(model_rf$results$RMSE) # Mean RMSE from cross-validation
            mae_test <- mean(model_rf$results$MAE) # Mean MAE from cross-validation
            r2_test <- mean(model_rf$results$R2) # Mean R² from cross-validation
            mape_test <- mean(model_rf$results$MAPE) # Assuming you've implemented MAPE in the model
            bias_test <- mean(model_rf$results$Bias) # Assuming you've implemented Bias in the model

            # Make predictions on the validation set
            cat("     Making predictions on the validation set...\n")
            predictions_validation <- predict(model_rf, newdata = validation_X)

            # Calculate metrics on the validation set
            rmse_validation <- rmse(validation_y, predictions_validation)
            mae_validation <- mae(validation_y, predictions_validation)
            r2_validation <- r_squared(validation_y, predictions_validation)
            mape_validation <- mape(validation_y, predictions_validation)
            bias_validation <- bias(validation_y, predictions_validation)

            # Results DF
            results <- rbind(results, data.frame(
                Method = method,
                Block_Size = block_size,
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
    }

    cat("BLOCKED CV Parallel processing completed and cluster stopped ", range, dataset, scene, "\n")

    return(results)
}

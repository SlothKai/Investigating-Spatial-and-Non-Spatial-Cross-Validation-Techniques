library(ranger)
library(caret)
library(readr)
library(dplyr)
library(densratio)
library(doParallel)

source("general_functions.r")

evaluate_with_ImportanceWeightedCV <- function(spatial_X, spatial_y, validation_X, validation_y, tune_grid, range, dataset, fold_array, scene) {
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

    on.exit(stopCluster(cl))

    for (fold in fold_array) {
        cat("   Evaluating with", fold, "folds...\n")

        # Set up cross-validation partitioning
        cv_split <- createFolds(spatial_y, k = fold, returnTrain = TRUE)

        # Temporary data frame to store metrics for this particular fold
        fold_results <- data.frame(
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


        for (i in seq_along(cv_split)) {
            cat("     Processing fold", i, "...\n")
            train_idx <- cv_split[[i]]
            test_idx <- setdiff(seq_len(nrow(spatial_X)), train_idx)

            train_X <- spatial_X[train_idx, ]
            train_y <- spatial_y[train_idx]
            test_X <- spatial_X[test_idx, ]
            test_y <- spatial_y[test_idx]

            # Estimate importance weights using RuLSIF
            densratio_result <- densratio(as.matrix(train_X), as.matrix(test_X), method = "RuLSIF", alpha = 0, kernel_num = 50)
            importance_weights <- densratio_result$compute_density_ratio(as.matrix(train_X))

            # Train the random forest model with importance weights
            train_control <- trainControl(method = "cv", number = fold, summaryFunction = custom_metrics)

            cat("     Training model...\n")
            model_rf <- train(
                x = train_X,
                y = train_y,
                method = "ranger",
                trControl = train_control,
                metric = "RMSE",
                tuneGrid = tune_grid,
                weights = importance_weights
            )

            # Predictions and metrics for test and validation sets
            predictions_test <- predict(model_rf, newdata = test_X)
            predictions_validation <- predict(model_rf, newdata = validation_X)

            # Calculate metrics
            rmse_test <- rmse(test_y, predictions_test)
            mae_test <- mae(test_y, predictions_test)
            r2_test <- r_squared(test_y, predictions_test)
            mape_test <- mape(test_y, predictions_test)
            bias_test <- bias(test_y, predictions_test)

            rmse_validation <- rmse(validation_y, predictions_validation)
            mae_validation <- mae(validation_y, predictions_validation)
            r2_validation <- r_squared(validation_y, predictions_validation)
            mape_validation <- mape(validation_y, predictions_validation)
            bias_validation <- bias(validation_y, predictions_validation)

            # Append the results for this fold to the fold-specific results data frame
            fold_results <- rbind(fold_results, data.frame(
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

        # After processing all iterations within a fold, compute the mean for each fold
        mean_fold_results <- fold_results %>%
            summarise(
                RMSE_Test = mean(RMSE_Test, na.rm = TRUE),
                MAE_Test = mean(MAE_Test, na.rm = TRUE),
                R_Squared_Test = mean(R_Squared_Test, na.rm = TRUE),
                MAPE_Test = mean(MAPE_Test, na.rm = TRUE),
                BIAS_Test = mean(BIAS_Test, na.rm = TRUE),
                RMSE_Validation = mean(RMSE_Validation, na.rm = TRUE),
                MAE_Validation = mean(MAE_Validation, na.rm = TRUE),
                R_Squared_Validation = mean(R_Squared_Validation, na.rm = TRUE),
                MAPE_Validation = mean(MAPE_Validation, na.rm = TRUE),
                BIAS_Validation = mean(BIAS_Validation, na.rm = TRUE)
            ) %>%
            mutate(
                Range = range,
                Dataset = dataset,
                Fold = fold
            )

        # Append the mean results for this fold to the overall results data frame
        results <- rbind(results, mean_fold_results)
    }

    # Stop the parallel backend
    # stopCluster(cl)
    cat("Parallel processing completed and cluster stopped.", range, dataset, scene, "\n")

    return(results)
}

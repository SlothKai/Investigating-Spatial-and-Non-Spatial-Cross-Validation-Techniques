# Load required libraries
library(sp)
library(gstat)
library(ggplot2)
library(akima)
library(readr)

source("general_functions.r")

# Define number of datasets to generate
number_of_Dataset <- 50

# Number of points for training and test sets
n_train <- 2000 # Number of training points
n_test <- 500 # Number of test points

# Define Gaussian variogram model for spatial dependence
sill <- 1.0 # Total variance (sill = partial sill + nugget)
nugget <- 0.1 # Nugget effect (small-scale random variance)
range_array <- c(1, 6, 12) # Range of spatial correlation
nugget_noise <- 0.1

# Define a small tolerance level for mean
tolerance <- 1e-10

results <- data.frame(
    Range = numeric(0),
    Dataset = numeric(0),
    Predictor = character(0), # Empty character vector for Predictor
    KS_Statistic = numeric(0), # Empty numeric vector for KS_Statistic
    P_Value = numeric(0) # Empty numeric vector for P_Value
)

for (range in range_array) {
    cat("Range: ", range, "...\n")
    # Create the variogram model
    vgm_model_predict <- vgm(psill = sill - nugget, model = "Gau", range = range, nugget = nugget)
    vgm_model_noise <- vgm(psill = sill - nugget, model = "Gau", range = range, nugget = nugget_noise)

    for (dataset in 1:number_of_Dataset) {
        cat("   Preparing Dataset: ", dataset, "...\n")
        complete_dataset <- paste("simulated_data/spatial/si/si_complete_", range, "_", dataset, ".csv", sep = "")

        # Combine into a data frame
        landscape_data <- invisible(read.csv(complete_dataset))

        # Split the data into training and testing sets
        train_data <- landscape_data[1:n_train, ] # Training data based on first 2000 points
        test_data <- landscape_data[(n_train + 1):nrow(landscape_data), ] # Remaining 500 points for test data

        # Check for CS and induce CS if necessary
        # Assume train_data and test_data are data frames with the same features
        features <- setdiff(colnames(train_data), c("lat", "lon", "X", "set", "z", "x4", "x5", "x6"))
        # print(features)

        # Adjust significance level for multiple comparisons (number of predictors) using Bonferroni correction
        alpha <- 0.05 / length(features)

        for (feature_name in features) {
            ks_test <- ks.test(train_data[[feature_name]], test_data[[feature_name]])
            # print(ks_test)
            # Continue adjusting until KS test p-value is below alpha
            while (ks_test$p.value > alpha) {
                # Example adjustments: you can modify this logic as needed
                if (feature_name == "x1") {
                    test_data[[feature_name]] <- test_data[[feature_name]] + 1 # Shift mean
                } else if (feature_name == "x2") {
                    test_data[[feature_name]] <- test_data[[feature_name]] * sqrt(2) # Increase variance
                } else if (feature_name == "x3") {
                    test_data[[feature_name]] <- (test_data[[feature_name]] + 0.3) * sqrt(1.5) # Shift mean and increase variance
                }
                # Re-run KS test
                ks_test <- ks.test(train_data[[feature_name]], test_data[[feature_name]])
            }

            new_row <- data.frame(
                Range = range,
                Dataset = dataset,
                Predictor = feature_name,
                KS_Statistic = ks_test$statistic,
                P_Value = ks_test$p.value,
                stringsAsFactors = FALSE
            )

            results <- rbind(results, new_row)
        }
        error_sd <- 0.1
        # Recalculate z in the test data after covariate shift, adding noise
        test_data$z <- with(test_data, x1 + x2 + x3 + (x2 * x3) + rnorm(nrow(test_data), mean = 0, sd = error_sd))

        combined_data <- rbind(train_data, test_data)

        # Define file paths for saving the data
        spatial_file <- paste("simulated_data/spatial/sics/complete_", range, "_", dataset, ".csv", sep = "")
        non_spatial_train <- paste("simulated_data/non_spatial/sics/train_", range, "_", dataset, ".csv", sep = "")
        non_spatial_test <- paste("simulated_data/non_spatial/sics/test_", range, "_", dataset, ".csv", sep = "")

        # Write full dataset into CSV file
        write.csv(combined_data, spatial_file, row.names = FALSE)

        # Write train_data to a CSV file
        write.csv(train_data, non_spatial_train, row.names = FALSE)

        # Write test_data to a CSV file
        write.csv(test_data, non_spatial_test, row.names = FALSE)
    }
}

# Output result of KS test to ensure covariate shift
write.csv(results, "SICS_KS_Test.csv", row.names = FALSE)

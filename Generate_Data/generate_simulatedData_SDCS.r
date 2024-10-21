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
n_test <- 500 # Number of test points in region 1

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
        cat("   Simulating Dataset: ", dataset, "...\n")

        # Simulate training data points (lat_train = 0-50, lon_train = 0-25)
        lat_train <- runif(n_train, min = 0, max = 50)
        lat_test <- runif(n_test, min = 0, max = 50)

        if (range == 1) {
            # No Separation when range = 1
            lon_train <- runif(n_train, min = 0, max = 25)
            lon_test <- runif(n_test, min = 25, max = 50)
        } else if (range == 6) {
            # Separation of 5 units, 1 unit of correlation
            lon_train <- runif(n_train, min = 0, max = 22.5)
            lon_test <- runif(n_test, min = 27.5, max = 50)
        } else {
            # Separation of 11 units, 1 unit of correlation
            lon_train <- runif(n_train, min = 0, max = 19.5)
            lon_test <- runif(n_test, min = 30.5, max = 50)
        }

        # Combine lat and lon for train and test data
        lat <- c(lat_train, lat_test)
        lon <- c(lon_train, lon_test)

        # Create a spatial points data frame
        coords <- data.frame(lat = lat, lon = lon)
        coordinates(coords) <- ~ lon + lat

        # Simulate spatially correlated variables x1 - x6
        g <- gstat(formula = x1 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_predict, nmax = 30)
        trashOutput <- capture.output(simulated_x1 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x1 <- center_variable(simulated_x1, tolerance)

        g <- gstat(formula = x2 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_predict, nmax = 30)
        trashOutput <- capture.output(simulated_x2 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x2 <- center_variable(simulated_x2, tolerance)

        g <- gstat(formula = x3 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_predict, nmax = 30)
        trashOutput <- capture.output(simulated_x3 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x3 <- center_variable(simulated_x3, tolerance)

        g <- gstat(formula = x4 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_noise, nmax = 30)
        trashOutput <- capture.output(simulated_x4 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x4 <- center_variable(simulated_x4, tolerance)

        g <- gstat(formula = x5 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_noise, nmax = 30)
        trashOutput <- capture.output(simulated_x5 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x5 <- center_variable(simulated_x5, tolerance)

        g <- gstat(formula = x6 ~ 1, locations = ~ lon + lat, dummy = TRUE, beta = 0, model = vgm_model_noise, nmax = 30)
        trashOutput <- capture.output(simulated_x6 <- predict(g, newdata = coords, nsim = 1)$sim1)
        simulated_x6 <- center_variable(simulated_x6, tolerance)

        # Combine into a data frame
        landscape_data <- data.frame(lat = coords$lat, lon = coords$lon, x1 = simulated_x1, x2 = simulated_x2, x3 = simulated_x3, x4 = simulated_x4, x5 = simulated_x5, x6 = simulated_x6)


        # Define the standard deviation for the error term
        error_sd <- 0.1 # You can adjust this value as needed

        # Generate the target variable z based on the equation: y = X1 + X2 + X3 + (X2 * X3) + error
        landscape_data$z <- with(landscape_data, x1 + x2 + x3 + (x2 * x3) + rnorm(nrow(landscape_data), mean = 0, sd = error_sd))

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

        # Recalculate z in the test data after covariate shift, adding noise
        test_data$z <- with(test_data, x1 + x2 + x3 + (x2 * x3) + rnorm(nrow(test_data), mean = 0, sd = error_sd))

        combined_data <- rbind(train_data, test_data)

        # Define file paths for saving the data
        spatial_file <- paste("simulated_data/spatial/sdcs/complete_", range, "_", dataset, ".csv", sep = "")
        non_spatial_train <- paste("simulated_data/non_spatial/sdcs/train_", range, "_", dataset, ".csv", sep = "")
        non_spatial_test <- paste("simulated_data/non_spatial/sdcs/test_", range, "_", dataset, ".csv", sep = "")

        # Write full dataset into CSV file
        write.csv(combined_data, spatial_file, row.names = FALSE)

        # Write train_data to a CSV file
        write.csv(train_data, non_spatial_train, row.names = FALSE)

        # Write test_data to a CSV file
        write.csv(test_data, non_spatial_test, row.names = FALSE)
    }
}

write.csv(results, "SDCD_KS_Test.csv", row.names = FALSE)

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

for (range in range_array) {
    cat("Range: ", range, "...\n")
    # Create the variogram model
    vgm_model_predict <- vgm(psill = sill - nugget, model = "Gau", range = range, nugget = nugget)
    vgm_model_noise <- vgm(psill = sill - nugget, model = "Gau", range = range, nugget = nugget_noise)

    for (dataset in 1:number_of_Dataset) {
        cat("   Simulating Dataset: ", dataset, "...\n")

        # Simulate training data points (lat_train = 0-50, lon_train = 0-25)
        lat_train <- runif(n_train, min = 0, max = 50)
        lat_test <- runif(n_test, min = 70, max = 120)

        lon_train <- runif(n_train, min = 0, max = 50)
        lon_test <- runif(n_test, min = 70, max = 120)

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

        # Define file paths for saving the data
        spatial_file <- paste("simulated_data/spatial/si/complete_", range, "_", dataset, ".csv", sep = "")
        non_spatial_train <- paste("simulated_data/non_spatial/si/train_", range, "_", dataset, ".csv", sep = "")
        non_spatial_test <- paste("simulated_data/non_spatial/si/test_", range, "_", dataset, ".csv", sep = "")

        # Write full dataset into CSV file
        write.csv(landscape_data, spatial_file, row.names = FALSE)

        # Write train_data to a CSV file
        write.csv(train_data, non_spatial_train, row.names = FALSE)

        # Write test_data to a CSV file
        write.csv(test_data, non_spatial_test, row.names = FALSE)
    }
}

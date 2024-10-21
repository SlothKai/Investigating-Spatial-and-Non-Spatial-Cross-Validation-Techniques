# Define a custom summary function to include MAE, RMSE, and R²
custom_metrics <- function(data, lev = NULL, model = NULL) {
  # Extract observed and predicted values
  obs <- data$obs # Actual values
  pred <- data$pred # Predicted values
  # Calculate RMSE
  rmse_val <- sqrt(mean((obs - pred)^2))
  # Calculate MAE
  mae_val <- mean(abs(obs - pred))
  # Calculate R² (coefficient of determination)
  r2_val <- 1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2))
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape_val <- mean(abs((obs - pred) / obs)) * 100
  # Calculate Bias (mean error)
  bias_val <- mean(pred - obs)
  # Return the custom metrics as a named list
  metrics <- c(RMSE = rmse_val, MAE = mae_val, R2 = r2_val, MAPE = mape_val, Bias = bias_val)
  return(metrics)
}

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# Custom R^2 function
r_squared <- function(actual, predicted) {
  ss_res <- sum((actual - predicted)^2) # Residual sum of squares
  ss_tot <- sum((actual - mean(actual))^2) # Total sum of squares
  r2 <- 1 - (ss_res / ss_tot) # Coefficient of determination
  return(r2)
}

# Custom MAE function
mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}

# Custom Bias function
bias <- function(actual, predicted) {
  mean(actual - predicted)
}

# Custom MAPE function
mape <- function(actual, predicted) {
  mean(abs((actual - predicted) / actual)) * 100 # Multiply by 100 for percentage
}

# Centering function, ensures mean of variable = 0 or close to 0
center_variable <- function(variable, tolerance) {
  mean_value <- mean(variable)

  # Center until the mean is below the tolerance
  while (abs(mean_value) > tolerance) {
    variable <- variable - mean_value # Center the variable
    mean_value <- mean(variable) # Recalculate the mean
  }

  return(variable)
}

# Function to calculate the distance between points
min_distance <- function(lon, lat) {
  dist_matrix <- rdist(cbind(lon, lat)) # Pairwise distances
  dist_matrix[upper.tri(dist_matrix)] # Return only upper triangular part
}

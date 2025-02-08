# Investigating-Spatial-and-Non-Spatial-Cross-Validation-Techniques

This repository contains R scripts for a Final Year Project (FYP) titled: Investigating Spatial and Non-Spatial Cross Validation techniques.

## Overview

This repository explores the impact of different cross-validation (CV) techniques on predictive modeling, particularly in spatial data contexts. The goal is to compare various spatial and non-spatial CV methods to assess their effectiveness in estimating model performance.

## Features

- **Multiple Cross-Validation Techniques:**
  - Blocked Cross-Validation
  - Bootstrap Cross-Validation
  - Buffered Cross-Validation
  - Importance-Weighted Cross-Validation
  - Random K-Fold Cross-Validation
  - Spatial K-Fold Cross-Validation
- **Data Generation Scripts:**
  - Simulated data with different spatial characteristics
  - Evaluates the effects of covariate shift on model validation
- **Performance Metrics:**
  - RMSE, MAE, MAPE, and R-squared
- **Implemented in R with Parallel Processing Support**

## Repository Structure

```
FYP/
│── main_script.r                  # Main entry point for running experiments
│── general_functions.r            # Common utility functions
│
├── Cross_Validation_Technique/    # Various cross-validation methods
│   ├── BlockedCV.r
│   ├── BootstrapCV.r
│   ├── BufferedCV.r
│   ├── Importance_Weighted_CV.r
│   ├── Random_KFoldCV.r
│   ├── Spatial_KFoldCV.r
│
├── Generate_Data/                 # Scripts for generating simulated datasets
│   ├── generate_simulatedData_SD.r
│   ├── generate_simulatedData_SDCS.r
│   ├── generate_simulatedData_SI.r
│   ├── generate_simulatedData_SICS.r
```

## Setup Instructions

### Prerequisites

Ensure you have R installed with the following packages:

```r
install.packages(c("ranger", "caret", "spatialsample", "parallel", "doParallel","readr", "dplyr", "sf", "sp","akima","ggplot2","densratio","gstat"))
```

### Running the Code

1. Clone this repository:
   ```sh
   git clone <https://github.com/SlothKai/Investigating-Spatial-and-Non-Spatial-Cross-Validation-Techniques.git>
   cd Investigating-Spatial-and-Non-Spatial-Cross-Validation-Techniques
   ```
2. Run the main script:
   ```r
   source("main_script.r")
   ```

## Usage

- Modify `Generate_Data/` scripts to create custom datasets.
- Choose different cross-validation techniques by modifying `main_script.r`.
- Analyze model performance metrics to evaluate the effectiveness of each CV method.

## Potential Applications

- Improving predictive modeling in spatial datasets.
- Evaluating the impact of covariate shift on model performance.
- Comparing spatial and non-spatial cross-validation approaches for real-world applications.

## Future Improvements

- Implement additional spatial CV techniques.
- Extend the analysis to different machine learning models.
- Optimize parallel processing for faster execution.

## Contact

For any questions or contributions, feel free to reach out!

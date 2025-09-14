# tcv: Determining the Number of Factors in Poisson Factor Models via Thinning Cross-Validation

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
## Overview

The `tcv` package provides a robust method for a crucial model selection problem: **choosing the number of factors in Poisson factor models**. This is particularly important for high-dimensional count data, where traditional cross-validation methods often fail. The package's core contribution is the implementation of **Thinning Cross-Validation (TCV)**, a technique specifically designed for count data that preserves the underlying data structure. This makes it a valuable tool for researchers in fields like genomics, text mining, and computational social science.

The package is currently under review and will be available on [CRAN](https://cran.r-project.org/) soon.

## Installation

### From CRAN (stable version - coming soon)
Once accepted, you will be able to install the package from CRAN:
```r
install.packages("tcv")
# Load the package
library(tcv)
```

### From GitHub (development version)
You can install the development version of `tcv` from GitHub:
```r
# Install devtools if you haven't already
# install.packages("devtools")

# Install tcv from GitHub
devtools::install_github("Wangzhijingwzj/tcv")

# Load the package
library(tcv)
```

## Key Features

- **Thinning Cross-Validation**: Implements the TCV method for robust and reliable model selection with count data.
- **Specialized for Poisson Factor Models**: Tailored specifically to determine the latent dimensionality in Poisson factor models.
- **Automated Factor Selection**: Automates the process of testing a range of factor numbers to find the optimal one.
- **Robust Implementation**: Uses the `GFM` package for model fitting with a C++ based fallback mechanism for stability.

## Functions Overview

### Core Function

- `multiDT()`: The main function for performing K-fold Thinning Cross-Validation to select the optimal number of factors.

### Supporting Functions

- `chooseFacNumber_ratio()`: An alternative, fast method for factor number estimation based on singular value ratios.
- `add_identifiability()`: An internal function to apply identifiability constraints, ensuring unique and stable model solutions.

## Quick Start: Finding the Optimal Number of Factors

Here is a complete example of how to use the `multiDT` function to find the best number of factors for a simulated dataset.

### 1. Load the Package and Helper Functions
```r
library(tcv)
set.seed(123) # for reproducibility
```

### 2. Data Generation
In a real-world scenario, `x` would be your own count data matrix (e.g., a gene expression matrix). Here, we simulate data for demonstration.
```r
# Parameters for data generation
n <- 50 # Number of samples
p <- 30 # Number of features
true_q <- 2  # True number of factors

# Factor matrix (scores)
FF <- matrix(rnorm(n * true_q), nrow = n, ncol = true_q)
# Loading matrix
BB <- matrix(runif(p * true_q, min = -1, max = 1), nrow = p, ncol = true_q)
# Intercept term
a <- runif(p, min = 0, max = 1)

# Enforce identifiability for a unique generating model
# Note: add_identifiability is an exported helper function in the package
ident <- add_identifiability(FF, BB, a)
FF0 <- ident$H
BB0 <- ident$B
alpha <- ident$mu

# Calculate the mean matrix (lambda) with some noise
lambda <- exp(FF0 %*% t(BB0) + rep(1, n) %*% t(alpha) + matrix(rnorm(n*p, 0, 0.5), n, p))

# Generate the final count data matrix 'x'
x <- matrix(rpois(n * p, lambda = as.vector(lambda)), nrow = n, ncol = p)
```

### 3. Run Thinning Cross-Validation
We use `multiDT()` to test a range of factor numbers (from 1 to 4) using 2-fold cross-validation. In practice, `K=5` is a common choice.
```r
# Run multiDT to find the best number of factors
cv_results <- multiDT(x, K = 2, rmax = 4)

# The output contains the total cross-validation error for each tested factor number
print(cv_results$TCV)
```

### 4. Select the Optimal Number of Factors
The best number of factors is the one that minimizes the total cross-validation error (TCV).
```r
best_r <- which.min(cv_results$TCV)
cat("The optimal number of factors is:", best_r, "\n")
```

## Parameters for `multiDT`

- **`x`**: An n x p numeric matrix of count data.
- **`K`**: The number of folds for cross-validation (e.g., 5).
- **`rmax`**: The maximum number of factors to test (e.g., 8).

## Methodological Details

The `tcv` package algorithm works as follows:
1.  **Data Thinning**: The input count matrix `x` is probabilistically partitioned into `K` folds using the `countsplit` package.
2.  **Iterative Model Fitting**: For each fold `k`, a Poisson factor model is fit on the training data (all folds except `k`).
3.  **Error Calculation**: The predictive performance (negative log-likelihood) is calculated on the hold-out test data (fold `k`).
4.  **Aggregation**: The errors are summed across all `K` folds for each candidate number of factors. The number of factors with the lowest total error is chosen as optimal.

## Dependencies

- R (≥ 3.5.0)
- `GFM`
- `countsplit`
- `irlba`
- `stats`

## Bug Reports and Issues

Please report any bugs or issues on the [GitHub Issues page](https://github.com/Wangzhijingwzj/tcv/issues).

## References

The data thinning methodology is detailed in:
- Dharamshi, A. et al. (2025) <doi:10.1080/01621459.2024.2353948>
- Wang, Z. et al. (2025) <doi:10.1080/01621459.2025.2546577>

## License

This package is licensed under GPL-3.

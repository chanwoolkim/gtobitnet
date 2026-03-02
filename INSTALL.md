## Prerequisites

You need:
- R (version 4.0.0 or higher)
- A C++ compiler (for Rcpp)
  - **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
  - **Linux**: `sudo apt-get install build-essential` (or equivalent)
  - **Windows**: Install Rtools from https://cran.r-project.org/bin/windows/Rtools/

## Installation Steps

### 1. Install Required R Packages

Open R and run:

```r
install.packages(c("Rcpp", "devtools", "AER", "survival", "glmnet"))
```

### 2. Install gtobitnet from Local Source

Navigate to the package directory in R:

```r
# Set working directory to the package root
setwd("/path/to/gtobitnet")

# Generate Rcpp exports (required for C++ code)
library(Rcpp)
compileAttributes()

# Install the package
library(devtools)
install()
```

### 3. Load and Use

```r
library(gtobitnet)

# Test with two-sided censoring
set.seed(123)
n <- 100
p <- 10
x <- matrix(rnorm(n * p), n, p)
y <- pmax(pmin(x[,1] + rnorm(n), 1), 0)  # Censored at 0 and 1

# Fit with two-sided censoring
fit <- gtobitnet(x, y, left = 0, right = 1, lambda1 = 0.01)
print(fit)
```

### After making code changes:
If you do modify R or C++ code later:
1. Run `Rcpp::compileAttributes()` again
2. Run `devtools::install()` again

## Using the Package

Once installed, you can use it from any R session:

```r
library(gtobitnet)

# Two-sided censoring
fit <- gtobitnet(x, y, left = 0, right = 1, ...)

# One-sided (backward compatible)
fit <- gtobitnet(x, y, left = 0, right = 1e10, ...)
```


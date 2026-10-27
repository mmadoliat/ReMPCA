# Regular Data Class

The `rdClass` class denotes an element of regular data.

- Data objects are constructed with observations in rows and variables
  in columns.

- The smoothing parameter is automatically set to zero (no smoothing).

- Users can assign sparsity parameters only!

- No grid points are involved in `rdClass` objects.

## Usage

``` r
rdClass(data, Sparsity_parameter = 0)
```

## Arguments

- data:

  A matrix representing the data, with rows indicating observations and
  columns representing variables.

- Sparsity_parameter:

  - A fixed number representing the level of sparsity for columns, or a
    vector of numerical values that will undergo cross-validation (CV)
    to determine the optimal value.

  - For no sparsity, set it to 0.

  - If `NULL`, the sparsity parameter will be tuned automatically.

## Examples

``` r
# Example for Regular Data (rd)
rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
rd_object <- rdClass(data = rd_data,
                     Sparsity_parameter = 3)  # Custom sparsity parameter

# Display the created rd object
print(rd_object)
#> Regular Data (rdClass) Object
#> -----------------------------------
#> Dimensions: 10 x 10
#> Sparsity Parameter: 3 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1         V2          V3          V4          V5
#> [1,] -0.3572989 -1.0905736  0.07730312 -2.52850069 -0.73732753
#> [2,] -1.1468141  2.2846593 -0.29686864 -0.93590256  0.29066665
#> [3,] -0.5174205 -0.8856176 -1.18324224 -0.96723946 -0.88484957
#> [4,] -0.3621238  0.1111064  0.01129269  0.04748859  0.20800648
#> [5,]  2.3505543  3.8102767  0.99160104 -0.40373679 -0.04773017
print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter
#> [1] 3
```

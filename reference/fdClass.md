# Functional Data Class

The `fdClass` class represents an element of functional data.

- Data objects are constructed using matrices, where columns represent
  grid points (user-defined or NULL) and rows represent observations.

- Users can assign both smoothing and sparsity parameters.

- If the smoothing parameter is set to zero, no smoothing is applied.
  Otherwise, users can specify a fixed smoothing value or provide a
  vector of values, which will be optimized using generalized
  cross-validation (GCV).

## Usage

``` r
fdClass(data, argval = NULL, Smoothing_parameter = 0, Sparsity_parameter = 0)
```

## Arguments

- data:

  A matrix representing the data, with rows indicating observations and
  columns representing grid points.

- argval:

  A vector of grid points, it assigns grid points to the columns, with a
  length equal to the number of columns in the data. If `NULL`, grid
  points are automatically assigned from 0 to 1.

- Smoothing_parameter:

  : Smoothing parameter for columns. It can be:

  - A fixed number representing the smoothing parameter.

  - A vector of numerical values, which will undergo generalized
    cross-validation (GCV) to determine the optimal value.

  - Set to 0 for no smoothing.

  - If `NULL`, it analyzes a sequence of
    `2^seq(-30, 5, length.out = 10)` and attempts to tune it.

- Sparsity_parameter:

  - A fixed number representing the level of sparsity for columns, or a
    vector of numerical values that will undergo cross-validation (CV)
    to determine the optimal value.

  - For no sparsity, set it to 0.

  - If `NULL`, the sparsity parameter will be tuned automatically.

## Examples

``` r
# Example for Functional Data (fd)
fd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
fd_object <- fdClass(data = fd_data,
                     argval = seq(0, 1, length.out = 10),  # Grid points for columns
                     Smoothing_parameter = 0.5,  # Custom smoothing parameter
                     Sparsity_parameter = 2)  # Custom sparsity parameter

# Display the created fd object
print(fd_object)
#> Functional Data (fdClass) Object
#> -----------------------------------
#> Dimensions: 10 x 10
#> Smoothing Parameter: 0.5 
#> Sparsity Parameter: 2 
#> GridPoints_v: 0 0.1111111 0.2222222 ... 0.8888889 1 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1         V2        V3         V4         V5
#> [1,]  1.5468813  1.0842790 0.9790561 -0.5505134  0.1603770
#> [2,]  0.1789210  0.1228316 0.3294698  0.9212599 -3.2532201
#> [3,] -0.2825466 -0.4942494 1.8796357  2.1894560 -0.3299305
#> [4,] -0.7672988  0.6695882 2.0428990  0.9639948  0.4570898
#> [5,] -0.5764042 -0.3514789 1.3285764  0.7418185 -1.3041038
print(attr(fd_object, "GridPoints_v"))  # Display grid points for columns
#>  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
#>  [8] 0.7777778 0.8888889 1.0000000
print(attr(fd_object, "Smoothing_parameter"))  # Display smoothing parameter
#> [1] 0.5
```

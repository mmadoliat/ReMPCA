# Coerce an object of class 'rdClass' or 'imgClass' to class 'fdClass'.

Coerce an object of class 'rdClass' or 'imgClass' to class 'fdClass'.

## Usage

``` r
as.fdClass(
  x,
  Sparsity_parameter = NULL,
  Smoothing_parameter = NULL,
  argval = NULL
)
```

## Arguments

- x:

  An object of class 'rdClass' or 'imgClass'.

- Sparsity_parameter:

  Optional sparsity parameter to assign. If `NULL`, defaults are used.

- Smoothing_parameter:

  Optional smoothing parameter to assign. If `NULL`, defaults are used.

- argval:

  Optional vector of grid points for columns. If `NULL`, defaults are
  used.

## Value

An object of class 'fdClass' with user-specified or inherited
regularization parameters.

## Examples

``` r
img_object <- imgClass(
  image = list(matrix(rnorm(100), nrow = 10), matrix(rnorm(100), nrow = 10)),
  argval = NULL,
  Smoothing_parameter = NULL,
  Sparsity_parameter = 0
)
newfd <- as.fdClass(img_object)
newfd
#> Functional Data (fdClass) Object
#> -----------------------------------
#> Dimensions: 2 x 100
#> Smoothing Parameter: 9.313226e-10 1.379661e-08 2.043829e-07 ... 2.160119 32 
#> Sparsity Parameter: 0 
#> GridPoints_v: 0.01 0.02 0.03 ... 0.99 1 
#> -----------------------------------
#> First few rows and columns of the data:
#>             V1         V2         V3         V4          V5
#> [1,] -1.086909 0.35339850 -0.1740555 -0.5443194  2.00571858
#> [2,] -1.541403 0.09535401  0.7372137  0.8140203 -0.04163922
```

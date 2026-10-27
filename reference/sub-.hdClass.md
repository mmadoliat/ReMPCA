# Indexing Operator for hdClass

Enables subsetting of an `hdClass` object by variables.

## Usage

``` r
# S3 method for class 'hdClass'
x[k = NULL]
```

## Arguments

- x:

  An object of class `hdClass`.

- k:

  Variables indices. If `NULL`, all variables are included.

## Value

A new `hdClass` object with subsetted data and inherited attributes.

## Examples

``` r
hd <- hdClass(list(fdClass(matrix(rnorm(100), 10)), rdClass(matrix(rnorm(30), 10))))
hd[2]
#> Hybrid Data (hdClass) Object
#> ===================================
#> Dimensions           : 10 rows x 3 columns
#> Number of Variables  :  1 
#> Columns per Variable :  3 
#> -----------------------------------
#> Grid Points (Rows)   : 0.1, 0.2, 0.3 , ...,  0.9, 1 
#> Grid Points (Columns):
#>   - Variable 1: NULL
#> Smoothing Parameter (Row): 0 
#> Smoothing Parameters (Col):
#>   - Variable 1: 0 
#> Sparsity Parameter (Row): 0 
#> Sparsity Parameters (Col):
#>   - Variable 1: 0 
#> ===================================
#> First few rows and columns of the data:
#>               V1           V2         V3
#> [1,]  0.77208630  1.781524475  0.2418959
#> [2,] -0.14083939  0.711213964 -1.1327594
#> [3,]  0.39309393 -0.337691156  1.4899074
#> [4,]  0.22421857 -0.009148952 -0.2482471
#> [5,]  0.02354199 -0.125309208  0.1835837

```

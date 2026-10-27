# Indexing operator for rdClass

Enables subsetting of an `rdClass` object by rows and columns.

## Usage

``` r
# S3 method for class 'rdClass'
x[i = NULL, j = NULL]
```

## Arguments

- x:

  An object of class `rdClass`.

- i:

  Row indices (observations). If `NULL`, all rows are included.

- j:

  Column indices. If `NULL`, all columns are included.

## Value

A new `rdClass` object with subsetted data and inherited attributes.

## Examples

``` r
rd <- rdClass(matrix(rnorm(50), 10, 5))
rd[1:3, 2:4]
#> Regular Data (rdClass) Object
#> -----------------------------------
#> Dimensions: 3 x 3
#> Sparsity Parameter: 0 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1         V2         V3
#> [1,] -1.2083178 -0.9810758  2.0737836
#> [2,] -0.4393227 -0.1392214 -0.7788302
#> [3,] -0.3755807  2.1854379  1.1315325
```

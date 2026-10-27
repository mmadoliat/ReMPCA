# Indexing operator for fdClass

Enables subsetting of an `fdClass` object by rows and columns.

## Usage

``` r
# S3 method for class 'fdClass'
x[i = NULL, j = NULL]
```

## Arguments

- x:

  An object of class `fdClass`.

- i:

  Row indices (observations). If `NULL`, all rows are included.

- j:

  Column indices (grid points). If `NULL`, all columns are included.

## Value

A new `fdClass` object with subsetted data and inherited attributes.

## Examples

``` r
fd <- fdClass(matrix(rnorm(200), 10, 20))
fd[1:5, 1:10]
#> Functional Data (fdClass) Object
#> -----------------------------------
#> Dimensions: 5 x 10
#> Smoothing Parameter: 0 
#> Sparsity Parameter: 0 
#> GridPoints_v: 0.05 0.1 0.15 ... 0.45 0.5 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1          V2          V3         V4          V5
#> [1,] -0.7598457 -0.08698633 -1.08690882  0.3533985 -0.17405549
#> [2,]  1.1489591  1.38228401 -1.82608301  0.7268137  0.96129056
#> [3,] -0.8424763  0.16849016  0.99528181  0.6682610  0.29382666
#> [4,]  0.3914133  0.82319095 -0.01186178 -2.4243173  0.08099936
#> [5,]  0.8913772 -0.22089460 -0.59962839 -0.2353574  0.18366184
```

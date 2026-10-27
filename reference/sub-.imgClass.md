# Indexing operator for imgClass

Enables subsetting of an `imgClass` object by rows.

## Usage

``` r
# S3 method for class 'imgClass'
x[i = NULL]
```

## Arguments

- x:

  An object of class `imgClass`.

- i:

  Row indices (Images). If `NULL`, all images are included.

## Value

A new `imgClass` object with subsetted data and inherited attributes.

## Examples

``` r
imgs <- lapply(1:3, function(k) matrix(rnorm(16), 4, 4))
img <- imgClass(imgs)
img[1:2]
#> Image Data (imgClass) Object
#> -----------------------------------
#> Dimensions: 2 x 16
#> Smoothing Parameter: 0 
#> Sparsity Parameter: 0 
#> GridPoints_v: None
#> Each image size: 4 x 4
#> -----------------------------------
#> First few rows and columns of the data:
#>               V1        V2         V3         V4         V5
#> [1,]  0.07383053 -1.402906 -0.1048521 -0.7185115  0.7059456
#> [2,] -0.21538448  1.542326  0.1798052  1.3663086 -2.4719617
```

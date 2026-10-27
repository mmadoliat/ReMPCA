# Image Data Class

The `imgClass` class is designed to handle images:

- Takes a list of images (matrices of equal dimension); a single image
  must be wrapped in [`list()`](https://rdrr.io/r/base/list.html).

- All images are vectorized (row-wise) and combined into one matrix,
  with one row per image.

- The resulting object behaves like a `fdClass` if `Smoothing_parameter`
  is not zero, otherwise as a `rdClass`.

- Grid points (optional) are assigned to the columns (i.e., pixel
  locations).

## Usage

``` r
imgClass(image, Smoothing_parameter = 0, Sparsity_parameter = 0, argval = NULL)
```

## Arguments

- image:

  A list of matrices, one per image (observation). All images must have
  the same dimension. Wrap a single image as `list(img)`.

- Smoothing_parameter:

  - A numeric scalar or vector controlling the level of smoothing (as in
    `fdClass`).

    - Set to 0 for no smoothing (produces an `rdClass` object).

    - If NULL, smoothing is tuned from a default sequence.

- Sparsity_parameter:

  - A fixed number or vector controlling the level of sparsity on pixel
    columns.

    - If NULL, it is tuned automatically based on image size.

- argval:

  Optional numeric vector for grid points over the pixels (columns).

## Value

An object of class `fdClass` or `rdClass`, depending on the
`Smoothing_parameter`.

## Examples

``` r
img1 <- matrix(rnorm(64), 8, 8)
img2 <- matrix(rnorm(64), 8, 8)

# Multiple images
rd_img <- imgClass(list(img1, img2), Smoothing_parameter = 0, Sparsity_parameter = NULL)
```

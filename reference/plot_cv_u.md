# Plot CV Error with 1-SE Rule for u

For each principal component (PC), plots cross-validation (CV) error
versus the sparsity tuning parameter (`gamma`) for the u-direction. The
red dashed line indicates the 1-SE rule threshold, and the selected
optimal gamma is highlighted. Optionally, standard error bars (+/- 1 SE)
can be displayed.

## Usage

``` r
plot_cv_u(ReMPCA_obj, show_se = TRUE, ...)
```

## Arguments

- ReMPCA_obj:

  Output list from the ReMPCA routine.

- show_se:

  Logical. If TRUE, adds standard error bars to each point. Default is
  TRUE.

- ...:

  Additional plotting options passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

No return value, called for its side effect (plotting base R plots
side-by-side for each PC).

## Examples

``` r
set.seed(1)
u <- sin(seq(0, pi, length.out = 30))
X1 <- outer(u, cos(seq(0, pi, length.out = 12))) + matrix(rnorm(360, sd = 0.2), 30)
X2 <- outer(u, c(1, -1, 0, 0.5)) + matrix(rnorm(120, sd = 0.2), 30)
hd <- hdClass(list(fdClass(X1, Smoothing_parameter = c(0.01, 0.1, 1)),
                   rdClass(X2, Sparsity_parameter = 0:2)),
              Smoothing_parameter = c(0.001, 0.01), Sparsity_parameter = 0:2)
fit <- ReMPCA(hd, num_pcs = 1, nfolds_u = 2, nfolds_v = c(2, 2),
              cv.pick = "min", verbose = FALSE)
plot_cv_u(fit)
```

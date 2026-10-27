# Plot GCV Scores for U

This function visualizes the GCV scores for the left singular vectors
(U) across different regularization values.

## Usage

``` r
plot_gcv_u(ReMPCA_obj, ...)
```

## Arguments

- ReMPCA_obj:

  Output list from ReMPCA.

- ...:

  Additional plotting parameters.

## Value

No return value, called for its side effect (plotting).

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
plot_gcv_u(fit)
```

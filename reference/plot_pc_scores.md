# Plot Principal Component Scores from ReMPCA Object

This function visualizes the principal component scores from a ReMPCA
object. Each principal component (i.e., each column in `PCScores`) is
plotted in a separate dot plot, displaying the score values across
observations.

## Usage

``` r
plot_pc_scores(ReMPCA_obj, ...)
```

## Arguments

- ReMPCA_obj:

  A list-like ReMPCA object, expected to contain a component named
  `PCScores` which is a matrix or data frame where each column
  corresponds to the scores of one principal component.

- ...:

  Additional plotting options.

## Value

No return value. The function produces a series of dot plots, one for
each principal component.

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
plot_pc_scores(fit)
```

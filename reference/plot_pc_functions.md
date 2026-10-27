# Plot Principal Component Functions

This function visualizes the estimated principal component functions (v)
from a ReMPCA object. Each component is plotted across its variables
using either a line plot (for functional variables) or a dot plot (for
regular variables).

## Usage

``` r
plot_pc_functions(ReMPCA_obj, ...)
```

## Arguments

- ReMPCA_obj:

  A ReMPCA object containing:

  `PCFunctions`

  :   A list of length equal to number of PCs. Each element is a list of
      PC functions for each variable.

  `variable_types`

  :   A character vector indicating type of each variable: either `"hd"`
      (functional) or `"rd"` (regular).

- ...:

  Additional plotting options.

## Value

No return value. This function is called for its side effect of
plotting.

## Details

The function arranges plots in a matrix layout with rows corresponding
to variables and columns to principal components. A light gray
horizontal line at 0 is added for reference unless the minimum value in
the plot is at least 5.

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
plot_pc_functions(fit)
```

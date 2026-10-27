# Changelog

## ReMPCA 0.1.0

First CRAN release.

### Bug fixes

- Row-sparsity cross-validation now predicts the held-out columns with
  the least-squares loadings `v = X'u / ||u||^2`, as in the manuscript.
  Previously the scores were not normalized, which mis-scaled the CV
  error and led to over-sparse row scores whenever several row-sparsity
  levels were compared.
- [`ReMPCA()`](https://mobinapourmoshir.github.io/ReMPCA/reference/ReMPCA.md)
  with `num_pcs > 1` now deflates the data by the estimated regularized
  component (Algorithm 1 of the accompanying paper) instead of by the
  leading unregularized singular component. Unpenalized fits are
  unchanged.
- `VarianceExplained` is now computed from the estimated (regularized)
  loadings rather than from the unregularized singular vectors.
- The one-standard-error rule for loading sparsity now uses the standard
  error of the K-fold CV mean (it previously used the standard deviation
  of the fold errors), and the row-sparsity tuning now honours `cv.pick`
  for both tuning orders (it was hard-coded to `"min"` when
  `tuning_order = "Sparsity"`).
- Cross-validation no longer resets the user’s random number generator
  ([`set.seed()`](https://rdrr.io/r/base/Random.html) inside the
  package); folds still use a fixed internal seed, so results are
  unchanged and reproducible. Empty folds (too many folds for the data
  size) no longer produce `NA` indices.
- `$` methods for all classes gave an error inside the error message
  ([`sprintf()`](https://rdrr.io/r/base/sprintf.html) was missing its
  argument).
- `[.hdClass` returned a data frame labelled as `hdClass`; it now
  returns a valid `hdClass` matrix. `[.imgClass` collapsed all selected
  images into one row; `*.imgClass` returned a list-matrix;
  [`plot.imgClass()`](https://mobinapourmoshir.github.io/ReMPCA/reference/plot.imgClass.md)
  displayed images scrambled (they are vectorized row-wise but were
  rebuilt column-wise) and no longer waits for input in non-interactive
  sessions. `[.fdClass` and `[.rdClass` keep matrix structure and subset
  the grid points with the columns.
- Internal fixes to roughness-matrix bookkeeping (`opt_alpha_Omega_*`)
  that did not affect fitted values.
- Unknown `sparse_tuning_type`, `smoothness_type`, `tuning_order`, or
  `cv.pick` values now give an informative error.

### Other changes

- New argument `verbose` in
  [`ReMPCA()`](https://mobinapourmoshir.github.io/ReMPCA/reference/ReMPCA.md);
  progress messages can be suppressed.
- Plot functions restore the user’s graphical parameters on exit.
- With `parallel = TRUE`, the internal cluster uses at most
  `min(6, detectCores())` workers.
- The `energydata` data set now contains an object named `energydata`
  (previously `appliance_dat`) and is xz-compressed;
  [`read.data()`](https://mobinapourmoshir.github.io/ReMPCA/reference/read.data.md)
  can load it.
- Documentation: return values and runnable examples for all exported
  functions, references, and a test suite.
- Reproducibility scripts for the bike-sharing and appliances-energy
  applications in `inst/reproducibility` now use the packaged data, can
  be run from any working directory, and save the manuscript figures.

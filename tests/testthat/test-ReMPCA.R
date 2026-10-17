make_hd <- function(seed = 1, sm = 0, sp = 0) {
  set.seed(seed)
  n <- 30
  u <- sin(seq(0, pi, length.out = n))
  X1 <- outer(u, cos(seq(0, pi, length.out = 12))) + matrix(rnorm(n * 12, sd = 0.3), n)
  X2 <- outer(u, c(1, -1, 0, 0.5)) + matrix(rnorm(n * 4, sd = 0.3), n)
  hdClass(list(fdClass(X1, Smoothing_parameter = sm), rdClass(X2, Sparsity_parameter = sp)))
}

test_that("unpenalized ReMPCA reproduces the leading singular vectors", {
  hd <- make_hd()
  fit <- ReMPCA(hd, num_pcs = 2, weights = 0, thresh = 1e-14, maxit = 5000,
                verbose = FALSE)
  X <- scale(as.matrix(hd$matrix), center = TRUE, scale = FALSE)
  s <- svd(X)
  v1 <- unlist(fit$PCFunctions[[1]]); v2 <- unlist(fit$PCFunctions[[2]])
  expect_equal(abs(sum(v1 * s$v[, 1])), 1, tolerance = 1e-6)
  expect_equal(abs(sum(v2 * s$v[, 2])), 1, tolerance = 1e-4)
  expect_equal(fit$VarianceExplained, s$d[1:2]^2 / sum(s$d^2), tolerance = 1e-6)
})

test_that("ReMPCA does not change the user's random number stream", {
  hd <- make_hd(sp = 0:2)
  set.seed(42); a <- runif(3)
  set.seed(42); invisible(ReMPCA(hd, num_pcs = 1, nfolds_u = 2, nfolds_v = c(2, 2),
                                 verbose = FALSE)); b <- runif(3)
  expect_identical(a, b)
})

test_that("ReMPCA is reproducible and returns the documented components", {
  hd <- make_hd(sm = c(0.01, 0.1, 1), sp = 0:2)
  f1 <- ReMPCA(hd, num_pcs = 2, nfolds_u = 2, nfolds_v = c(2, 2), thresh = 1e-6,
               maxit = 500, verbose = FALSE)
  f2 <- ReMPCA(hd, num_pcs = 2, nfolds_u = 2, nfolds_v = c(2, 2), thresh = 1e-6,
               maxit = 500, verbose = FALSE)
  expect_identical(f1, f2)
  expect_named(f1, c("ReconstructedData", "PCFunctions", "PCScores", "OptimalAlphaV",
                     "OptimalAlphaU", "OptimalGammaV", "OptimalGammaU", "GCVResultsV",
                     "GCVResultsU", "CVResultsV", "CVResultsU", "variable_types",
                     "VarianceExplained"))
  expect_equal(dim(f1$PCScores), c(30, 2))
  expect_length(f1$VarianceExplained, 2)
})

test_that("invalid options are rejected", {
  hd <- make_hd()
  expect_error(ReMPCA(hd, sparse_tuning_type = "lasso", verbose = FALSE))
  expect_error(ReMPCA(hd, tuning_order = "both", verbose = FALSE))
  expect_error(ReMPCA(hd, cv.pick = "max", verbose = FALSE))
  expect_error(ReMPCA(unclass(hd), verbose = FALSE), "hdClass")
})

test_that("verbose = FALSE is silent", {
  hd <- make_hd()
  expect_silent(ReMPCA(hd, num_pcs = 1, verbose = FALSE))
})

test_that("$ methods give an informative error for unknown fields", {
  fd <- fdClass(matrix(rnorm(20), 4))
  rd <- rdClass(matrix(rnorm(20), 4))
  hd <- hdClass(list(fd, rd))
  img <- imgClass(list(matrix(1:4, 2), matrix(5:8, 2)))
  expect_error(fd$foo, "Unknown field 'foo'")
  expect_error(rd$foo, "Unknown field 'foo'")
  expect_error(hd$foo, "Unknown field 'foo'")
  expect_error(img$foo, "Unknown field 'foo'")
  expect_s3_class(hd$matrix, "data.frame")
})

test_that("[.hdClass returns a valid hdClass object", {
  X1 <- matrix(rnorm(50), 10); X2 <- matrix(rnorm(30), 10)
  hd <- hdClass(list(fdClass(X1), rdClass(X2)))
  sub <- hd[2]
  expect_s3_class(sub, "hdClass")
  expect_true(is.matrix(unclass(sub)))
  expect_equal(dim(sub), c(10, 3))
  expect_equal(attr(sub, "n_var"), 1)
  expect_equal(unname(unclass(sub)[, 1:3]), unname(X2))
  expect_equal(attr(sub, "variable_types"), "rd")
})

test_that("[.fdClass subsets grid points together with columns", {
  fd <- fdClass(matrix(rnorm(200), 10, 20), argval = 1:20)
  sub <- expect_silent(fd[1:5, 3:8])
  expect_equal(dim(sub), c(5, 6))
  expect_equal(attr(sub, "GridPoints_v"), 3:8)
})

test_that("imgClass vectorizes row-wise and subsetting/scaling keep the structure", {
  img1 <- matrix(1:6, 2, 3); img2 <- matrix(7:12, 2, 3)
  img <- imgClass(list(img1, img2))
  expect_equal(as.numeric(unclass(img)[1, ]), c(t(img1)))
  sub <- img[2]
  expect_s3_class(sub, "imgClass")
  expect_equal(dim(sub), c(1, 6))
  expect_equal(attr(sub, "nrow"), 2)
  dbl <- 2 * img
  expect_s3_class(dbl, "imgClass")
  expect_true(is.numeric(unclass(dbl)))
  expect_equal(as.numeric(unclass(dbl)), 2 * as.numeric(unclass(img)))
})

test_that("read.data loads all packaged data sets", {
  expect_s3_class(read.data("bike_day"), "data.frame")
  expect_s3_class(read.data("bike_hour"), "data.frame")
  expect_equal(dim(read.data("energydata")), c(19735L, 29L))
})

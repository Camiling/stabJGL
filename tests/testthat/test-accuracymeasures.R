library(stabJGL)

context("test-accuracymeasures.R")


test_that("Test confusion.matrix gives valid results", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  confmat <- stabJGL::confusion.matrix(adj.mat, adj.mat.2)

  # Test that the results are valid
  expect_true(confmat[1, 1] >= 0) # Value larger or equal to 0.
  expect_true(confmat[1, 2] >= 0) # Value larger or equal to 0.
  expect_true(confmat[2, 1] >= 0) # Value larger or equal to 0.
  expect_true(confmat[2, 2] >= 0) # Value larger or equal to 0.
})

test_that("confusion.matrix results in a 2 times 2 matrix", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  confmat <- stabJGL::confusion.matrix(adj.mat, adj.mat.2)

  # Tests -----------
  testthat::expect_equal(dim(confmat), c(2, 2)) # Check that result has correct dimension.
})

test_that("confusion.matrix gives right answer for identical matrices", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  confmat <- stabJGL::confusion.matrix(adj.mat, adj.mat.2)

  # Test that results are correct for identical precision matrices.
  expect_equal(stabJGL::confusion.matrix(adj.mat, adj.mat), matrix(c((sum(adj.mat != 0) - p) / 2, 0, 0, sum(adj.mat == 0) / 2), ncol = 2, byrow = T))
})

test_that("confusion.matrix sums up rows and cols correctly", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  confmat <- stabJGL::confusion.matrix(adj.mat, adj.mat.2)

  # Test that sum of columns and rows are correct.
  expect_equal(confmat[1, 1] + confmat[2, 1], (sum(adj.mat != 0) - p) / 2) # Number of edges in adj.mat
  expect_equal(confmat[1, 1] + confmat[1, 2], (sum(adj.mat.2 != 0) - p) / 2) # Number of edges in adj.mat.2
  expect_equal(confmat[1, 2] + confmat[2, 2], sum(adj.mat == 0) / 2) # Number of 'no edges' in adj.mat
  expect_equal(confmat[2, 1] + confmat[2, 2], sum(adj.mat.2 == 0) / 2) # Number of 'no edges' in adj.mat
})

test_that("confusion.matrix throws error", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.

  # Test errors
  expect_error(stabJGL::confusion.matrix(adj.mat, prec.mat.2)) # Not providing adjacency matric.
  expect_error(stabJGL::confusion.matrix(adj.mat, adj.mat[1:(p - 1), 1:(p - 1)])) # Different dimensions.
})

test_that("Test precision", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  res <- stabJGL::precision(adj.mat, adj.mat.2)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that the results are valid
  expect_true(res >= 0) # Value larger or equal to 0.
  expect_true(res <= 1) # Value smaller or equal to 0.

  # Test special cases.
  expect_equal(stabJGL::precision(adj.mat, adj.mat), 1) # Test that results are correct for identical adjacency matrices.
  expect_equal(stabJGL::precision(adj.mat, matrix(0, p, p)), 1) # Test that an an empty predicting graph gives precision 1.
  expect_equal(stabJGL::precision(matrix(0, p, p), adj.mat), 1) # Test that an an empty graph to predict gives precision 1.



  # Test errors
  expect_error(stabJGL::precision(adj.mat, prec.mat2)) # Not providing adjacency matric.
  expect_error(stabJGL::precision(adj.mat, adj.mat[1:(p - 1), 1:(p - 1)])) # Different dimensions.
})

test_that("Test recall", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  adj.mat.2 <- abs(prec.mat.2) >= 1e-8 # Avoid rounding errors.
  res <- stabJGL::recall(adj.mat, adj.mat.2)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that the results are valid
  expect_true(res >= 0) # Value larger or equal to 0.
  expect_true(res <= 1) # Value smaller or equal to 0.

  # Test special cases.
  expect_equal(stabJGL::recall(adj.mat, adj.mat), 1) # Test that results are correct for identical adjacency matrices.
  expect_equal(stabJGL::recall(adj.mat, matrix(0, p, p)), 0) # Test that an an empty predicting graph gives recall 0.
  expect_equal(stabJGL::recall(matrix(0, p, p), adj.mat), 1) # Test that an an empty graph to predict gives recall 1.

  # Test errors
  expect_error(stabJGL::recall(adj.mat, prec.mat2)) # Not providing adjacency matric.
  expect_error(stabJGL::recall(adj.mat, adj.mat[1:(p - 1), 1:(p - 1)])) # Different dimensions.
})

test_that("Test matrix.distance", {

  # Example: Two unrelated graphs
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  dat.2 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.2 <- dat.2$omega # true precision matrix
  res <- stabJGL::matrix.distance(prec.mat, prec.mat.2)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that the results are valid
  expect_true(res >= 0) # Value larger or equal to 0.
  expect_true(res <= 1) # Value smaller or equal to 0.

  # Test special cases.
  expect_equal(stabJGL::matrix.distance(prec.mat, prec.mat), 0) # Test that results are correct for identical precision matrices.
  empty.graph <- matrix(0, p, p)
  diag(empty.graph) <- 1
  expect_equal(stabJGL::matrix.distance(prec.mat, empty.graph), 1) # One empty graph gives maximum matrix distance.
  expect_equal(stabJGL::matrix.distance(empty.graph, empty.graph), 0) # Empty graph has matrix distance 0 to itself.

  # Test errors
  expect_error(stabJGL::matrix.distance(prec.mat, matrix(0, p, p))) # One precision matrix with diagonal elements
  expect_error(stabJGL::matrix.distance(prec.mat, prec.mat[1:(p - 1), 1:(p - 1)])) # Different dimensions.
})

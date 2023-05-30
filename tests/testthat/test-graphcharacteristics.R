library(stabJGL)

context("test-graphcharacteristics.R")

test_that("Test sparsity", {

  # Example
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  adj.mat <- abs(prec.mat) >= 1e-8 # Avoid rounding errors.
  res <- stabJGL::sparsity(adj.mat)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that the results are valid
  expect_true(res >= 0) # Value larger or equal to 0.
  expect_true(res <= 1) # Value smaller or equal to 0.

  # Test special cases.
  expect_equal(stabJGL::sparsity(matrix(0, p, p)), 0) # Test that results are correct for empty graphs
  expect_equal(stabJGL::sparsity(matrix(1, p, p)), 1) # Test that results are correct for full graphs

  # Test errors
  expect_error(stabJGL::sparsity(adj.mat[1:p, 1:(p - 1)])) # Different numbers of columns and rows.
})

test_that("Test gaussianloglik", {

  # Example
  set.seed(123)
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  cov.mat <- cov(dat$data)
  prec.mat <- dat$omega # true precision matrix
  res <- stabJGL::gaussianloglik(cov.mat, prec.mat, n)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that log likelihood is better for true precision matrix than a wrong one.
  dat.new <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.wrong <- dat.new$omega
  expect_true(res > stabJGL::gaussianloglik(cov.mat, prec.mat.wrong, n))

  # Test errors
  expect_error(stabJGL::gaussianloglik(cov.mat, prec.mat[1:(p - 1), 1:(p - 1)], n)) # Different dimensions
  expect_error(stabJGL::gaussianloglik(cov.mat, prec.mat, n = 0)) # Unvalid number of observations
  cov.unsym <- cov.mat
  cov.unsym[5, 8] <- 0.3
  expect_error(stabJGL::gaussianloglik(cov.unsym, prec.mat, n)) # Unsymmetric sample covariance matrix.
  prec.mat.notpos <- prec.mat
  prec.mat.notpos[which(abs(prec.mat.notpos) < 1e-7)] <- 1.2 # No zero elements
  expect_error(stabJGL::gaussianloglik(cov.mat, prec.mat.notpos, n)) # Precision matric not positive definite.
})

test_that("Test gaussianAIC", {

  # Example
  set.seed(123)
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  cov.mat <- cov(dat$data)
  prec.mat <- dat$omega # true precision matrix
  res <- stabJGL::gaussianAIC(cov.mat, prec.mat, n)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that AIC is better for true precision matrix than a wrong one.
  dat.new <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.wrong <- dat.new$omega
  expect_true(res < stabJGL::gaussianAIC(cov.mat, prec.mat.wrong, n))

  # Test errors
  expect_error(stabJGL::gaussianAIC(cov.mat, prec.mat[1:(p - 1), 1:(p - 1)], n)) # Different dimensions
  expect_error(stabJGL::gaussianAIC(cov.mat, prec.mat, n = 0)) # Unvalid number of observations
  cov.unsym <- cov.mat
  cov.unsym[5, 8] <- 0.3
  expect_error(stabJGL::gaussianAIC(cov.unsym, prec.mat, n)) # Unsymmetric sample covariance matrix.
  prec.mat.notpos <- prec.mat
  prec.mat.notpos[which(abs(prec.mat.notpos) < 1e-7)] <- 1.2 # No zero elements
  expect_error(stabJGL::gaussianAIC(cov.mat, prec.mat.notpos, n)) # Precision matric not positive definite.
})

test_that("Test eBIC", {

  # Example
  set.seed(123)
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  cov.mat <- cov(dat$data)
  prec.mat <- dat$omega # true precision matrix
  res <- stabJGL::eBIC(cov.mat, prec.mat, n)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that eBIC is better for true precision matrix than a wrong one.
  dat.new <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.wrong <- dat.new$omega
  expect_true(res < stabJGL::eBIC(cov.mat, prec.mat.wrong, n))

  # Test default argument
  expect_equal(res, stabJGL::eBIC(cov.mat, prec.mat, n, gamma = 0))

  # Test that larger gamma gives larger eBIC score
  expect_true(res < stabJGL::eBIC(cov.mat, prec.mat, n, gamma = 0.5))

  # Test errors
  expect_error(stabJGL::eBIC(cov.mat, prec.mat[1:(p - 1), 1:(p - 1)], n)) # Different dimensions
  expect_error(stabJGL::eBIC(cov.mat, prec.mat, n = 0)) # Unvalid number of observations
  expect_error(stabJGL::eBIC(cov.mat, prec.mat, n, gamma = -1)) # Negative gamma
  cov.unsym <- cov.mat
  cov.unsym[5, 8] <- 0.3
  expect_error(stabJGL::eBIC(cov.unsym, prec.mat, n)) # Unsymmetric sample covariance matrix.
  prec.mat.notpos <- prec.mat
  prec.mat.notpos[which(abs(prec.mat.notpos) < 1e-7)] <- 1.2 # No zero elements
  expect_error(stabJGL::eBIC(cov.mat, prec.mat.notpos, n)) # Precision matrix not positive definite.
})

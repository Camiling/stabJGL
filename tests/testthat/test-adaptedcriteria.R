library(stabJGL)

context("test-adaptedcriteria.R")


test_that("Test eBIC_adapted", {

  # Example
  set.seed(123)
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  cov.list = list(cov(dat$data), cov(dat$data))
  prec.mat=dat$omega
  cov.mat = cov(dat$data)
  prec.mat.list <- list(dat$omega, dat$omega) # true precision matrix
  n.list = c(n,n)
  res <- stabJGL::eBIC_adapted(prec.mat.list, cov.list, n.list)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that adapted eBIC is better for true precision matrix than a wrong one.
  dat.new <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.wrong <- dat.new$omega
  expect_true(res < stabJGL::eBIC_adapted(list(prec.mat.wrong,prec.mat.wrong),cov.list, n.list))

  # Test default argument
  expect_equal(res, stabJGL::eBIC_adapted(prec.mat.list,cov.list, n.list, gamma = 0.2))

  # Test that larger gamma gives larger eBIC score
  expect_true(res < stabJGL::eBIC_adapted(prec.mat.list,cov.list, n.list, gamma = 0.8))

  # Test errors
  prec.diff = prec.mat[1:(p - 1), 1:(p - 1)]
  expect_error(stabJGL::eBIC_adapted(list(prec.diff,prec.diff) , cov.list, n.list)) # Different dimensions
  expect_error(stabJGL::eBIC_adapted(prec.mat.list, cov.list, n.vals = c(0,0))) # Unvalid number of observations
  expect_error(stabJGL::eBIC_adapted(prec.mat.list, cov.list, n.list, gamma = -1)) # Negative gamma
  cov.unsym <- cov.mat
  cov.unsym[5, 8] <- 0.3
  expect_error(stabJGL::eBIC_adapted(prec.mat.list, list(cov.unsym, cov.unsym), n.list)) # Unsymmetric sample covariance matrix.
  prec.mat.notpos <- prec.mat
  prec.mat.notpos[which(abs(prec.mat.notpos) < 1e-7)] <- 1.2 # No zero elements
  expect_error(stabJGL::eBIC_adapted(list(prec.mat.notpos, prec.mat.notpos),cov.list, n.list)) # Precision matrix not positive definite.
  expect_error(stabJGL::eBIC_adapted(list(prec.mat,prec.mat,prec.mat) , cov.list, n.list)) # Different number of elements in cov list and prec mat list

})


test_that("Test AIC_adapted", {

  # Example
  set.seed(123)
  n <- 30
  p <- 10
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  cov.list = list(cov(dat$data), cov(dat$data))
  prec.mat=dat$omega
  cov.mat = cov(dat$data)
  prec.mat.list <- list(dat$omega, dat$omega) # true precision matrix
  n.list = c(n,n)
  res <- stabJGL::AIC_adapted(prec.mat.list, cov.list, n.list)

  # Tests -----------
  expect_equal(class(res), "numeric") # Test that a numeric is returned.

  # Test that adapted eBIC is better for true precision matrix than a wrong one.
  dat.new <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat.wrong <- dat.new$omega
  expect_true(res < stabJGL::AIC_adapted(list(prec.mat.wrong,prec.mat.wrong),cov.list, n.list))

  # Test errors
  prec.diff = prec.mat[1:(p - 1), 1:(p - 1)]
  expect_error(stabJGL::AIC_adapted(list(prec.diff,prec.diff) , cov.list, n.list)) # Different dimensions
  expect_error(stabJGL::AIC_adapted(prec.mat.list, cov.list, n.vals = c(0,0))) # Unvalid number of observations
  expect_error(stabJGL::AIC_adapted(prec.mat.list, cov.list, n.list, gamma = -1)) # Negative gamma
  cov.unsym <- cov.mat
  cov.unsym[5, 8] <- 0.3
  expect_error(stabJGL::AIC_adapted(prec.mat.list, list(cov.unsym, cov.unsym), n.list)) # Unsymmetric sample covariance matrix.
  prec.mat.notpos <- prec.mat
  prec.mat.notpos[which(abs(prec.mat.notpos) < 1e-7)] <- 1.2 # No zero elements
  expect_error(stabJGL::AIC_adapted(list(prec.mat.notpos, prec.mat.notpos),cov.list, n.list)) # Precision matrix not positive definite.
  expect_error(stabJGL::AIC_adapted(list(prec.mat,prec.mat,prec.mat) , cov.list, n.list)) # Different number of elements in cov list and prec mat list

})

library(JoStARS)

context("test-JoStARS.R")

test_that("JoStARS returns objects of correct classes", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  res <- JoStARS(dat.list,verbose = F)

  # Tests -----------
  # Test that the classes of results are correct.
  list.names <- c('opt.fit','opt.ebic', 'opt.sparsities', 'opt.lambda1', 'opt.lambda2','lambda1s', 'lambda2s', 'ebic.vals', 'opt.fit.lambda1', 'opt.sparsities.lambda1','total.variability', 'variability')
  expect_equal(class(res), "list")
  expect_true(mean(names(res) %in% list.names) == 1)
  expect_equal(class(res$opt.fit), "list")
  expect_equal(class(res$opt.ebic), "numeric")
  expect_equal(class(res$opt.sparsities), "numeric")
  expect_equal(class(res$opt.lambda1), "numeric")
  expect_equal(class(res$opt.lambda2), "numeric")
  expect_equal(class(res$lambda1s), "numeric")
  expect_equal(class(res$lambda2s), "numeric")
  expect_equal(class(res$ebic.vals), "numeric")
  expect_equal(class(res$opt.fit.lambda1), "list")
  expect_equal(class(res$opt.sparsities.lambda1), "numeric")
  expect_equal(class(res$total.variability), "numeric")
})

test_that("JoStARS throws error", {
  # Example: Highly informative prior.
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)

  # Test errors
  expect_error(JoStARS(list(dat$data, dat$data[1:n,1:(p-1)]),verbose=F)) # Different dimensions of data matrices
  expect_error(JoStARS(dat.list, ebic.gamma = -1, verbose = F)) # Negative ebic.gamma
  expect_error(JoStARS(dat.list, subsample.ratio = 1.2, verbose = F)) # Too large subsample ratio
  expect_error(JoStARS(dat.list, nCores = -1, verbose = F)) # Negative nCores
  expect_error(JoStARS(dat.list, nlambda1=1, verbose = F)) # Too few lambda1 values
  expect_error(JoStARS(dat.list, nlambda2=1, verbose = F)) # Too few lambda2 values
  expect_error(JoStARS(dat.list, lambda1.min=0.2, lambda1.max = 0.1, verbose = F)) # lambda1 max smaller than min
  expect_error(JoStARS(dat.list, lambda2.min=0.2, lambda2.max = 0.1, verbose = F)) # lambda2 max smaller than min
  expect_error(JoStARS(dat.list, var.thresh = 1.2 , verbose = F)) # Too large var.thresh
  expect_error(JoStARS(dat.list, rho=0 , verbose = F)) # Non-positive step size
  expect_error(JoStARS(dat.list, rep.num=0 , verbose = F)) # Non-positive number of subsamplings
  expect_error(JoStARS(dat.list, weights = 'blabla',verbose = F)) # Not valid weights argument
})

test_that("JoStARS default arguments", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  res <- JoStARS(dat.list,verbose = F)

  # Test defalt arguments
  expect_equal(res$opt.ebic, JoStARS(dat.list, ebic.gamma = 0.2, verbose = F)$opt.ebic) # ebic.gamma
  expect_equal(res$opt.ebic, JoStARS(dat.list, var.thresh = 0.1, verbose = F)$opt.ebic) # var.thresh
  expect_equal(res$opt.ebic, JoStARS(dat.list, scale = TRUE, verbose = F)$opt.ebic) # scaling
  expect_equal(res$opt.ebic, JoStARS(dat.list, penalize.diagonal = FALSE, verbose = F)$opt.ebic) # penalize diagonal
  expect_equal(res$opt.ebic, JoStARS(dat.list, subsample.ratio = NULL, verbose = F)$opt.ebic) # subsample ratio
  expect_equal(res$opt.ebic, JoStARS(dat.list, rep.num=20, verbose = F)$opt.ebic) # number of subsamplings
  expect_equal(res$opt.ebic, JoStARS(dat.list, retune.lambda1 = FALSE, verbose = F)$opt.ebic) # retuning lambda
  expect_equal(res$opt.ebic, JoStARS(dat.list, nlambda1 = 20, verbose = F)$opt.ebic) # subsample ratio
  expect_equal(res$opt.ebic, JoStARS(dat.list, nlambda2 = 20, verbose = F)$opt.ebic) # subsample ratio
  expect_equal(res$opt.ebic, JoStARS(dat.list, weights='equal', verbose = F)$opt.ebic) # equal weighing of classes
})

test_that("JoStARS works given optional parameters", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  expect_equal(class(JoStARS(dat.list, scale=FALSE,penalize.diagonal = TRUE, subsample.ratio=0.7,verbose = F)$opt.lambda1),'numeric')
  # One case with more observations
  n=150
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  expect_equal(class(JoStARS(dat.list, scale=FALSE,penalize.diagonal = TRUE,verbose = F)$opt.lambda1),'numeric')
})

test_that("JoStARS arguments control what they should", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  res <- JoStARS(dat.list,verbose = F)

  # Test that arguments control what they are supposed to.
  expect_true(mean(res$opt.sparsities) >= mean(JoStARS(dat.list, ebic.gamma = 1, verbose = F)$opt.sparsities)) # Test argument ebic.gamma controls sparsity
  expect_true(mean(res$opt.sparsities) >= mean(JoStARS(dat.list, var.thresh = 0.01, verbose = F)$opt.sparsities)) # Test argument var.thresh controls sparsity
  expect_output(JoStARS(dat.list, verbose = T)) # Test that verbose=T gives output
  expect_true(res$opt.ebic!= JoStARS(dat.list, verbose = FALSE, scale = FALSE)$opt.ebic) # Test that not scaling gives different results.
  expect_true(res$opt.ebic== JoStARS(dat.list, verbose = TRUE, parallelize = FALSE)$opt.ebic) # Test that parallelizing does not affect final estimate
  expect_true(any(res$total.variability != JoStARS(dat.list, retune.lambda1 = TRUE ,verbose = FALSE)$total.variability)) # Test that retuning changes the resulting variability of lambda1
})


test_that("JoStARS results are valid", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  res <- JoStARS(dat.list,verbose = F)

  # Test that the results are valid
  expect_true(det(res$opt.fit[[1]]) > 0) # Positive definite.
  expect_true(det(res$opt.fit[[2]]) > 0) # Positive definite.
  expect_true(det(res$opt.fit.lambda1[[1]]) > 0) # Positive definite.
  expect_true(det(res$opt.fit.lambda1[[2]]) > 0) # Positive definite.
  expect_true(res$opt.sparsities[1] >= 0 & res$opt.sparsities[1] <= 1) # Valid optimal sparsity.
  expect_true(res$opt.sparsities[2] >= 0 & res$opt.sparsities[2] <= 1) # Valid optimal sparsity.
  expect_true(res$opt.sparsities.lambda1[1] >= 0 & res$opt.sparsities.lambda1[1] <= 1) # Valid optimal sparsity.
  expect_true(res$opt.sparsities.lambda1[2] >= 0 & res$opt.sparsities.lambda1[2] <= 1) # Valid optimal sparsity.
  expect_true(res$opt.sparsities[2] >= 0 & res$opt.sparsities[2] <= 1) # Valid optimal sparsity.
  expect_equal(res$opt.ebic, min(res$ebic.vals)) # Check that eBIC is minimized.
  expect_true(mean(res$variability>= 0)==1) # Only positive variability values.
  expect_true(mean(res$total.variability>= 0)==1) # Only positive variability values.
  expect_true(res$opt.lambda1>= 0) # Positive lambda1.
  expect_true(res$opt.lambda2>= 0) # Positive lambda2.
  expect_true(mean(res$lambda1s>= 0)==1) # Positive lambda1.
  expect_true(mean(res$lambda2s>= 0)==1) # Positive lambda2.
})

test_that("JoStARS results have correct lengths", {
  # Example: The two data sets are identical
  set.seed(123)
  n <- 30
  p <- 5
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  dat.list = list(dat$data, dat$data)
  res <- JoStARS(dat.list,verbose = F)

  # Check that results have correct lengths.
  expect_equal(nrow(res$variability), length(res$lambda1s)) # Check length of results
  expect_equal(length(res$ebic.vals), length(res$lambda2s)) # Check length of results
  expect_equal(dim(res$opt.fit[[1]] ), c(p, p)) # Check dimension of array.
  expect_equal(dim(res$opt.fit[[2]] ), c(p, p)) # Check dimension of array.
})

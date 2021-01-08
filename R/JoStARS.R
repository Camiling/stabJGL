#' Perform JoStARS
#'
#' @description  Implements JoStARS for penalty selection in joint network reconstruction of multiple graphs. JoStARS performs penalty parameter selection in the joint graphical lasso, selecting both the sparsity- and the similarity controlling penalty parameters.
#'
#' @details The objective is to borrow strength across simialar classes to increase statistical power, while ensuring that the joint modelling may not decrease the accuracy of the resulting inferred graphs. The method takes a set of data matrices for which graphs are to be inferred with the joint graphical lasso. The method takes a list \code{Y} of \eqn{K} data matrices for which separate graphs are to be inferred, selects the sparsity controlling penalty parameter \eqn{\lambda_1} and similarity controlling penalty parameter \eqn{\lambda_2}, and performs the joint graphical lasso, resulting in \eqn{K} precision matrix estimates. To increase computational efficiency, the code can be run in parallel with \code{nCores} threads.
#'
#' @aliases jostars Jostars JoStars
#'
#' @param Y A list of \eqn{K} data matrices, each of dimension \eqn{n_k} by \eqn{p} where \eqn{n_k} is the sample size of class \eqn{K} and \eqn{p} is the dimension.
#'
#' @param scale If \code{scale=TRUE}, all variables will be scaled. Default value is \code{TRUE}.
#'
#' @param penalize.diagonal Should the diagonal elements of the precision matrices are to be penalized with \eqn{\lambda_1}? Default value is \code{FALSE}.
#'
#' @param var.thresh The variability threshold to use in the selection of \eqn{\lambda_1}. The default value is \eqn{0.1}.
#'
#' @param subsample.ratio The subsampling ratio to use when selecting \eqn{\lambda_1}. The default value is \eqn{10*\sqrt(n)/n} when \eqn{n>144} and \eqn{0.8} when \eqn{n\leq 144}, where \eqn{n} is the sample size.
#'
#' @param rep.num The number of subsamplings to use when selecting \eqn{\lambda_1}. The default value is \eqn{20}.
#'
#' @param nlambda1 The number of \eqn{\lambda_1} values to consider. The default value is \eqn{20}.
#'
#' @param lambda1.min The smallest value of \eqn{\lambda_1} to consider. The default value is \eqn{0.01}.
#'
#' @param lambda1.max The largest value of \eqn{\lambda_1} to consider. The default value is \eqn{0.01}.
#'
#' @param nlambda2 The number of \eqn{\lambda_2} values to consider. The default value is \eqn{20}.
#'
#' @param lambda2.min The smallest value of \eqn{\lambda_2} to consider. The default value is \eqn{0}.
#'
#' @param lambda2.max The largest value of \eqn{\lambda_2} to consider. The default value is \eqn{0.1}.
#'
#' @param lambda2.init The initial value to fix \eqn{\lambda_2} to when selecting \eqn{\lambda_1}. The default value is \eqn{0.01}.
#'
#' @param ebic.gamma The value of \eqn{\gamma} to use in the adapted extended BIC (adapted eBIC) selection criterion. Negative values are not valid. The default value is \eqn{0}.
#'
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#'
#' @param retune.lambda1 Should the sparsity controlling parameter \eqn{\lambda_1} be re-tuned after \eqn{\lambda_2} has been selected? The default value is \code{FALSE}.
#'
#' @param parallelize Should the code be parallelized? The default value is \code{TRUE}.
#'
#' @param nCores If \code{parallelize=TRUE}, the number of threads to initialize. The default value is 2.
#'
#' @param rho A step size parameter to use in the joint graphical lasso. Large values decrease step size. Default value is 1.
#'
#' @param weights Determines the putative sample size of each class's data in the joint graphical lasso. Allowed values; a vector with length equal to the number of classes; "\code{equal}", giving each class weight 1; "\code{sample.size}", giving each class weight corresponding to its sample size. The default value is "\code{equal}".
#'
#'
#' @return Object of class \code{"list"}. Contains the following items:
#' \describe{
#'   \item{opt.fit}{The JoStARS precision matrix estimates. A list of length \eqn{K} precision matrices, each of dimension \eqn{p} by \eqn{p}.}
#'   \item{opt.ebic}{The adapted eBIC value of the set of inferred graphs.}
#'   \item{opt.sparsities}{The sparsities of the inferred graphs. A \eqn{K} dimensional vector.}
#'   \item{opt.lambda1}{The selected value of \eqn{\lambda_1}.}
#'   \item{opt.lambda2}{The selected value of \eqn{\lambda_2}.}
#'   \item{lambda1s}{The sequence of \eqn{\lambda_1} values considered in the selection.}
#'   \item{lambda2s}{The sequence of \eqn{\lambda_2} values considered in the selection.}
#'   \item{ebic.vals}{The adapted eBIC scores of the models corresponding to the different values of \eqn{\lambda_2}.}
#'   \item{opt.fit.lambda1}{The precision matrix estimates found after selecting \eqn{\lambda_1} while \eqn{\lambda_2} is fixed to its initial value. A list of length \eqn{K} precision matrices, each of dimension \eqn{p} by \eqn{p}.}
#'   \item{opt.sparsities.lambda1}{The sparsities of the graphs found after selecting \eqn{\lambda_1} while \eqn{\lambda_2} is fixed to its initial value. A vector of length \eqn{K}.}
#'   \item{total.variability}{The total variability along the subsampling path when selecting \eqn{\lambda_1}. A vector of length \code{nlambda1}.}
#'   \item{variability}{The variability of each class along the subsampling path when selecting \eqn{\lambda_1}. A matrix of dimension \code{nlambda1} by \eqn{K}.}
#'   }
#'
#'
#' @importFrom foreach %dopar%
#'
#' @seealso \code{\link[JGL]{JGL}}
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#'
#' @examples
#'
#' # example 1: simple example where the data sets are from
#' #            from completely different distributions
#' #            and have no edges in their graph structure
#' # generate data with independent variables
#' set.seed(1234)
#' x1 <- matrix(rnorm(5 * 20), ncol = 5)
#' x2 <- matrix(rnorm(5 * 20), ncol = 5)
#' Y = list(x1,x2)
#' # perform JoStARS with variability threshold 0.1
#' res <- JoStARS(Y)
#' res$opt.fit # the list of estimated JoStARS precision matrices
#' res$opt.lambda1 # the optimal selected value of lambda1
#' res$opt.lambda2 # the optimal selected value of lambda2
#' res$opt.sparsities # the sparsity of the estimated precision matrices
#'
#' # example 2: scaling the data
#' set.seed(123)
#' res <- JoStARS(Y, scale=TRUE)
#'
#' # example 3: scale-free data where where the data sets are
#' #            from identical distributions
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' x1 = MASS::mvrnorm(n, mu=rep(0,p),Sigma=dat$sigma)
#' x2 = MASS::mvrnorm(n, mu=rep(0,p),Sigma=dat$sigma)
#' Y=list(x1, x2)
#' res <- JoStARS(Y, lambda2.max=0.3)
#' res$opt.lambda1 # the optimal selected value of lambda1
#' res$opt.lambda2 # the optimal selected value of lambda2
#' res$opt.sparsities # the sparsity of the estimated precision matrices
#' # Look at precision of inferred graphs
#' precision(abs(dat$omega) > 1e-7, res$opt.fit[[1]] != 0)
#' precision(abs(dat$omega) > 1e-7, res$opt.fit[[2]] != 0)
#'
#' # example 4: scale-free data where where the data sets are
#' #            from completely unrelated distributions
#' # Create a completely unrelated prior data set
#' set.seed(123)
#' n1 <- 80
#' n2 <- 60
#' p <- 20
#' dat1 <- huge::huge.generator(n = n1, d = p, graph = "scale-free")
#' dat2 <- huge::huge.generator(n = n2, d = p, graph = "scale-free")
#' x1 = MASS::mvrnorm(n1, mu=rep(0,p),Sigma=dat1$sigma)
#' x2 = MASS::mvrnorm(n2, mu=rep(0,p),Sigma=dat2$sigma)
#' Y = list(x1, x2)
#' res <- JoStARS(Y, scale=TRUE,lambda2.max=0.3)
#' res$opt.lambda1 # the optimal selected value of lambda1
#' res$opt.lambda2 # the optimal selected value of lambda2
#' res$opt.sparsities # the sparsity of the estimated precision matrices
#' # Look at precision of inferred graphs
#' precision(abs(dat1$omega) > 1e-7, res$opt.fit[[1]] != 0)
#' precision(abs(dat2$omega) > 1e-7, res$opt.fit[[2]] != 0)
#'
#'
JoStARS = function(Y,scale=T,penalize.diagonal=FALSE,var.thresh = 0.1, subsample.ratio = NULL,
                   rep.num = 20,  nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                   ebic.gamma=0.2,verbose=T, retune.lambda1=F,parallelize=T,nCores=2,rho=1,weights="equal"){

  if (length(unique(unlist(lapply(Y, ncol)))) !=1) {
    stop("dimension (no of variables) not equal for all data sets.")
  }
  if (nlambda1<2) {
    stop("more values of lambda1 should be considered. Try nlambda1=20. \n")
  }
  if (lambda1.min > lambda1.max) {
    stop("lambda1.min must be smaller than lambda1.max. \n")
  }
  if (nlambda2<2) {
    stop("more values of lambda2 should be considered. Try nlambda2=20. \n")
  }
  if (lambda2.min > lambda2.max) {
    stop("lambda2.min must be smaller than lambda2.max. \n")
  }
  if(!is.null(subsample.ratio)){
    if(subsample.ratio < 0 | subsample.ratio > 1) stop("subsample ratio must be between 0 and 1. \n")
  }
  if (var.thresh < 0 | var.thresh > 1) {
    stop("variability threshold should be between 0 and 1. \n")
  }
  if (ebic.gamma < 0) {
    stop("ebic.gamma cannot have a negative value. \n")
  }
  if (parallelize) {
    if(nCores < 2) stop("if method is to be run in parallel, at least two threads should be initiated. Try nCores=2. \n")
  }
  if (rho <= 0) {
    stop("step size rho must be positive. \n")
  }
  if (rep.num < 1) {
    stop("Number of subsamplings must be positive. \n")
  }

  est = list()
  if (scale) Y=lapply(Y,scale)
  # Start by selecting lambda1 while fixing lambda2 to its initial value
  est.lambda1 = JoStARS_select_lambda1(Y,rho=rho,weights=weights, penalize.diagonal=penalize.diagonal, stars.thresh=var.thresh,
                                   stars.subsample.ratio=subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                   lambda1.max=lambda1.max, lambda2=lambda2.init,verbose=verbose,parallelize=parallelize,nCores=nCores)
  # Select lambda2
  est.lambda2 = JoStARS_select_lambda2_eBIC(Y,rho=rho,weights=weights,penalize.diagonal=penalize.diagonal,
                                        nlambda2=nlambda2,lambda2.min=lambda2.min,lambda2.max=lambda2.max,
                                        lambda1=est.lambda1$opt.lambda1,gamma=ebic.gamma,verbose=verbose,parallelize=parallelize,
                                        nCores=nCores)
  est$opt.ebic = est.lambda2$opt.ebic
  est$ebic.vals = est.lambda2$ebic.vals
  est$lambda2s = est.lambda2$lambda2s
  est$opt.lambda2 = est.lambda2$opt.lambda2
  est$opt.fit = est.lambda2$opt.fit
  est$opt.sparsities = est.lambda2$opt.sparsities
  # If lambda1 should be retuned, select it while fixing lambda2
  if(retune.lambda1){
    est.lambda1 = JoStARS_select_lambda1(Y,rho=rho,weights=weights, penalize.diagonal=penalize.diagonal, stars.thresh=var.thresh,
                                     stars.subsample.ratio=subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                     lambda1.max=lambda1.max, lambda2=est.lambda2$opt.lambda2,verbose=verbose,parallelize=parallelize,nCores=nCores)
    est$opt.fit = est.lambda1$opt.fit
    est$opt.fit.lambda2 = est.lambda2$opt.fit
    est$opt.sparsities = est.lambda1$opt.sparsities
    est$opt.sparsities.lambda2 = est.lambda2$opt.sparsities
  }
  else {
    est$opt.fit.lambda1  = est.lambda1$opt.fit
    est$opt.sparsities.lambda1 = est.lambda1$opt.sparsities
  }
  est$lambda1s = est.lambda1$lambda1s
  est$total.variability = est.lambda1$total.variability
  est$variability = est.lambda1$variability
  est$opt.lambda1 = est.lambda1$opt.lambda1
  return(est)
}

#' Internal function for selecting lambda1 in JoStARS
#'
#' @keywords internal
JoStARS_select_lambda1 = function(Y,rho=1,weights="equal",penalize.diagonal=FALSE,stars.thresh = 0.1, stars.subsample.ratio = NULL,rep.num = 20,
                                  nlambda1=20,lambda1.min,lambda1.max, lambda2,verbose,parallelize=F,nCores=4){
  K = length(Y)
  n.vals = unlist(lapply(Y,nrow))
  p = ncol(Y[[1]])
  stars.subsample.ratios = rep(0,K)
  # If code is run with parallelization, to make results reproducible we must define the seeds for the drawing of each subsample
  seeds=sample(1:1000,rep.num)
  est=list()
  # Make vector of lambda1 values to consider
  lambda1s=seq(lambda1.max,lambda1.min,length.out=nlambda1)
  if(verbose) cat('Tuning lambda1... \n ')
  # Find the subsample ratio for each class
  if(is.null(stars.subsample.ratio))
  {
    for(k in 1:K){
      if(n.vals[k]>144) stars.subsample.ratios[k] = 10*sqrt(n.vals[k])/n.vals[k]
      if(n.vals[k]<=144) stars.subsample.ratios[k] = 0.8
    }
  }
  # If a single subsample ratio is provided, make a vector of the value
  else {
    if(length(stars.subsample.ratio)<=1){
      stars.subsample.ratios = rep(stars.subsample.ratio,K)
    }
  }
  est = list()
  est$merge=list()
  est$lambda1s=lambda1s
  # For each lambda, make a Kxpxp array to store the cumulative adjacency matrices in
  for(i in 1:nlambda1)
  {
    est$merge[[i]] = array(0,dim=c(K,p,p))
  }
  # Unparallelized analysis, drawing and analysing each subsample sequentially
  if(!parallelize){
    # Perform rep.num samplings
    for(i in 1:rep.num)
    {
      Y.sample = list()
      # Sample from each of the K classes
      for(k in 1:K){
        # Sample indices from class k
        ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE)
        # The data subsample from class k
        Y.sample[[k]] = Y[[k]][ind.sample, ]
      }
      # For each lambda1 value, fit a JGL model using these samples
      for(j in 1:nlambda1)
      {
        lambda = lambda1s[j]
        tmp = JGL::JGL(Y.sample,penalty='fused',rho=rho,lambda1=lambda,lambda2=lambda2,return.whole.theta=T,penalize.diagonal = penalize.diagonal)$theta
        # tmp is a K-length list of estimated precision matrices

        # Add up estimated adjacency matrices to find out how many of the graphs that agree on each edge.
        for (k in 1:K){
          est$merge[[j]][k,,] = est$merge[[j]][k,,] + (tmp[[k]]!=0)
        }
      }
      # Clear memory
      rm(ind.sample,Y.sample,tmp)
      gc()

      if (verbose) {
        # Print how far we have come when percentage is dividable by 5.
        done <- round(100 * i / rep.num)
        done.next <- round(100 * (i + 1) / rep.num)
        if (i == rep.num | (done %% 5) == 0 & (done.next %% 5) != 0) cat('Tuning lambda1: ', done, ' % done \n')
      }
    }
  }
  # Parallelized analysis, drawing and analysing each subsample using threads
  else{
    #doParallel::registerDoParallel(nCores)
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    res.list = foreach::foreach(i=1:rep.num, .packages = 'JoStARS') %dopar% {
      JoStARS_select_lambda1_parallel(Y,rep.num=rep.num,rho=rho,n.vals=n.vals,stars.subsample.ratios=stars.subsample.ratios,
                                      lambda1s=lambda1s,lambda2=lambda2,penalize.diagonal = penalize.diagonal,
                                      seed=seeds[i], array.list=est$merge)
    }
    parallel::stopCluster(cl)
    #doParallel::stopImplicitCluster()
    #foreach::registerDoSEQ()
    # Collapse results, adding up estimated adjacency matrices to find out how many of the graphs that agree on each edge.
    for(j in 1:length(lambda1s)){
      for(k in 1:K){
        est$merge[[j]][k,,] = Reduce('+',lapply(res.list,function(mat) mat[[j]][k,,]))
      }
    }
    rm(res.list)
  }

  # Estimate edge and graph stability
  est$variability = matrix(0,nlambda1,K)
  for(i in 1:nlambda1){
    for(k in 1:K){
      # Fraction of times subsample graphs agree on the different edges in class k
      est$merge[[i]][k,,] = est$merge[[i]][k,,]/rep.num
      # 2 x Bernoulli indicator for edge (i,j)
      est$variability[i,k] = 4*sum(est$merge[[i]][k,,]*(1-est$merge[[i]][k,,]))/(p*(p-1))
    }
  }
  # The total variability is the mean across all K classes
  est$total.variability = rowMeans(est$variability)
  # Select the smallest lambda that does not exceed variability threshold
  est$opt.index = max(which.max(est$total.variability >= stars.thresh)[1]-1,1)
  est$opt.lambda1 = est$lambda1s[est$opt.index]
  # Use selected lambda1 to get optimal fit
  est$opt.fit = JGL::JGL(Y,penalty='fused',rho=rho,lambda1=est$opt.lambda1,lambda2=lambda2,return.whole.theta=T,penalize.diagonal=penalize.diagonal)$theta # Optimal fit.
  est$opt.sparsities = unlist(lapply(est$opt.fit,sparsity))
  return(est)

}


#' Internal function used when selecting lambda1 with several threads in JoStARS
#'
#' @keywords internal
JoStARS_select_lambda1_parallel = function(Y,rep.num,rho,n.vals,stars.subsample.ratios,
                                           lambda1s,lambda2,penalize.diagonal,
                                           seed,array.list){
  # This function draws one subsample and performs JGL on it with each lambda1 value to consider
  set.seed(seed)
  Y.sample = list()
  K=length(Y)
  # Sample from each of the K classes
  for(k in 1:K){
    # Sample indices of observations to use from class k
    ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE)
    # The subsample from class k
    Y.sample[[k]] = Y[[k]][ind.sample, ]
  }
  # For each lambda1 value, fit a JGL model using this subsample
  for(j in 1:length(lambda1s))
  {
    lambda = lambda1s[j]
    tmp = JGL::JGL(Y.sample,penalty='fused',rho=rho,lambda1=lambda,lambda2=lambda2,return.whole.theta=T,penalize.diagonal = penalize.diagonal)$theta

    # Save adjacency matrix
    for (k in 1:K){
      array.list[[j]][k,,] = (tmp[[k]]!=0)
    }
  }
  # Clear memory
  rm(ind.sample,Y.sample,tmp)
  gc()
  return(array.list)
}

#' Internal function for selecting lambda2 in JoStARS
#'
#'
#' @keywords internal
JoStARS_select_lambda2_eBIC = function(Y,rho=1,weights="equal",penalize.diagonal=FALSE,
                                       nlambda2=30,lambda2.min,lambda2.max, lambda1=NULL,gamma=NULL,verbose=F,
                                       parallelize=F,nCores=4){
  # Select lambda_2 by eBIC
  ebic.vals = rep(0,nlambda2)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out= nlambda2)
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,FUN = function(s) stats::cov(s))
  mods.lam2=list()
  if(verbose) cat('Tuning lambda2...\n')
  # For each lambda2 value, fit a JGL model
  for (i in 1:nlambda2){
    mods.lam2[[i]] = JGL::JGL(Y,penalty='fused',rho=rho,lambda1=lambda1,lambda2 = lambda2.vals[i],penalize.diagonal = penalize.diagonal,
                              return.whole.theta = T,weights=weights)$theta
    # Find the eBIC score for the resultign precision matrix estimate.
    ebic.vals[i] = eBIC_adapted(mods.lam2[[i]],sample.cov=sample.cov,n.vals=n.vals,gamma=gamma)

    if (verbose) {
      # Print how far we have come when percentage is dividable by 5.
      done <- round(100 * i / nlambda2)
      done.next <- round(100 * (i + 1) / nlambda2)
      if (i == nlambda2| (done %% 5) == 0 & (done.next %% 5) != 0) cat('Tuning lambda2: ', done, ' % done \n')
    }
  }
  opt.ind = which.min(ebic.vals)
  # Resturn list where mod.opt is the optimal JGL object
  # res$mod.opt$theta gives the list of the K optimal precision matrices
  res=list(opt.fit=mods.lam2[[opt.ind]],opt.lambda2 = lambda2.vals[opt.ind], opt.index = opt.ind,opt.ebic=ebic.vals[opt.ind],
           ebic.vals=ebic.vals, lambda2s=lambda2.vals, mods=mods.lam2, opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
}

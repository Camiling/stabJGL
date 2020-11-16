JoStARS = function(Y,penalty="fused",rho=1,weights="equal",scale=T,penalize.diagonal=FALSE,stars.thresh = 0.1, stars.subsample.ratio = NULL,
                   rep.num = 20,  nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                   ebic.gamma=0.2,verbose=T, retune.lambda1=F,parallelize=F,nCores=4){
  # Function for tuning both lambda1 and lambda2 in JGL using a modified version of StARS for several graphs - simultaneous sparsity tuning
  # Y : A list of length K with the observed n_k x p data matrices.
  # rho: Step size parameter
  # stars.subsample.ratio: The fraction of all observations we sample in each subsampling
  # rep.num : The number of subsamplings 
  
  # Start by tuning lambda1, fixing lambda2
  est = list()
  if (scale) Y=lapply(Y,scale)
  est.lambda1 = JoStARS_select_lambda1(Y,penalty=penalty,rho=rho,weights=weights, penalize.diagonal=penalize.diagonal, stars.thresh=stars.thresh, 
                                   stars.subsample.ratio=stars.subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                   lambda1.max=lambda1.max, lambda2=lambda2.init,verbose=verbose,parallelize=parallelize,nCores=nCores)
  est.lambda2 = JoStARS_select_lambda2_eBIC(Y,penalty=penalty,rho=rho,weights=weights,penalize.diagonal=penalize.diagonal,
                                        nlambda2=nlambda2,lambda2.min=lambda2.min,lambda2.max=lambda2.max, 
                                        lambda1=est.lambda1$opt.lambda1,gamma=ebic.gamma,verbose=verbose,parallelize=parallelize,
                                        nCores=nCores)
  est$opt.ebic = est.lambda2$opt.ebic
  est$ebic.vals = est.lambda2$ebic.vals
  est$lambda2s = est.lambda2$lambda2s
  est$opt.index.2 = est.lambda2$opt.index
  est$opt.lambda2 = est.lambda2$opt.lambda2
  est$opt.fit = est.lambda2$opt.fit # The final fit
  est$opt.sparsities = est.lambda2$opt.sparsities # Final optimal sparsity
  if(retune.lambda1){
    est.lambda1 = JoStARS_select_lambda1(Y,penalty=penalty,rho=rho,weights=weights, penalize.diagonal=penalize.diagonal, stars.thresh=stars.thresh, 
                                     stars.subsample.ratio=stars.subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                     lambda1.max=lambda1.max, lambda2=est.lambda2$opt.lambda2,verbose=verbose,parallelize=parallelize,nCores=nCores)
    est$opt.fit = est.lambda1$opt.fit # The final fit
    est$opt.fit.lambda2 = est.lambda2$opt.fit # The final fit
    est$opt.sparsities = est.lambda1$opt.sparsities # Final optimal sparsity
    est$opt.sparsities.lambda2 = est.lambda2$opt.sparsities # Final optimal sparsity
      }
  else {
    est$opt.fit.lambda1  = est.lambda1$opt.fit
    est$opt.sparsities.lambda1 = est.lambda1$opt.sparsities # Optimal sparsity after initial tuning. 
  }
  est$lambda1s = est.lambda1$lambda1s
  est$opt.index.1 = est.lambda1$opt.index
  est$total.variability.1 = est.lambda1$total.variability
  est$variability.1 = est.lambda1$variability
  est$opt.lambda1 = est.lambda1$opt.lambda1

  return(est)
}

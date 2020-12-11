JoStARS_select_lambda2_eBIC = function(Y,penalty="fused",rho=1,weights="equal",penalize.diagonal=FALSE,
                                  nlambda2=30,lambda2.min,lambda2.max, lambda1=NULL,gamma=ebic.gamma,verbose,
                                  parallelize=F,nCores=4){
  # Select lambda_2 by eBIC
  ebic.vals = rep(0,nlambda2)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out= nlambda2)
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,cov)
  mods.lam2=list()
  if(verbose) cat('Tuning lambda2...\n')
  for (i in 1:nlambda2){
    mods.lam2[[i]] = JGL(Y,penalty=penalty,rho=rho,lambda1=lambda1,lambda2 = lambda2.vals[i],penalize.diagonal = penalize.diagonal,
                         return.whole.theta = T,weights=weights)$theta
    ebic.vals[i] = eBIC_adapted(mods.lam2[[i]],sample.cov=sample.cov,n.vals=n.vals,gamma=gamma)
    if(verbose){
      done=round(100*i/nlambda2) # Print how far we have come when percentage is dividable by 5.
      if((done%%5) == 0) cat('Tuning lambda2: ', done, ' % done \n')      
    }
    
    flush.console()
  }
  opt.ind = which.min(ebic.vals)
  # Resturns optimal jgl object as mod.opt.
  # res$mod.opt$theta gives list of optimal thetas.  
  res=list(opt.fit=mods.lam2[[opt.ind]],opt.lambda2 = lambda2.vals[opt.ind], opt.index = opt.ind,opt.ebic=ebic.vals[opt.ind],
           ebic.vals=ebic.vals, lambda2s=lambda2.vals, mods=mods.lam2, opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
}

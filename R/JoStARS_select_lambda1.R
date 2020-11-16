JoStARS_select_lambda1 = function(Y,penalty="fused",rho=1,weights="equal",penalize.diagonal=FALSE,stars.thresh = 0.1, stars.subsample.ratio = NULL,rep.num = 20,
                              nlambda1=20,lambda1.min,lambda1.max, lambda2,verbose,parallelize=F,nCores=4,seeds=NULL){
  # Tune lambda1
  # lambda2: value to fix lambda2 to. 
  K = length(Y) 
  n.vals = unlist(lapply(Y,nrow))
  p = ncol(Y[[1]]) # p is the same for all
  stars.subsample.ratios = rep(0,K) # Note that unless provided, the subsample ratio will vary slightly between classes. 
  if(is.null(seeds)) seeds=sample(1:1000,rep.num)
  est=list() 
  lambda1s=seq(lambda1.max,lambda1.min,length.out=nlambda1) # From largest to smallest. 
  if(verbose) cat('Tuning lambda1... \n ')
  # Find the subsample ratio for each class
  if(is.null(stars.subsample.ratio))
  {
    for(k in 1:K){
      if(n.vals[k]>144) stars.subsample.ratios[k] = 10*sqrt(n.vals[k])/n.vals[k]
      if(n.vals[k]<=144) stars.subsample.ratios[k] = 0.8
    }
  }
  
  est = list()
  est$merge=list()
  est$lambda1s=lambda1s
  # For each lambda, make a Kxpxp array to store the cumulative adacency matrices in
  for(i in 1:nlambda1) 
  { 
    est$merge[[i]] = array(0,dim=c(K,p,p)) 
  }
  if(!parallelize){
    # Perform rep.num samplings
    for(i in 1:rep.num)
    { 
      Y.sample = list() 
      # Sample from each of the K classes
      for(k in 1:K){
        ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE) # Sample indices from class k
        Y.sample[[k]] = Y[[k]][ind.sample, ]  # The subsample from class k
      }
      # For each lambda, fit a JGL model using these samples
      for(j in 1:nlambda1)
      { 
        lambda = lambda1s[j]
        tmp = JGL(Y.sample,penalty,rho=rho,lambda1=lambda,lambda2=lambda2,return.whole.theta=T,penalize.diagonal = penalize.diagonal)$theta# A K-length list of prec. matrices
        
        for (k in 1:K){
          est$merge[[j]][k,,] = est$merge[[j]][k,,] + (tmp[[k]]!=0) # add to find how many graphs agree on each edge. 
        }
      }
      rm(ind.sample,Y.sample,tmp)
      gc() # Clear memory
      
      if(verbose){
        done=round(100*i/rep.num) # Print how far we have come when percentage is dividable by 5.
        if((done%%5) == 0) cat('Tuning lambda1: ', done, ' % done \n')      
      }
      flush.console()
    }
  }
  # Paralellized
  else{
    registerDoParallel(nCores)
    res.list = foreach (i=1:rep.num) %dopar% {
      JoStARS_select_lambda1_parallel(Y,rep.num=rep.num,rho=rho,n.vals=n.vals,stars.subsample.ratios=stars.subsample.ratios,
                                  penalty=penalty,lambda1s=lambda1s,lambda2=lambda2,penalize.diagonal = penalize.diagonal,
                                  seed=seeds[i], array.list=est$merge)
    }
    registerDoSEQ()
    # Collapse results. 
    for(j in 1:length(lambda1s)){
      for(k in 1:K){
        est$merge[[j]][k,,] = Reduce('+',lapply(res.list,function(mat) mat[[j]][k,,])) # sum up
      }
    }
    rm(res.list)
  }
  
  # Estimate edge and graph stability
  est$variability = matrix(0,nlambda1,K)
  for(i in 1:nlambda1){
    for(k in 1:K){
      est$merge[[i]][k,,] = est$merge[[i]][k,,]/rep.num # Frac of times edges agree in graphs of class k
      est$variability[i,k] = 4*sum(est$merge[[i]][k,,]*(1-est$merge[[i]][k,,]))/(p*(p-1)) # 2x Bernoulli indicator for edge (i,j)
    }
  }
  est$total.variability = rowMeans(est$variability) # Mean across all K classes. Length nlambda
  est$opt.index = max(which.max(est$total.variability >= stars.thresh)[1]-1,1) # Smallest lambda that does not exceed variability threshold
  est$opt.lambda1 = est$lambda1s[est$opt.index]
  est$opt.fit = JGL(Y,penalty,rho=rho,lambda1=est$opt.lambda1,lambda2=lambda2,return.whole.theta=T,penalize.diagonal=penalize.diagonal)$theta # Optimal fit.
  est$opt.sparsities = unlist(lapply(est$opt.fit,sparsity))
  return(est)
  
}

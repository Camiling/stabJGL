JoStARS_select_lambda1_parallel = function(Y,rep.num,rho,n.vals,stars.subsample.ratios,
                            penalty,lambda1s,lambda2,penalize.diagonal,
                            seed,array.list)
  {
    set.seed(seed)
    Y.sample = list() 
    K=length(Y)
    # Sample from each of the K classes
    for(k in 1:K){
      ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE) # Sample indices from class k
      Y.sample[[k]] = Y[[k]][ind.sample, ]  # The subsample from class k
    }
    # For each lambda, fit a JGL model using these samples
    for(j in 1:length(lambda1s))
    { 
      lambda = lambda1s[j]
      tmp = JGL(Y.sample,penalty,rho=rho,lambda1=lambda,lambda2=lambda2,return.whole.theta=T,penalize.diagonal = penalize.diagonal)$theta# A K-length list of prec. matrices
  
      for (k in 1:K){
        array.list[[j]][k,,] = (tmp[[k]]!=0)
      }
    }
    rm(ind.sample,Y.sample,tmp)
    gc() # Clear memory
    return(array.list)
}

# JoStARS simulation functiona ----------------------------------------

JoStARS_simulation = function(K,n.vals,p,N,frac.disagreement,ebic.gamma=0,scale=T,penalize.diagonal=FALSE,stars.thresh = 0.1, stars.subsample.ratio = NULL,
                              rep.num = 20, nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                              verbose=T,seeds=sample(1:1000,N),nCores=3, method='symmetric'){
  # Function for performing a simulation study assessing the performance of JoStARS, comparing its performance to that of the AIC and glasso. 
  # K: The number of classes to consider. 
  # n.vals: A vector of the K different n values. 
  # p: The number of variables to consider. 
  # N: The number of simulations to perform. 
  # frac.disagreement: The fraction of disagreement we want the graphs to have. 
  # ebic.gamma: Value of gamma to use in the eBIC. 
  # method: How should the similarity between the prec matrices be? Symmetric by default, meaning all are equally different. If not, one will stand out as less similar.  
  
  res=list()
  # JoStARS results
  res$opt.lambda1s = rep(0,N)
  res$opt.lambda2s = rep(0,N)
  res$opt.sparsities = matrix(0,N,K)
  res$precisions =  matrix(0,N,K)
  res$specificities =  matrix(0,N,K)
  res$recalls =  matrix(0,N,K)
  res$matrix.distances =  matrix(0,N,K)
  # AIC results
  res$opt.lambda1s.aic = rep(0,N)
  res$opt.lambda2s.aic = rep(0,N)
  res$opt.sparsities.aic = matrix(0,N,K)
  res$precisions.aic =  matrix(0,N,K)
  res$specificities.aic =  matrix(0,N,K)
  res$recalls.aic =  matrix(0,N,K)
  res$matrix.distances.aic =  matrix(0,N,K)
  # Glasso results
  res$opt.lambdas.glasso = matrix(0,N,K) # Must be a matrix
  res$opt.sparsities.glasso = matrix(0,N,K)
  res$precisions.glasso =  matrix(0,N,K)
  res$specificities.glasso =  matrix(0,N,K)
  res$recalls.glasso =  matrix(0,N,K)
  res$matrix.distances.glasso =  matrix(0,N,K)
  
  # Start by generating the precision matrices 
  cov.matrices = list()
  prec.matrices = list()
  # Start by generating the first prec matrix
  huge.init = huge.generator(n.vals[1],p,graph='scale-free',verbose = F)
  theta.init = huge.init$omega
  theta.init[which(abs(theta.init)<1e-5)] = 0 
  spars.init = huge.init$sparsity
  cov.matrices[[1]] = huge.init$sigma
  prec.matrices[[1]] = theta.init
  if(method=='symmetric'){
    for(k in 2:K){
      huge.tmp = mutate.graph(huge.init,frac.disagreement)
      cov.matrices[[k]] = huge.tmp$cov.mat
      prec.matrices[[k]] = huge.tmp$prec.mat
    }
  }
  else{ # First K-1 graphs are similar
    for(k in 2:(K-1)){
      huge.tmp = mutate.graph(huge.init,frac.disagreement)
      cov.matrices[[k]] = huge.tmp$cov.mat
      prec.matrices[[k]] = huge.tmp$prec.mat
    }
    # Last graphs completely different
    huge.tmp = mutate.graph(huge.init,fraction =1)
    cov.matrices[[K]] = huge.tmp$cov.mat
    prec.matrices[[K]] = huge.tmp$prec.mat
  }
  registerDoParallel(nCores)
  res.list = foreach (i=1:N) %dopar% {
    JoStARS_simulation_one_iteration(n.vals=n.vals,cov.matrices=cov.matrices,prec.matrices=prec.matrices,stars.thresh = stars.thresh,
                                     stars.subsample.ratio = stars.subsample.ratio,rep.num = rep.num,scale=scale,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                     lambda1.max=lambda1.max, nlambda2=nlambda2,lambda2.min=lambda2.min,lambda2.max = lambda2.max, lambda2.init = lambda2.init, 
                                     ebic.gamma=ebic.gamma,penalize.diagonal=penalize.diagonal,seed=seeds[i]);
  }
  registerDoSEQ()
  for(i in 1:N){
    est.tmp = res.list[[i]]
    # Results from JoStARS
    res$opt.lambda1s[i] = est.tmp$opt.lambda1 
    res$opt.lambda2s[i] = est.tmp$opt.lambda2 
    res$opt.sparsities[i,] = est.tmp$opt.sparsities 
    res$matrix.distances[i,] = est.tmp$matrix.distances
    res$precisions[i,] = est.tmp$precisions
    res$recalls[i,] = est.tmp$recalls
    res$specificities[i,] =  est.tmp$specificities
    # Results from AIC
    res$opt.lambda1s.aic[i] = est.tmp$opt.lambda1.aic
    res$opt.lambda2s.aic[i] = est.tmp$opt.lambda2.aic 
    res$opt.sparsities.aic[i,] = est.tmp$opt.sparsities.aic
    res$matrix.distances.aic[i,] = est.tmp$matrix.distances.aic
    res$precisions.aic[i,] = est.tmp$precisions.aic
    res$recalls.aic[i,] = est.tmp$recalls.aic
    res$specificities.aic[i,] =  est.tmp$specificities.aic
    # Results from glasso
    res$matrix.distances.glasso[i,] = est.tmp$matrix.distances.glasso
    res$precisions.glasso[i,] = est.tmp$precisions.glasso
    res$recalls.glasso[i,] = est.tmp$recalls.glasso
    res$specificities.glasso[i,] =  est.tmp$specificities.glasso
    res$opt.sparsities.glasso[i,] = est.tmp$opt.sparsities.glasso
    res$opt.lambdas.glasso[i,] = est.tmp$opt.lambdas.glasso
  }
  # Mean results from JoStARS
  res$mean.opt.lambda1= mean(res$opt.lambda1s)
  res$mean.opt.lambda2= mean(res$opt.lambda2s)
  res$mean.opt.sparsities= colMeans(res$opt.sparsities)
  res$mean.precisions= colMeans(res$precisions)
  res$mean.recalls = colMeans(res$recalls)
  res$mean.specificities =  colMeans(res$specificities)
  res$mean.matrix.distances = colMeans(res$matrix.distances)
  
  # Mean results from AIC
  res$mean.opt.lambda1.aic= mean(res$opt.lambda1s.aic)
  res$mean.opt.lambda2.aic = mean(res$opt.lambda2s.aic)
  res$mean.opt.sparsities.aic = colMeans(res$opt.sparsities.aic)
  res$mean.precisions.aic = colMeans(res$precisions.aic)
  res$mean.recalls.aic = colMeans(res$recalls.aic)
  res$mean.specificities.aic =  colMeans(res$specificities.aic)
  res$mean.matrix.distances.aic = colMeans(res$matrix.distances.aic)
  
  # Mean results from glasso
  res$mean.opt.lambdas.glasso= colMeans(res$opt.lambdas.glasso)
  res$mean.opt.sparsities.glasso= colMeans(res$opt.sparsities.glasso)
  res$mean.precisions.glasso= colMeans(res$precisions.glasso)
  res$mean.recalls.glasso = colMeans(res$recalls.glasso)
  res$mean.specificities.glasso =  colMeans(res$specificities.glasso)
  res$mean.matrix.distances.glasso = colMeans(res$matrix.distances.glasso)
  
  res$true.sparsity = spars.init
  return(res)
}

# Function for performing one iteration -----------------------------------------

JoStARS_simulation_one_iteration = function(n.vals,cov.matrices,prec.matrices,stars.thresh,stars.subsample.ratio,rep.num,scale,nlambda1,lambda1.min,lambda1.max, 
                                            nlambda2,lambda2.min,lambda2.max, lambda2.init,ebic.gamma,penalize.diagonal,seed) {
  y = list()
  K=length(n.vals)
  p=ncol(prec.matrices[[1]])
  glasso.res=list()
  set.seed(seed)
  # Generate data. 
  for(k in 1:K){
    y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
    if (scale) y[[k]]=scale(y[[k]])
    # Use the ordinary graphical lasso for comparison
    glasso.tmp = huge(y[[k]],method='glasso',verbose = F)
    glasso.res[[k]] = huge.select(glasso.tmp,criterion='stars', stars.thresh=stars.thresh,verbose = F)
  }
  jgl.tmp = JoStARS(Y=y,penalty='fused',stars.thresh = stars.thresh,stars.subsample.ratio = stars.subsample.ratio,rep.num = rep.num,scale=F,nlambda1=nlambda1,
                    lambda1.min=lambda1.min,lambda1.max=lambda1.max, nlambda2=nlambda2,lambda2.min=lambda2.min,lambda2.max = lambda2.max, lambda2.init = lambda2.init,
                    ebic.gamma=ebic.gamma,verbose=F,penalize.diagonal=penalize.diagonal)
  jgl.aic.tmp = JGL_select_AIC(Y=y,penalty='fused',nlambda1=nlambda1,lambda1.min=lambda1.min,lambda1.max=lambda1.max,nlambda2=nlambda2,lambda2.min=lambda2.min,
                               lambda2.max=lambda2.max,lambda2.init = lambda2.init,penalize.diagonal=penalize.diagonal)
  est=list()
  # Results from JoStARS
  est$opt.lambda1 = jgl.tmp$opt.lambda1
  est$opt.lambda2 = jgl.tmp$opt.lambda2 
  est$opt.sparsities  = jgl.tmp$opt.sparsities
  est$matrix.distances = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], jgl.tmp$opt.fit[[k]]))
  est$precisions =  sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
  est$recalls =  sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
  est$specificities =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
  
  # Results from AIC
  est$opt.lambda1.aic = jgl.aic.tmp$opt.lambda1
  est$opt.lambda2.aic = jgl.aic.tmp$opt.lambda2 
  est$opt.sparsities.aic = jgl.aic.tmp$opt.sparsities
  est$matrix.distances.aic = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], jgl.aic.tmp$opt.fit[[k]]))
  est$precisions.aic = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jgl.aic.tmp$opt.fit[[k]]!=0))
  est$recalls.aic = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jgl.aic.tmp$opt.fit[[k]]!=0))
  est$specificities.aic =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, jgl.aic.tmp$opt.fit[[k]]!=0))
  
  # Results from glasso
  est$matrix.distances.glasso = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], glasso.res[[k]]$opt.icov))
  est$precisions.glasso = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$recalls.glasso = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$specificities.glasso =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$opt.sparsities.glasso = unlist(lapply(glasso.res,FUN=function(l) l$opt.sparsity))
  est$opt.lambdas.glasso = unlist(lapply(glasso.res,FUN=function(l) l$opt.lambda))
  return(est)
}

# Adapted AIC crierion (Danaher et al.) --------------------------------------------------

JGL_select_AIC = function(Y,penalty='fused',nlambda1,lambda1.min,lambda1.max,nlambda2,lambda2.min,lambda2.max, 
                          lambda2.init,penalize.diagonal){
  K=length(Y)
  p=ncol(Y[[1]])
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,cov)
  lambda1.vals = seq(lambda1.min,lambda1.max,length.out=nlambda1)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out=nlambda2)
  mods.lam1=list()
  aic.lam1=rep(0,length(lambda1.vals)) 
  for (i in 1:length(lambda1.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1.vals[i],lambda2 = lambda2.init,return.whole.theta = T,penalize.diagonal=penalize.diagonal)$theta
    mods.lam1[[i]] = mod.temp
    aic.lam1[i] = AIC.adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind.lam1 = which.min(aic.lam1)
  lambda1=lambda1.vals[opt.ind.lam1]
  
  mods.lam2=list()
  aic.lam2=rep(0,length(lambda2.vals))
  for (i in 1:length(lambda2.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1,lambda2 = lambda2.vals[i],penalize.diagonal = penalize.diagonal,
                   return.whole.theta = T)$theta
    mods.lam2[[i]] = mod.temp
    aic.lam2[i] = AIC.adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind = which.min(aic.lam2) 
  res=list(opt.fit=mods.lam2[[opt.ind]],opt.lambda1 = lambda1,opt.lambda2 = lambda2.vals[opt.ind],
           opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
  return(res)
}


print_results_JGL = function(obj.list,fracs.mutated,show.distance=F, show.sd = F, show.specificity=F){
  if(show.sd==T){
    print_results_JGL_show_SD(obj.list,fracs.mutated,show.distance, show.specificity)
    return()
  }
  # obj is a list of objects returned by a JoStARS_simulation.
  # show distance: should the matrix distance be printed?
  # Function for printing mean lambda2, lambda2, sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  K=length(obj.list[[1]]$mean.opt.sparsities)
  for (i in 1:length(obj.list)){
    obj=obj.list[[i]]
    cat(fracs.mutated[i],' & Glasso  & ',  round(mean(obj$mean.opt.lambdas.glasso),3), ' & - ')
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities.glasso[k],3), '[',paste(round(quantile(obj$opt.sparsities.glasso[,k],probs=c(.025,.975)),3),collapse=','),']',
          ' & ',
          round(obj$mean.precisions.glasso[k],2), '[',paste(round(quantile(obj$precisions.glasso[,k],probs=c(.025,.975)),2),collapse=','),'] & ',
          round(obj$mean.recalls.glasso[k],2), '[',paste(round(quantile(obj$recalls.glasso[,k],probs=c(.025,.975)),2),collapse=','),']')
      if(show.specificity)cat('&',round(obj$mean.specificities.glasso[k],2), '[',paste(round(quantile(obj$specificities.glasso[,k],probs=c(.025,.975)),2),collapse=','),']')
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.glasso[k],3))
    }    
    cat(' \\\\ \n')
    cat('  & JGL (AIC)  & ',  round(obj$mean.opt.lambda1.aic,3),' & ', round(obj$mean.opt.lambda2.aic,3))  
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities.aic[k],3), '[',paste(round(quantile(obj$opt.sparsities.aic[,k],probs=c(.025,.975)),3),collapse=','),'] & ',
          round(obj$mean.precisions.aic[k],2),'[',paste(round(quantile(obj$precisions.aic[,k],probs=c(.025,.975)),2),collapse=','),'] & ',
          round(obj$mean.recalls.aic[k],2),'[',paste(round(quantile(obj$recalls.aic[,k],probs=c(.025,.975)),2),collapse=','),']')
      if(show.specificity) cat('&',round(obj$mean.specificities.aic[k],2), '[',paste(round(quantile(obj$specificities.aic[,k],probs=c(.025,.975)),2),collapse=','),']') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.aic[k],3))
    }  
    cat(' \\\\ \n')
    cat('  & JoStARS  & ',  round(obj$mean.opt.lambda1,3),' & ', round(obj$mean.opt.lambda2,3))
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities[k],3), '[',paste(round(quantile(obj$opt.sparsities[,k],probs=c(.025,.975)),3),collapse=','),'] & ',
          round(obj$mean.precisions[k],2),'[',paste(round(quantile(obj$precisions[,k],probs=c(.025,.975)),2),collapse=','),'] & ',
          round(obj$mean.recalls[k],2),'[',paste(round(quantile(obj$recalls[,k],probs=c(.025,.975)),2),collapse=','),']')
      if(show.specificity)cat('&',round(obj$mean.specificities[k],2),'[',paste(round(quantile(obj$specificities[,k],probs=c(.025,.975)),2),collapse=','),']') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances[k],3))
    }  
    cat(' \\\\ \n \\hline \n')
  }
}

print_results_JGL_show_SD = function(obj.list,fracs.mutated,show.distance=F,show.specificity=F){
  # obj is a list of objects returned by a JoStARS_simulation.
  # show distance: should the matrix distance be printed?
  # Function for printing mean lambda2, lambda2, sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  K=length(obj.list[[1]]$mean.opt.sparsities)
  for (i in 1:length(obj.list)){
    obj=obj.list[[i]]
    cat(fracs.mutated[i],' & Glasso  & ',  round(mean(obj$mean.opt.lambdas.glasso),3),'(',round(sd(obj$mean.opt.lambdas.glasso),3),')', ' & - ')
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities.glasso[k],3), '(',round(sd(obj$opt.sparsities.glasso[,k]),3),')',' & ',
          round(obj$mean.precisions.glasso[k],2),'(',round(sd(obj$precisions.glasso[,k]),2),')',' & ',
          round(obj$mean.recalls.glasso[k],2), '(',round(sd(obj$recalls.glasso[,k]),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities.glasso[k],2), '(',round(sd(obj$specificities.glasso[,k]),2),')')
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.glasso[k],3))
    }    
    cat(' \\\\ \n')
    cat('  & JGL (AIC)  & ',  round(obj$mean.opt.lambda1.aic,3), '(',round(sd(obj$opt.lambda1s.aic),3),')',' & ',
        round(obj$mean.opt.lambda2.aic,3),'(',round(sd(obj$opt.lambda2s.aic),3),')')  
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities.aic[k],3), '(',round(sd(obj$opt.sparsities.aic[,k]),3),')',' & ',
          round(obj$mean.precisions.aic[k],2),'(',round(sd(obj$precisions.aic[,k]),2),')',' & ',
          round(obj$mean.recalls.aic[k],2),'(',round(sd(obj$recalls.aic[,k]),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities.aic[k],2), '(',round(sd(obj$specificities.aic[,k]),2),')') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.aic[k],3))
    }  
    cat(' \\\\ \n')
    cat('  & JoStARS  & ',  round(obj$mean.opt.lambda1,3),'(',round(sd(obj$opt.lambda1s),3),')', ' & ',
        round(obj$mean.opt.lambda2,3),'(',round(sd(obj$opt.lambda2s),3),')')
    for(k in 1:K){
      cat(' && ',round(obj$mean.opt.sparsities[k],3), '(',round(sd(obj$opt.sparsities[,k]),3),')',' & ',
          round(obj$mean.precisions[k],2),'(',round(sd(obj$precisions[,k]),2),')',' & ',
          round(obj$mean.recalls[k],2),'(',round(sd(obj$recalls[,k]),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities[k],2), '(',round(sd(obj$specificities[,k]),2),')') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances[k],3))
    }  
    cat(' \\\\ \n \\hline \n')
  }
}







AIC.adapted = function(theta, sample.cov,n.vals) {
  # theta is the list of estimated prec matrices
  # sample cov the list of sample cov matrices
  p = dim(sample.cov[[1]])[2]
  K = length(sample.cov)
  aic.vals = rep(0,K)
  for(k in 1:length(theta))
  {
    theta.temp=theta[[k]]
    diag(theta.temp) = rep(0,p) # Do not count diagonal
    d=sum(theta.temp!=0)/2 # Number of edges
    aic.vals[k]  =  2*d - n.vals[k]*log(det(theta[[k]])) + n.vals[k]*sum(diag(sample.cov[[k]]%*%theta[[k]]))
  }
  return(sum(aic.vals))
}

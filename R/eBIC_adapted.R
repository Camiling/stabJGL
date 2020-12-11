eBIC_adapted = function(theta, sample.cov,n.vals,gamma=0.2) {
  # theta: list of estimated precision matrices
  # sample cov: list of sample cov matrices
  # gamma: additional edge penalty parameter
  p = dim(sample.cov[[1]])[2]
  K = length(sample.cov)
  ebic.vals = rep(0,K)
  for(k in 1:length(theta))
  {
    theta.temp=theta[[k]]
    diag(theta.temp) = rep(0,p) # Do not count diagonal
    d=sum(theta.temp!=0)/2 # Number of edges
    ebic.vals[k]=d*log(n.vals[k])+4*d*gamma*log(p)-n.vals[k]*log(det(theta[[k]]))+n.vals[k]*sum(diag(sample.cov[[k]]%*%theta[[k]]))
  }
  return(sum(ebic.vals))
}
